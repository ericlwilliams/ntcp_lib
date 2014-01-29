classdef classDVHAnalysis
    properties
        DVHorg % original DVHs
        Beta2Alpha % beta to alpha ration in the LQ model
        FxNum % fraction number
        DateSurvive % the survival dates (months)
        flgCensor % censorship
        DoseStep % for (di,vj) analysis, the dose step size
        VolStep % for (di,jv) analysis, the volume step size
        DoseBins % in case the bins are not regular, it can be overwriten
        VolBins % in case the volumes are not regular, it can be overwriten
        CoxMinSize % minimim size in the cox model
        LogRankMinSize % minimum group size in the log rank test
        
        DVH_LQ % corrected DVHs by Linear quardratic model

        CoxStatsFx % statistic results of cox model for Fx
        CoxStatsDx % statistic results of cox model for Dx
        CoxStatsVx % statistic results of cox model for Vx
        
        pDxmat % p-value matrix of the pair (di,vj) in terms of Dx
        flgCorDx % flag for negative correlation
        pDxdata % details of results (n1,censor1, n2,censor2)
        pVxmat % p-value matrix of (di,vj) in terms of Vx
        flgCorVx % flag for negative correlation
        pVxdata % details of results (n1,censor1, n2,cnesor2)
    end

    methods
        function DAobj = classDVHAnalysis()
            DAobj.CoxMinSize=2;
            DAobj.LogRankMinSize=2;
        end
        function DAobj = CalculateDoseBins(DAobj)
            dv=cellfun(@(x) x(end,1), DAobj.DVH_LQ(:,2)); dmax=max(dv);
            DAobj.DoseBins=0:DAobj.DoseStep:dmax;
        end
        function DAobj = CalculateVolBins(DAobj)
            f=cellfun(@(x) x(1,3),DAobj.DVH_LQ(:,2)); vmax=max(f); % maximum volume
            DAobj.VolBins=0:DAobj.VolStep:vmax;
        end
        function DAobj = LinearQuardraticCorrection(DAobj)
            DAobj.DVH_LQ=DAobj.DVHorg;
            if DAobj.Beta2Alpha>0
                for m=1:length(DAobj.DVH_LQ)
                    DAobj.DVH_LQ{m,2}(:,1) = DAobj.DVH_LQ{m,2}(:,1) .* ( 1+DAobj.Beta2Alpha/DAobj.FxNum(m)*(DAobj.DVH_LQ{m,2}(:,1)));
                    DAobj.DVH_LQ{m,3}(:,1) = DAobj.DVH_LQ{m,3}(:,1) .* ( 1+DAobj.Beta2Alpha/DAobj.FxNum(m)*(DAobj.DVH_LQ{m,3}(:,1)));
                end
            end
        end
        function DAobj = CoxModel(DAobj)
            warning('off','MATLAB:singularMatrix');
            warning('off','stats:coxphfit:FitWarning');
            % dmax
                dv=cellfun(@(x) x(end,1), DAobj.DVH_LQ(:,2)); % maximum physical doses
                [~,logl,h,stats]=coxphfit(dv,DAobj.DateSurvive,'baseline',0,'censoring',DAobj.flgCensor);
                stats.logl=logl; stats.h=h;
                ptnum=length(dv); % number of patients
            % fx
                fxnum=unique(DAobj.FxNum); % all diferent fractions
                flgfx = length(fxnum)>1; % if there is more than one fraction numbers, do the cox model
                if flgfx
                    [~,logl,h,stats]=coxphfit(DAobj.FxNum,DAobj.DateSurvive,'baseline',0,'censoring',DAobj.flgCensor);
                    stats.logl=logl; stats.h=h;
                else
                    stats.beta=-inf; stats.logl=-inf; stats.h=-inf;
                end
                DAobj.CoxStatsFx=stats;
                
            % Dx
                volnum=length(DAobj.VolBins);
                DAobj.CoxStatsDx=repmat(stats,[volnum,2]); % 2 columns, 1 - for Dx, 2 - for Dx and fx
                
                % check for Dx one by one
                for v=1:volnum
                    % doses under v
                    dv(:)=0;
                    for m=1:ptnum
                        f=find(DAobj.DVH_LQ{m,2}(:,3)>=DAobj.VolBins(v));
                        if isempty(f) % no volume large enough, skip it (it is zero by default)
                            continue;
                        end
                        f=f(end);
                        if DAobj.DVH_LQ{m,2}(f,3)==DAobj.VolBins(v)
                            dv(m)=DAobj.DVH_LQ{m,2}(f,1);
                        else
                            dv(m) = interp1( [DAobj.DVH_LQ{m,2}(f,3);DAobj.DVH_LQ{m,2}(f+1,3)], [DAobj.DVH_LQ{m,2}(f,1);DAobj.DVH_LQ{m,2}(f+1,1)], DAobj.VolBins(v) );
                        end
                    end
                    g=find(dv); % patients with non-zero doses
                    
                    % check if sample size is too small
                    if length(g)<DAobj.CoxMinSize
                        stats.beta=-inf; stats.logl=-inf; stats.h=-inf;
                        DAobj.CoxStatsDx(v,1)=stats; DAobj.CoxStatsDx(v,2)=stats;
                        continue;
                    end
                    
                    % otherwise calcualte the cox model
                    % Dx
                    [~,logl,h,stats]=coxphfit(dv(g),DAobj.DateSurvive(g),'baseline',0,'censoring',DAobj.flgCensor(g));
                    stats.logl=logl; stats.h=h;
                    DAobj.CoxStatsDx(v,1)=stats;
                    % Dx + fx
                    if flgfx
                        [~,logl,h,stats]=coxphfit([dv(g),DAobj.FxNum(g)],DAobj.DateSurvive(g),'baseline',0,'censoring',DAobj.flgCensor(g));
                        stats.logl=logl; stats.h=h;
                    else
                        stats.beta=-inf; stats.logl=-inf; stats.h=-inf;
                    end
                    DAobj.CoxStatsDx(v,2)=stats;
                end
            % Vx
                dosenum=length(DAobj.DoseBins);
                DAobj.CoxStatsVx=repmat(stats,[dosenum,2]); % 2 columns, 1 - for Dx, 2 - for Dx and fx
                vd=zeros(ptnum,1);
                % check for Vx one by one
                for d=1:dosenum
                    % volumes under d
                    vd(:)=0;
                    for m=1:ptnum
                        f=find(DAobj.DVH_LQ{m,2}(:,1)<=DAobj.DoseBins(d)); f=f(end); % f won't be empty by definition
                        if DAobj.DVH_LQ{m,2}(f,1)<DAobj.DVH_LQ{m,2}(end,1) % not the last element, interpolate to get the best estimation of volume, for the last element, it is zero by default
                            vd(m) = interp1( [DAobj.DVH_LQ{m,2}(f,1);DAobj.DVH_LQ{m,2}(f+1,1)], [DAobj.DVH_LQ{m,2}(f,3);DAobj.DVH_LQ{m,2}(f+1,3)], DAobj.DoseBins(d) );
                        end
                    end
                    g=find(vd); % non-zeros dose cases
                    
                    % check if sample size is too small
                    if length(g)<DAobj.CoxMinSize
                        stats.beta=-inf; stats.logl=-inf; stats.h=-inf;
                        DAobj.CoxStatsVx(d,1)=stats; DAobj.CoxStatsVx(d,2)=stats;
                        continue;
                    end
                    
                    % otherwise calcualte the cox model
                    % Vx
                    [~,logl,h,stats]=coxphfit(vd(g),DAobj.DateSurvive(g),'baseline',0,'censoring',DAobj.flgCensor(g));
                    stats.logl=logl; stats.h=h;
                    DAobj.CoxStatsVx(d,1)=stats;
                    % Vx + fx
                    if flgfx
                        [~,logl,h,stats]=coxphfit([vd(g),DAobj.FxNum(g)],DAobj.DateSurvive(g),'baseline',0,'censoring',DAobj.flgCensor(g));
                        stats.logl=logl; stats.h=h;
                    else
                        stats.beta=-inf; stats.logl=-inf; stats.h=-inf;
                    end
                    DAobj.CoxStatsVx(d,2)=stats;
                end
                warning('on','MATLAB:singularMatrix');
                warning('on','stats:coxphfit:FitWarning');
        end
        function DAobj = LogRankTest(DAobj)
            % prepare
                volnum=length(DAobj.VolBins);
                dosenum=length(DAobj.DoseBins);
                ptnum=length(DAobj.DVH_LQ);
            
                DAobj.pDxmat=ones(dosenum,volnum);
                DAobj.flgCorDx=false(dosenum,volnum);
                DAobj.pDxdata=zeros(dosenum,volnum,4);
                
                DAobj.pVxmat=ones(dosenum,volnum);
                DAobj.flgCorVx=false(dosenum,volnum);
                DAobj.pVxdata=zeros(dosenum,volnum,4);
               
                sa=classSurvivalAnalysis(); % initialize a survivalanalysis obj
                
            % Dx
                dv=zeros(ptnum,1);
                for v=1:volnum
                    % doses under v
                    dv(:)=0;
                    for m=1:ptnum
                        f=find(DAobj.DVH_LQ{m,2}(:,3)>=DAobj.VolBins(v));
                        if isempty(f) % no volume large enough, skip it (it is zero by default)
                            continue;
                        end
                        f=f(end);
                        if DAobj.DVH_LQ{m,2}(f,3)==DAobj.VolBins(v)
                            dv(m)=DAobj.DVH_LQ{m,2}(f,1);
                        else
                            dv(m) = interp1( [DAobj.DVH_LQ{m,2}(f,3);DAobj.DVH_LQ{m,2}(f+1,3)], [DAobj.DVH_LQ{m,2}(f,1);DAobj.DVH_LQ{m,2}(f+1,1)], DAobj.VolBins(v) );
                        end
                    end
                    g=find(dv); % patients with non-zero doses
                    
                    % (di,vj)
                    flg_dosebelow2=-1;
                    numstart=DAobj.LogRankMinSize;
                    numend=length(g)-numstart;
                    for d=1:dosenum
                        % check the sample size
                        flg_dosebelow1=dv(g)<=DAobj.DoseBins(d); f=length(find(flg_dosebelow1));
                        if f<numstart || f>numend % one group has too few patients, skip it
                            continue;
                        end
                        % check the change of grouping
                        if isequal(flg_dosebelow1,flg_dosebelow2) % if it is the same grouping, skip the computation and save the result directly
                            DAobj.pDxmat(d,v)=DAobj.pDxmat(d-1,v);
                            DAobj.flgCorDx(d,v)=DAobj.flgCorDx(d-1,v); % the group with lower volume had worse survival curve, record it
                            DAobj.pDxdata(d,v,:)=DAobj.pDxdata(d-1,v,:);
                            continue;
                        end
                        flg_dosebelow2=flg_dosebelow1; % keep the current grouping
                        
                        % assign properties of object sa
                        survivedate={DAobj.DateSurvive(g(flg_dosebelow1)); DAobj.DateSurvive(g(~flg_dosebelow1))}; % survive time of each group
                        fcensor={DAobj.flgCensor(g(flg_dosebelow1));DAobj.flgCensor(g(~flg_dosebelow1))}; % censor flag for each group
                        sa.SurvivalTime=survivedate;
                        sa.FlagCensor=fcensor;
                        % compute survival curves and compare them
                        sa=sa.CalculateSurvivalCurve();
                        sa=sa.CombineSurvivalTime();
                        sa=sa.CompareSurvivalByLogrank();
                        DAobj.pDxmat(d,v)=sa.pValue;
                        DAobj.flgCorDx(d,v)=sa.CurveArea(1)<sa.CurveArea(2); % the group with lower volume had worse survival curve, record it
                        DAobj.pDxdata(d,v,:)=[length(survivedate{1,1}),sum(fcensor{1,1}), length(survivedate{2,1}),sum(fcensor{2,1})];
                    end
                end
            
            % Vx
                vd=zeros(ptnum,1); % volume v at dose d
                for d=1:dosenum
                    % volume under d
                    vd(:)=0;
                    for m=1:ptnum
                        f=find(DAobj.DVH_LQ{m,2}(:,1)<=DAobj.DoseBins(d)); f=f(end); % f won't be empty by definition
                        if DAobj.DVH_LQ{m,2}(f,1)<DAobj.DVH_LQ{m,2}(end,1) % not the last element, interpolate to get the best estimation of volume, for the last element, it is zero by default
                            vd(m) = interp1( [DAobj.DVH_LQ{m,2}(f,1);DAobj.DVH_LQ{m,2}(f+1,1)], [DAobj.DVH_LQ{m,2}(f,3);DAobj.DVH_LQ{m,2}(f+1,3)], DAobj.DoseBins(d) );
                        end
                    end
                    g=find(vd); % non-zeros volume cases
                    
                    % (di,vj)
                    flg_volbelow2=-1;
                    numend=length(g)-numstart;
                    for v=1:volnum
                        % check smaple size
                        flg_volbelow1=vd(g)<=DAobj.VolBins(v); f=length(find(flg_volbelow1)); % group DVHs by (d,v)
                        if f<numstart || f>numend % one group has too less patients, or the volume at the dose is too small, skip it
                            continue;
                        end
                        % check the change of grouping
                        if isequal(flg_volbelow1,flg_volbelow2) % if it is the same grouping, skip the computation and save the result directly
                            DAobj.pVxmat(d,v)=DAobj.pVxmat(d,v-1);
                            DAobj.flgCorVx(d,v)=DAobj.flgCorVx(d,v-1); % the group with lower volume had worse survival curve, record it
                            DAobj.pVxdata(d,v,:)=DAobj.pVxdata(d,v-1,:);
                            continue;
                        end
                        flg_volbelow2=flg_volbelow1; % keep the current grouping

                        % assign properties of object sa
                        survivedate={DAobj.DateSurvive(g(flg_volbelow1));DAobj.DateSurvive(g(~flg_volbelow1))}; % survive time of each group
                        fcensor={DAobj.flgCensor(g(flg_volbelow1));DAobj.flgCensor(g(~flg_volbelow1))}; % censor flag for each group
                        sa.SurvivalTime=survivedate;
                        sa.FlagCensor=fcensor;
                        % compute survival curves and compare them
                        sa=sa.CalculateSurvivalCurve();
                        sa=sa.CombineSurvivalTime();
                        sa=sa.CompareSurvivalByLogrank();
                        DAobj.pVxmat(d,v)=sa.pValue;
                        DAobj.flgCorVx(d,v)=sa.CurveArea(1)<sa.CurveArea(2); % the group with lower volume had worse survival curve, record it
                        DAobj.pVxdata(d,v,:)=[length(survivedate{1,1}),sum(fcensor{1,1}), length(survivedate{2,1}),sum(fcensor{2,1})];
                    end
                end
        end
    end
end