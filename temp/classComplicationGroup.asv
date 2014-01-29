classdef classComplicationGroup
    properties
        % general
        ptGrp = classComplicationIndividual.empty(1,0); % patient group of classComplicationIndividual objects
        numGrp = 0; % number of patients in the group

        Beta2Alpha = 0; % beta to alpha ratio
        CoxMinSize = 2; % minimum sample size for Cox model
        LogRankMinSize =2; % minimum sample size for log rank test

        % DVH
        TimeStep_DVH = 3; % time step size in month for dynamic analysis
        TimeBins_DVH % the time bins can be overwritten in case irregular bins steps were used
        
        DoseStep_DVH = 1; % dose step size for computation based on DVH (in Gy)
        VolStep_DVH = 1; % volume step size for computation based on DVH (in cc)
        DoseBins_DVH % the dose bins can be generated from DoseStep_DVH, or overwritten in case irregular bin steps were used
        VolBins_DVH % similar to DoseBins_DVH

        PatientTotal_DVH % matrix for total patients computed from DVH
        PatientComp_DVH % matrix for patients with complications computed from DVH
        PatientComplicationExact_DVH % matrix for patient survival curve from DVH
        PatientComplicationExactMat_DVH % matrix for patient survival curve at specific time point, from DVH
        PatientComplicationBin_DVH % matrix of patient survival from binned atlas.
        PatientComplicationCurveOverall_DVH = classSurvivalAnalysis.empty(1,0); % overall complication curve

        PatientSurvivalOverall_DVH = classSurvivalAnalysis.empty(1,0); % overall survival curve
        PatientSurvivalRelapseFree_DVH = classSurvivalAnalysis.empty(1,0); % relapse free survival curve
        
        CoxParameters_DVH % cox model pairs (cell 1 -- label of the Cox model fitting, cell 2 -- fitting result
        LogRank_DVH % log rank test for (Di,Vj) (cell 1 -- label of the test (Dx, Vx, DVx), cell 2 -- n1, censor1, n2, censor2, p-value for each grid point

        % EUD
%         lnn = (-1:0.1:1)'; % log10(n) for EUD computation, default is from -1 to 1 with step of 0.1
        lgn = (-1:0.1:1)';
        DoseStep_EUD = 1; % dose step size for computation based on EUD (in Gy)
        DoseBins_EUD; % the dose bins can be generated from DoseStep_EUD, or overwritten in case irregular bin steps were used
        
        PatientTotal_EUD % matrix for total patients computed from EUD
        PatientComp_EUD % matrix for patients with complications computed from EUD

        BetaCumulativeMat_EUD % based on EUD, the probability that the complication rate is larger than a specific number BetaCumulativeThreshold_EUD
        BetaCumulativeThreshold_EUD = 0.2;
        BetaInverseMat_EUD % based on EUD, The chi-square where the probability is BetaInverseThreshold_EUD
        BetaInverseThreshold_EUD = 0.16;
        
        LogisticRegressionMatExact_EUD % logistic regression model using exact EUD
        LogisticRegressionMatBin_EUD % logistic regression model using binned EUD
        LymanModelExact_EUD % Lyman model parameters
        LymanModelBin_EUD
        LymanModelGridExact_EUD % Lyman results (likelihood/probability) for each grid point
        LymanModelGridBin_EUD
    end

    methods % set member vaulues
        function CGobj = classComplicationGroup()
%             CGobj.ptGrp = classComplicationIndividual.empty(1,0);
        end
        function CGobj = set.Beta2Alpha(CGobj,beta2alpha)
            if any (beta2alpha < 0)
                error('Beta to Alpha ratio is negative, in set.Beta2Alpha of classComplicationGroup');
            end
            if length(beta2alpha)>1
                warning('classComplicationGroup:beta2alpha','Parameter for Beta2Alpha not a scalar, choose the first element');
            end
            for m = 1:CGobj.numGrp
                CGobj.ptGrp(m).Beta2Alpha = beta2alpha;
            end
            CGobj.Beta2Alpha = beta2alpha;
        end
        function CGobj = set.TimeBins_DVH(CGobj,timebins)
            if any( diff( timebins ) <=0 )
                disp('time bins is not strictly monotonic increasing when assign to TimeBins_DVH in classComplicationGroup');
                disp('thus they are sorted');
                timebins = sort(timebins);
                timebins(diff(timebins)==0)=[];
            end
            CGobj.TimeBins_DVH = timebins(:);
        end
        function CGobj = set.DoseBins_DVH(CGobj,dosebins)
            if any(diff(dosebins)<=0)
                disp('dose bins is not monotonic increasing when assigned to DoseBins_DVH in classComplicationGroup');
                disp('thus they are sorted');
                dosebins = sort( dosebins );
                dosebins(diff(dosebins)==0)=[];
            end
            CGobj.DoseBins_DVH = dosebins(:);
        end
        function CGobj = set.VolBins_DVH(CGobj,volbins)
            if any(diff(volbins)<=0)
                disp('volume bins is not monotonic increasing when assigned to VolBins_DVH in classComplicationGroup');
                disp('thus they are sorted');
                volbins = sort( volbins );
                volbins(diff(bolbins)==0)=[];
            end
            CGobj.VolBins_DVH = volbins(:);
        end
        function CGobj = set.lgn(CGobj,lgn)
            CGobj.lgn = lgn(:);
        end
        function CGobj = set.DoseBins_EUD(CGobj,dosebins)
            if isempty(dosebins)
                CGobj.DoseBins_EUD = dosebins;
                return;
            end
            
            if any(diff(dosebins)<=0)
                disp('dose bins is not monotonic increasing when assigned to DoseBins_EUD in classComplicationGroup');
                disp('thus they are sorted');
                dosebins = sort( dosebins );
                dosebins(diff(dosebins)==0)=[];
            end
            euds = [CGobj.ptGrp.EUD]; dmax = max(euds(:));
            if dmax > dosebins(end)
                disp('assigned dose bins can not cover the maximum EUD');
            end
            CGobj.DoseBins_EUD = dosebins(:);
        end
        function CGobj = set.BetaCumulativeThreshold_EUD(CGobj, betacumthreshold)
            CGobj.BetaCumulativeThreshold_EUD = betacumthreshold(:);
        end
        function CGobj = set.BetaInverseThreshold_EUD(CGobj,betainvthreshold)
            CGobj.BetaInverseThreshold_EUD = betainvthreshold(:);
        end
    end

    methods % patient operations like adding and removing
        function CGobj = AddPatient(CGobj,ptobjs)
            % check if the new patient info matches the class definition
            if ~all(isa(ptobjs,'classComplicationIndividual'))
                disp('Not an instance of patient individual class when adding a new patient');
                return;
            end
            
            % make sure the log10(n) and the beta2alpha are the same for all patients
            for m = 1:length(ptobjs)
                ptobjs(m).lgn = CGobj.lgn;
                ptobjs(m).Beta2Alpha = CGobj.Beta2Alpha;
            end

            % for empty patient object, add the patient directly
            if isempty([CGobj.ptGrp.PatientID])
                CGobj.ptGrp = ptobjs(:);
                CGobj.numGrp = size(CGobj.ptGrp,1);
                return;
            end
            
            % update existing patients in the group
            [fg,g] = ismember( {ptobjs.PatientID}, {CGobj.ptGrp.PatientID} );
            f=find(g);
            for m = 1:length(f)
                CGobj.ptGrp(g(f(m))) = ptobjs(f(m));
            end
            
            % add new patient
            CGobj.numGrp = CGobj.numGrp + sum(~fg);
            CGobj.ptGrp(end+1:CGobj.numGrp) = ptobjs(~fg);
        end
        function CGobj = RemovePatient(CGobj,idx)
            if ~exist('idx','var') % if idx is not passed on, remove all patients (default)
                idx = true(CGobj.numGrp,1);
            end
            CGobj.ptGrp(idx) = [];
            CGobj.numGrp = length(CGobj.ptGrp);
        end
        function CGobj = LinearQuartraticCorrection(CGobj)
            for k = 1:CGobj.numGrp
                CGobj.ptGrp(k) = CGobj.ptGrp(k).LinearQuartraticCorrection();
            end
        end
        function flg = PatientsWithComplicationData(CGobj)
            pt = CGobj.ptGrp;
            f1 = cellfun(@(x) isempty(x)|isinf(x)|isnan(x),{pt.BaselineDate}); % patients with no baseline date
            f2 = cellfun(@(x) isempty(x)|isinf(x)|isnan(x),{pt.CompOccurDate}); % patients with no complication date
            flg = find(~f2);  % patients with complciation date were not censored, so this information is restated
            for k = 1:length(flg)
                pt(flg(k)).flgCensor = 0;
            end
            
            f3 = cellfun(@(x) isempty(x)|isinf(x)|isnan(x),{pt.LastFollowupDate}); % patients with no last follow up date
            f4 = cellfun(@(x) isempty(x)|isinf(x)|isnan(x),{pt.flgCensor}); % patients with no censor info
            flg = f1 | (f2&f3) | f4; % patients missing at least one data
            flg = ~flg; % patients with all data
        end
    end

    methods % DVH operations
        function CGobj = CalculateDoseBins_DVH(CGobj)
            doses = cellfun(@(x) x(end), {CGobj.ptGrp.DoseBins_LQ}); dmax = max(doses);
            CGobj.DoseBins_DVH = ( 0 : CGobj.DoseStep_DVH : dmax+CGobj.DoseStep_DVH )';
        end
        function CGobj = CalculateDoseBinsLog_DVH(CGobj)
            doses = cellfun(@(x) x(end), {CGobj.ptGrp.DoseBins_LQ}); dmax = log10(max(doses));
            doses = cellfun(@(x) x(1), {CGobj.ptGrp.DoseBins_LQ}); dmin = min(doses); dmin = log10(max(1,dmin));
            dosebinslog = dmin : CGobj.DoseStep_DVH : (dmax+CGobj.DoseStep_DVH);
            CGobj.DoseBins_DVH = [0, 10.^dosebinslog]';
        end
        function CGobj = CalculateVolBins_DVH(CGobj)
            vols = cellfun(@(x) x(1), {CGobj.ptGrp.VolCum}); vmax = max(vols);
            CGobj.VolBins_DVH = ( 0 : CGobj.VolStep_DVH : vmax+CGobj.VolStep_DVH )';
        end
        function CGobj = CalculateVolBinsLog_DVH(CGobj)
            vols = cellfun(@(x) x(1), {CGobj.ptGrp.VolCum}); vmax = log10(max(vols));
            vols = cellfun(@(x) min(x(x>0)), {CGobj.ptGrp.VolCum}); vmin = min(vols); vmin = log10(vmin);
            vmin = fix(vmin/CGobj.VolStep_DVH)*CGobj.VolStep_DVH;
            volbinslog = vmin : CGobj.VolStep_DVH : (vmax + CGobj.VolStep_DVH);
            CGobj.VolBins_DVH = [0, 10.^volbinslog]';
        end
        function CGobj = CalculateTimeBins_DVH(CGobj)
            timebins = ( [CGobj.ptGrp.LastFollowupDate] - [CGobj.ptGrp.BaselineDate] ) / 30; tmax = max(timebins);
            CGobj.TimeBins_DVH = ( 0 : CGobj.TimeStep_DVH : tmax+CGobj.TimeStep_DVH )';
        end
        function CGobj = AtlasAlongSampleTime_DVH(CGobj)
            % prepare
            CGobj.PatientTotal_DVH = zeros( length(CGobj.DoseBins_DVH), length(CGobj.VolBins_DVH), length(CGobj.TimeBins_DVH) );
            CGobj.PatientComp_DVH = zeros( length(CGobj.DoseBins_DVH), length(CGobj.VolBins_DVH), length(CGobj.TimeBins_DVH) );
            
            % select patients with data
            f = CGobj.PatientsWithComplicationData();
            pt = CGobj.ptGrp(f);
            num = sum(f); % number of patients

            % complication time
            f2 = ~cellfun(@(x) isempty(x)|isinf(x)|isnan(x),{pt.CompOccurDate}); % patients with complication date
            f3 = ~cellfun(@(x) isempty(x)|isinf(x)|isnan(x),{pt.LastFollowupDate}); % patients with last follow up date
            complicationdays = inf(num,1);
            lastfollowup = inf(num,1);
            complicationdays(f2) = ([pt(f2).CompOccurDate] - [pt(f2).BaselineDate])' / 30;
            lastfollowup(f3) = ([pt(f3).LastFollowupDate] - [pt(f3).BaselineDate])' / 30;
            lastfollowup = min( lastfollowup, complicationdays );

            % Vx
            vols = zeros( num, 1 );
            for ii = 1:length(CGobj.DoseBins_DVH)
                % volume of each patient at current dose (x)
                vols(:) = 0;
                for jj = 1:num
                    vols(jj) = pt(jj).VolAtDose( CGobj.DoseBins_DVH(ii) );
                end
                vols(vols==0) = -1; % exclude zero volume patients
                
                % matrix at each (Di, Vj, Tk)
                for jj = 1:length(CGobj.VolBins_DVH) % for each volume point under dose x
                    f = find( vols >= CGobj.VolBins_DVH(jj) ); % patient at the grid point
                    % for each time point
                    for kk = 1:length(CGobj.TimeBins_DVH)
                        % patients with last followup or complications before the time point
                        f1 = find( lastfollowup(f) <= CGobj.TimeBins_DVH(kk) );
                        % total patients
                        CGobj.PatientTotal_DVH(ii,jj,kk) = length(f1);
                        % patients with complications
                        f1 = find( complicationdays(f) <= CGobj.TimeBins_DVH(kk) );
                        CGobj.PatientComp_DVH(ii,jj,kk) = length(f1);
                    end
                end
            end
        end
        function CGobj = ComplicationCurves_DVH(CGobj)
            % select patients with data
            f = CGobj.PatientsWithComplicationData();
            pt = CGobj.ptGrp(f);
            num = sum(f); % number of patients

            % complication time
            f2 = ~cellfun(@(x) isempty(x)|isinf(x)|isnan(x),{pt.CompOccurDate}); % patients with complication date
            f3 = ~cellfun(@(x) isempty(x)|isinf(x)|isnan(x),{pt.LastFollowupDate}); % patients with last follow up date
            complicationdays = inf(num,1);
            lastfollowup = inf(num,1);
            complicationdays(f2) = ([pt(f2).CompOccurDate] - [pt(f2).BaselineDate])' / 30;
            lastfollowup(f3) = ([pt(f3).LastFollowupDate] - [pt(f3).BaselineDate])' / 30;
            lastfollowup = min( lastfollowup, complicationdays );
            flgcensor = [pt.flgCensor]';

            % overall survival
            sa = classSurvivalAnalysis();
            sa.SurvivalTime = {lastfollowup};
            sa.flgCensor = {flgcensor};
            sa = sa.CalculateSurvivalCurve();

            CGobj.PatientComplicationCurveOverall_DVH = sa;
        end
        function CGobj = SurvivalCurves_DVH(CGobj)
            % select patients with data
            pt = CGobj.ptGrp;
            f1 = cellfun(@(x) isempty(x)|isinf(x)|isnan(x),{pt.BaselineDate}); % patients with no baseline date
            f2 = cellfun(@(x) isempty(x)|isinf(x)|isnan(x),{pt.DeathDate}); % patients with no death date
            f3 = cellfun(@(x) isempty(x)|isinf(x)|isnan(x),{pt.LastFollowupDate}); % patients with no last follow-up date
            f4 = cellfun(@(x) isempty(x)|isinf(x)|isnan(x),{pt.CompOccurDate}); % patients with no complication date
            flg = f1 | (f2 & f3 & f4); % patient with out enough data
            pt(flg) = [];

            num = sum(~flg); % number of patients

            % survival time
            f2 = ~cellfun(@(x) isempty(x)|isinf(x)|isnan(x),{pt.DeathDate}); % patients with death date
            f3 = ~cellfun(@(x) isempty(x)|isinf(x)|isnan(x),{pt.LastFollowupDate}); % patients with last follow up date
            f4 = ~cellfun(@(x) isempty(x)|isinf(x)|isnan(x),{pt.CompOccurDate}); % patients with complication date
            survivaltime = inf(num,1);
            survivaltime(f4) = ([pt(f4).CompOccurDate] - [pt(f4).BaselineDate])' / 30;
            survivaltime(f3) = ([pt(f3).LastFollowupDate] - [pt(f3).BaselineDate])' / 30;
            survivaltime(f2) = ([pt(f2).DeathDate] - [pt(f2).BaselineDate])' / 30;
            flgcensor = ~f2';

            % overall survival
            sa = classSurvivalAnalysis();
            sa.SurvivalTime = {survivaltime};
            sa.flgCensor = {flgcensor};
            sa = sa.CalculateSurvivalCurve();

            CGobj.PatientSurvivalOverall_DVH = sa;
            
            % free of relapse survival
            f2 = ~cellfun(@(x) isempty(x)|isinf(x)|isnan(x),{pt.RelapseDate}); % patients with relapse date
            flgcensor(:) = 1;
            flgcensor(f2) = 0;
            survivaltime(f2) = ([pt(f2).RelapseDate] - [pt(f2).BaselineDate])'/30;
            sa.SurvivalTime = {survivaltime};
            sa.flgCensor = {flgcensor};
            sa = sa.CalculateSurvivalCurve();

            CGobj.PatientSurvivalRelapseFree_DVH = sa;
        end
        function CGobj = ComplicationActuary_DVH(CGobj)
            % select patients with data
            f = CGobj.PatientsWithComplicationData();
            CG = CGobj.RemovePatient(~f);

            % complication time
            f2 = ~cellfun('isempty',{CG.ptGrp.CompOccurDate}); % patients with no complication date
            f3 = ~cellfun('isempty',{CG.ptGrp.LastFollowupDate}); % patients with no last follow up date
            complicationdays = inf(CG.numGrp,1);
            lastfollowup = inf(CG.numGrp,1);
            complicationdays(f2) = ([CG.ptGrp(f2).CompOccurDate] - [CG.ptGrp(f2).BaselineDate])' / 30;
            lastfollowup(f3) = ([CG.ptGrp(f3).LastFollowupDate] - [CG.ptGrp(f3).BaselineDate])' / 30;
            lastfollowup = min( lastfollowup, complicationdays );
            flgcensor = [CG.ptGrp.flgCensor]';

            % populate to all grid points
            sa = classSurvivalAnalysis(); % a new clean object
            pse = repmat(sa, [length(CG.DoseBins_DVH), length(CG.VolBins_DVH)]); % patient survival exact

            % Vx
            vols = zeros( CG.numGrp, 1 );
            for ii = 1:length(CG.DoseBins_DVH)
                % volume of each patient at current dose (x)
                vols(:) = 0;
                for jj = 1:CG.numGrp
                    vols(jj) = CG.ptGrp(jj).VolAtDose( CG.DoseBins_DVH(ii) );
                end
                vols(vols==0) = -1; % exclude zero volume patients
                
                % matrix at each (Di, Vj)
                for jj = 1:length(CG.VolBins_DVH) % for each volume point under dose x
                    f = vols >= CG.VolBins_DVH(jj); % patient at the grid point
                    if ~any(f) % no patient for the grid, clean its contents
                        continue;
                    end

                    % use clasSurvivalAnalysis to compute the Kelplan Meier curve
                    pse(ii,jj).SurvivalTime = {lastfollowup(f)};
                    pse(ii,jj).flgCensor = {flgcensor(f)};
                    pse(ii,jj) = pse(ii,jj).CalculateSurvivalCurve();
                end
            end
            CGobj.PatientComplicationExact_DVH = pse;
        end
        function CGobj = ComplicationAtSampleTimeActuary_DVH(CGobj)
            % prepare
            psem = zeros( length(CGobj.DoseBins_DVH), length(CGobj.VolBins_DVH), length(CGobj.TimeBins_DVH) );
            pse = CGobj.PatientComplicationExact_DVH(:); pse = reshape(pse,size(CGobj.PatientComplicationExact_DVH));
            tb = CGobj.TimeBins_DVH(:);
            db = CGobj.DoseBins_DVH(:);
            vb = CGobj.VolBins_DVH(:);

            % Vx
            for ii = 1:length(db)
                % matrix at each (Di, Vj, Tk)
                for jj = 1:length(vb) % for each volume point under dose x
                    sa = pse(ii,jj);
                    if isempty(sa.SurvivalTime) % no patient for the grid, skip it
                        continue;
                    end
                    
                    % sample the survival curve for each time point
                    f = find( tb <= sa.SurvivalTimeSorted{1}(end) );
                    for kk = 1:length(f)
                        f1 = find( tb(kk) >= sa.SurvivalTimeSorted{1} );
                        psem(ii,jj,kk) = sa.SurvivalCurve{1}(f1(end));
                    end
                end
            end
            CGobj.PatientComplicationExactMat_DVH = psem;
        end
        function CGobj = ComplicationAtSampleTimeBins_DVH(CGobj)
            % prepare
            CGobj.PatientComplicationBin_DVH = ones( length(CGobj.DoseBins_DVH), length(CGobj.VolBins_DVH), length(CGobj.TimeBins_DVH) );
            
            % survival along the time bins
            for kk = 2:length(CGobj.TimeBins_DVH)
                % patients at risk
                ptrisk = CGobj.PatientTotal_DVH(:,:,end) - CGobj.PatientTotal_DVH(:,:,kk-1);
                % patient with complication
                ptcomp = CGobj.PatientComp_DVH(:,:,end) - CGobj.PatientComp_DVH(:,:,kk-1);
                % rate of patients without complication
                CGobj.PatientComplicationBin_DVH(:,:,kk) = CGobj.PatientComplicationBin_DVH(:,:,kk-1).*(1-ptcomp./ptrisk);
            end
        end
        function CGobj = CoxModelExact_DVH(CGobj)
            % select patients with data
            f = CGobj.PatientsWithComplicationData();
            CG = CGobj.RemovePatient(~f);
            
            % survival/complication time
            f2 = ~cellfun('isempty',{CG.ptGrp.CompOccurDate}); % patients with no complication date
            f3 = ~cellfun('isempty',{CG.ptGrp.LastFollowupDate}); % patients with no last follow up date
            compdate = inf(CG.numGrp,1);
            lastfollowup = inf(CG.numGrp,1);
            compdate(f2) = ([CG.ptGrp(f2).CompOccurDate] - [CG.ptGrp(f2).BaselineDate])' / 30;
            lastfollowup(f3) = ([CG.ptGrp(f3).LastFollowupDate] - [CG.ptGrp(f3).BaselineDate])' / 30;
            compdate = min( lastfollowup, compdate );
            flgcensor = [CG.ptGrp.flgCensor]';

            warning('off');
            % dmax
                dv = cellfun(@(x) x(end), {CG.ptGrp.DoseBins_LQ}');

                [~,logl,h,stats]=coxphfit(dv,compdate,'baseline',0,'censoring',flgcensor);
                stats.logl=logl; stats.h=h;
                stats.data_exposure = dv; stats.data_hazard = compdate;

            % fx
                fxnum=[CG.ptGrp.FxNum]'; % all fractions
                flgfx = length(unique(fxnum))>1; % if there is more than one fraction numbers, do the cox model
                if flgfx
                    [~,logl,h,stats]=coxphfit([CG.ptGrp.FxNum]',compdate,'baseline',0,'censoring',flgcensor);
                    stats.logl=logl; stats.h=h;
                    stats.data_exposure = [CG.ptGrp.FxNum]'; stats.data_hazard = compdate;
                    
                    if isempty(CGobj.CoxParameters_DVH)
                            CGobj.CoxParameters_DVH{end+1,1}='Fx';
                            CGobj.CoxParameters_DVH{end,2}=stats;
                    else
                        f = cellfun(@(x) strcmpi('Fx',x),CGobj.CoxParameters_DVH(:,1));
                        if any(f)
                            CGobj.CoxParameters_DVH{f,2}=stats;
                        else
                            CGobj.CoxParameters_DVH{end+1,1}='Fx';
                            CGobj.CoxParameters_DVH{end,2}=stats;
                        end
                    end
                end

            % prescription dose
                dv = cell2mat([CG.ptGrp.DosePrescription]');
                [~,logl,h,stats]=coxphfit(dv,compdate,'baseline',0,'censoring',flgcensor);
                stats.logl=logl; stats.h=h;
                stats.data_exposure = dv; stats.data_hazard = compdate;

                    if isempty(CGobj.CoxParameters_DVH)
                            CGobj.CoxParameters_DVH{end+1,1}='Dp'; % dose prescription
                            CGobj.CoxParameters_DVH{end,2}=stats;
                    else
                        f = cellfun(@(x) strcmpi('Dp',x),CGobj.CoxParameters_DVH(:,1));
                        if any(f)
                            CGobj.CoxParameters_DVH{f,2}=stats;
                        else
                            CGobj.CoxParameters_DVH{end+1,1}='Dp';
                            CGobj.CoxParameters_DVH{end,2}=stats;
                        end
                    end

                % fx + prescription dose
                if flgfx
                    [~,logl,h,stats]=coxphfit([dv,fxnum],compdate,'baseline',0,'censoring',flgcensor);
                    stats.logl=logl; stats.h=h;
                    stats.data_exposure = [dv,fxnum]; stats.data_hazard = compdate;
                    
                        f = cellfun(@(x) strcmpi('DpFx',x),CGobj.CoxParameters_DVH(:,1));
                        if any(f)
                            CGobj.CoxParameters_DVH{f,2}=stats;
                        else
                            CGobj.CoxParameters_DVH{end+1,1}='DpFx';
                            CGobj.CoxParameters_DVH{end,2}=stats;
                        end
                end
                
            % Dx
                volnum=length(CG.VolBins_DVH);
                CoxDx=repmat(stats,[volnum,1]);
                CoxDVx=repmat(stats,[volnum,1]);
                CoxDxFx=repmat(stats,[volnum,1]);
                CoxDVxFx=repmat(stats,[volnum,1]);
                
                % check for Dx one by one
                dv = zeros(CG.numGrp,1);
                for v=1:volnum
                    % doses under v
                    dv(:) = 0;
                    for k = 1:CG.numGrp
                        dv(k) = CG.ptGrp(k).DoseAtVol( CG.VolBins_DVH(v) );
                    end
                    
                    % Cox model for all patients
                    % Dx
                    try
                        if length(unique(dv))>1
                            [~,logl,h,stats]=coxphfit(dv,compdate,'baseline',0,'censoring',flgcensor);
                            stats.logl=logl; stats.h=h;
                        else
                            stats.beta=-inf; stats.logl=-inf; stats.h=-inf;
                        end
                    catch
                        stats.beta=-inf; stats.logl=-inf; stats.h=-inf;
                    end
                    stats.data_exposure = dv; stats.data_hazard = compdate;
                    CoxDVx(v)=stats;
                    % Dx + fx
                    if flgfx
                        try
                            [~,logl,h,stats]=coxphfit([dv,fxnum],compdate,'baseline',0,'censoring',flgcensor);
                            stats.logl=logl; stats.h=h;
                        catch
                            stats.beta=-inf; stats.logl=-inf; stats.h=-inf;
                        end
                        stats.data_exposure = [dv,fxnum]; stats.data_hazard = compdate;
                        CoxDVxFx(v)=stats;
                    end
                    
                    % Cox model for patients with non-zero data
                    g=find(dv); % patients with non-zero doses
                    % check if sample size is too small
                    if length(g)<CG.CoxMinSize
                        stats.beta=-inf; stats.logl=-inf; stats.h=-inf;
                        CoxDx(v)=stats; CoxDxFx(v)=stats;
                        continue;
                    end
                    
                    % otherwise calcualte the cox model
                    % Dx
                    [~,logl,h,stats]=coxphfit(dv(g),compdate(g),'baseline',0,'censoring',flgcensor(g));
                    stats.logl=logl; stats.h=h;
                    stats.data_exposure = dv(g); stats.data_hazard = compdate(g);
                    CoxDx(v)=stats;
                    
                    % Dx + fx
                    if flgfx
                        [~,logl,h,stats]=coxphfit([dv(g),fxnum(g)],compdate(g),'baseline',0,'censoring',flgcensor(g));
                        stats.logl=logl; stats.h=h;
                        stats.data_exposure = [dv(g),fxnum(g)]; stats.data_hazard = compdate(g);
                        CoxDxFx(v)=stats;
                    end
                end
                
                % save results
                if isempty(CGobj.CoxParameters_DVH)
                    CGobj.CoxParameters_DVH{end+1,1}='Dx';
                    CGobj.CoxParameters_DVH{end,2}=CoxDx;
                else
                    f = cellfun(@(x) strcmpi('Dx',x),CGobj.CoxParameters_DVH(:,1));
                    if any(f)
                        CGobj.CoxParameters_DVH{f,2}=CoxDx;
                    else
                        CGobj.CoxParameters_DVH{end+1,1}='Dx';
                        CGobj.CoxParameters_DVH{end,2}=CoxDx;
                    end
                end
                f = cellfun(@(x) strcmpi('DVx',x),CGobj.CoxParameters_DVH(:,1));
                if any(f)
                    CGobj.CoxParameters_DVH{f,2}=CoxDVx;
                else
                    CGobj.CoxParameters_DVH{end+1,1}='DVx';
                    CGobj.CoxParameters_DVH{end,2}=CoxDVx;
                end
                if flgfx
                    f = cellfun(@(x) strcmpi('DVxFx',x),CGobj.CoxParameters_DVH(:,1));
                    if any(f)
                        CGobj.CoxParameters_DVH{f,2}=CoxDVxFx;
                    else
                        CGobj.CoxParameters_DVH{end+1,1}='DVxFx';
                        CGobj.CoxParameters_DVH{end,2}=CoxDVxFx;
                    end
                    f = cellfun(@(x) strcmpi('DxFx',x),CGobj.CoxParameters_DVH(:,1));
                    if any(f)
                        CGobj.CoxParameters_DVH{f,2}=CoxDxFx;
                    else
                        CGobj.CoxParameters_DVH{end+1,1}='DxFx';
                        CGobj.CoxParameters_DVH{end,2}=CoxDxFx;
                    end
                end

            % Vx
                dosenum=length(CG.DoseBins_DVH);
                CoxVx=repmat(stats,[dosenum,1]);
                CoxVDx=repmat(stats,[dosenum,1]);
                CoxVxFx=repmat(stats,[dosenum,1]);
                CoxVDxFx=repmat(stats,[dosenum,1]);
                vd=zeros(CG.numGrp,1);
                % check for Vx one by one
                for d=1:dosenum
                    % volumes under d
                    vd(:)=0;
                    for k=1:CG.numGrp
                        vd(k) = CG.ptGrp(k).VolAtDose( CG.DoseBins_DVH(d) );
                    end

                    % Cox model for all patients
                    % Vx
                    try
                        if length(unique(vd))>1
                            [~,logl,h,stats]=coxphfit(vd,compdate,'baseline',0,'censoring',flgcensor);
                            stats.logl=logl; stats.h=h;
                        else
                            stats.beta=-inf; stats.logl=-inf; stats.h=-inf;
                        end
                    catch
                        stats.beta=-inf; stats.logl=-inf; stats.h=-inf;
                    end
                    stats.data_exposure = vd; stats.data_hazard = compdate;
                    CoxVDx(d)=stats;
                    
                    % Vx + fx
                    if flgfx
                        try
                            [~,logl,h,stats]=coxphfit([vd,fxnum],compdate,'baseline',0,'censoring',flgcensor);
                            stats.logl=logl; stats.h=h;
                        catch
                            stats.beta=-inf; stats.logl=-inf; stats.h=-inf;
                        end
                        stats.data_exposure = [vd,fxnum]; stats.data_hazard = compdate;
                        CoxVDxFx(d)=stats;
                    end
                    
                    % Cox model for patients with non-zero data
                    g=find(vd); % non-zeros dose cases
%                     vd(vd==0) = -1; % exclude zero volume patients
                    % check if sample size is too small
                    if length(g)<CG.CoxMinSize
                        stats.beta=-inf; stats.logl=-inf; stats.h=-inf;
                        CoxVx(d)=stats; CoxVxFx(d)=stats;
                        continue;
                    end
                    
                    % otherwise calcualte the cox model
                    % Vx
                    [~,logl,h,stats]=coxphfit(vd(g),compdate(g),'baseline',0,'censoring',flgcensor(g));
                    stats.logl=logl; stats.h=h;
                    stats.data_exposure = vd(g); stats.data_hazard = compdate(g);
                    CoxVx(d)=stats;
                    
                    % Vx + fx
                    if flgfx
                        [~,logl,h,stats]=coxphfit([vd(g),fxnum(g)],compdate(g),'baseline',0,'censoring',flgcensor(g));
                        stats.logl=logl; stats.h=h;
                        stats.data_exposure = [vd(g),fxnum(g)]; stats.data_hazard = compdate(g);
                        CoxVxFx(d)=stats;
                    end
                end
                
                % save results
                f = cellfun(@(x) strcmpi('VDx',x),CGobj.CoxParameters_DVH(:,1));
                if any(f)
                    CGobj.CoxParameters_DVH{f,2}=CoxVDx;
                else
                    CGobj.CoxParameters_DVH{end+1,1}='VDx';
                    CGobj.CoxParameters_DVH{end,2}=CoxVDx;
                end
                f = cellfun(@(x) strcmpi('Vx',x),CGobj.CoxParameters_DVH(:,1));
                if any(f)
                    CGobj.CoxParameters_DVH{f,2}=CoxVx;
                else
                    CGobj.CoxParameters_DVH{end+1,1}='Vx';
                    CGobj.CoxParameters_DVH{end,2}=CoxVx;
                end
                if flgfx
                    f = cellfun(@(x) strcmpi('VDxFx',x),CGobj.CoxParameters_DVH(:,1));
                    if any(f)
                        CGobj.CoxParameters_DVH{f,2}=CoxVDxFx;
                    else
                        CGobj.CoxParameters_DVH{end+1,1}='VDxFx';
                        CGobj.CoxParameters_DVH{end,2}=CoxVDxFx;
                    end
                    f = cellfun(@(x) strcmpi('VxFx',x),CGobj.CoxParameters_DVH(:,1));
                    if any(f)
                        CGobj.CoxParameters_DVH{f,2}=CoxVxFx;
                    else
                        CGobj.CoxParameters_DVH{end+1,1}='VxFx';
                        CGobj.CoxParameters_DVH{end,2}=CoxVxFx;
                    end
                end
                warning('on');
        end
        function CGobj = CoxModelAtlas_DVH(CGobj)
            if isempty(CGobj.PatientTotal_DVH)
                disp('Atlas does not exist');
                disp('in CoxModelAtlas_DVH');
                return;
            end

            warning('off');
%             warning('off','MATLAB:singularMatrix');
%             warning('off','stats:coxphfit:FitWarning');
%             warning('off','stats:coxphfit:RankDeficient');
%             warning('off','stats:coxphfit:IterOrEvalLimit');
            
            % prepare
            flgcensor = false(CGobj.PatientTotal_DVH(1,1,end),1);
            vx = zeros(size(flgcensor));
            dx = zeros(size(flgcensor));
            compdate = zeros(size(flgcensor));
            
            [dimd,dimv,dimt] = size( CGobj.PatientTotal_DVH ); % dimensions of dose, volume, and time
            
            dosebins = CGobj.DoseBins_DVH; dosebins(1:end-1) = (dosebins(1:end-1) + dosebins(2:end)) / 2; % dose bins
            volbins = CGobj.VolBins_DVH; volbins(1:end-1) = (volbins(1:end-1) + volbins(2:end)) / 2; % volume bins
            timebins = CGobj.TimeBins_DVH; timebins(2:end) = (timebins(1:end-1) + timebins(2:end)) / 2; % time bins
            
            % Cox model
            [~,logl,h,stats]=coxphfit([0; 1],[0; 1],'baseline',0,'censoring',[0; 1]);
            stats.logl=logl; stats.h=h;
            stats.data_exposure = [0; 1]; stats.data_hazard = [0; 1];
            % Vx
            CoxVxAtlas = repmat(stats,[dimd,1]);
            for d = 1:dimd % Cox model for each dose (x)
                % extract the matrix corresponding to the dose
                vc = squeeze( CGobj.PatientComp_DVH(d,:,:) ); % complication at different volume and time points
                vt = squeeze( CGobj.PatientTotal_DVH(d,:,:) ); % complication and censored at different time points
                vt = vt-vc; % cumulative censored info 
                
                vt(:,2:end) = diff(vt,1,2); % censored info at different time points
                vt(1:end-1,:) = abs(diff(vt,1,1)); % censored info at different volume bins
                vc(:,2:end) = diff(vc,1,2); % complication info at different time points
                vc(1:end-1,:) = abs(diff(vc,1,1)); % complication info at different volume bins
                
                % extract the volume, complication date, and censor info
                n = 0;
                vx(:) = 0; flgcensor(:) = 0; compdate(:) = 0;
                % censor info
                [v,t] = find(vt);
                for m = 1:length(v)
                    vx(n+1 : n+vt(v(m),t(m))) = volbins(v(m));
                    compdate(n+1 : n+vt(v(m),t(m))) = timebins(t(m));
                    flgcensor(n+1 : n+vt(v(m),t(m))) = true;
                    n = n + vt(v(m),t(m));
                end
                % complication info
                [v,t] = find(vc);
                for m = 1:length(v)
                    vx(n+1 : n+vc(v(m),t(m))) = volbins(v(m));
                    compdate(n+1 : n+vc(v(m),t(m))) = timebins(t(m));
                    n = n + vc(v(m),t(m));
                end
                
                
%                 for t = 1:dimt % extract patient info for each time point
%                     % censored patient
%                     f = find(vt(:,t)); % patient with censored info
%                     for v = 1:length(f)
%                         vx(n+1:n+vt(f(v),t)) = volbins(f(v));
%                         compdate(n+1:n+vt(f(v),t)) = timebins(t);
%                         flgcensor(n+1:n+vt(f(v),t)) = true;
%                         n = n+vt(f(v),t);
%                     end
%                     % complication
%                     f = find(vc(:,t)); % patient with complication info
%                     for v = 1:length(f)
%                         vx(n+1:n+vc(f(v),t)) = volbins(f(v));
%                         compdate(n+1:n+vc(f(v),t)) = timebins(t);
%                         n = n+vc(f(m),t);
%                     end
%                 end
                % cox model
                try
                    if length(unique(vx(1:n)))>1
                        [~,logl,h,stats]=coxphfit(vx(1:n),compdate(1:n),'baseline',0,'censoring',flgcensor(1:n));
                        stats.logl=logl; stats.h=h;
                    else
                        stats.beta=-inf; stats.logl=-inf; stats.h=-inf;
                    end
                catch
                    stats.beta=-inf; stats.logl=-inf; stats.h=-inf;
                end
                stats.data_exposure = vx(1:n); stats.data_hazard = compdate(1:n);
                CoxVxAtlas(d)=stats;
            end
            
            % Dx
            CoxDxAtlas = repmat(stats,[dimv,1]);
            for v = 1:dimv % Cox model for each volume (x)
                % extract the matrix corresponding to the dose
                dc = squeeze( CGobj.PatientComp_DVH(:,v,:) ); % complication at different volume and time points
                dt = squeeze( CGobj.PatientTotal_DVH(:,v,:) ); % complication and censored at different time points
                dt = dt-dc; % cumulative censored info 
                
                dt(:,2:end) = diff(dt,1,2); % censored info at different time points
                dt(1:end-1,:) = abs(diff(dt,1,1));
                dc(:,2:end) = diff(dc,1,2); % complicated info at different time points
                dc(1:end-1,:) = abs(diff(dc,1,1));
                
                % extract the volume, complication date, and censor info
                n = 0;
                dx(:) = 0; flgcensor(:) = 0; compdate(:) = 0;
                % censor info
                [d,t] = find(dt);
                for m = 1:length(d)
                    dx(n+1 : n+dt(d(m),t(m))) = dosebins(d(m));
                    compdate(n+1 : n+dt(d(m),t(m))) = timebins(t(m));
                    flgcensor(n+1 : n+dt(d(m),t(m))) = true;
                    n = n + dt(d(m),t(m));
                end
                % complication info
                [d,t] = find(dc);
                for m = 1:length(d)
                    dx(n+1 : n+dc(d(m),t(m))) = dosebins(d(m));
                    compdate(n+1 : n+dc(d(m),t(m))) = timebins(t(m));
                    n = n + dc(d(m),t(m));
                end

                % cox model
                try
                    if length(unique(dx(1:n)))>1
                        [~,logl,h,stats]=coxphfit(dx(1:n),compdate(1:n),'baseline',0,'censoring',flgcensor(1:n));
                        stats.logl=logl; stats.h=h;
                    else
                        stats.beta=-inf; stats.logl=-inf; stats.h=-inf;
                    end
                catch
                    stats.beta=-inf; stats.logl=-inf; stats.h=-inf;
                end
                stats.data_exposure = dx(1:n); stats.data_hazard = compdate(1:n);
                CoxDxAtlas(v)=stats;
            end
            
            % save result
            f = cellfun(@(x) strcmpi('VxAtlas',x),CGobj.CoxParameters_DVH(:,1));
            if any(f)
                CGobj.CoxParameters_DVH{f,2}=CoxVxAtlas;
            else
                CGobj.CoxParameters_DVH{end+1,1}='VxAtlas';
                CGobj.CoxParameters_DVH{end,2}=CoxVxAtlas;
            end
            f = cellfun(@(x) strcmpi('DxAtlas',x),CGobj.CoxParameters_DVH(:,1));
            if any(f)
                CGobj.CoxParameters_DVH{f,2}=CoxDxAtlas;
            else
                CGobj.CoxParameters_DVH{end+1,1}='DxAtlas';
                CGobj.CoxParameters_DVH{end,2}=CoxDxAtlas;
            end
            
            warning('on');
%             warning('on','MATLAB:singularMatrix');
%             warning('on','stats:coxphfit:FitWarning');
%             warning('on','stats:coxphfit:RankDeficient');
%             warning('on','stats:coxphfit:IterOrEvalLimit');
        end
        function CGobj = LogRankTestDxExact_DVH(CGobj)
            % select patients with data
            f = CGobj.PatientsWithComplicationData();
            CG = CGobj.RemovePatient(~f);
            
            % survival/complication time
            f2 = ~cellfun('isempty',{CG.ptGrp.CompOccurDate}); % patients with no complication date
            f3 = ~cellfun('isempty',{CG.ptGrp.LastFollowupDate}); % patients with no last follow up date
            compdate = inf(CG.numGrp,1);
            lastfollowup = inf(CG.numGrp,1);
            compdate(f2) = ([CG.ptGrp(f2).CompOccurDate] - [CG.ptGrp(f2).BaselineDate])' / 30;
            lastfollowup(f3) = ([CG.ptGrp(f3).LastFollowupDate] - [CG.ptGrp(f3).BaselineDate])' / 30;
            compdate = min( lastfollowup, compdate );
            flgcensor = [CG.ptGrp.flgCensor]';

            % prepare
                volnum=length(CG.VolBins_DVH);
                dosenum=length(CG.DoseBins_DVH);
                Dxmat = ones(dosenum,volnum,6); % n1,c1,n2,c2,p,flg (0 -- positive corelation, 1 -- negative corelation, 2 -- not available)
                Dxmat(:,:,6) = 2; % default is not available
                
                sa=classSurvivalAnalysis(); % initialize a survivalanalysis obj
                
            % Dx
                dv=zeros(CG.numGrp,1);
                numstart=CG.LogRankMinSize;
                for v=1:volnum
                    % doses under v
                    dv(:) = 0;
                    for k = 1:CG.numGrp
                        dv(k) = CG.ptGrp(k).DoseAtVol( CG.VolBins_DVH(v) );
                    end
                    
                    % Dx at (di,vj)
                    g=find(dv); % patients with non-zero doses
                    flg_dosebelow2=-1;
                    numend=length(g)-numstart;
                    for d=1:dosenum
                        % check the sample size
                        flg_dosebelow1=dv(g)<=CG.DoseBins_DVH(d); f=length(find(flg_dosebelow1));
                        if f<numstart || f>numend % one group has too few patients, skip it
                            continue;
                        end
                        % check the change of grouping
                        if isequal(flg_dosebelow1,flg_dosebelow2) % if it is the same grouping, skip the computation and save the result directly
                            Dxmat(d,v,:)=Dxmat(d-1,v,:);
                            continue;
                        end
                        flg_dosebelow2=flg_dosebelow1; % keep the current grouping
                        
                        % assign properties of object sa
                        survivedate={compdate(g(flg_dosebelow1)); compdate(g(~flg_dosebelow1))}; % survive time of each group
                        fcensor={flgcensor(g(flg_dosebelow1));flgcensor(g(~flg_dosebelow1))}; % censor flag for each group
                        sa.SurvivalTime=survivedate;
                        sa.flgCensor=fcensor;
                        % compute survival curves and compare them
                        sa=sa.CalculateSurvivalCurve();
                        sa=sa.CombineSurvivalTime();
                        sa=sa.CompareSurvivalByLogrank();
                        Dxmat(d,v,1:5)=[length(survivedate{1,1}),sum(fcensor{1,1}), length(survivedate{2,1}),sum(fcensor{2,1}),sa.pValue];
                        Dxmat(d,v,6)=sa.CurveArea(1)<sa.CurveArea(2); % the group with lower volume had worse survival curve, record it
                    end
                end
            % save reslt
                if isempty(CGobj.LogRank_DVH)
                        CGobj.LogRank_DVH{end+1,1}='Dx';
                        CGobj.LogRank_DVH{end,2}=Dxmat;
                else
                    f = cellfun(@(x) strcmpi('Dx',x),CGobj.LogRank_DVH(:,1));
                    if any(f)
                        CGobj.LogRank_DVH{f,2}=Dxmat;
                    else
                        CGobj.LogRank_DVH{end+1,1}='Dx';
                        CGobj.LogRank_DVH{end,2}=Dxmat;
                    end
                end
        end
        function CGobj = LogRankTestVxExact_DVH(CGobj)
            % select patients with data
            f = CGobj.PatientsWithComplicationData();
            CG = CGobj.RemovePatient(~f);
            
            % survival/complication time
            f2 = ~cellfun('isempty',{CG.ptGrp.CompOccurDate}); % patients with no complication date
            f3 = ~cellfun('isempty',{CG.ptGrp.LastFollowupDate}); % patients with no last follow up date
            compdate = inf(CG.numGrp,1);
            lastfollowup = inf(CG.numGrp,1);
            compdate(f2) = ([CG.ptGrp(f2).CompOccurDate] - [CG.ptGrp(f2).BaselineDate])' / 30;
            lastfollowup(f3) = ([CG.ptGrp(f3).LastFollowupDate] - [CG.ptGrp(f3).BaselineDate])' / 30;
            compdate = min( lastfollowup, compdate );
            flgcensor = [CG.ptGrp.flgCensor]';

            % prepare
                volnum=length(CG.VolBins_DVH);
                dosenum=length(CG.DoseBins_DVH);
                Vxmat = ones(dosenum,volnum,6); % (n1,c1,n2,c2,p,flg (0 -- positive corelation, 1 -- negative corelation, 2 -- not available))
                Vxmat(:,:,6) = 2;
                
                sa=classSurvivalAnalysis(); % initialize a survivalanalysis obj
                
            % Vx
                vd=zeros(CG.numGrp,1); % volume v at dose d
                numstart=CG.LogRankMinSize;
                for d=1:dosenum
                    % volume under d
                    vd(:)=0;
                    for k=1:CG.numGrp
                        vd(k) = CG.ptGrp(k).VolAtDose( CG.DoseBins_DVH(d) );
                    end
                    g=find(vd); % non-zeros volume cases
                    
                    % (di,vj)
                    flg_volbelow2=-1;
                    numend=length(g)-numstart;
                    for v=1:volnum
                        % check smaple size
                        flg_volbelow1=vd(g)<=CG.VolBins_DVH(v); f=length(find(flg_volbelow1)); % group DVHs by (d,v)
                        if f<numstart || f>numend % one group has too less patients, or the volume at the dose is too small, skip it
                            continue;
                        end
                        % check the change of grouping
                        if isequal(flg_volbelow1,flg_volbelow2) % if it is the same grouping, skip the computation and save the result directly
                            Vxmat(d,v,:)=Vxmat(d,v-1,:);
                            continue;
                        end
                        flg_volbelow2=flg_volbelow1; % keep the current grouping

                        % assign properties of object sa
                        survivedate={compdate(g(flg_volbelow1)); compdate(g(~flg_volbelow1))}; % survive time of each group
                        fcensor={flgcensor(g(flg_volbelow1)); flgcensor(g(~flg_volbelow1))}; % censor flag for each group
                        sa.SurvivalTime=survivedate;
                        sa.flgCensor=fcensor;
                        % compute survival curves and compare them
                        sa=sa.CalculateSurvivalCurve();
                        sa=sa.CombineSurvivalTime();
                        sa=sa.CompareSurvivalByLogrank();
                        Vxmat(d,v,1:5)=[length(survivedate{1,1}),sum(fcensor{1,1}), length(survivedate{2,1}),sum(fcensor{2,1}),sa.pValue];
                        Vxmat(d,v,6)=sa.CurveArea(1)<sa.CurveArea(2); % the group with lower volume had worse survival curve, record it
                    end
                end
                
            % save reslt
                if isempty(CGobj.LogRank_DVH)
                        CGobj.LogRank_DVH{end+1,1}='Vx';
                        CGobj.LogRank_DVH{end,2}=Vxmat;
                else
                    f = cellfun(@(x) strcmpi('Vx',x),CGobj.LogRank_DVH(:,1));
                    if any(f)
                        CGobj.LogRank_DVH{f,2}=Vxmat;
                    else
                        CGobj.LogRank_DVH{end+1,1}='Vx';
                        CGobj.LogRank_DVH{end,2}=Vxmat;
                    end
                end
        end
        function CGobj = LogRankTestDVxExact_DVH(CGobj)
            % select patients with data
            f = CGobj.PatientsWithComplicationData();
            CG = CGobj.RemovePatient(~f);
            
            % survival/complication time
            f2 = ~cellfun('isempty',{CG.ptGrp.CompOccurDate}); % patients with no complication date
            f3 = ~cellfun('isempty',{CG.ptGrp.LastFollowupDate}); % patients with no last follow up date
            compdate = inf(CG.numGrp,1);
            lastfollowup = inf(CG.numGrp,1);
            compdate(f2) = ([CG.ptGrp(f2).CompOccurDate] - [CG.ptGrp(f2).BaselineDate])' / 30;
            lastfollowup(f3) = ([CG.ptGrp(f3).LastFollowupDate] - [CG.ptGrp(f3).BaselineDate])' / 30;
            compdate = min( lastfollowup, compdate );
            flgcensor = [CG.ptGrp.flgCensor]';

            % prepare
                volnum=length(CG.VolBins_DVH);
                dosenum=length(CG.DoseBins_DVH);
                DVxmat = ones(dosenum,volnum,6); % n1,c1,n2,c2,p,flg (0 -- positive corelation, 1 -- negative corelation, 2 -- not available)
                DVxmat(:,:,6) = 2; % default is not available
                
                sa=classSurvivalAnalysis(); % initialize a survivalanalysis obj
                
            % Dx
                dv=zeros(CG.numGrp,1);
                numstart=CG.LogRankMinSize;
                numend=CG.numGrp-numstart;
                for v=1:volnum
                    % doses under v
                    dv(:) = 0;
                    for k = 1:CG.numGrp
                        dv(k) = CG.ptGrp(k).DoseAtVol( CG.VolBins_DVH(v) );
                    end
                    
                    % DVx at (dj,vj)
                    flg_dosebelow2=-1;
                    for d=1:dosenum
                        % check the sample size
                        flg_dosebelow1=dv<=CG.DoseBins_DVH(d); f=length(find(flg_dosebelow1));
                        if f<numstart || f>numend % one group has too few patients, skip it
                            continue;
                        end
                        % check the change of grouping
                        if isequal(flg_dosebelow1,flg_dosebelow2) % if it is the same grouping, skip the computation and save the result directly
                            DVxmat(d,v,:)=DVxmat(d-1,v,:);
                            continue;
                        end
                        flg_dosebelow2=flg_dosebelow1; % keep the current grouping
                        
                        % assign properties of object sa
                        survivedate={compdate(flg_dosebelow1); compdate(~flg_dosebelow1)}; % survive time of each group
                        fcensor={flgcensor(flg_dosebelow1);flgcensor(~flg_dosebelow1)}; % censor flag for each group
                        sa.SurvivalTime=survivedate;
                        sa.flgCensor=fcensor;
                        % compute survival curves and compare them
                        sa=sa.CalculateSurvivalCurve();
                        sa=sa.CombineSurvivalTime();
                        sa=sa.CompareSurvivalByLogrank();
                        DVxmat(d,v,1:5)=[length(survivedate{1,1}),sum(fcensor{1,1}), length(survivedate{2,1}),sum(fcensor{2,1}),sa.pValue];
                        DVxmat(d,v,6)=sa.CurveArea(1)<sa.CurveArea(2); % the group with lower volume had worse survival curve, record it
                    end
                end
            
            % save reslt
                if isempty(CGobj.LogRank_DVH)
                        CGobj.LogRank_DVH{end+1,1}='DVx';
                        CGobj.LogRank_DVH{end,2}=DVxmat;
                else
                    f = cellfun(@(x) strcmpi('DVx',x),CGobj.LogRank_DVH(:,1));
                    if any(f)
                        CGobj.LogRank_DVH{f,2}=DVxmat;
                    else
                        CGobj.LogRank_DVH{end+1,1}='DVx';
                        CGobj.LogRank_DVH{end,2}=DVxmat;
                    end
                end
        end
        function [allCox,flgCox,flgAnti] = CoxPar_DVH(CGobj,strCoxVx)
            f = cellfun(@(x) strcmpi(strCoxVx,x),CGobj.CoxParameters_DVH(:,1)); % search the label
            if any(f)
                allCox = CGobj.CoxParameters_DVH{f,2}; % extract Cox model result
                flgCox = ~arrayfun( @(y) any(structfun(@(x) any(isempty(x(:)))|any(isinf(x(:))), y)), allCox); % some fields are empty or infinite, indicating no data for those values
                % correlation
                f = [allCox.beta]';
                flgAnti = f<0;
            else
                allCox = []; flgCox = []; flgAnti = [];
            end
        end
    end

    methods % DVH plot and atlas writing
        function CoxFig_DVH(CGobj,t)
            % this is only a demo of how to use the cox model to plot response function

            % Vx
            % search the best Cox model
            [allCox,flgCox,flganti] = CoxPar_DVH(CGobj,'VDx'); % find availabe Cox models
            flgCox(flganti)=false; % anti-correlations were not be considered
            logl = [allCox.logl]'; logl(~flgCox) = -inf; % log likelihood of Cox model, anti-correlation points not counted
            [~,doseloc]=max(logl); % the best fitting of Cox model

%             num = cellfun(@(x) size(x,1),{allCox.data_exposure});
            figure(1); clf reset; plot(CGobj.DoseBins_DVH(flgCox), [allCox(flgCox).logl],'.-'); grid on;
            xlabel('Dose (Gy)','fontsize',18); ylabel('log likelihood','fontsize',18);
            figure(2); clf reset; semilogy(CGobj.DoseBins_DVH(flgCox),[allCox(flgCox).p],'.-'); grid on;
            xlabel('Dose (Gy)','fontsize',18); ylabel('p-value','fontsize',18);
            p = ones(size(allCox));
            p(flgCox) = [allCox(flgCox).p]';
            [m,ploc] = min(p);
            disp(['best fit by log likelihood and p-value: ',num2str([doseloc,ploc])]);
            
%             doseloc = 31;
            allCox = allCox(doseloc);
            disp([CGobj.DoseBins_DVH(doseloc),m,allCox.p]); % display the selected Cox model
            
            % median patient complication time
            f = CGobj.PatientsWithComplicationData(); % select patients with data
            CG = CGobj.RemovePatient(~f);
            f2 = ~cellfun('isempty',{CG.ptGrp.CompOccurDate}); % patients with no complication date
            compdate = inf(CG.numGrp,1);
            compdate(f2) = ([CG.ptGrp(f2).CompOccurDate] - [CG.ptGrp(f2).BaselineDate])' / 30;
            t = median(compdate(isfinite(compdate)));
%             t = 50;
            flg = compdate > t;

%             f3 = ~cellfun('isempty',{CG.ptGrp.LastFollowupDate}); % patients with no last follow up date
%             lastfollowup = inf(CG.numGrp,1);
%             lastfollowup(f3) = ([CG.ptGrp(f3).LastFollowupDate] - [CG.ptGrp(f3).BaselineDate])' / 30;
%             compdate = min( lastfollowup, compdate );

            % observed data at time t
            % volumes of patients
            vd=zeros(CG.numGrp,1);
            d = CG.DoseBins_DVH(doseloc);
            for k=1:CG.numGrp
                vd(k) = CG.ptGrp(k).VolAtDose( d );
            end
            % quatiles of patients
            numintv = 4;
            [sortQ,indxQ,indxorg] = EqualIntervals(vd,numintv);
            meanvol = zeros(numintv,1);
            prob = zeros(numintv,1);
            stdprob = zeros(numintv,1);
            erroru = zeros(numintv,1);
            errorl = zeros(numintv,1);
            f = ~flg(indxorg);
            for m = 1 : numintv
                meanvol(m) = mean(sortQ(indxQ==m));
                numcomp = sum(f(indxQ==m)); numtotal = length(find(indxQ==m));
                prob(m) = numcomp/numtotal;
                stdprob(m) = sqrt(numtotal*prob(m)*(1-prob(m)))/numtotal;
                erroru(m) = betainv( .84, numcomp+1, numtotal - numcomp + 1 );
                errorl(m) = betainv( .16, numcomp+1, numtotal - numcomp + 1 );
            end
            figure(3); clf reset; hold on;
            errorbar(meanvol,prob,max(0,prob-errorl),max(0,erroru-prob),'r*','linewidth',1,'markersize',12);

            % plot the cox model at time point t
            % h(t)
            f = find(diff(allCox.h(:,1))==0); % find duplicate time values of h(t)
            allCox.h(f+1,1) = allCox.h(f,1)+1e-6; % adjust it a bit to avoid ambiguius
            h = interp1(allCox.h(:,1),allCox.h(:,2),t,'linear','extrap');
            % cox curve
            vol = 0:meanvol(end)+10;
            yfit = h * exp(allCox.beta * vol);
%             yfit = log(yfit);
            plot(vol,yfit); grid on;
            yfitu = h * exp((allCox.beta+1.96*allCox.se) * vol); % upper CI
            plot(vol,yfitu,'--');
            yfitl = h * exp((allCox.beta-1.96*allCox.se) * vol); % lower CI
            plot(vol,yfitl,'--');
        end
        function DVHCurvesSummary_DVH(CGobj)
            % prepare
            f = [CGobj.ptGrp.flgCensor]; % censor info
            dosebins = cat(1,CGobj.ptGrp.DoseBins_LQ); % all doses
            dosebins = (0:max(dosebins))'; % dose bins for complication patients
            vol_censor = zeros(length(dosebins),5); % 5 columns for -95%, -68%, median, 68%, 95% lines
            vol_comp = vol_censor;

            % volume computation
            vol = -inf(CGobj.numGrp,1);
            for kk = 1:length(dosebins)
                % volumes of each patient at dose kk
                vol(:) = -inf;
                for mm = 1:CGobj.numGrp
                    vol(mm) = CGobj.ptGrp(mm).VolAtDose( dosebins(kk) );
                end
                % volumes of censored patients
                v = sort(vol(f)); v(v==0) = []; % censored patients, but patients with no volumes at current dose shall be excluded
                if ~isempty(v)
                    vl = length(v);
                    num95 = round((vl-0.95*vl)/2);
                    num68 = round((vl-0.68*vl)/2);
                    if num95>0
                        vol_censor(kk,1) = v(num95);
                        vol_censor(kk,5) = v(vl-num95);
                    else
                        vol_censor(kk,1) = -inf;
                        vol_censor(kk,5) = -inf;
                    end
                    if num68>0
                        vol_censor(kk,2) = v(num68);
                        vol_censor(kk,4) = v(vl-num68);
                    else
                        vol_censor(kk,2) = -inf;
                        vol_censor(kk,4) = -inf;
                    end
                    vol_censor(kk,3) = median(v);
                end
                % volumes of complication patients
                v = sort(vol(~f)); v(v==0) = []; % complication patients, but patients with no volumes at current dose shall be excluded
                if ~isempty(v)
                    vl = length(v);
                    num95 = round((vl-0.95*vl)/2);
                    num68 = round((vl-0.68*vl)/2);
                    if num95>0
                        vol_comp(kk,1) = v(num95);
                        vol_comp(kk,5) = v(vl-num95);
                    else
                        vol_comp(kk,1) = -inf;
                        vol_comp(kk,5) = -inf;
                    end
                    if num68>0
                        vol_comp(kk,2) = v(num68);
                        vol_comp(kk,4) = v(vl-num68);
                    else
                        vol_comp(kk,2) = -inf;
                        vol_comp(kk,4) = -inf;
                    end
                    vol_comp(kk,3) = median(v);
                end
            end
            
            % plot the curves
            hold on;
            plot(dosebins,vol_censor(:,3),'b.-');
            plot(dosebins,vol_censor(:,[2,4]),'b');
            plot(dosebins,vol_censor(:,[1,5]),'b--');
            plot(dosebins,vol_comp(:,3),'r.-');
            plot(dosebins,vol_comp(:,[2,4]),'r');
            plot(dosebins,vol_comp(:,[1,5]),'r--');
            set(gca,'xminortick','on','yminortick','on');
            set(gca,'box','on');
            xlabel('Dose (Gy)'); ylabel('Volume (cc)');
        end
    end

    methods % EUD operations
        function CGobj = CalculateDoseBins_EUD(CGobj)
            euds = [CGobj.ptGrp.EUD]'; dmax = max(euds(:));
            CGobj.DoseBins_EUD = (0:CGobj.DoseStep_EUD:dmax+CGobj.DoseStep_EUD)';
        end
        function CGobj = CalculateDoseBinsLog_EUD(CGobj)
            euds = [CGobj.ptGrp.EUD]'; dmax = log10(max(euds(:)));
%             f = euds(:)>0; dmin = log10(min(euds(f)));
%             dosebinslog = ((dmin-CGobj.DoseStep_EUD): CGobj.DoseStep_EUD:(dmax+CGobj.DoseStep_EUD))';
            dosebinslog = 0 : CGobj.DoseStep_EUD : (dmax+CGobj.DoseStep_EUD)';
            CGobj.DoseBins_EUD = 10.^dosebinslog;
        end
        function CGobj = CalculateEUD(CGobj)
            for k = 1:CGobj.numGrp
                CGobj.ptGrp(k) = CGobj.ptGrp(k).CalculateEUD();
            end
        end
        function CGobj = CrudeAtlas_EUD(CGobj)
            % extract all EUDs from all patients
                euds = [CGobj.ptGrp.EUD]';
%                 
%                 euds = zeros( CGobj.numGrp, size(CGobj.lgn,1) );
%                 for k = 1:size(CGobj.lgn,1)
%                     euds(k,:) = CGobj.ptGrp(k).EUD';
%                 end
%                 dmax = max(euds(:));
                
            % generate dose bins if DoseBins_EUD is not specifically assigned
                if isempty(CGobj.DoseBins_EUD)
                    error('DoseBins_EUD not determined before method "CrudeAtlas_EUD" is called (classComplicationGroup)');
                end
                numdosebins=size(CGobj.DoseBins_EUD,1);
                
            % censor and complication info
                flgcensor = [CGobj.ptGrp.flgCensor]; flgcomp = ~flgcensor; % a patient either was censored or had complication
                
            % for each log10(n) and each dose step, compute the total patients and their complications
                CGobj.PatientTotal_EUD = zeros( numdosebins, size(CGobj.lgn,1) );
                CGobj.PatientComp_EUD = CGobj.PatientTotal_EUD;
                for n = 1:size(CGobj.lgn,1)
                    for m = 1:numdosebins
                        f = find( euds(:,n) >= CGobj.DoseBins_EUD(m) );
                        g = find( flgcomp(f) );
                        CGobj.PatientTotal_EUD(m,n) = length(f);
                        CGobj.PatientComp_EUD(m,n) = length(g);
                    end
                end
        end
        function CGobj = BetaCumulativeProbability_EUD(CGobj)
            CGobj.BetaCumulativeMat_EUD = zeros( [size(CGobj.PatientTotal_EUD), size(CGobj.BetaCumulativeThreshold_EUD,1)] );
            for k = 1:size(CGobj.BetaCumulativeThreshold_EUD,1)
                CGobj.BetaCumulativeMat_EUD(:,:,k) = betacdf( CGobj.BetaCumulativeThreshold_EUD(k), CGobj.PatientComp_EUD+1, CGobj.PatientTotal_EUD - CGobj.PatientComp_EUD + 1 );
            end
        end
        function CGobj = BetaInverseProbability_EUD(CGobj)
            CGobj.BetaInverseMat_EUD = zeros( [size(CGobj.PatientTotal_EUD), size(CGobj.BetaInverseThreshold_EUD,1)] );
            for k=1:length(CGobj.BetaInverseThreshold_EUD)
                CGobj.BetaInverseMat_EUD(:,:,k) = betainv( CGobj.BetaInverseThreshold_EUD(k), CGobj.PatientComp_EUD+1, CGobj.PatientTotal_EUD - CGobj.PatientComp_EUD + 1 );
            end
        end
        function CGobj = LogisticRegressionAnalysisExact_EUD(CGobj)
            CGobj.LogisticRegressionMatExact_EUD = repmat( struct('b',[],'dev',[],'stats',[]), [size(CGobj.lgn,1),1] );
            
            % using exact EUD
            euds = [CGobj.ptGrp.EUD]';
            pttotal = ones(CGobj.numGrp,1);
            ptcomp = ones(CGobj.numGrp,1); ptcomp([CGobj.ptGrp.flgCensor])=0;
            for k=1:size(CGobj.lgn,1)
                doses=euds(:,k);
                % regression using exact EUD
                [b,dev,s]=glmfit(doses,[ptcomp pttotal],'binomial','link','logit');
                CGobj.LogisticRegressionMatExact_EUD(k).b=b;
                CGobj.LogisticRegressionMatExact_EUD(k).dev=dev;
                CGobj.LogisticRegressionMatExact_EUD(k).stats=s;
            end
            
        end
        function CGobj = LogisticRegressionAnalysisBin_EUD(CGobj)
            CGobj.LogisticRegressionMatBin_EUD = repmat( struct('b',[],'dev',[],'stats',[]), [size(CGobj.lgn,1),1] );
            
            % using bins
            if isempty(CGobj.DoseBins_EUD) || isempty(CGobj.PatientTotal_EUD)
                error('DoseBins_EUD or PatientTotal_EUD not determined before method "LogisticRegressionAnalysisBin_EUD" is called (classComplicationGroup)');
            end
            
            dosebins = (CGobj.DoseBins_EUD(1:end-1)+CGobj.DoseBins_EUD(2:end))/2; % dose bins are at the middle of the intervals
            pttotal = ones( CGobj.PatientTotal_EUD(1,1), 1 ); % each patient has his own row
            ptcomp = true( CGobj.PatientTotal_EUD(1,1), 1 ); % allocate space for complication of each patient
            doseval = -inf(CGobj.numGrp,1); % allocate space for dose of each patient.

            % patient complication and censor info
            pta = CGobj.PatientTotal_EUD; % patient total from atlas
            pca = CGobj.PatientComp_EUD; % patient complication from atlas
            pta = pta - pca; % censored info of patients
            
            pta(1:end-1,:) = abs(diff(pta,1,1)); pta(end,:) = 0; % locations of patients with censor info
            pca(1:end-1,:) = abs(diff(pca,1,1)); pca(end,:) = 0; % locations of patient complications
            
            warning('off','stats:glmfit:IterationLimit');
            for k=1:length(CGobj.lgn) % Logistic Regression for each ln(n)
                % determine the dose for each censored patient
                n = 0;
                f = find( pta(:,k) );
                for m = 1:length(f)
                    doseval(n+1:n+pta(f(m),k)) = dosebins(f(m));
                    ptcomp(n+1:n+pta(f(m),k)) = false;
                    n = n + pta(f(m),k);
                end
                % determine the dose for each complicated patient
                f = find( pca(:,k) );
                for m = 1:length(f)
                    doseval(n+1:n+pca(f(m),k)) = dosebins(f(m));
                    n = n + pca(f(m),k);
                end
                
                % logistic regression
                [b,dev,s]=glmfit(doseval,[ptcomp pttotal],'binomial','link','logit');
                CGobj.LogisticRegressionMatBin_EUD(k).b=b;
                CGobj.LogisticRegressionMatBin_EUD(k).dev=dev;
                CGobj.LogisticRegressionMatBin_EUD(k).stats=s;
            end
            warning('on','stats:glmfit:IterationLimit');



%             CGobj.LogisticRegressionMatBin_EUD = repmat( struct('b',[],'dev',[],'stats',[]), [size(CGobj.lgn,1),1] );
%             
%             % using bins
%             if isempty(CGobj.DoseBins_EUD)
%                 error('DoseBins_EUD not determined before method "LogisticRegressionAnalysis_EUD" is called (classComplicationGroup)');
%             end
%             dosebins = (CGobj.DoseBins_EUD(1:end-1)+CGobj.DoseBins_EUD(2:end))/2; % dose bins are at the middle of the intervals
%             pttotal = ones( CGobj.numGrp, 1 ); % each patient has his own row
%             ptcomp = ones( CGobj.numGrp, 1 ); ptcomp(CGobj.ptGrp.flgCensor)=0; % allocate space for complication of each patient
%             doseval = -inf(CGobj.numGrp,1); % allocate space for dose of each patient.
%             euds = [CGobj.ptGrp.EUD]';
%             warning('off','stats:glmfit:IterationLimit');
%             for k=1:length(CGobj.lgn)
%                 % determine the dose for each patient
%                 for m = 1:length(CGobj.DoseBins_EUD)-1
%                     f = euds(:,k)>=CGobj.DoseBins_EUD(m);
%                     doseval(f) = dosebins(m);
%                 end
%                 f = euds(:,k)>=CGobj.DoseBins_EUD(end); % patient with dose outside the required dose bins should be removed
%                 doseval(f) = min(CGobj.DoseBins_EUD)-1;
%                 
%                 f = doseval>=min(CGobj.DoseBins_EUD); % pick up the patients whose dose fall into the dose bins
%                 % logistic regression
%                 [b,dev,s]=glmfit(doseval(f),[ptcomp(f) pttotal(f)],'binomial','link','logit');
%                 CGobj.LogisticRegressionMatBin_EUD(k).b=b;
%                 CGobj.LogisticRegressionMatBin_EUD(k).dev=dev;
%                 CGobj.LogisticRegressionMatBin_EUD(k).stats=s;
%             end
%             warning('on','stats:glmfit:IterationLimit');
        end
        function CGobj = LymanModelAnalysisGridExact_EUD(CGobj)
            % preparation
            euds = [CGobj.ptGrp.EUD]';
            flg = [CGobj.ptGrp.flgCensor]';
            
            TD50 = [0:0.1:CGobj.DoseBins_EUD(end)]';
            m = [0:0.01:10]'; m(1) = 0.001;
%             m = [0:0.01:1, 1.1:0.1:2, 3:1:10]; m(1) = 0.001;
%             lgm = -1:0.01:1;
%             m = 10.^lgm';

            % for each lgn,TD50, and m, compute the log likelihood
            loglikelihood = -inf(length(TD50),length(m),length(CGobj.lgn));
            for kk = 1:length(CGobj.lgn)
                for jj = 1:length(m)
                    for ii = 1:length(TD50)
                        pr = normcdf((euds(:,kk)-TD50(ii))/(m(jj)*TD50(ii)),0,1); % Lyman probability
                        pr(flg) = 1-pr(flg); % non-complication patients
                        pr = log10(pr); % log likelihood of each patients
                        loglikelihood(ii,jj,kk) = sum(pr); % loglikelihood of all
                    end
                end
            end
            CGobj.LymanModelGridExact_EUD = struct('TD50',TD50,'m',m,'loglikelihood',loglikelihood);
        end
        function CGobj = LymanModelAnalysisGridBin_EUD(CGobj)
            % preparation
            dosebins = (CGobj.DoseBins_EUD(1:end-1)+CGobj.DoseBins_EUD(2:end))/2; % dose bins are at the middle of the intervals
            ptcomp = true( CGobj.PatientTotal_EUD(1,1), 1 ); % allocate space for complication of each patient
            doseval = -inf(CGobj.numGrp,1); % allocate space for dose of each patient.

            pta = CGobj.PatientTotal_EUD; % patient total from atlas
            pca = CGobj.PatientComp_EUD; % patient complication from atlas
            pta = pta - pca; % censored info of patients
            pta(1:end-1,:) = abs(diff(pta,1,1)); pta(end,:) = 0; % locations of patients with censor info
            pca(1:end-1,:) = abs(diff(pca,1,1)); pca(end,:) = 0; % locations of patient complications

            TD50 = [0:0.1:CGobj.DoseBins_EUD(end)]';
            m = 0:0.01:10; m(1) = 0.001;

            % for each lgn, TD50, and m, compute the log likelihood
            loglikelihood = -inf(length(TD50),length(m),length(CGobj.lgn));
            for kk = 1:length(CGobj.lgn)
                % parse patient dose and complication info
                % determine the dose for each censored patient
                ii = 0;
                f = find( pta(:,kk) );
                for jj = 1:length(f)
                    doseval(ii+1:ii+pta(f(jj),kk)) = dosebins(f(jj));
                    ptcomp(ii+1:ii+pta(f(jj),kk)) = false;
                    ii = ii + pta(f(jj),kk);
                end
                % determine the dose for each complicated patient
                f = find( pca(:,kk) );
                for jj = 1:length(f)
                    doseval(ii+1:ii+pca(f(jj),kk)) = dosebins(f(jj));
                    ii = ii + pca(f(jj),kk);
                end
                flg = ~ptcomp; % complication flag becomes censor flag

                % compute log likelihood
                for jj = 1:length(m)
                    for ii = 1:length(TD50)
                        pr = normcdf((doseval-TD50(ii))/(m(jj)*TD50(ii)),0,1); % Lyman probability
                        pr(flg) = 1-pr(flg); % non-complication patients
                        pr = log10(pr); % log likelihood of each patients
                        loglikelihood(ii,jj,kk) = sum(pr); % loglikelihood of all
                    end
                end
            end
            CGobj.LymanModelGridBin_EUD = struct('TD50',TD50,'m',m,'loglikelihood',loglikelihood);
        end
    end

    methods % EUD plot and atlas writing
        function WriteXls_EUD(CGobj) % write results from EUD into a spread sheet
            warning('off','MATLAB:xlswrite:AddSheet');
            if CGobj.flgColOutput
                % EUD
                xlswrite(CGobj.xlsFile_output,CGobj.log10n',strcat(CGobj.xlsSheet,'_EUD'),'B1');
                if ~isempty(CGobj.PatientCode)
                    xlswrite(CGobj.xlsFile_output,CGobj.PatientCode(CGobj.PatientRows),strcat(CGobj.xlsSheet,'_EUD'),'A2');
                end
                if ~isempty(CGobj.EUD)
                    xlswrite(CGobj.xlsFile_output,CGobj.EUD(CGobj.PatientRows,:),strcat(CGobj.xlsSheet,'_EUD'),'B2');
                end
                
                % total patients at log10n and dose
                xlswrite(CGobj.xlsFile_output,CGobj.log10n', strcat(CGobj.xlsSheet,'_Total'),'B1');
                if ~isempty(CGobj.PatientTotal_EUD);
%                     doses = (0:size(CGobj.PatientTotal_EUD,1)-1)*CGobj.GyStep;
                    xlswrite(CGobj.xlsFile_output,CGobj.AtlasBins,strcat(CGobj.xlsSheet,'_Total'),'A2');
                    xlswrite(CGobj.xlsFile_output,CGobj.PatientTotal_EUD,strcat(CGobj.xlsSheet,'_Total'),'B2');
                end
                
                % complication patients at log10n and dose
                xlswrite(CGobj.xlsFile_output,CGobj.log10n', strcat(CGobj.xlsSheet,'_Comp'),'B1');
                if ~isempty(CGobj.PatientComp_EUD)
%                     doses = (0:size(CGobj.PatientComp_EUD,1)-1)*CGobj.GyStep;
                    xlswrite(CGobj.xlsFile_output,CGobj.AtlasBins,strcat(CGobj.xlsSheet,'_Comp'),'A2');
                    xlswrite(CGobj.xlsFile_output,CGobj.PatientComp_EUD,strcat(CGobj.xlsSheet,'_Comp'),'B2');
                end
                
                % Beta probability
                if ~isempty(CGobj.BetaCumulativeMat)
%                     doses = (0:size(CGobj.PatientComp_EUD,1)-1)*CGobj.GyStep;
                    for k=1:length(CGobj.BetaCumulativeThreshold_EUD)
                        xlswrite(CGobj.xlsFile_output, CGobj.log10n', strcat(CGobj.xlsSheet,'_prob_',num2str(CGobj.BetaCumulativeThreshold_EUD(k))),'B1');
                        xlswrite(CGobj.xlsFile_output,CGobj.AtlasBins,strcat(CGobj.xlsSheet,'_prob_',num2str(CGobj.BetaCumulativeThreshold_EUD(k))),'A2');
                        xlswrite(CGobj.xlsFile_output,CGobj.BetaCumulativeMat(:,:,k),strcat(CGobj.xlsSheet,'_prob_',num2str(CGobj.BetaCumulativeThreshold_EUD(k))),'B2');
                    end
                end
                if ~isempty(CGobj.BetaInverseMat)
%                     doses = (0:size(CGobj.PatientComp_EUD,1)-1)*CGobj.GyStep;
                    for k=1:length(CGobj.BetaInverseThreshold_EUD)
                        xlswrite(CGobj.xlsFile_output,CGobj.log10n', strcat(CGobj.xlsSheet,'_Low_',num2str(CGobj.BetaInverseThreshold_EUD(k))),'B1');
                        xlswrite(CGobj.xlsFile_output,CGobj.AtlasBins,strcat(CGobj.xlsSheet,'_Low_',num2str(CGobj.BetaInverseThreshold_EUD(k))),'A2');
                        xlswrite(CGobj.xlsFile_output,CGobj.BetaCumulativeMat(:,:,k),strcat(CGobj.xlsSheet,'_Low_',num2str(CGobj.BetaInverseThreshold_EUD(k))),'B2');
                    end
                end
            else
                % EUD
                xlswrite(CGobj.xlsFile_output, CGobj.log10n, strcat(CGobj.xlsSheet,'_EUD'),'A2');
                if ~isempty(CGobj.PatientCode)
                    xlswrite(CGobj.xlsFile_output,CGobj.PatientCode(CGobj.PatientRows)',strcat(CGobj.xlsSheet,'_EUD'),'B1');
                end
                if ~isempty(CGobj.EUD)
                    xlswrite(CGobj.xlsFile_output,CGobj.EUD(CGobj.PatientRows,:)',strcat(CGobj.xlsSheet,'_EUD'),'B2');
                end
                
                % total patients at log10n and dose
                xlswrite(CGobj.xlsFile_output, CGobj.log10n, strcat(CGobj.xlsSheet,'_Total'),'A2');
                if ~isempty(CGobj.PatientTotal_EUD);
%                     doses = (0:size(CGobj.PatientTotal_EUD,1)-1)*CGobj.GyStep;
                    xlswrite(CGobj.xlsFile_output,CGobj.AtlasBins',strcat(CGobj.xlsSheet,'_Total'),'B1');
                    xlswrite(CGobj.xlsFile_output,CGobj.PatientTotal_EUD',strcat(CGobj.xlsSheet,'_Total'),'B2');
                end
                
                % complication patients at log10n and dose
                xlswrite(CGobj.xlsFile_output, CGobj.log10n, strcat(CGobj.xlsSheet,'_Comp'),'A2');
                if ~isempty(CGobj.PatientComp_EUD)
%                     doses = (0:size(CGobj.PatientComp_EUD,1)-1)*CGobj.GyStep;
                    xlswrite(CGobj.xlsFile_output,CGobj.AtlasBins',strcat(CGobj.xlsSheet,'_Comp'),'B1');
                    xlswrite(CGobj.xlsFile_output,CGobj.PatientComp_EUD',strcat(CGobj.xlsSheet,'_Comp'),'B2');
                end
                
                % Beta probability
                if ~isempty(CGobj.BetaCumulativeMat)
%                     doses = (0:size(CGobj.PatientComp_EUD,1)-1)*CGobj.GyStep;
                    for k=1:length(CGobj.BetaCumulativeThreshold_EUD)
                        xlswrite(CGobj.xlsFile_output, CGobj.log10n, strcat(CGobj.xlsSheet,'_prob_',num2str(CGobj.BetaCumulativeThreshold_EUD(k))),'A2');
                        xlswrite(CGobj.xlsFile_output,CGobj.AtlasBins',strcat(CGobj.xlsSheet,'_prob_',num2str(CGobj.BetaCumulativeThreshold_EUD(k))),'B1');
                        xlswrite(CGobj.xlsFile_output,CGobj.BetaCumulativeMat(:,:,k)',strcat(CGobj.xlsSheet,'_prob_',num2str(CGobj.BetaCumulativeThreshold_EUD(k))),'B2');
                    end
                end
                if ~isempty(CGobj.BetaInverseMat)
%                     doses = (0:size(CGobj.PatientComp_EUD,1)-1)*CGobj.GyStep;
                    for k=1:length(CGobj.BetaInverseThreshold_EUD)
                        xlswrite(CGobj.xlsFile_output, CGobj.log10n, strcat(CGobj.xlsSheet,'_Low_',num2str(CGobj.BetaInverseThreshold_EUD(k))),'A2');
                        xlswrite(CGobj.xlsFile_output,CGobj.AtlasBins',strcat(CGobj.xlsSheet,'_Low_',num2str(CGobj.BetaInverseThreshold_EUD(k))),'B1');
                        xlswrite(CGobj.xlsFile_output,CGobj.BetaInverseMat(:,:,k)',strcat(CGobj.xlsSheet,'_Low_',num2str(CGobj.BetaInverseThreshold_EUD(k))),'B2');
                    end
                end
            end
            warning('on','MATLAB:xlswrite:AddSheet');
        end
        function EUDCurvesFig_EUD(CGobj)
            f=[CGobj.ptGrp.flgCensor]; g=find(f);
            figure(1); hold on;
            for m = 1:length(g)
                CGobj.ptGrp(g(m)).EUDCurve_a_fig();
            end
            g=find(~f);
            for m = 1:length(g)
                CGobj.ptGrp(g(m)).EUDCurve_a_fig('r');
            end
            set(gca,'YTickLabel',num2str(CGobj.lgn(end:-2:1)));
            set(gca,'xminortick','on','yminortick','on');
            set(gca,'box','on');
            xlabel('EUD'); ylabel('log_1_0(a)');
        end
        function AtlasFig_EUD(CGobj)
            if isempty(CGobj.PatientTotal_EUD)
                disp('Empty member "PatientTotal_EUD", cannot display its figure.'); return;
            end
            dosestep = 5;
            doses=CGobj.DoseBins_EUD; x=mod(doses,dosestep)==0;
            [xx,yy]=ndgrid(doses(x),1:length(CGobj.lgn)); xx=num2cell(xx); yy=num2cell(yy);
            strComp=(CGobj.PatientComp_EUD(x,:)); strTotal=(CGobj.PatientTotal_EUD(x,:)); strAtlas=arrayfun(@(a,b) strcat(num2str(a),'/',num2str(b)),strComp,strTotal,'UniformOutpu',false);
            cellfun(@(a,b,c) text(a,b,c,'fontsize',16),xx,yy,strAtlas);
            set(gca,'XLim',[xx{1,1}-1,xx{end,1}+2]); set(gca,'YLim',[1-1,length(CGobj.lgn)+1]);
            set(gca,'YTick',1:2:length(CGobj.lgn)); set(gca,'YTickLabel',CGobj.lgn(1:2:end));
            pos=get(gcf,'Position'); set(gcf,'Position',[pos(1:2), 900/3*4, 900]);
            set(gca,'xminortick','on','yminortick','on');
            set(gca,'box','on');
%             title([CGobj.xlsSheet,', number of patients and numbers with severe pneumonitis treated to at least a given EUD'],'fontsize',16);
            xlabel('gEUD doses (Gy)'); ylabel('log_1_0a'); % set(gca,'fontsize',16);
        end
        function AtlasRotatedFig_EUD(CGobj,fonts,ticks)
            if ~exist('fonts','var')
                fonts = 10;
            end
            if ~exist('ticks','var')
                ticks = 10;
            end
            
            if isempty(CGobj.PatientTotal_EUD)
                disp('Empty member "PatientTotal_EUD", cannot display its figure.'); return;
            end
            % check columns with informaiton
            % complication part
            f = diff(CGobj.PatientComp_EUD);
            for k = 1:size(CGobj.PatientComp_EUD,1) % upper side
                if any(f(k,:))
                    colcompu = k;
                    break;
                end
            end
            for k = size(CGobj.PatientComp_EUD,1)-1 : -1 : 1 % lower side
                if any(f(k,:))
                    colcompl = k+1;
                    break;
                end
            end
            % total part
            f = diff(CGobj.PatientTotal_EUD);
            for k = 1:size(CGobj.PatientTotal_EUD,1) % upper side
                if any(f(k,:))
                    coltotalu = k;
                    break;
                end
            end
            for k = size(CGobj.PatientTotal_EUD,1)-1 : -1 : 1 % right side
                if any(f(k,:))
                    coltotall = k+1;
                    break;
                end
            end
            % combination of the comp and total
            colu = min( colcompu, coltotalu );
            coll = max( colcompl, coltotall );
            % prepare to write the table in a figure
%             doses = zeros(size(CGobj.DoseBins_EUD));
            doses=CGobj.DoseBins_EUD(colu:coll);
            [xx,yy]=ndgrid(doses,1:length(CGobj.lgn)); xx=num2cell(xx); yy=num2cell(yy);
            strComp=CGobj.PatientComp_EUD(colu:coll,:); strTotal=CGobj.PatientTotal_EUD(colu:coll,:);
            strAtlas=arrayfun(@(a,b) strcat(num2str(a),'/',num2str(b)),strComp,strTotal,'UniformOutpu',false);
%             figure(1);
            clf reset;
            cellfun(@(a,b,c) text(b,a,c,'fontsize',fonts),xx,yy,strAtlas);
            set(gca,'YLim',[xx{1,1}-1,xx{end,1}+2]);
            set(gca,'XLim',[1-1,length(CGobj.lgn)+1]);
            set(gca,'XTick',1:2:length(CGobj.lgn)); set(gca,'XTickLabel',CGobj.lgn(1:2:end));
            set(gca,'fontsize',ticks);
            set(gca,'xminortick','on','yminortick','on');
            set(gca,'Position',[0.04,0.04,0.955,0.955]);
            set(gca,'box','on');
%             pos=get(gcf,'Position'); set(gcf,'Position',[pos(1:2), 900/3*4, 900]);
%             title([CGobj.xlsSheet,', number of patients and numbers with severe pneumonitis treated to at least a given EUD'],'fontsize',16);
            xlabel('gEUD doses (Gy)'); ylabel('log_1_0a'); % set(gca,'fontsize',16);
        end
        function AtlasCompactFig_EUD(CGobj,fonts,ticks)
            if ~exist('fonts','var')
                fonts = 10;
            end
            if ~exist('ticks','var')
                ticks = 10;
            end
            
            if isempty(CGobj.PatientTotal_EUD)
                disp('Empty member "PatientTotal_EUD", cannot display its figure.'); return;
            end
            
            % check duplicated columns for each n to determine the talbe size
            dupcol = logical(diff(CGobj.PatientComp_EUD));
            dupcol = dupcol | logical(diff(CGobj.PatientTotal_EUD));
            shiftpos = zeros(2,length(CGobj.lgn));
            for n = 1:length(CGobj.lgn)
                f = find(dupcol(:,n));
                shiftpos(1,n) = f(1);
                shiftpos(2,n) = f(end);
            end
            f = diff(shiftpos); fy = max(f)+1;

            % prepare data
            strComp = zeros(fy+1,length(CGobj.lgn));
            strTotal = zeros(fy+1,length(CGobj.lgn));
            for n = 1:length(CGobj.lgn)
                strComp(:,n) = CGobj.PatientComp_EUD(shiftpos(1,n):shiftpos(1,n)+fy,n);
                strTotal(:,n) = CGobj.PatientTotal_EUD(shiftpos(1,n):shiftpos(1,n)+fy,n);
            end
            strAtlas = arrayfun(@(a,b) strcat(num2str(a),'/',num2str(b)),strComp,strTotal,'UniformOutput',false);
%             strshift = num2cell(shiftpos(1,:))';
            strshift = arrayfun(@(a) num2str(CGobj.DoseBins_EUD(a)),shiftpos(1,:),'UniformOutput',false);

            % plot the table
            clf reset;
            [x,y] = ndgrid(0:fy,1:length(CGobj.lgn));
            xx = num2cell(x); yy = num2cell(y);
            cellfun(@(a,b,c) text(b,a,c,'fontsize',fonts),xx,yy,strAtlas);
            set(gca,'YLim',[xx{1,1}-1,xx{end,1}+2]);
            set(gca,'XLim',[1-1,length(CGobj.lgn)+1]);
            set(gca,'XTick',1:1:length(CGobj.lgn)); set(gca,'XTickLabel',CGobj.lgn(end:-1:1));
            set(gca,'fontsize',ticks);

%             axis manual;
            cellfun(@(b,x) text(b,-4,x,'fontsize',fonts),yy(1,:),strshift);
            set(gca,'Position',[0.04,0.1,0.956,0.89]);
            set(gca,'xminortick','off','yminortick','on');
            set(gca,'box','on');
            xlabel('log_1_0a'); ylabel('gEUD doses (Gy)'); % set(gca,'fontsize',16);
        end
        function LogisticRegressionLikelyhoodExactFig_a_EUD(CGobj,strMarker,lw)
            if ~exist('strMarker','var')
                strMarker = 'b';
            elseif ~ischar(strMarker)
                warning('classComplicationGroup:ComplicationObservedFig_EUD:strMarker','strMarker is not a string, omit');
                strMarker = 'b';
            end
            if ~exist('lw','var')
                lw = 1;
            elseif ~isnumeric(lw)
                warning('classComplicationGroup:ComplicationObservedFig_EUD:linewidth','lw is not a number, omit');
                lw = 1;
            end
            st = [CGobj.LogisticRegressionMatExact_EUD];
            dpf = [st.dev]; % deviations
            st =[st.stats];
            df = [st.dfe]; % degree of freedom
            dpf = dpf./df; % deviations per degree of freedom
            loglikelyhood = -0.5*dpf;
            [mx,loc] = max(loglikelyhood); % the maximum loglikelyhood
            loglikelyhood68 = repmat(mx-0.5/df(loc),size(CGobj.lgn));
            loglikelyhood95 = repmat(mx-0.5*1.96/df(loc),size(CGobj.lgn));
            hold on;
            plot(-CGobj.lgn, loglikelyhood,strMarker,'LineWidth',lw);
            plot(-CGobj.lgn, loglikelyhood68,strcat(strMarker(1),'--'),'LineWidth',1);
            plot(-CGobj.lgn, loglikelyhood95,strcat(strMarker(1),'-'),'LineWidth',1);
            plot(-CGobj.lgn(loc), loglikelyhood(loc),strMarker,'LineWidth',lw+2);
            hold off;
            set(gca,'xminortick','on','yminortick','on');
            set(gca,'box','on');
            xlabel('log_1_0a'); ylabel('loglikely per degree of freedom');
            disp(['best log10(a) is: ',num2str(-CGobj.lgn(loc))]);
            disp('the corresponding coefficients beta0 and beta1, 95% CI are:');
            disp(num2str([st(loc).beta, st(loc).beta-1.96*st(loc).se, st(loc).beta+1.96*st(loc).se]));
        end
        function LogisticRegressionLikelyhoodBinFig_a_EUD(CGobj,strMarker,lw)
            if ~exist('strMarker','var')
                strMarker = 'b';
            elseif ~ischar(strMarker)
                warning('classComplicationGroup:ComplicationObservedFig_EUD:strMarker','strMarker is not a string, omit');
                strMarker = 'b';
            end
            if ~exist('lw','var')
                lw = 1;
            elseif ~isnumeric(lw)
                warning('classComplicationGroup:ComplicationObservedFig_EUD:linewidth','lw is not a number, omit');
                lw = 1;
            end
            st = [CGobj.LogisticRegressionMatBin_EUD];
            dpf = [st.dev]; % deviations
            st =[st.stats];
            df = [st.dfe]; % degree of freedom
            dpf = dpf./df; % deviations per degree of freedom
            loglikelyhood = -0.5*dpf;
            plot(-CGobj.lgn, loglikelyhood,strMarker,'LineWidth',lw);
            set(gca,'xminortick','on','yminortick','on');
            set(gca,'box','on');
            xlabel('log_1_0a'); ylabel('loglikely per degree of freedom');
        end
        function LogisticRegressionPvalueExactFig_a_EUD(CGobj,strMarker,lw)
            if ~exist('strMarker','var')
                strMarker = 'b';
            elseif ~ischar(strMarker)
                warning('classComplicationGroup:ComplicationObservedFig_EUD:strMarker','strMarker is not a string, omit');
                strMarker = 'b';
            end
            if ~exist('lw','var')
                lw = 1;
            elseif ~isnumeric(lw)
                warning('classComplicationGroup:ComplicationObservedFig_EUD:linewidth','lw is not a number, omit');
                lw = 1;
            end
            st = [CGobj.LogisticRegressionMatExact_EUD];
            st =[st.stats];
            pvalue = [st.p];
            pvalue = pvalue(2,:); % the p-value corresponding to gEUD
            [~,loc] = min(pvalue); % the location of mininum p-value
            semilogy(-CGobj.lgn, pvalue,strMarker,'LineWidth',lw);
            hold on;
            semilogy(-CGobj.lgn(loc), pvalue(loc),strMarker,'LineWidth',lw+2);
            semilogy(-CGobj.lgn, repmat(0.05,size(CGobj.lgn)), strcat(strMarker(1),'--'),'LineWidth',1);
            hold off;
            set(gca,'xminortick','on','yminortick','on');
            set(gca,'box','on');
            xlabel('log_1_0a'); ylabel('p-value');
        end
        function LogisticRegressionPvalueBinFig_a_EUD(CGobj,strMarker,lw)
            if ~exist('strMarker','var')
                strMarker = 'b';
            elseif ~ischar(strMarker)
                warning('classComplicationGroup:ComplicationObservedFig_EUD:strMarker','strMarker is not a string, omit');
                strMarker = 'b';
            end
            if ~exist('lw','var')
                lw = 1;
            elseif ~isnumeric(lw)
                warning('classComplicationGroup:ComplicationObservedFig_EUD:linewidth','lw is not a number, omit');
                lw = 1;
            end
            st = [CGobj.LogisticRegressionMatBin_EUD];
            st =[st.stats];
            pvalue = [st.p];
            pvalue = pvalue(2,:); % the p-value corresponding to gEUD
            semilogy(-CGobj.lgn, pvalue,strMarker,'LineWidth',lw);
            hold on;
            semilogy(-CGobj.lgn, repmat(0.05,size(CGobj.lgn)), strcat(strMarker(1),'--'),'LineWidth',1);
            hold off;
            set(gca,'xminortick','on','yminortick','on');
            set(gca,'box','on');
            xlabel('log_1_0a'); ylabel('p-value');
        end
        function LogisticRegressionRespondingCurveExactFig_a_EUD(CGobj,loga,strMarker,lw)
            % parameters
            if ~exist('strMarker','var')
                strMarker = 'b';
            elseif ~ischar(strMarker)
                warning('classComplicationGroup:ComplicationObservedFig_EUD:strMarker','strMarker is not a string, omit');
                strMarker = 'b';
            end
            if ~exist('lw','var')
                lw = 1;
            elseif ~isnumeric(lw)
                warning('classComplicationGroup:ComplicationObservedFig_EUD:linewidth','lw is not a number, omit');
                lw = 1;
            end
            if ~exist('loga','var')
                st = [CGobj.LogisticRegressionMatExact_EUD];
                dpf = [st.dev]; % deviations
                [~,loc] = min(dpf); % the minimum deviation corresponds to the maximum likelihood
                loga = -CGobj.lgn(loc);
            end

            % responding curve
            n = find(CGobj.lgn == -loga); % the n whose corresponding responding function will be ploted
            st = CGobj.LogisticRegressionMatExact_EUD(n); % the fitting result of that n
            euds = [CGobj.ptGrp.EUD]; euds = euds(n,:); % the gEUDs of that n
            doses = (0:max(euds))'; % doses (gEUD) whose RP probability will be computed
            [rpb,rplo,rphi] = glmval(st.b, doses,'logit',st.stats); % the responding function values at doses

            % plot
            hold on;
            plot(doses,rpb,strMarker,'LineWidth',lw); % responding function
            plot(doses,rpb-rplo,strMarker,'LineWidth',lw); % low CI curve
            plot(doses,rpb+rphi,strMarker,'LineWidth',lw); % high CI curve
            hold off;
            set(gca,'xminortick','on','yminortick','on');
            set(gca,'box','on');
            xlabel('gEUD'); ylabel('RP probability');
        end
        function LogisticRegressionRespondingCurveBinFig_a_EUD(CGobj,loga,strMarker,lw)
            % parameters
            if ~exist('strMarker','var')
                strMarker = 'b';
            elseif ~ischar(strMarker)
                warning('classComplicationGroup:ComplicationObservedFig_EUD:strMarker','strMarker is not a string, omit');
                strMarker = 'b';
            end
            if ~exist('lw','var')
                lw = 1;
            elseif ~isnumeric(lw)
                warning('classComplicationGroup:ComplicationObservedFig_EUD:linewidth','lw is not a number, omit');
                lw = 1;
            end
            if ~exist('loga','var')
                st = [CGobj.LogisticRegressionMatBin_EUD];
                dpf = [st.dev]; % deviations
                [~,loc] = min(dpf); % the minimum deviation corresponds to the maximum likelihood
                loga = -CGobj.lgn(loc);
            end

            % responding curve
            n = find(CGobj.lgn == -loga); % the n whose corresponding responding function will be ploted
            st = CGobj.LogisticRegressionMatBin_EUD(n); % the fitting result of that n
            euds = [CGobj.ptGrp.EUD]; euds = euds(n,:); % the gEUDs of that n
            doses = (0:max(euds))'; % doses (gEUD) whose RP probability will be computed
            [rpb,rplo,rphi] = glmval(st.b, doses,'logit',st.stats); % the responding function values at doses

            % plot
            hold on;
            plot(doses,rpb,strMarker,'LineWidth',lw); % responding function
            plot(doses,rpb-rplo,strMarker,'LineWidth',lw); % low CI curve
            plot(doses,rpb+rphi,strMarker,'LineWidth',lw); % high CI curve
            hold off;
            set(gca,'xminortick','on','yminortick','on');
            set(gca,'box','on');
            xlabel('gEUD'); ylabel('RP probability');
        end
        
        function LymanModelGridExactFig_TD50_m_EUD(CGobj,zdepth)
%             if ~exist('strMarker','var')
%                 strMarker = 'b';
%             elseif ~ischar(strMarker)
%                 warning('classComplicationGroup:ComplicationObservedFig_EUD:strMarker','strMarker is not a string, omit');
%                 strMarker = 'b';
%             end
%             if ~exist('lw','var')
%                 lw = 1;
%             elseif ~isnumeric(lw)
%                 warning('classComplicationGroup:ComplicationObservedFig_EUD:linewidth','lw is not a number, omit');
%                 lw = 1;
%             end

            if ~exist('zdepth','var')
                zdepth = 10;
            end
            
            % find the best "a"
            [mx,loc] = max(CGobj.LymanModelGridExact_EUD.loglikelihood(:));
            [~,~,loc] = ind2sub(size(CGobj.LymanModelGridExact_EUD.loglikelihood),loc);
%             st = [CGobj.LogisticRegressionMatExact_EUD];
%             dpf = [st.dev]; % deviations
%             st =[st.stats];
%             df = [st.dfe]; % degree of freedom
%             dpf = dpf./df; % deviations per degree of freedom
%             loglikelyhood = -0.5*dpf;
%             [~,loc] = max(loglikelyhood); % the maximum loglikelyhood
            disp(['the "log10(a)" for the maximum likelihood in Lyman model is: ',num2str(-CGobj.lgn(loc))]);
            
            % TD50 and gamma of the best "log10(a)"
%             lymangrid = CGobj.LymanModelGridExact_EUD;
            ll = CGobj.LymanModelGridExact_EUD.loglikelihood(:,:,loc);
            [dd,mm] = find(ll == mx,1);
            disp(['TD50 & m are: ',num2str([CGobj.LymanModelGridExact_EUD.TD50(dd),CGobj.LymanModelGridExact_EUD.m(mm)])]);
            
            % map of the grid of the best "a"
            ll( ll < mx-zdepth ) = -inf;
%             surf(CGobj.LymanModelGridExact_EUD.m,CGobj.LymanModelGridExact_EUD.TD50,ll,'EdgeColor','none');
%             ylabel('TD50 (Gy)'); xlabel('m'); zlabel('Log likelihood');
            contour(CGobj.LymanModelGridExact_EUD.TD50,CGobj.LymanModelGridExact_EUD.m,ll');
            hold on; plot(CGobj.LymanModelGridExact_EUD.TD50(dd),CGobj.LymanModelGridExact_EUD.m(mm),'k*'); hold off;
            xlabel('TD50 (Gy)'); ylabel('m');
        end
        function LymanModelGridBinFig_TD50_m_EUD(CGobj,zdepth)
            if ~exist('zdepth','var')
                zdepth = 10;
            end
            
            % find the best "a"
            [mx,loc] = max(CGobj.LymanModelGridBin_EUD.loglikelihood(:));
            [~,~,loc] = ind2sub(size(CGobj.LymanModelGridBin_EUD.loglikelihood),loc);
            disp(['the "log10(a)" for the maximum likelihood in Lyman model is: ',num2str(-CGobj.lgn(loc))]);
            
            % TD50 and gamma of the best "log10(a)"
            ll = CGobj.LymanModelGridBin_EUD.loglikelihood(:,:,loc);
            [dd,mm] = find(ll == mx,1);
            disp(['TD50 & m are: ',num2str([CGobj.LymanModelGridBin_EUD.TD50(dd),CGobj.LymanModelGridBin_EUD.m(mm)])]);
            
            % map of the grid of the best "a"
            ll( ll < mx-zdepth ) = -inf;
%             surf(CGobj.LymanModelGridBin_EUD.m,CGobj.LymanModelGridBin_EUD.TD50,ll,'EdgeColor','none');
%             ylabel('TD50 (Gy)'); xlabel('m'); zlabel('Log likelihood');
            contour(CGobj.LymanModelGridBin_EUD.TD50,CGobj.LymanModelGridBin_EUD.m,ll');
            hold on; plot(CGobj.LymanModelGridBin_EUD.TD50(dd),CGobj.LymanModelGridBin_EUD.m(mm),'k*'); hold off;
            xlabel('TD50 (Gy)'); ylabel('m');
        end
        function LymanModelGridExactFig_TD50_a_EUD(CGobj,zdepth)
            if ~exist('zdepth','var')
                zdepth = 10;
            end
            
            % find the best "m"
            [mx,loc] = max(CGobj.LymanModelGridExact_EUD.loglikelihood(:));
            [~,loc,~] = ind2sub(size(CGobj.LymanModelGridExact_EUD.loglikelihood),loc);
            disp(['the "m" for the maximum likelihood in Lyman model is: ',num2str(CGobj.LymanModelGridExact_EUD.m(loc))]);
            
            % TD50 and log10(a) corresponding to the best m
            ll = squeeze(CGobj.LymanModelGridExact_EUD.loglikelihood(:,loc,:));
            [dd,aa] = find(ll == mx,1);
            disp(['TD50 & a are: ',num2str([CGobj.LymanModelGridExact_EUD.TD50(dd),-CGobj.lgn(aa)])]);
            
            % map of the grid of the best "m"
            ll( ll < mx-zdepth ) = -inf;
%             surf(-CGobj.lgn,CGobj.LymanModelGridExact_EUD.TD50,ll,'EdgeColor','none');
%             ylabel('TD50 (Gy)'); xlabel('log_1_0(a)'); zlabel('Log likelihood');
            contour(-CGobj.lgn,CGobj.LymanModelGridExact_EUD.TD50,ll);
            hold on; plot(-CGobj.lgn(aa),CGobj.LymanModelGridExact_EUD.TD50(dd),'k*'); hold off;
            xlabel('log_1_0(a)'); ylabel('TD50 (Gy)');
        end
        function LymanModelGridBinFig_TD50_a_EUD(CGobj,zdepth)
            if ~exist('zdepth','var')
                zdepth = 10;
            end
            
            % find the best "m"
            [mx,loc] = max(CGobj.LymanModelGridBin_EUD.loglikelihood(:));
            [~,loc,~] = ind2sub(size(CGobj.LymanModelGridBin_EUD.loglikelihood),loc);
            disp(['the "m" for the maximum likelihood in Lyman model is: ',num2str(CGobj.LymanModelGridBin_EUD.m(loc))]);
            
            % TD50 and log10(a) corresponding to the best m
            ll = squeeze(CGobj.LymanModelGridBin_EUD.loglikelihood(:,loc,:));
            [dd,aa] = find(ll == mx,1);
            disp(['TD50 & a are: ',num2str([CGobj.LymanModelGridBin_EUD.TD50(dd),-CGobj.lgn(aa)])]);
            
            % map of the grid of the best "m"
            ll( ll < mx-zdepth ) = -inf;
%             surf(-CGobj.lgn,CGobj.LymanModelGridBin_EUD.TD50,ll,'EdgeColor','none');
%             ylabel('TD50 (Gy)'); xlabel('log_1_0(a)'); zlabel('Log likelihood');
            contour(-CGobj.lgn,CGobj.LymanModelGridBin_EUD.TD50,ll);
            hold on; plot(-CGobj.lgn(aa),CGobj.LymanModelGridBin_EUD.TD50(dd),'k*'); hold off;
            xlabel('log_1_0(a)'); ylabel('TD50 (Gy)');
        end
        function LymanModelGridExactFig_m_a_EUD(CGobj,zdepth)
            if ~exist('zdepth','var')
                zdepth = 10;
            end
            
            % find the best "m"
            [mx,loc] = max(CGobj.LymanModelGridExact_EUD.loglikelihood(:));
            [loc,~,~] = ind2sub(size(CGobj.LymanModelGridExact_EUD.loglikelihood),loc);
            disp(['the "TD50" for the maximum likelihood in Lyman model is: ',num2str(CGobj.LymanModelGridExact_EUD.TD50(loc))]);
            
            % TD50 and log10(a) corresponding to the best m
            ll = squeeze(CGobj.LymanModelGridExact_EUD.loglikelihood(loc,:,:));
            [mm,aa] = find(ll == mx,1);
            disp(['m & a are: ',num2str([CGobj.LymanModelGridExact_EUD.m(mm),-CGobj.lgn(aa)])]);
            
            % map of the grid of the best "m"
            ll( ll < mx-zdepth ) = -inf;
%             surf(-CGobj.lgn,CGobj.LymanModelGridExact_EUD.m,ll,'EdgeColor','none');
%             ylabel('m'); xlabel('log_1_0(a)'); zlabel('Log likelihood');
            contour(-CGobj.lgn,CGobj.LymanModelGridExact_EUD.m,ll);
            hold on; plot(-CGobj.lgn(aa),CGobj.LymanModelGridExact_EUD.m(mm),'k*'); hold off;
            xlabel('log_1_0(a)'); ylabel('m');
        end
        function LymanModelGridBinFig_m_a_EUD(CGobj,zdepth)
            if ~exist('zdepth','var')
                zdepth = 10;
            end
            
            % find the best "m"
            [mx,loc] = max(CGobj.LymanModelGridBin_EUD.loglikelihood(:));
            [loc,~,~] = ind2sub(size(CGobj.LymanModelGridBin_EUD.loglikelihood),loc);
            disp(['the "TD50" for the maximum likelihood in Lyman model is: ',num2str(CGobj.LymanModelGridBin_EUD.TD50(loc))]);
            
            % TD50 and log10(a) corresponding to the best m
            ll = squeeze(CGobj.LymanModelGridBin_EUD.loglikelihood(loc,:,:));
            [mm,aa] = find(ll == mx,1);
            disp(['m & a are: ',num2str([CGobj.LymanModelGridBin_EUD.m(mm),-CGobj.lgn(aa)])]);
            
            % map of the grid of the best "m"
            ll( ll < mx-zdepth ) = -inf;
%             surf(-CGobj.lgn,CGobj.LymanModelGridBin_EUD.m,ll,'EdgeColor','none');
%             ylabel('m'); xlabel('log_1_0(a)'); zlabel('Log likelihood');
            contour(-CGobj.lgn,CGobj.LymanModelGridBin_EUD.m,ll);
            hold on; plot(-CGobj.lgn(aa),CGobj.LymanModelGridBin_EUD.m(mm),'k*'); hold off;
            xlabel('log_1_0(a)'); ylabel('m');
        end
        function LymanModelGridExactFig_a_loglikelihood(CGobj,strMarker,lw)
            % parameters
            if ~exist('strMarker','var')
                strMarker = 'b-';
            elseif ~ischar(strMarker)
                warning('classComplicationGroup:ComplicationObservedFig_EUD:strMarker','strMarker is not a string, omit');
                strMarker = 'b-';
            end
            if ~exist('lw','var')
                lw = 1;
            elseif ~isnumeric(lw)
                warning('classComplicationGroup:ComplicationObservedFig_EUD:linewidth','lw is not a number, omit');
                lw = 1;
            end

            % find the best "a"
            [mx,loc] = max(CGobj.LymanModelGridExact_EUD.loglikelihood(:));
            [~,~,loc] = ind2sub(size(CGobj.LymanModelGridExact_EUD.loglikelihood),loc);
            disp(['the maximum log likelihood and its corresponding "log10(a)" in Lyman model are: ',num2str([mx, -CGobj.lgn(loc)])]);

            % plot the maximum log likelihood curve w.r.t. log10(a)
            llmx = zeros(size(CGobj.lgn)); %log likelihood's max for each log10(a)
            for kk = 1:length(llmx)
                ll = CGobj.LymanModelGridExact_EUD.loglikelihood(:,:,kk);
                llmx(kk) = max(ll(:));
            end
            llmx = llmx / (CGobj.numGrp-1);
            mx = mx / (CGobj.numGrp-1);
            loglikelyhood68 = repmat(mx-0.5/(CGobj.numGrp-1),size(llmx));
            loglikelyhood95 = repmat(mx-0.5*1.96/(CGobj.numGrp-1),size(llmx));
            hold on;
            plot(-CGobj.lgn,llmx,strMarker,'LineWidth',lw);
            plot(-CGobj.lgn, loglikelyhood68,strcat(strMarker(1),'--'),'LineWidth',lw);
            plot(-CGobj.lgn, loglikelyhood95,strcat(strMarker(1),'-'),'LineWidth',lw);
            if length(strMarker)==1
                strMarker = [strMarker,'*'];
            elseif strMarker(2) == '-' || strMarker(2) == '.'
                strMarker(2) = '*';
            end
            plot(-CGobj.lgn(loc),mx,strMarker,'LineWidth',lw+2);
            hold off;
            set(gca,'xminortick','on','yminortick','on');
            set(gca,'box','on');
            xlabel('log_1_0(a)'); ylabel('log likelihood per degree of freedom');
        end
        function LymanModelGridBinFig_a_loglikelihood(CGobj,strMarker,lw)
            % parameters
            if ~exist('strMarker','var')
                strMarker = 'b';
            elseif ~ischar(strMarker)
                warning('classComplicationGroup:ComplicationObservedFig_EUD:strMarker','strMarker is not a string, omit');
                strMarker = 'b';
            end
            if ~exist('lw','var')
                lw = 1;
            elseif ~isnumeric(lw)
                warning('classComplicationGroup:ComplicationObservedFig_EUD:linewidth','lw is not a number, omit');
                lw = 1;
            end

            % find the best "a"
            [mx,loc] = max(CGobj.LymanModelGridBin_EUD.loglikelihood(:));
            [~,~,loc] = ind2sub(size(CGobj.LymanModelGridBin_EUD.loglikelihood),loc);
            disp(['the maximum log likelihood and its corresponding "log10(a)" in Lyman model are: ',num2str([mx, -CGobj.lgn(loc)])]);

            % plot the maximum log likelihood curve w.r.t. log10(a)
            llmx = zeros(size(CGobj.lgn)); %log likelihood's max for each log10(a)
            for kk = 1:length(llmx)
                ll = CGobj.LymanModelGridBin_EUD.loglikelihood(:,:,kk);
                llmx(kk) = max(ll(:));
            end
            llmx = llmx / (CGobj.numGrp-1);
            mx = mx / (CGobj.numGrp-1);
            hold on;
            plot(-CGobj.lgn,llmx,strMarker,'LineWidth',lw);
            plot(-CGobj.lgn(loc),mx,strMarker,'LineWidth',lw+2);
            hold off;
            set(gca,'xminortick','on','yminortick','on');
            set(gca,'box','on');
            xlabel('log_1_0(a)'); ylabel('log likelihood per degree of freedom');
        end
        function LymanModelGridExactFig_TD50_loglikelihood(CGobj,strMarker,lw)
            % parameters
            if ~exist('strMarker','var')
                strMarker = 'b';
            elseif ~ischar(strMarker)
                warning('classComplicationGroup:ComplicationObservedFig_EUD:strMarker','strMarker is not a string, omit');
                strMarker = 'b';
            end
            if ~exist('lw','var')
                lw = 1;
            elseif ~isnumeric(lw)
                warning('classComplicationGroup:ComplicationObservedFig_EUD:linewidth','lw is not a number, omit');
                lw = 1;
            end

            % find the best "TD50"
            [mx,loc] = max(CGobj.LymanModelGridExact_EUD.loglikelihood(:));
            [loc,~,~] = ind2sub(size(CGobj.LymanModelGridExact_EUD.loglikelihood),loc);
            disp(['the maximum log likelihood and its corresponding "TD50" in Lyman model are: ',num2str([mx, CGobj.LymanModelGridExact_EUD.TD50(loc)])]);

            % plot the maximum log likelihood curve w.r.t. TD50
            llmx = zeros(size(CGobj.LymanModelGridExact_EUD.TD50)); %log likelihood's max for each log10(a)
            for kk = 1:length(llmx)
                ll = CGobj.LymanModelGridExact_EUD.loglikelihood(kk,:,:);
                llmx(kk) = max(ll(:));
            end
            llmx = llmx / (CGobj.numGrp-1);
            mx = mx / (CGobj.numGrp-1);
            loglikelyhood68 = repmat(mx-0.5/(CGobj.numGrp-1),size(llmx));
            loglikelyhood95 = repmat(mx-0.5*1.96/(CGobj.numGrp-1),size(llmx));
            hold on;
            plot(CGobj.LymanModelGridExact_EUD.TD50,llmx,strMarker,'LineWidth',lw);
            plot(CGobj.LymanModelGridExact_EUD.TD50, loglikelyhood68,strcat(strMarker(1),'--'),'LineWidth',lw);
            plot(CGobj.LymanModelGridExact_EUD.TD50, loglikelyhood95,strcat(strMarker(1),'-'),'LineWidth',lw);
            if length(strMarker)==1
                strMarker = [strMarker,'*'];
            elseif strMarker(2) == '-' || strMarker(2) == '.'
                strMarker(2) = '*';
            end
            plot(CGobj.LymanModelGridExact_EUD.TD50(loc),mx,strMarker,'LineWidth',lw+2);
            hold off;
            set(gca,'xminortick','on','yminortick','on');
            set(gca,'box','on');
            xlabel('TD50'); ylabel('log likelihood per degree of freedom');
        end
        function LymanModelGridBinFig_TD50_loglikelihood(CGobj,strMarker,lw)
            % parameters
            if ~exist('strMarker','var')
                strMarker = 'b';
            elseif ~ischar(strMarker)
                warning('classComplicationGroup:ComplicationObservedFig_EUD:strMarker','strMarker is not a string, omit');
                strMarker = 'b';
            end
            if ~exist('lw','var')
                lw = 1;
            elseif ~isnumeric(lw)
                warning('classComplicationGroup:ComplicationObservedFig_EUD:linewidth','lw is not a number, omit');
                lw = 1;
            end

            % find the best "TD50"
            [mx,loc] = max(CGobj.LymanModelGridBin_EUD.loglikelihood(:));
            [loc,~,~] = ind2sub(size(CGobj.LymanModelGridBin_EUD.loglikelihood),loc);
            disp(['the maximum log likelihood and its corresponding "TD50" in Lyman model are: ',num2str([mx, CGobj.LymanModelGridBin_EUD.TD50(loc)])]);

            % plot the maximum log likelihood curve w.r.t. TD50
            llmx = zeros(size(CGobj.LymanModelGridBin_EUD.TD50)); %log likelihood's max for each log10(a)
            for kk = 1:length(llmx)
                ll = CGobj.LymanModelGridBin_EUD.loglikelihood(kk,:,:);
                llmx(kk) = max(ll(:));
            end
            llmx = llmx / (CGobj.numGrp-1);
            mx = mx / (CGobj.numGrp-1);
            hold on;
            plot(CGobj.LymanModelGridBin_EUD.TD50,llmx,strMarker,'LineWidth',lw);
            plot(CGobj.LymanModelGridBin_EUD.TD50(loc),mx,strMarker,'LineWidth',lw+2);
            hold off;
            set(gca,'xminortick','on','yminortick','on');
            set(gca,'box','on');
            xlabel('TD50'); ylabel('log likelihood per degree of freedom');
        end
        function LymanModelGridExactFig_m_loglikelihood(CGobj,strMarker,lw)
            % parameters
            if ~exist('strMarker','var')
                strMarker = 'b';
            elseif ~ischar(strMarker)
                warning('classComplicationGroup:ComplicationObservedFig_EUD:strMarker','strMarker is not a string, omit');
                strMarker = 'b';
            end
            if ~exist('lw','var')
                lw = 1;
            elseif ~isnumeric(lw)
                warning('classComplicationGroup:ComplicationObservedFig_EUD:linewidth','lw is not a number, omit');
                lw = 1;
            end

            % find the best "TD50"
            [mx,loc] = max(CGobj.LymanModelGridExact_EUD.loglikelihood(:));
            [~,loc,~] = ind2sub(size(CGobj.LymanModelGridExact_EUD.loglikelihood),loc);
            disp(['the maximum log likelihood and its corresponding "m" in Lyman model are: ',num2str([mx, CGobj.LymanModelGridExact_EUD.m(loc)])]);

            % plot the maximum log likelihood curve w.r.t. TD50
            llmx = zeros(size(CGobj.LymanModelGridExact_EUD.m)); %log likelihood's max for each log10(a)
            for kk = 1:length(llmx)
                ll = CGobj.LymanModelGridExact_EUD.loglikelihood(:,kk,:);
                llmx(kk) = max(ll(:));
            end
            llmx = llmx / (CGobj.numGrp-1);
            mx = mx / (CGobj.numGrp-1);
            loglikelyhood68 = repmat(mx-0.5/(CGobj.numGrp-1),size(llmx));
            loglikelyhood95 = repmat(mx-0.5*1.96/(CGobj.numGrp-1),size(llmx));
            hold on;
            plot(CGobj.LymanModelGridExact_EUD.m,llmx,strMarker,'LineWidth',lw);
            plot(CGobj.LymanModelGridExact_EUD.m, loglikelyhood68,strcat(strMarker(1),'--'),'LineWidth',lw);
            plot(CGobj.LymanModelGridExact_EUD.m, loglikelyhood95,strcat(strMarker(1),'-'),'LineWidth',lw);
            if length(strMarker)==1
                strMarker = [strMarker,'*'];
            elseif strMarker(2) == '-' || strMarker(2) == '.'
                strMarker(2) = '*';
            end
            plot(CGobj.LymanModelGridExact_EUD.m(loc),mx,strMarker,'LineWidth',lw+2);
            hold off;
            set(gca,'xminortick','on','yminortick','on');
            set(gca,'box','on');
            xlabel('m'); ylabel('log likelihood per degree of freedom');
        end
        function LymanModelGridBinFig_m_loglikelihood(CGobj,strMarker,lw)
            % parameters
            if ~exist('strMarker','var')
                strMarker = 'b';
            elseif ~ischar(strMarker)
                warning('classComplicationGroup:ComplicationObservedFig_EUD:strMarker','strMarker is not a string, omit');
                strMarker = 'b';
            end
            if ~exist('lw','var')
                lw = 1;
            elseif ~isnumeric(lw)
                warning('classComplicationGroup:ComplicationObservedFig_EUD:linewidth','lw is not a number, omit');
                lw = 1;
            end

            % find the best "TD50"
            [mx,loc] = max(CGobj.LymanModelGridBin_EUD.loglikelihood(:));
            [~,loc,~] = ind2sub(size(CGobj.LymanModelGridBin_EUD.loglikelihood),loc);
            disp(['the maximum log likelihood and its corresponding "m" in Lyman model are: ',num2str([mx, CGobj.LymanModelGridBin_EUD.m(loc)])]);

            % plot the maximum log likelihood curve w.r.t. TD50
            llmx = zeros(size(CGobj.LymanModelGridBin_EUD.m)); %log likelihood's max for each log10(a)
            for kk = 1:length(llmx)
                ll = CGobj.LymanModelGridBin_EUD.loglikelihood(:,kk,:);
                llmx(kk) = max(ll(:));
            end
            llmx = llmx / (CGobj.numGrp-1);
            mx = mx / (CGobj.numGrp-1);
            hold on;
            plot(CGobj.LymanModelGridBin_EUD.m,llmx,strMarker,'LineWidth',lw);
            plot(CGobj.LymanModelGridBin_EUD.m(loc),mx,strMarker,'LineWidth',lw+2);
            hold off;
            set(gca,'xminortick','on','yminortick','on');
            set(gca,'box','on');
            xlabel('m'); ylabel('log likelihood per degree of freedom');
        end
        
        function ComplicationObservedFig_EUD(CGobj,loga,numIntervals,strMarker,lw)
            % parameters
            if ~exist('strMarker','var')
                strMarker = 'b';
            elseif ~ischar(strMarker)
                warning('classComplicationGroup:ComplicationObservedFig_EUD:strMarker','strMarker is not a string, omit');
                strMarker = 'b';
            end
            if ~exist('lw','var')
                lw = 1;
            elseif ~isnumeric(lw)
                warning('classComplicationGroup:ComplicationObservedFig_EUD:linewidth','lw is not a number, omit');
                lw = 1;
            end
            if ~exist('loga','var')
                st = [CGobj.LogisticRegressionMatBin_EUD];
                dpf = [st.dev]; % deviations
                [~,loc] = min(dpf); % the minimum deviation corresponds to the maximum likelihood
                loga = -CGobj.lgn(loc);
            end
            if ~exist('numIntervals','var')
                numIntervals = 4;
            end
            
            % grouping patients
            n = (CGobj.lgn == -loga); % the n whose gEUDs will be ploted
            euds = [CGobj.ptGrp.EUD]; euds = euds(n,:); % the gEUDs of that n
            [sortQ,indxQ,indxorg] = EqualIntervals(euds,numIntervals); % group patients according to their gEUDs
            flg=[CGobj.ptGrp.flgCensor]; % censor flags of patients
            meaneud = zeros(numIntervals,1);
            prob = zeros(numIntervals,1);
            stdprob = zeros(numIntervals,1);
            erroru = zeros(numIntervals,1);
            errorl = zeros(numIntervals,1);
            flg = ~flg(indxorg);
            for m = 1 : numIntervals
                meaneud(m) = mean(sortQ(indxQ==m));
                numcomp = sum(flg(indxQ==m)); numtotal = length(find(indxQ==m));
                prob(m) = numcomp/numtotal;
                stdprob(m) = sqrt(numtotal*prob(m)*(1-prob(m)))/numtotal;
                erroru(m) = betainv( .84, numcomp+1, numtotal - numcomp + 1 );
                errorl(m) = betainv( .16, numcomp+1, numtotal - numcomp + 1 );
            end
            
            % plot
            errorbar(meaneud,prob,max(0,prob-errorl),max(0,erroru-prob),strMarker,'LineWidth',lw);
            set(gca,'xminortick','on','yminortick','on');
            set(gca,'box','on');
            xlabel('gEUD'); ylabel('RP rate observed');
        end
        function ProbabilityFig_EUD(CGobj)
            % map of the probability of have complication rate of at least value (e.g. 20%)
            if isempty(CGobj.BetaCumulativeMat_EUD)
                disp('Empty member (BetaCumulativeMat), can not display its figure.'); return;
            end

            % prepare
            cm = colormap(jet(300)); cm=cm(1:256,:); %cm(end,:) = 0.5;
            imgmsk = CGobj.PatientTotal_EUD > 0; imgmsk=imgmsk';
            eudmx = CGobj.DoseBins_EUD(end); eudmn = CGobj.DoseBins_EUD(1);

            % image data 
            img=1-CGobj.BetaCumulativeMat_EUD; img=img';
            mx = ceil(max(img(imgmsk))*256);
            img(~imgmsk) = NaN;
            colormap(cm(1:mx,:));
            contourf(img); axis xy;
            colorbar;

%             set(gca,'Position',[0.05,0.05,0.8,0.9]);
            
            set(gca,'YTick',1:length(CGobj.lgn)); set(gca,'YTickLabel',-CGobj.lgn);
            xlim = get(gca,'XLim'); % x limit
            xticklabel = 0:5:eudmx; xtick = diff(xlim)/(eudmx-eudmn)*xticklabel+xlim(1);
            set(gca,'XTick',xtick);
            set(gca,'XTickLabel',xticklabel);
            
            set(gca,'xminortick','on','yminortick','on');
            set(gca,'box','on');
            xlabel('EUD doses (Gy)'); ylabel('log a');
            
%         xlim=get(gca,'XLim');
%         xtick=linspace(xlim(1),xlim(2),xticknum);
%         set(gca,'XTick',xtick);
%         set(gca,'XTickLabel',doseticks);
%         xlabel('EUD dose (Gy)');
% 
%         set(gca,'YTick',1:length(CGobjs(k).lgn)); set(gca,'YTickLabel',CGobjs(k).lgn);

%             img=1-CGobj.BetaCumulativeMat'; img(1,end)=1; img(2,end)=0;
%             contourf(img); %contourf(flipud(rot90(CGobj.BetaCumulativeMat,1)));
%             xsteps=0:CGobj.GyStep:max(CGobj.EUD(:)); xtickvec=10:10:xsteps(end); xtickstep=(size(CGobj.BetaCumulativeMat,1)/xsteps(end)*10);
%             set(gca,'XTick',xtickstep:xtickstep:size(CGobj.BetaCumulativeMat,1)); set(gca,'XTickLabel',xtickvec);
%             ytickvec=CGobj.log10n; ytickstep=size(CGobj.BetaCumulativeMat,2)/length(ytickvec);
%             set(gca,'YTick',ytickstep:ytickstep:size(CGobj.BetaCumulativeMat,2)); set(gca,'YTickLabel',ytickvec);
%             %             pmin=min(CGobj.BetaCumulativeMat(:)); pmax=max(CGobj.BetaCumulativeMat(:)); %bstep=(pmax-pmin)/10*256;
%             colorbar; %colorbar('YTickLabel',round([pmin:(pmax-pmin)/9:pmax]*10)/10);
%             title([CGobj.xlsSheet,', the probability that observed complications arise from true rate > 20%']);
        end
        function Low68pctConfidenceFig_EUD(CGobj)
            % map of lower 68% confidence limit on complication probability
            if isempty(CGobj.BetaCumulativeMat_EUD)
                disp('Empty member (BetaCumulativeMat), can not display its figure.'); return;
            end

            % prepare
            cm = colormap(jet(300)); cm=cm(1:256,:); %cm(end,:) = 0.5;
            imgmsk = CGobj.PatientTotal_EUD > 0; imgmsk=imgmsk';
            eudmx = CGobj.DoseBins_EUD(end); eudmn = CGobj.DoseBins_EUD(1);

            img=CGobj.BetaInverseMat_EUD';
            mx = ceil(max(img(imgmsk))*256);
            img(~imgmsk) = NaN;
            colormap(cm(1:mx,:));
            colormap(cm);
            contourf(img); axis xy;
            colorbar;

%             set(gca,'Position',[0.05,0.05,0.8,0.9]);
           
            set(gca,'YTick',1:length(CGobj.lgn)); set(gca,'YTickLabel',-CGobj.lgn);
            xlim = get(gca,'XLim'); % x limit
            xticklabel = 0:5:eudmx; xtick = diff(xlim)/(eudmx-eudmn)*xticklabel+xlim(1);
            set(gca,'XTick',xtick);
            set(gca,'XTickLabel',xticklabel);
            
            set(gca,'xminortick','on','yminortick','on');
            set(gca,'box','on');
            xlabel('EUD doses (Gy)'); ylabel('log a');
            
%             if isempty(CGobj.BetaInverseMat)
%                 disp('Empty member (BetaInverseMat), can not display its figure.'); return;
%             end
%             img=CGobj.BetaInverseMat'; img(1,end)=1; img(2,end)=0;
%             contourf(img); %contourf(flipud(rot90(CGobj.BetaInverseMat,1)));
%             xsteps=0:CGobj.GyStep:max(CGobj.EUD(:)); xtickvec=10:10:xsteps(end); xtickstep=(size(CGobj.BetaInverseMat,1)/xsteps(end)*10);
%             set(gca,'XTick',xtickstep:xtickstep:size(CGobj.BetaInverseMat,1)); set(gca,'XTickLabel',xtickvec);
%             ytickvec=CGobj.log10n; ytickstep=size(CGobj.BetaInverseMat,2)/length(ytickvec);
%             set(gca,'YTick',ytickstep:ytickstep:size(CGobj.BetaInverseMat,2)); set(gca,'YTickLabel',ytickvec);
%             colorbar;
%             title([CGobj.xlsSheet,', lower 68% confidence limit on complication probability']);
        end
    end
    
    methods(Static)
        function CGobj = CombineAtlas_EUD(CGobj1,CGobj2)
            % integrity check
            % lgn
            if any(abs(CGobj1.lgn-CGobj2.lgn)>1e-6)
                error('The log10(n) is not equal in the two CGobjs so they can not be combined');
            end
            % dose bins in atlas computation
            if max(CGobj1.DoseBins_EUD) >= max(CGobj2.DoseBins_EUD)
                dosebins1 = CGobj1.DoseBins_EUD;
                dosebins2 = CGobj2.DoseBins_EUD;
            else
                dosebins1 = CGobj2.DoseBins_EUD;
                dosebins2 = CGobj1.DoseBins_EUD;
            end
            d = abs( dosebins1(1:size(dosebins2,1)) - dosebins2 );
            if any(d>1e-6)
                error('The dose bins in atlas computation were not consistant so they are not combined');
            end
            
            % combination
            CGobj = CGobj1;
            CGobj = CGobj.AddPatient(CGobj2.ptGrp);
            
            CGobj.DoseBins_EUD = dosebins1;
            size1 = size(CGobj1.PatientTotal_EUD); size2 = size(CGobj2.PatientTotal_EUD);
            CGobj.PatientTotal_EUD = zeros(max(size1,size2));
            CGobj.PatientTotal_EUD(1:size1(1),1:size1(2)) = CGobj1.PatientTotal_EUD;
            CGobj.PatientTotal_EUD(1:size2(1),1:size2(2)) = CGobj.PatientTotal_EUD(1:size2(1),1:size2(2))+CGobj2.PatientTotal_EUD;
            CGobj.PatientComp_EUD = zeros(max(size1,size2));
            CGobj.PatientComp_EUD(1:size1(1),1:size1(2)) = CGobj1.PatientComp_EUD;
            CGobj.PatientComp_EUD(1:size2(1),1:size2(2)) = CGobj.PatientComp_EUD(1:size2(1),1:size2(2))+CGobj2.PatientComp_EUD;
            
%             % compute probabilities
%             CGobj=CGobj.BetaCumulativeProbability_EUD();
%             CGobj=CGobj.BetaInverseProbability_EUD();
%             CGobj=CGobj.LogisticRegressionAnalysis_EUD();
        end
    end
    
    methods % unused
        function CGobj = LymanModelAnalysisFittingExact_EUD(CGobj)
            % preparation
            fttp = fittype('normcdf((eud-TD50)/(m*TD50),0,1)', 'independent','eud', 'coefficients',{'TD50','m'});
            fttp = fittype('normcdf(eud,TD50,m)', 'independent','eud', 'coefficients',{'TD50','m'});
            ftopt = fitoptions('method','LinearLeastSquares','algorith','Levenberg-Marquardt','Display','on');
            euds = [CGobj.ptGrp.EUD]';
            flg = double([CGobj.ptGrp.flgCensor]');

%             warning('off','MATLAB:singularMatrix');
            % for each lgn, compute the parameters of Lyman model, TD50 and m
            kk = 1; % the first fit
            [dose,indx] = sort(euds(:,kk),'ascend');
            comp = ~flg(indx); f = find(comp); % rearrange the complication info so it corresponds to the dose
            comp = cumsum(comp)/CGobj.numGrp; % cumulative complication for dose
            % found the best grid as start point
                llh = CGobj.LymanModelGridExact_EUD.loglikelihood(:,:,kk); %log likeli-hood from grid computation
                mn = max(llh(:)); [dd,mm]=find(llh==mn);
            [fittedobj,goodness,output] = fit(dose,comp,fttp,fitoptions(ftopt,'StartPoint',[CGobj.LymanModelGridExact_EUD.TD50(dd(1)),CGobj.LymanModelGridExact_EUD.m(mm(1))]+10));
            [fittedobj,goodness,output] = fit(euds(:,:,kk),flg,fttp,fitoptions(ftopt,'StartPoint',[CGobj.LymanModelGridExact_EUD.TD50(dd(1)),CGobj.LymanModelGridExact_EUD.m(mm(1))]+10));
            [fittedobj,goodness,output] = fit(euds(:,:,kk),flg,fttp,'StartPoint',[CGobj.LymanModelGridExact_EUD.TD50(dd(1)),CGobj.LymanModelGridExact_EUD.m(mm(1))]+10);
            fitobj = cell(size(CGobj.lgn));
            fitobj{kk} = fittedobj;
            goodness = repmat(goodness,size(CGobj.lgn));
            output = repmat(output,size(CGobj.lgn));
            for kk = 2:size(CGobj.lgn,1)
                llh = CGobj.LymanModelGridExact_EUD.loglikelihood(:,:,kk); %log likeli-hood from grid computation
                mn = max(llh(:)); [dd,mm]=find(llh==mn);

                [fittedobj,goodness(kk),output(kk)] = fit(euds(:,kk),flg,fttp,fitoptions(ftopt,'StartPoint',[CGobj.LymanModelGridExact_EUD.TD50(dd(1)),CGobj.LymanModelGridExact_EUD.m(mm(1))]));
                fittedobj
                fitobj{kk} = fittedobj;
            end
%             warning('on','MATLAB:singularMatrix');
            CGobj.LymanModelExact_EUD = struct('FitObj',fitobj,'Goodness',goodness,'Output',output);
        end
    end
end