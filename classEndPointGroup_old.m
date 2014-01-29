classdef classEndPointGroup
    properties
        % general
        mGrp = classEndPointIndividual.empty(1,0); % patient group of classEndPointIndividual objects
        mNumInGrp = 0; % number of patients in the group

        mBeta2Alpha = 0; % beta to alpha ratio
        mLymanN

        mStepX % dose step size for computation based on DVH (in Gy)
        mStepY % volume step size for computation based on DVH (in cc)
        mStepZ % time step size in month for dynamic analysis
        mBinsX % the dose bins can be generated from mStepX, or overwritten in case irregular bin steps were used
        mBinsY % similar to mBinsX
        mBinsZ % the time bins can be overwritten in case irregular bins steps were used

        mAtlasTotal % matrix for total patients computed from DVH
        mAtlasComp % matrix for patients with complications computed from DVH
        mBetaCumMat % based on EUD, the probability that the complication rate is larger than a specific number mBetaCumTh
        mBetaCumTh = 0.2;
        mBetaInvMat % based on EUD, The chi-square where the probability is mBetaInvTh
        mBetaInvTh = 0.16;

        mKaplanMeierCompMat % matrix for patient complication curve from DVH
        mKaplanMeierCompSample % matrix for patient complication curve at specific time point, from DVH
        mKaplanMeierCompFromAtlas % matrix of patient complication curve from binned atlas.
        mKaplanMeierCompOverall = classSurvivalAnalysis.empty(1,0); % overall complication curve

        mKaplanMeierSurvivalOverall = classSurvivalAnalysis.empty(1,0); % overall survival curve
        mKaplanMeierRelapseFree = classSurvivalAnalysis.empty(1,0); % relapse free survival curve
        
        mCoxMinSize = 2; % minimum sample size for Cox model
        mLogRankMinSize =2; % minimum sample size for log rank test
        mCoxPar % cox model pairs (cell 1 -- label of the Cox model fitting, cell 2 -- fitting result
        mLogRank % log rank test for (Di,Vj) (cell 1 -- label of the test (Dx, Vx, DVx), cell 2 -- n1, censor1, n2, censor2, p-value for each grid point

        mLogitMat % logistic regression model using exact EUD
        mLogitMatBin % logistic regression model using binned EUD

        mLogitHosmerLemeshow = struct('n',0,'Chi2',[],'df',[],'p_value',[]); % goodness of fit for logistic regression (binary only) using Hosmer Lemshow test, ...
        % structure: n (number of risk groups), Chi2 (test result which is chi-square), df (degree of freedom, usually n-2), p-value (p-value of the chi-square)
        mLogitHosmerLemeshowBin = struct('n',0,'Chi2',[],'df',[],'p_value',[]);
        mLogitGTest = struct('n',0,'Chi2',[],'df',[],'p_value',[]); % goodness of fit for logistic regression (binary only) using G-test, ...
        % structure: n (number of risk groups), Chi2 (test result which is chi-square), df (degree of freedom, usually n-2), p-value (p-value of the chi-square)
        mLogitGTestBin = struct('n',0,'Chi2',[],'df',[],'p_value',[]); % goodness of fit for logistic regression (binary only), ...
        mLogitPearson = struct('n',0,'Chi2',[],'df',[],'p_value',[]); % goodness of fit for logistic regression (binary only) using Pearson test, ...
        % structure: n (number of risk groups), Chi2 (test result which is chi-square), df (degree of freedom, usually n-2), p-value (p-value of the chi-square)
        mLogitPearsonBin = struct('n',0,'Chi2',[],'df',[],'p_value',[]); % goodness of fit for logistic regression (binary only), ...
        mLogitGridBetaRange % beta ranges for grid searching
        mLogitGrid % logistic regression model using exact EUD for grid
        mLogitGridBin % logistic regression model using binned EUD for grid
        mLogitGOFSim = struct('SSRSim',[],'SSRObserve',[],'p_value',[]);  % goodness of fit using simulations
        mLogitGOFSimBin = struct('SSRSim',[],'SSRObserve',[],'p_value',[]); % goodness of fit using simulations

        mLyman % Lyman model parameters
        mLymanBin
        mLymanGridTD50Range % range of TD50 for grid searching
        mLymanGridMRange % range of m for grid searching
        mLymanGrid % Lyman results (likelihood/probability) for each grid point
        mLymanGridBin
        mLymanHosmerLemeshow = struct('n',0,'Chi2',[],'df',[],'p_value',[]); % goodness of fit for logistic regression (binary only) using Hosmer Lemeshow test, ...
        % structure: n (number of risk groups), Chi2 (test result which is chi-square), df (degree of freedom, usually n-2), p-value (p-value of the chi-square)
        mLymanHosmerLemeshowBin = struct('n',0,'Chi2',[],'df',[],'p_value',[]);
        mLymanGTest = struct('n',0,'Chi2',[],'df',[],'p_value',[]); % goodness of fit for logistic regression (binary only) using G-test, ...
        % structure: n (number of risk groups), Chi2 (test result which is chi-square), df (degree of freedom, usually n-2), p-value (p-value of the chi-square)
        mLymanGTestBin = struct('n',0,'Chi2',[],'df',[],'p_value',[]); % goodness of fit for logistic regression (binary only), ...
        mLymanPearson = struct('n',0,'Chi2',[],'df',[],'p_value',[]); % goodness of fit for logistic regression (binary only) using Pearson test, ...
        % structure: n (number of risk groups), Chi2 (test result which is chi-square), df (degree of freedom, usually n-2), p-value (p-value of the chi-square)
        mLymanPearsonBin = struct('n',0,'Chi2',[],'df',[],'p_value',[]); % goodness of fit for logistic regression (binary only), ...
        mLymanGOFSim = struct('SSRSim',[],'SSRObserve',[],'p_value',[]);  % goodness of fit using simulations
        mLymanGOFSimBin = struct('SSRSim',[],'SSRObserve',[],'p_value',[]); % goodness of fit using simulations
    end

    methods % set member vaulues
        function EPGobj = classEndPointGroup()
%             EPGobj.mGrp = classEndPointIndividual.empty(1,0);
        end
        function EPGobj = set.mBeta2Alpha(EPGobj,mBeta2Alpha)
            if any (mBeta2Alpha < 0)
                error('Beta to Alpha ratio is negative, in set.mBeta2Alpha of classComplicationGroup');
            end
            if length(mBeta2Alpha)>1
                warning('classComplicationGroup:mBeta2Alpha','Parameter for mBeta2Alpha not a scalar, choose the first element');
            end
            for m = 1:EPGobj.mNumInGrp
                EPGobj.mGrp(m).mBeta2Alpha = mBeta2Alpha;
            end
            EPGobj.mBeta2Alpha = mBeta2Alpha;
        end
        function EPGobj = set.mBinsZ(EPGobj,timebins)
            if any( diff( timebins ) <=0 )
                disp('time bins is not strictly monotonic increasing when assign to mBinsZ in classComplicationGroup');
                disp('thus they are sorted');
                timebins = sort(timebins);
                timebins(diff(timebins)==0)=[];
            end
            EPGobj.mBinsZ = timebins(:);
        end
        function EPGobj = set.mBinsX(EPGobj,dosebins)
            if any(diff(dosebins)<=0)
                disp('dose bins is not monotonic increasing when assigned to mBinsX in classComplicationGroup');
                disp('thus they are sorted');
                dosebins = sort( dosebins );
                dosebins(diff(dosebins)==0)=[];
            end
            EPGobj.mBinsX = dosebins(:);
        end
        function EPGobj = set.mBinsY(EPGobj,volbins)
            if any(diff(volbins)<=0)
                disp('volume bins is not monotonic increasing when assigned to mBinsY in classComplicationGroup');
                disp('thus they are sorted');
                volbins = sort( volbins );
                volbins(diff(bolbins)==0)=[];
            end
            EPGobj.mBinsY = volbins(:);
        end
        function EPGobj = set.mLgN(EPGobj,mLgN)
            EPGobj.mLgN = mLgN(:);
        end
        function EPGobj = set.mEUD_DoseBins(EPGobj,dosebins)
            if isempty(dosebins)
                EPGobj.mEUD_DoseBins = dosebins;
                return;
            end
            
            if any(diff(dosebins)<=0)
                disp('dose bins is not monotonic increasing when assigned to mEUD_DoseBins in classComplicationGroup');
                disp('thus they are sorted');
                dosebins = sort( dosebins );
                dosebins(diff(dosebins)==0)=[];
            end
            euds = [EPGobj.mGrp.EUD]; dmax = max(euds(:));
            if dmax > dosebins(end)
                disp('assigned dose bins can not cover the maximum EUD');
            end
            EPGobj.mEUD_DoseBins = dosebins(:);
        end
        function EPGobj = set.mBetaCumTh(EPGobj, betacumthreshold)
            EPGobj.mBetaCumTh = betacumthreshold(:);
        end
        function EPGobj = set.mBetaInvTh(EPGobj,betainvthreshold)
            EPGobj.mBetaInvTh = betainvthreshold(:);
        end
    end

    methods % patient operations like adding and removing
        function EPGobj = fAddPatient(EPGobj,ptobjs)
            % check if the new patient info matches the class definition
            if ~all(isa(ptobjs,'classEndPointIndividual'))
                disp('Not an instance of patient individual class when adding a new patient');
                return;
            end
            
            % make sure the log10(n) and the mBeta2Alpha are the same for all patients
            for m = 1:length(ptobjs)
                ptobjs(m).mLymanN = EPGobj.mLymanN;
                ptobjs(m).mBeta2Alpha = EPGobj.mBeta2Alpha;
            end

            % for empty patient object, add the patient directly
            if isempty([EPGobj.mGrp.mID])
                EPGobj.mGrp = ptobjs(:);
                EPGobj.mNumInGrp = size(EPGobj.mGrp,1);
                return;
            end
            
            % update existing patients in the group
            [fg,g] = ismember( {ptobjs.mID}, {EPGobj.mGrp.mID} );
            f=find(g);
            for m = 1:length(f)
                EPGobj.mGrp(g(f(m))) = ptobjs(f(m));
            end
            
            % add new patient
            EPGobj.mNumInGrp = EPGobj.mNumInGrp + sum(~fg);
            EPGobj.mGrp(end+1:EPGobj.mNumInGrp) = ptobjs(~fg);
        end
        function EPGobj = fRemovePatient(EPGobj,idx)
            if ~exist('idx','var') % if idx is not passed on, remove all patients (default)
                idx = true(EPGobj.mNumInGrp,1);
            end
            EPGobj.mGrp(idx) = [];
            EPGobj.mNumInGrp = length(EPGobj.mGrp);
        end
        function EPGobj = fLinearQuartraticCorrection(EPGobj)
            for k = 1:EPGobj.mNumInGrp
                EPGobj.mGrp(k) = EPGobj.mGrp(k).LinearQuartraticCorrection();
            end
        end
        function flg = fPatientsWithComplicationData(EPGobj)
            pt = EPGobj.mGrp;
            f1 = cellfun(@(x) isempty(x)|isinf(x)|isnan(x),{pt.BaselineDate}); % patients with no baseline date
            f2 = cellfun(@(x) isempty(x)|isinf(x)|isnan(x),{pt.CompOccurDate}); % patients with no complication date
            flg = find(~f2);  % patients with complciation date were not censored, so this information is restated
            for k = 1:length(flg)
                pt(flg(k)).mFlgCensor = 0;
            end
            
            f3 = cellfun(@(x) isempty(x)|isinf(x)|isnan(x),{pt.LastFollowupDate}); % patients with no last follow up date
            f4 = cellfun(@(x) isempty(x)|isinf(x)|isnan(x),{pt.mFlgCensor}); % patients with no censor info
            flg = f1 | (f2&f3) | f4; % patients missing at least one data
            flg = ~flg; % patients with all data
        end
    end

    methods % DVH operations
        function EPGobj = fCalculateDoseBins_DVH(EPGobj)
            doses = cellfun(@(x) x(end), {EPGobj.mGrp.DoseBins_LQ}); dmax = max(doses);
            EPGobj.mBinsX = ( 0 : EPGobj.mStepX : dmax+EPGobj.mStepX )';
        end
        function EPGobj = fCalculateDoseBinsLog_DVH(EPGobj)
            doses = cellfun(@(x) x(end), {EPGobj.mGrp.DoseBins_LQ}); dmax = log10(max(doses));
            doses = cellfun(@(x) x(1), {EPGobj.mGrp.DoseBins_LQ}); dmin = min(doses); dmin = log10(max(1,dmin));
            dosebinslog = dmin : EPGobj.mStepX : (dmax+EPGobj.mStepX);
            EPGobj.mBinsX = [0, 10.^dosebinslog]';
        end
        function EPGobj = fCalculateVolBins_DVH(EPGobj)
            vols = cellfun(@(x) x(1), {EPGobj.mGrp.VolCum}); vmax = max(vols);
            EPGobj.mBinsY = ( 0 : EPGobj.mStepY : vmax+EPGobj.mStepY )';
        end
        function EPGobj = fCalculateVolBinsLog_DVH(EPGobj)
            vols = cellfun(@(x) x(1), {EPGobj.mGrp.VolCum}); vmax = log10(max(vols));
            vols = cellfun(@(x) min(x(x>0)), {EPGobj.mGrp.VolCum}); vmin = min(vols); vmin = log10(vmin);
            vmin = fix(vmin/EPGobj.mStepY)*EPGobj.mStepY;
            volbinslog = vmin : EPGobj.mStepY : (vmax + EPGobj.mStepY);
            EPGobj.mBinsY = [0, 10.^volbinslog]';
        end
        function EPGobj = fCalculateTimeBins_DVH(EPGobj)
            timebins = ( [EPGobj.mGrp.LastFollowupDate] - [EPGobj.mGrp.BaselineDate] ) / 30; tmax = max(timebins);
            EPGobj.mBinsZ = ( 0 : EPGobj.mStepZ : tmax+EPGobj.mStepZ )';
        end
        function EPGobj = fAtlasAlongSampleTime_DVH(EPGobj)
            % prepare
            EPGobj.mAtlasTotal = zeros( length(EPGobj.mBinsX), length(EPGobj.mBinsY), length(EPGobj.mBinsZ) );
            EPGobj.mAtlasComp = zeros( length(EPGobj.mBinsX), length(EPGobj.mBinsY), length(EPGobj.mBinsZ) );
            
            % select patients with data
            f = EPGobj.fPatientsWithComplicationData();
            pt = EPGobj.mGrp(f);
            num = sum(f); % number of patients

            % complication time
            f2 = ~cellfun(@(x) isempty(x)|isinf(x)|isnan(x),{pt.CompOccurDate}); % patients with complication date
            f3 = ~cellfun(@(x) isempty(x)|isinf(x)|isnan(x),{pt.LastFollowupDate}); % patients with last follow up date
            complicationdays = inf(num,1);
            lastfollowup = inf(num,1);
            complicationdays(f2) = ([pt(f2).CompOccurDate] - [pt(f2).BaselineDate])' / 30; % convert to months
            lastfollowup(f3) = ([pt(f3).LastFollowupDate] - [pt(f3).BaselineDate])' / 30;
            lastfollowup = min( lastfollowup, complicationdays );

            % Vx
            vols = zeros( num, 1 );
            for ii = 1:length(EPGobj.mBinsX)
                % volume of each patient at current dose (x)
                vols(:) = 0;
                for jj = 1:num
                    vols(jj) = pt(jj).VolAtDose( EPGobj.mBinsX(ii) );
                end
                vols(vols==0) = -1; % exclude zero volume patients
                
                % matrix at each (Di, Vj, Tk)
                for jj = 1:length(EPGobj.mBinsY) % for each volume point under dose x
                    f = find( vols >= EPGobj.mBinsY(jj) ); % patient at the grid point
                    % for each time point
                    for kk = 1:length(EPGobj.mBinsZ)
                        % patients with last followup or complications before the time point
                        f1 = find( lastfollowup(f) <= EPGobj.mBinsZ(kk) );
                        % total patients
                        EPGobj.mAtlasTotal(ii,jj,kk) = length(f1);
                        % patients with complications
                        f1 = find( complicationdays(f) <= EPGobj.mBinsZ(kk) );
                        EPGobj.mAtlasComp(ii,jj,kk) = length(f1);
                    end
                end
            end
        end
        function EPGobj = fComplicationCurves_DVH(EPGobj)
            % select patients with data
            f = EPGobj.fPatientsWithComplicationData();
            pt = EPGobj.mGrp(f);
            num = sum(f); % number of patients

            % complication time
            f2 = ~cellfun(@(x) isempty(x)|isinf(x)|isnan(x),{pt.CompOccurDate}); % patients with complication date
            f3 = ~cellfun(@(x) isempty(x)|isinf(x)|isnan(x),{pt.LastFollowupDate}); % patients with last follow up date
            complicationdays = inf(num,1);
            lastfollowup = inf(num,1);
            complicationdays(f2) = ([pt(f2).CompOccurDate] - [pt(f2).BaselineDate])' / 30;
            lastfollowup(f3) = ([pt(f3).LastFollowupDate] - [pt(f3).BaselineDate])' / 30;
            lastfollowup = min( lastfollowup, complicationdays );
            flgcensor = [pt.mFlgCensor]';

            % overall survival
            sa = classSurvivalAnalysis();
            sa.mSurvivalTime = {lastfollowup};
            sa.mFlgCensor = {flgcensor};
            sa = sa.fCalculateSurvivalCurve();

            EPGobj.mKaplanMeierCompOverall = sa;
        end
        function EPGobj = fSurvivalCurves_DVH(EPGobj)
            % select patients with data
            pt = EPGobj.mGrp;
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
            sa.mSurvivalTime = {survivaltime};
            sa.mFlgCensor = {flgcensor};
            sa = sa.fCalculateSurvivalCurve();

            EPGobj.mKaplanMeierSurvivalOverall = sa;
            
            % free of relapse survival
            f2 = ~cellfun(@(x) isempty(x)|isinf(x)|isnan(x),{pt.RelapseDate}); % patients with relapse date
            flgcensor(:) = 1;
            flgcensor(f2) = 0;
            survivaltime(f2) = ([pt(f2).RelapseDate] - [pt(f2).BaselineDate])'/30;
            sa.mSurvivalTime = {survivaltime};
            sa.mFlgCensor = {flgcensor};
            sa = sa.fCalculateSurvivalCurve();

            EPGobj.mKaplanMeierRelapseFree = sa;
        end
        function EPGobj = fComplicationActuary_DVH(EPGobj)
            % select patients with data
            f = EPGobj.fPatientsWithComplicationData();
            CG = EPGobj.fRemovePatient(~f);

            % complication time
            f2 = ~cellfun('isempty',{CG.mGrp.CompOccurDate}); % patients with no complication date
            f3 = ~cellfun('isempty',{CG.mGrp.LastFollowupDate}); % patients with no last follow up date
            complicationdays = inf(CG.mNumInGrp,1);
            lastfollowup = inf(CG.mNumInGrp,1);
            complicationdays(f2) = ([CG.mGrp(f2).CompOccurDate] - [CG.mGrp(f2).BaselineDate])' / 30;
            lastfollowup(f3) = ([CG.mGrp(f3).LastFollowupDate] - [CG.mGrp(f3).BaselineDate])' / 30;
            lastfollowup = min( lastfollowup, complicationdays );
            flgcensor = [CG.mGrp.mFlgCensor]';

            % populate to all grid points
            sa = classSurvivalAnalysis(); % a new clean object
            pse = repmat(sa, [length(CG.mBinsX), length(CG.mBinsY)]); % patient survival exact

            % Vx
            vols = zeros( CG.mNumInGrp, 1 );
            for ii = 1:length(CG.mBinsX)
                % volume of each patient at current dose (x)
                vols(:) = 0;
                for jj = 1:CG.mNumInGrp
                    vols(jj) = CG.mGrp(jj).VolAtDose( CG.mBinsX(ii) );
                end
                vols(vols==0) = -1; % exclude zero volume patients
                
                % matrix at each (Di, Vj)
                for jj = 1:length(CG.mBinsY) % for each volume point under dose x
                    f = vols >= CG.mBinsY(jj); % patient at the grid point
                    if ~any(f) % no patient for the grid, clean its contents
                        continue;
                    end

                    % use clasSurvivalAnalysis to compute the Kelplan Meier curve
                    pse(ii,jj).mSurvivalTime = {lastfollowup(f)};
                    pse(ii,jj).mFlgCensor = {flgcensor(f)};
                    pse(ii,jj) = pse(ii,jj).fCalculateSurvivalCurve();
                end
            end
            EPGobj.mKaplanMeierCompMat = pse;
        end
        function EPGobj = fComplicationAtSampleTimeActuary_DVH(EPGobj)
            % prepare
            psem = zeros( length(EPGobj.mBinsX), length(EPGobj.mBinsY), length(EPGobj.mBinsZ) );
            pse = EPGobj.mKaplanMeierCompMat(:); pse = reshape(pse,size(EPGobj.mKaplanMeierCompMat));
            tb = EPGobj.mBinsZ(:);
            db = EPGobj.mBinsX(:);
            vb = EPGobj.mBinsY(:);

            % Vx
            for ii = 1:length(db)
                % matrix at each (Di, Vj, Tk)
                for jj = 1:length(vb) % for each volume point under dose x
                    sa = pse(ii,jj);
                    if isempty(sa.mSurvivalTime) % no patient for the grid, skip it
                        continue;
                    end
                    
                    % sample the survival curve for each time point
                    f = find( tb <= sa.mSurvivalTimeSorted{1}(end) );
                    for kk = 1:length(f)
                        f1 = find( tb(kk) >= sa.mSurvivalTimeSorted{1} );
                        psem(ii,jj,kk) = sa.mSurvivalCurve{1}(f1(end));
                    end
                end
            end
            EPGobj.mKaplanMeierCompSample = psem;
        end
        function EPGobj = fComplicationAtSampleTimeBins_DVH(EPGobj)
            % prepare
            EPGobj.mKaplanMeierCompFromAtlas = ones( length(EPGobj.mBinsX), length(EPGobj.mBinsY), length(EPGobj.mBinsZ) );
            
            % survival along the time bins
            for kk = 2:length(EPGobj.mBinsZ)
                % patients at risk
                ptrisk = EPGobj.mAtlasTotal(:,:,end) - EPGobj.mAtlasTotal(:,:,kk-1);
                % patient with complication
                ptcomp = EPGobj.mAtlasComp(:,:,end) - EPGobj.mAtlasComp(:,:,kk-1);
                % rate of patients without complication
                EPGobj.mKaplanMeierCompFromAtlas(:,:,kk) = EPGobj.mKaplanMeierCompFromAtlas(:,:,kk-1).*(1-ptcomp./ptrisk);
            end
        end
        function EPGobj = fCoxModelExact_DVH(EPGobj)
            % select patients with data
            f = EPGobj.fPatientsWithComplicationData();
            CG = EPGobj.fRemovePatient(~f);
            
            % survival/complication time
            f2 = ~cellfun('isempty',{CG.mGrp.CompOccurDate}); % patients with no complication date
            f3 = ~cellfun('isempty',{CG.mGrp.LastFollowupDate}); % patients with no last follow up date
            compdate = inf(CG.mNumInGrp,1);
            lastfollowup = inf(CG.mNumInGrp,1);
            compdate(f2) = ([CG.mGrp(f2).CompOccurDate] - [CG.mGrp(f2).BaselineDate])' / 30;
            lastfollowup(f3) = ([CG.mGrp(f3).LastFollowupDate] - [CG.mGrp(f3).BaselineDate])' / 30;
            compdate = min( lastfollowup, compdate );
            flgcensor = [CG.mGrp.mFlgCensor]';

            warning('off');
            % dmax
                dv = cellfun(@(x) x(end), {CG.mGrp.DoseBins_LQ}');

                [~,logl,h,stats]=coxphfit(dv,compdate,'baseline',0,'censoring',flgcensor);
                stats.logl=logl; stats.h=h;
                stats.data_exposure = dv; stats.data_hazard = compdate;

            % fx
                fxnum=[CG.mGrp.FxNum]'; % all fractions
                flgfx = length(unique(fxnum))>1; % if there is more than one fraction numbers, do the cox model
                if flgfx
                    [~,logl,h,stats]=coxphfit([CG.mGrp.FxNum]',compdate,'baseline',0,'censoring',flgcensor);
                    stats.logl=logl; stats.h=h;
                    stats.data_exposure = [CG.mGrp.FxNum]'; stats.data_hazard = compdate;
                    
                    if isempty(EPGobj.mCoxPar)
                            EPGobj.mCoxPar{end+1,1}='Fx';
                            EPGobj.mCoxPar{end,2}=stats;
                    else
                        f = cellfun(@(x) strcmpi('Fx',x),EPGobj.mCoxPar(:,1));
                        if any(f)
                            EPGobj.mCoxPar{f,2}=stats;
                        else
                            EPGobj.mCoxPar{end+1,1}='Fx';
                            EPGobj.mCoxPar{end,2}=stats;
                        end
                    end
                end

            % prescription dose
                dv = cell2mat([CG.mGrp.DosePrescription]');
                [~,logl,h,stats]=coxphfit(dv,compdate,'baseline',0,'censoring',flgcensor);
                stats.logl=logl; stats.h=h;
                stats.data_exposure = dv; stats.data_hazard = compdate;

                    if isempty(EPGobj.mCoxPar)
                            EPGobj.mCoxPar{end+1,1}='Dp'; % dose prescription
                            EPGobj.mCoxPar{end,2}=stats;
                    else
                        f = cellfun(@(x) strcmpi('Dp',x),EPGobj.mCoxPar(:,1));
                        if any(f)
                            EPGobj.mCoxPar{f,2}=stats;
                        else
                            EPGobj.mCoxPar{end+1,1}='Dp';
                            EPGobj.mCoxPar{end,2}=stats;
                        end
                    end

                % fx + prescription dose
                if flgfx
                    [~,logl,h,stats]=coxphfit([dv,fxnum],compdate,'baseline',0,'censoring',flgcensor);
                    stats.logl=logl; stats.h=h;
                    stats.data_exposure = [dv,fxnum]; stats.data_hazard = compdate;
                    
                        f = cellfun(@(x) strcmpi('DpFx',x),EPGobj.mCoxPar(:,1));
                        if any(f)
                            EPGobj.mCoxPar{f,2}=stats;
                        else
                            EPGobj.mCoxPar{end+1,1}='DpFx';
                            EPGobj.mCoxPar{end,2}=stats;
                        end
                end
                
            % Dx
                volnum=length(CG.mBinsY);
                CoxDx=repmat(stats,[volnum,1]);
                CoxDVx=repmat(stats,[volnum,1]);
                CoxDxFx=repmat(stats,[volnum,1]);
                CoxDVxFx=repmat(stats,[volnum,1]);
                
                % check for Dx one by one
                dv = zeros(CG.mNumInGrp,1);
                for v=1:volnum
                    % doses under v
                    dv(:) = 0;
                    for k = 1:CG.mNumInGrp
                        dv(k) = CG.mGrp(k).DoseAtVol( CG.mBinsY(v) );
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
                    if length(g)<CG.mCoxMinSize
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
                if isempty(EPGobj.mCoxPar)
                    EPGobj.mCoxPar{end+1,1}='Dx';
                    EPGobj.mCoxPar{end,2}=CoxDx;
                else
                    f = cellfun(@(x) strcmpi('Dx',x),EPGobj.mCoxPar(:,1));
                    if any(f)
                        EPGobj.mCoxPar{f,2}=CoxDx;
                    else
                        EPGobj.mCoxPar{end+1,1}='Dx';
                        EPGobj.mCoxPar{end,2}=CoxDx;
                    end
                end
                f = cellfun(@(x) strcmpi('DVx',x),EPGobj.mCoxPar(:,1));
                if any(f)
                    EPGobj.mCoxPar{f,2}=CoxDVx;
                else
                    EPGobj.mCoxPar{end+1,1}='DVx';
                    EPGobj.mCoxPar{end,2}=CoxDVx;
                end
                if flgfx
                    f = cellfun(@(x) strcmpi('DVxFx',x),EPGobj.mCoxPar(:,1));
                    if any(f)
                        EPGobj.mCoxPar{f,2}=CoxDVxFx;
                    else
                        EPGobj.mCoxPar{end+1,1}='DVxFx';
                        EPGobj.mCoxPar{end,2}=CoxDVxFx;
                    end
                    f = cellfun(@(x) strcmpi('DxFx',x),EPGobj.mCoxPar(:,1));
                    if any(f)
                        EPGobj.mCoxPar{f,2}=CoxDxFx;
                    else
                        EPGobj.mCoxPar{end+1,1}='DxFx';
                        EPGobj.mCoxPar{end,2}=CoxDxFx;
                    end
                end

            % Vx
                dosenum=length(CG.mBinsX);
                CoxVx=repmat(stats,[dosenum,1]);
                CoxVDx=repmat(stats,[dosenum,1]);
                CoxVxFx=repmat(stats,[dosenum,1]);
                CoxVDxFx=repmat(stats,[dosenum,1]);
                vd=zeros(CG.mNumInGrp,1);
                % check for Vx one by one
                for d=1:dosenum
                    % volumes under d
                    vd(:)=0;
                    for k=1:CG.mNumInGrp
                        vd(k) = CG.mGrp(k).VolAtDose( CG.mBinsX(d) );
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
                    if length(g)<CG.mCoxMinSize
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
                f = cellfun(@(x) strcmpi('VDx',x),EPGobj.mCoxPar(:,1));
                if any(f)
                    EPGobj.mCoxPar{f,2}=CoxVDx;
                else
                    EPGobj.mCoxPar{end+1,1}='VDx';
                    EPGobj.mCoxPar{end,2}=CoxVDx;
                end
                f = cellfun(@(x) strcmpi('Vx',x),EPGobj.mCoxPar(:,1));
                if any(f)
                    EPGobj.mCoxPar{f,2}=CoxVx;
                else
                    EPGobj.mCoxPar{end+1,1}='Vx';
                    EPGobj.mCoxPar{end,2}=CoxVx;
                end
                if flgfx
                    f = cellfun(@(x) strcmpi('VDxFx',x),EPGobj.mCoxPar(:,1));
                    if any(f)
                        EPGobj.mCoxPar{f,2}=CoxVDxFx;
                    else
                        EPGobj.mCoxPar{end+1,1}='VDxFx';
                        EPGobj.mCoxPar{end,2}=CoxVDxFx;
                    end
                    f = cellfun(@(x) strcmpi('VxFx',x),EPGobj.mCoxPar(:,1));
                    if any(f)
                        EPGobj.mCoxPar{f,2}=CoxVxFx;
                    else
                        EPGobj.mCoxPar{end+1,1}='VxFx';
                        EPGobj.mCoxPar{end,2}=CoxVxFx;
                    end
                end
                warning('on');
        end
        function EPGobj = fCoxModelAtlas_DVH(EPGobj)
            if isempty(EPGobj.mAtlasTotal)
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
            flgcensor = false(EPGobj.mAtlasTotal(1,1,end),1);
            vx = zeros(size(flgcensor));
            dx = zeros(size(flgcensor));
            compdate = zeros(size(flgcensor));
            
            [dimd,dimv,dimt] = size( EPGobj.mAtlasTotal ); % dimensions of dose, volume, and time
            
            dosebins = EPGobj.mBinsX; dosebins(1:end-1) = (dosebins(1:end-1) + dosebins(2:end)) / 2; % dose bins
            volbins = EPGobj.mBinsY; volbins(1:end-1) = (volbins(1:end-1) + volbins(2:end)) / 2; % volume bins
            timebins = EPGobj.mBinsZ; timebins(2:end) = (timebins(1:end-1) + timebins(2:end)) / 2; % time bins
            
            % Cox model
            [~,logl,h,stats]=coxphfit([0; 1],[0; 1],'baseline',0,'censoring',[0; 1]);
            stats.logl=logl; stats.h=h;
            stats.data_exposure = [0; 1]; stats.data_hazard = [0; 1];
            % Vx
            CoxVxAtlas = repmat(stats,[dimd,1]);
            for d = 1:dimd % Cox model for each dose (x)
                % extract the matrix corresponding to the dose
                vc = squeeze( EPGobj.mAtlasComp(d,:,:) ); % complication at different volume and time points
                vt = squeeze( EPGobj.mAtlasTotal(d,:,:) ); % complication and censored at different time points
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
                dc = squeeze( EPGobj.mAtlasComp(:,v,:) ); % complication at different volume and time points
                dt = squeeze( EPGobj.mAtlasTotal(:,v,:) ); % complication and censored at different time points
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
            f = cellfun(@(x) strcmpi('VxAtlas',x),EPGobj.mCoxPar(:,1));
            if any(f)
                EPGobj.mCoxPar{f,2}=CoxVxAtlas;
            else
                EPGobj.mCoxPar{end+1,1}='VxAtlas';
                EPGobj.mCoxPar{end,2}=CoxVxAtlas;
            end
            f = cellfun(@(x) strcmpi('DxAtlas',x),EPGobj.mCoxPar(:,1));
            if any(f)
                EPGobj.mCoxPar{f,2}=CoxDxAtlas;
            else
                EPGobj.mCoxPar{end+1,1}='DxAtlas';
                EPGobj.mCoxPar{end,2}=CoxDxAtlas;
            end
            
            warning('on');
%             warning('on','MATLAB:singularMatrix');
%             warning('on','stats:coxphfit:FitWarning');
%             warning('on','stats:coxphfit:RankDeficient');
%             warning('on','stats:coxphfit:IterOrEvalLimit');
        end
        function EPGobj = fLogRankDxExact_DVH(EPGobj)
            % select patients with data
            f = EPGobj.fPatientsWithComplicationData();
            CG = EPGobj.fRemovePatient(~f);
            
            % survival/complication time
            f2 = ~cellfun('isempty',{CG.mGrp.CompOccurDate}); % patients with no complication date
            f3 = ~cellfun('isempty',{CG.mGrp.LastFollowupDate}); % patients with no last follow up date
            compdate = inf(CG.mNumInGrp,1);
            lastfollowup = inf(CG.mNumInGrp,1);
            compdate(f2) = ([CG.mGrp(f2).CompOccurDate] - [CG.mGrp(f2).BaselineDate])' / 30;
            lastfollowup(f3) = ([CG.mGrp(f3).LastFollowupDate] - [CG.mGrp(f3).BaselineDate])' / 30;
            compdate = min( lastfollowup, compdate );
            flgcensor = [CG.mGrp.mFlgCensor]';

            % prepare
                volnum=length(CG.mBinsY);
                dosenum=length(CG.mBinsX);
                Dxmat = ones(dosenum,volnum,6); % n1,c1,n2,c2,p,flg (0 -- positive corelation, 1 -- negative corelation, 2 -- not available)
                Dxmat(:,:,6) = 2; % default is not available
                
                sa=classSurvivalAnalysis(); % initialize a survivalanalysis obj
                
            % Dx
                dv=zeros(CG.mNumInGrp,1);
                numstart=CG.mLogRankMinSize;
                for v=1:volnum
                    % doses under v
                    dv(:) = 0;
                    for k = 1:CG.mNumInGrp
                        dv(k) = CG.mGrp(k).DoseAtVol( CG.mBinsY(v) );
                    end
                    
                    % Dx at (di,vj)
                    g=find(dv); % patients with non-zero doses
                    flg_dosebelow2=-1;
                    numend=length(g)-numstart;
                    for d=1:dosenum
                        % check the sample size
                        flg_dosebelow1=dv(g)<=CG.mBinsX(d); f=length(find(flg_dosebelow1));
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
                        sa.mSurvivalTime=survivedate;
                        sa.mFlgCensor=fcensor;
                        % compute survival curves and compare them
                        sa=sa.fCalculateSurvivalCurve();
                        sa=sa.fCombineSurvivalTime();
                        sa=sa.fCompareSurvivalByLogrank();
                        Dxmat(d,v,1:5)=[length(survivedate{1,1}),sum(fcensor{1,1}), length(survivedate{2,1}),sum(fcensor{2,1}),sa.mpValue];
                        Dxmat(d,v,6)=sa.mCurveArea(1)<sa.mCurveArea(2); % the group with lower volume had worse survival curve, record it
                    end
                end
            % save reslt
                if isempty(EPGobj.mLogRank)
                        EPGobj.mLogRank{end+1,1}='Dx';
                        EPGobj.mLogRank{end,2}=Dxmat;
                else
                    f = cellfun(@(x) strcmpi('Dx',x),EPGobj.mLogRank(:,1));
                    if any(f)
                        EPGobj.mLogRank{f,2}=Dxmat;
                    else
                        EPGobj.mLogRank{end+1,1}='Dx';
                        EPGobj.mLogRank{end,2}=Dxmat;
                    end
                end
        end
        function EPGobj = fLogRankVxExact_DVH(EPGobj)
            % select patients with data
            f = EPGobj.fPatientsWithComplicationData();
            CG = EPGobj.fRemovePatient(~f);
            
            % survival/complication time
            f2 = ~cellfun('isempty',{CG.mGrp.CompOccurDate}); % patients with no complication date
            f3 = ~cellfun('isempty',{CG.mGrp.LastFollowupDate}); % patients with no last follow up date
            compdate = inf(CG.mNumInGrp,1);
            lastfollowup = inf(CG.mNumInGrp,1);
            compdate(f2) = ([CG.mGrp(f2).CompOccurDate] - [CG.mGrp(f2).BaselineDate])' / 30;
            lastfollowup(f3) = ([CG.mGrp(f3).LastFollowupDate] - [CG.mGrp(f3).BaselineDate])' / 30;
            compdate = min( lastfollowup, compdate );
            flgcensor = [CG.mGrp.mFlgCensor]';

            % prepare
                volnum=length(CG.mBinsY);
                dosenum=length(CG.mBinsX);
                Vxmat = ones(dosenum,volnum,6); % (n1,c1,n2,c2,p,flg (0 -- positive corelation, 1 -- negative corelation, 2 -- not available))
                Vxmat(:,:,6) = 2;
                
                sa=classSurvivalAnalysis(); % initialize a survivalanalysis obj
                
            % Vx
                vd=zeros(CG.mNumInGrp,1); % volume v at dose d
                numstart=CG.mLogRankMinSize;
                for d=1:dosenum
                    % volume under d
                    vd(:)=0;
                    for k=1:CG.mNumInGrp
                        vd(k) = CG.mGrp(k).VolAtDose( CG.mBinsX(d) );
                    end
                    g=find(vd); % non-zeros volume cases
                    
                    % (di,vj)
                    flg_volbelow2=-1;
                    numend=length(g)-numstart;
                    for v=1:volnum
                        % check smaple size
                        flg_volbelow1=vd(g)<=CG.mBinsY(v); f=length(find(flg_volbelow1)); % group DVHs by (d,v)
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
                        sa.mSurvivalTime=survivedate;
                        sa.mFlgCensor=fcensor;
                        % compute survival curves and compare them
                        sa=sa.fCalculateSurvivalCurve();
                        sa=sa.fCombineSurvivalTime();
                        sa=sa.fCompareSurvivalByLogrank();
                        Vxmat(d,v,1:5)=[length(survivedate{1,1}),sum(fcensor{1,1}), length(survivedate{2,1}),sum(fcensor{2,1}),sa.mpValue];
                        Vxmat(d,v,6)=sa.mCurveArea(1)<sa.mCurveArea(2); % the group with lower volume had worse survival curve, record it
                    end
                end
                
            % save reslt
                if isempty(EPGobj.mLogRank)
                        EPGobj.mLogRank{end+1,1}='Vx';
                        EPGobj.mLogRank{end,2}=Vxmat;
                else
                    f = cellfun(@(x) strcmpi('Vx',x),EPGobj.mLogRank(:,1));
                    if any(f)
                        EPGobj.mLogRank{f,2}=Vxmat;
                    else
                        EPGobj.mLogRank{end+1,1}='Vx';
                        EPGobj.mLogRank{end,2}=Vxmat;
                    end
                end
        end
        function EPGobj = fLogRankDVxExact_DVH(EPGobj)
            % select patients with data
            f = EPGobj.fPatientsWithComplicationData();
            CG = EPGobj.fRemovePatient(~f);
            
            % survival/complication time
            f2 = ~cellfun('isempty',{CG.mGrp.CompOccurDate}); % patients with no complication date
            f3 = ~cellfun('isempty',{CG.mGrp.LastFollowupDate}); % patients with no last follow up date
            compdate = inf(CG.mNumInGrp,1);
            lastfollowup = inf(CG.mNumInGrp,1);
            compdate(f2) = ([CG.mGrp(f2).CompOccurDate] - [CG.mGrp(f2).BaselineDate])' / 30;
            lastfollowup(f3) = ([CG.mGrp(f3).LastFollowupDate] - [CG.mGrp(f3).BaselineDate])' / 30;
            compdate = min( lastfollowup, compdate );
            flgcensor = [CG.mGrp.mFlgCensor]';

            % prepare
                volnum=length(CG.mBinsY);
                dosenum=length(CG.mBinsX);
                DVxmat = ones(dosenum,volnum,6); % n1,c1,n2,c2,p,flg (0 -- positive corelation, 1 -- negative corelation, 2 -- not available)
                DVxmat(:,:,6) = 2; % default is not available
                
                sa=classSurvivalAnalysis(); % initialize a survivalanalysis obj
                
            % Dx
                dv=zeros(CG.mNumInGrp,1);
                numstart=CG.mLogRankMinSize;
                numend=CG.mNumInGrp-numstart;
                for v=1:volnum
                    % doses under v
                    dv(:) = 0;
                    for k = 1:CG.mNumInGrp
                        dv(k) = CG.mGrp(k).DoseAtVol( CG.mBinsY(v) );
                    end
                    
                    % DVx at (dj,vj)
                    flg_dosebelow2=-1;
                    for d=1:dosenum
                        % check the sample size
                        flg_dosebelow1=dv<=CG.mBinsX(d); f=length(find(flg_dosebelow1));
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
                        sa.mSurvivalTime=survivedate;
                        sa.mFlgCensor=fcensor;
                        % compute survival curves and compare them
                        sa=sa.fCalculateSurvivalCurve();
                        sa=sa.fCombineSurvivalTime();
                        sa=sa.fCompareSurvivalByLogrank();
                        DVxmat(d,v,1:5)=[length(survivedate{1,1}),sum(fcensor{1,1}), length(survivedate{2,1}),sum(fcensor{2,1}),sa.mpValue];
                        DVxmat(d,v,6)=sa.mCurveArea(1)<sa.mCurveArea(2); % the group with lower volume had worse survival curve, record it
                    end
                end
            
            % save reslt
                if isempty(EPGobj.mLogRank)
                        EPGobj.mLogRank{end+1,1}='DVx';
                        EPGobj.mLogRank{end,2}=DVxmat;
                else
                    f = cellfun(@(x) strcmpi('DVx',x),EPGobj.mLogRank(:,1));
                    if any(f)
                        EPGobj.mLogRank{f,2}=DVxmat;
                    else
                        EPGobj.mLogRank{end+1,1}='DVx';
                        EPGobj.mLogRank{end,2}=DVxmat;
                    end
                end
        end
        function [allCox,flgCox,flgAnti] = fCoxPar_DVH(EPGobj,strCoxVx)
            f = cellfun(@(x) strcmpi(strCoxVx,x),EPGobj.mCoxPar(:,1)); % search the label
            if any(f)
                allCox = EPGobj.mCoxPar{f,2}; % extract Cox model result
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
        function fCoxRiskVDxFig_DVH(EPGobj,d,t)
            numintv = 4; % group number, i.e., interval numbers
            % Vx
            % search the best Cox model
            [allCox,flgCox,flganti] = CoxPar_DVH(EPGobj,'VDx'); % find availabe Cox models
            flgCox(flganti)=false; % anti-correlations were not be considered
            logl = [allCox.logl]'; logl(~flgCox) = -inf; % log likelihood of Cox model, anti-correlation points not counted
            [~,doseloc]=max(logl); % the best fitting of Cox model
            disp(' ');
            disp(['Best Cox Model at dose :',num2str(EPGobj.mBinsX(doseloc))]);
            disp(allCox(doseloc));
            if exist('d','var')
                if d~=-1
                    doseloc = EPGobj.mBinsX == d;
                end
            end
            allCox = allCox(doseloc);
%             disp(' ');
            disp(['Cox model at dose: ',num2str(EPGobj.mBinsX(doseloc))]);
            disp(allCox);

            flgcensor = [EPGobj.mGrp.mFlgCensor]';% the overall patient complication info
            % specify event time
            comptime = ([EPGobj.mGrp.CompOccurDate]' - [EPGobj.mGrp.BaselineDate]')/30;
%             disp(' ');
            disp(['median complication time: ',num2str(median(comptime(~flgcensor)))]);
            if ~exist('t','var')
                t = median(comptime(~flgcensor));
            end
            disp(['Cox model analysis at time: ',num2str(t)]);

            % volumes of patients at best Vx
            Vx=zeros(EPGobj.mNumInGrp,1);
            d = EPGobj.mBinsX(doseloc);
            for k=1:EPGobj.mNumInGrp
                Vx(k) = EPGobj.mGrp(k).VolAtDose( d );
            end
            % observed volume at best Vx in groups
            flgcensor(comptime>t) = 1; % by the time t some patients might not develop complications so they shall be excluded.
%             for k = 1:EPGobj.mNumInGrp
%                 EPGobj.mGrp(k).mFlgCensor = flgcensor(k);
%             end
            [sortQ,indxQ,indxorg] = EqualIntervals(Vx,numintv);
            meanvol = zeros(numintv,1);
            prob = zeros(numintv,1);
            betainv84 = zeros(numintv,1);
            betainv16 = zeros(numintv,1);
            f = ~flgcensor(indxorg);
            for m = 1 : numintv
                meanvol(m) = median(sortQ(indxQ==m));
                numcomp = sum(f(indxQ==m)); numtotal = sum(indxQ==m);
                prob(m) = numcomp/numtotal;
                betainv84(m) = betainv( .84, numcomp+1, numtotal - numcomp + 1 );
                betainv16(m) = betainv( .16, numcomp+1, numtotal - numcomp + 1 );
            end
            
            % servivial curves for each group
%             str = 'rbkc';
            CoxComplicationTime = cell(numintv,1);
            CoxComplicationCurve = cell(numintv,1);
            CG = repmat(EPGobj,[numintv,1]);
            for m = 1:numintv
                % compute the complication of the group
                CG(m) = EPGobj.fRemovePatient(indxorg(indxQ ~= m));
                
%                 for k = 1:CG(m).mNumInGrp
%                     Vx(k) = CG(m).mGrp(k).VolAtDose(d);
%                 end
%                 hold on; plot(m,Vx(1:CG(m).mNumInGrp),[str(m),'*']); hold off;
                
                CG(m) = CG(m).ComplicationCurves_DVH();
                [CoxComplicationTime{m},CoxComplicationCurve{m}] = CoxRiskVDxComplicationFig_DVH(CG(m),allCox,meanvol(m));

                f = find(CG(m).mKaplanMeierCompOverall.mSurvivalTimeSorted{1} <= t);
                prob(m) = CG(m).mKaplanMeierCompOverall.mSurvivalCurve{1}(f(end));
            end
            prob = 1-prob;
            
            % observed risk ratio
            disp(' ');
            disp(['observed risk and the ratio of fourth and first quartile: ',num2str([prob', prob(end)/prob(1)])]);

            % computed risk ratio at time t
            % h(t) clean up
            f = find(diff(allCox.h(:,1))==0); % find duplicate time values of h(t)
            while ~isempty(f)
                allCox.h(f,1) = allCox.h(f,1)-eps*10; % adjust it a bit to avoid ambiguius
                f = find(diff(allCox.h(:,1))==0); % find duplicate time values of h(t)
            end
            h = interp1(allCox.h(:,1),allCox.h(:,2),t,'linear','extrap');
            % complication rate from Cox model at time t w.r.t volume
            vol = (0:max(Vx)*2)';
            os = exp(-h*exp(allCox.beta.*vol));
            osup68 = exp(-h*exp((allCox.beta+allCox.se).*vol));
            oslow68 = exp(-h*exp((allCox.beta-allCox.se).*vol));
            s = exp(-h*exp(allCox.beta.*meanvol));
            % highest risk from Cox model
            disp(['computed risk at quartile and its risk ratio: ',num2str([1-s',(1-s(end))/(1-s(1))])]);
            mxprobchangeloc = -log(h)/allCox.beta;
            mxprobchange = exp(-1)*allCox.beta;
            disp(['computed maximum probability change rate is at: ', num2str(mxprobchangeloc)]);
            disp(['computed maximum probability change rate is: ', num2str(mxprobchange)]);
            % plot curves
            plot(vol,1-[os,osup68,oslow68]);
%             hold on; errorbar(meanvol,prob,max(0,prob-betainv16),max(0,betainv84-prob),'r*','linewidth',1,'markersize',12); hold off;
            hold on; plot(meanvol,prob,'r*','linewidth',1,'markersize',12); hold off;



            % plot Cox survival curves at different volume resolutions
            % survival curve from Cox model at volume vol
            CoxCompTime_HighResolution = 0:0.1:44;
            h = interp1(allCox.h(:,1),allCox.h(:,2),CoxCompTime_HighResolution,'linear','extrap');
            expbetax = exp(allCox.beta*meanvol);
            CoxCompCurve_HighResolution = zeros(length(CoxCompTime_HighResolution),numintv);
            for m = 1:numintv
                CoxCompCurve_HighResolution(:,m) = exp( -h * expbetax(m) );
            end
            % plot curves
            str = 'rbkc'; clf reset;hold on;
            for m = 1:numintv
                plot(CoxCompTime_HighResolution,CoxCompCurve_HighResolution(:,m),[str(mod(m,length(str))+1),'-'],'LineWidth',1); % high resolution Cox curve
                CoxTime = [CoxComplicationTime{m}; CoxCompTime_HighResolution(end)];
                CoxCurve = [CoxComplicationCurve{m}; CoxComplicationCurve{m}(end)];
                stairs(CoxTime,CoxCurve,[str(mod(m,length(str))+1),'--']); % Cox curve at event time
            end
            hold off;



            % plot K-M survival curve and cox survival curve
            str = 'rbkc'; clf reset;hold on;
            for m = 1:numintv
%                 CoxTime = [CoxComplicationTime{m}; CG(m).mKaplanMeierCompOverall.mSurvivalTimeSorted{1}(end)];
%                 CoxCurve = [CoxComplicationCurve{m}; CoxComplicationCurve{m}(end)];
%                 stairs(CoxTime,CoxCurve,[str(mod(m,length(str))+1),'-']);
                plot(CoxCompTime_HighResolution,CoxCompCurve_HighResolution(:,m),[str(mod(m,length(str))+1),'--'],'LineWidth',1); % high resolution Cox curve
                stairs(CG(m).mKaplanMeierCompOverall.mSurvivalTimeSorted{1}, CG(m).mKaplanMeierCompOverall.mSurvivalCurve{1}, [str(mod(m,length(str))+1),'-']);
            end
            hold off;



            set(gca,'xminortick','on','yminortick','on');
            set(gca,'box','on');
            xlabel('Volume (cc)'); ylabel('Probability of Complication');
        end
        function [CoxComplicationTime,CoxComplicationCurve] = fCoxRiskVDxComplicationFig_DVH(EPGobj,CoxPar,vol) % plot K-M survival curves for patients in the EPGobj and that predicted by Cox model "CoxPar" at specified volume "vol"
%             % K-M survival curve
%             EPGobj = EPGobj.ComplicationCurves_DVH();

            % h(t) clean up
            f = find(diff(CoxPar.h(:,1))==0); % find duplicate time values of h(t)
            while ~isempty(f)
                CoxPar.h(f,1) = CoxPar.h(f,1)-eps*10; % adjust it a bit to avoid ambiguius
                f = find(diff(CoxPar.h(:,1))==0); % find duplicate time values of h(t)
            end

            % survival curve from Cox model at volume vol
            CoxComplicationTime = EPGobj.mKaplanMeierCompOverall.mSurvivalTimeSorted{1};
            CoxComplicationTime(EPGobj.mKaplanMeierCompOverall.CensorStatistics{1}(:,1)) = [];
            h = interp1(CoxPar.h(:,1),CoxPar.h(:,2),CoxComplicationTime,'linear','extrap');
            expbetax = exp(CoxPar.beta*vol);
            CoxComplicationCurve = exp( -h * expbetax );
        end
        function fCoxFig_DVH(EPGobj,t)
            % this is only a demo of how to use the cox model to plot response function

            % Vx
            % search the best Cox model
            [allCox,flgCox,flganti] = CoxPar_DVH(EPGobj,'VDx'); % find availabe Cox models
            flgCox(flganti)=false; % anti-correlations were not be considered
            logl = [allCox.logl]'; logl(~flgCox) = -inf; % log likelihood of Cox model, anti-correlation points not counted
            [mx,doseloc]=max(logl); % the best fitting of Cox model
            lowCI68 = mx - 0.5; % 68% confidence
            lowCI95 = mx - 1.96; % 95% confidence

%             num = cellfun(@(x) size(x,1),{allCox.data_exposure});
            figure(1); clf reset; plot(EPGobj.mBinsX(flgCox), [allCox(flgCox).logl],'.-');
            hold on; plot(EPGobj.mBinsX(flgCox),repmat(lowCI68,size(EPGobj.mBinsX(flgCox))),'r--'); hold off;
            hold on; plot(EPGobj.mBinsX(flgCox),repmat(lowCI95,size(EPGobj.mBinsX(flgCox))),'c--'); hold off;
            set(gca,'xminortick','on','yminortick','on');
            set(gca,'box','on');
            xlabel('Dose (Gy)','fontsize',18); ylabel('log likelihood','fontsize',18);
            figure(2); clf reset; semilogy(EPGobj.mBinsX(flgCox),[allCox(flgCox).p],'.-');
            hold on; semilogy(EPGobj.mBinsX(flgCox),repmat(0.05,size(EPGobj.mBinsX(flgCox))),'r--'); hold off;
            set(gca,'xminortick','on','yminortick','on');
            set(gca,'box','on');
            xlabel('Dose (Gy)','fontsize',18); ylabel('p-value','fontsize',18);
            p = ones(size(allCox));
            p(flgCox) = [allCox(flgCox).p]';
            [m,ploc] = min(p);
            disp(['best fit by log likelihood and p-value: ',num2str([doseloc,ploc])]);
            
%             doseloc = 31;
            allCox = allCox(doseloc);
            disp('Best Cox Model:');
            disp(allCox);
%             disp([EPGobj.mBinsX(doseloc),m,allCox.p]); % display the selected Cox model
            
%             % median patient complication time
%             f = EPGobj.fPatientsWithComplicationData(); % select patients with data
%             CG = EPGobj.fRemovePatient(~f);
%             f2 = ~cellfun('isempty',{CG.mGrp.CompOccurDate}); % patients with no complication date
%             compdate = inf(CG.mNumInGrp,1);
%             compdate(f2) = ([CG.mGrp(f2).CompOccurDate] - [CG.mGrp(f2).BaselineDate])' / 30;
%             t = median(compdate(isfinite(compdate)));
% %             t = 50;
%             flg = compdate > t;
% 
% %             f3 = ~cellfun('isempty',{CG.mGrp.LastFollowupDate}); % patients with no last follow up date
% %             lastfollowup = inf(CG.mNumInGrp,1);
% %             lastfollowup(f3) = ([CG.mGrp(f3).LastFollowupDate] - [CG.mGrp(f3).BaselineDate])' / 30;
% %             compdate = min( lastfollowup, compdate );
% 
%             % observed data at time t
%             % volumes of patients
%             vd=zeros(CG.mNumInGrp,1);
%             d = CG.mBinsX(doseloc);
%             for k=1:CG.mNumInGrp
%                 vd(k) = CG.mGrp(k).VolAtDose( d );
%             end
%             % quatiles of patients
%             numintv = 4;
%             [sortQ,indxQ,indxorg] = EqualIntervals(vd,numintv);
%             meanvol = zeros(numintv,1);
%             prob = zeros(numintv,1);
%             stdprob = zeros(numintv,1);
%             betainv84 = zeros(numintv,1);
%             betainv16 = zeros(numintv,1);
%             f = ~flg(indxorg);
%             for m = 1 : numintv
%                 meanvol(m) = mean(sortQ(indxQ==m));
%                 numcomp = sum(f(indxQ==m)); numtotal = length(find(indxQ==m));
%                 prob(m) = numcomp/numtotal;
%                 stdprob(m) = sqrt(numtotal*prob(m)*(1-prob(m)))/numtotal;
%                 betainv84(m) = betainv( .84, numcomp+1, numtotal - numcomp + 1 );
%                 betainv16(m) = betainv( .16, numcomp+1, numtotal - numcomp + 1 );
%             end
%             figure(3); clf reset; hold on;
%             errorbar(meanvol,prob,max(0,prob-betainv16),max(0,betainv84-prob),'r*','linewidth',1,'markersize',12);
% 
%             % plot the cox model at time point t
%             % h(t)
%             f = find(diff(allCox.h(:,1))==0); % find duplicate time values of h(t)
%             allCox.h(f+1,1) = allCox.h(f,1)+1e-6; % adjust it a bit to avoid ambiguius
%             h = interp1(allCox.h(:,1),allCox.h(:,2),t,'linear','extrap');
%             % cox curve
%             vol = 0:meanvol(end)+10;
%             yfit = h * exp(allCox.beta * vol);
% %             yfit = log(yfit);
%             plot(vol,yfit); grid on;
%             yfitu = h * exp((allCox.beta+1.96*allCox.se) * vol); % upper CI
%             plot(vol,yfitu,'--');
%             yfitl = h * exp((allCox.beta-1.96*allCox.se) * vol); % lower CI
%             plot(vol,yfitl,'--');
%             hold off;
%             set(gca,'xminortick','on','yminortick','on');
%             set(gca,'box','on');
%             xlabel('Volume (cc)'); ylabel('Probability of Complication');
        end
        function fDVHCurvesSummary_DVH(EPGobj)
            % prepare
            f = [EPGobj.mGrp.mFlgCensor]; % censor info
            dosebins = cat(1,EPGobj.mGrp.DoseBins_LQ); % all doses
            dosebins = (0:max(dosebins))'; % dose bins for complication patients
            vol_center = zeros(length(dosebins),5); % 5 columns for -95%, -68%, median, 68%, 95% lines
            vol_comp = vol_center;

            % volume computation
            vol = -inf(EPGobj.mNumInGrp,1);
            for kk = 1:length(dosebins)
                % volumes of each patient at dose kk
                vol(:) = -inf;
                for mm = 1:EPGobj.mNumInGrp
                    vol(mm) = EPGobj.mGrp(mm).VolAtDose( dosebins(kk) );
                end
                % volumes of censored patients
                v = sort(vol(f)); v(v==0) = []; % censored patients, but patients with no volumes at current dose shall be excluded
                if ~isempty(v)
                    vl = length(v);
                    num95 = round((vl-0.95*vl)/2);
                    num68 = round((vl-0.68*vl)/2);
                    if num95>0
                        vol_center(kk,1) = v(num95);
                        vol_center(kk,5) = v(vl-num95);
                    else
                        vol_center(kk,1) = -inf;
                        vol_center(kk,5) = -inf;
                    end
                    if num68>0
                        vol_center(kk,2) = v(num68);
                        vol_center(kk,4) = v(vl-num68);
                    else
                        vol_center(kk,2) = -inf;
                        vol_center(kk,4) = -inf;
                    end
                    vol_center(kk,3) = median(v);
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
            plot(dosebins,vol_center(:,3),'b.-');
            plot(dosebins,vol_center(:,[2,4]),'b');
            plot(dosebins,vol_center(:,[1,5]),'b--');
            plot(dosebins,vol_comp(:,3),'r.-');
            plot(dosebins,vol_comp(:,[2,4]),'r');
            plot(dosebins,vol_comp(:,[1,5]),'r--');
            set(gca,'xminortick','on','yminortick','on');
            set(gca,'box','on');
            xlabel('Dose (Gy)'); ylabel('Volume (cc)');
            hold off;
        end
    end

    methods % EUD operations
        function EPGobj = fCalculateDoseBins_EUD(EPGobj)
            euds = [EPGobj.mGrp.EUD]'; dmax = max(euds(:));
            EPGobj.mEUD_DoseBins = (0:EPGobj.mEUD_DoseStep:dmax+EPGobj.mEUD_DoseStep)';
        end
        function EPGobj = fCalculateDoseBinsLog_EUD(EPGobj)
            euds = [EPGobj.mGrp.EUD]'; dmax = log10(max(euds(:)));
%             f = euds(:)>0; dmin = log10(min(euds(f)));
%             dosebinslog = ((dmin-EPGobj.mEUD_DoseStep): EPGobj.mEUD_DoseStep:(dmax+EPGobj.mEUD_DoseStep))';
            dosebinslog = 0 : EPGobj.mEUD_DoseStep : (dmax+EPGobj.mEUD_DoseStep)';
            EPGobj.mEUD_DoseBins = 10.^dosebinslog;
        end
        function EPGobj = fCalculateEUD(EPGobj)
            for k = 1:EPGobj.mNumInGrp
                EPGobj.mGrp(k) = EPGobj.mGrp(k).CalculateEUD();
            end
        end
        function EPGobj = fCrudeAtlas_EUD(EPGobj)
            % extract all EUDs from all patients
                euds = [EPGobj.mGrp.EUD]';
%                 
%                 euds = zeros( EPGobj.mNumInGrp, size(EPGobj.mLymanN,1) );
%                 for k = 1:size(EPGobj.mLymanN,1)
%                     euds(k,:) = EPGobj.mGrp(k).EUD';
%                 end
%                 dmax = max(euds(:));
                
            % generate dose bins if mEUD_DoseBins is not specifically assigned
                if isempty(EPGobj.mEUD_DoseBins)
                    error('mEUD_DoseBins not determined before method "CrudeAtlas_EUD" is called (classComplicationGroup)');
                end
                numdosebins=size(EPGobj.mEUD_DoseBins,1);
                
            % censor and complication info
                flgcensor = [EPGobj.mGrp.mFlgCensor]; flgcomp = ~flgcensor; % a patient either was censored or had complication
                
            % for each log10(n) and each dose step, compute the total patients and their complications
                EPGobj.mEUD_AtlasTotal = zeros( numdosebins, size(EPGobj.mLgN,1) );
                EPGobj.mEUD_AtlasComp = EPGobj.mEUD_AtlasTotal;
                for n = 1:size(EPGobj.mLgN,1)
                    for m = 1:numdosebins
                        f = find( euds(:,n) >= EPGobj.mEUD_DoseBins(m) );
                        g = find( flgcomp(f) );
                        EPGobj.mEUD_AtlasTotal(m,n) = length(f);
                        EPGobj.mEUD_AtlasComp(m,n) = length(g);
                    end
                end
        end
        function EPGobj = fBetaCumProb_EUD(EPGobj)
            EPGobj.mBetaCumMat = zeros( [size(EPGobj.mEUD_AtlasTotal), size(EPGobj.mBetaCumTh,1)] );
            for k = 1:size(EPGobj.mBetaCumTh,1)
                EPGobj.mBetaCumMat(:,:,k) = betacdf( EPGobj.mBetaCumTh(k), EPGobj.mEUD_AtlasComp+1, EPGobj.mEUD_AtlasTotal - EPGobj.mEUD_AtlasComp + 1 );
            end
        end
        function EPGobj = fBetaInvProb_EUD(EPGobj)
            EPGobj.mBetaInvMat = zeros( [size(EPGobj.mEUD_AtlasTotal), size(EPGobj.mBetaInvTh,1)] );
            for k=1:length(EPGobj.mBetaInvTh)
                EPGobj.mBetaInvMat(:,:,k) = betainv( EPGobj.mBetaInvTh(k), EPGobj.mEUD_AtlasComp+1, EPGobj.mEUD_AtlasTotal - EPGobj.mEUD_AtlasComp + 1 );
            end
        end
        
        function EPGobj = fLogitAnalysisExact_EUD(EPGobj)
            EPGobj.mLogitMat = repmat( struct('b',[],'dev',[],'stats',[]), [size(EPGobj.mLgN,1),1] );
            
            % using exact EUD
            euds = [EPGobj.mGrp.EUD]';
            pttotal = ones(EPGobj.mNumInGrp,1);
            ptcomp = ones(EPGobj.mNumInGrp,1); ptcomp([EPGobj.mGrp.mFlgCensor])=0;
            for k=1:size(EPGobj.mLgN,1)
                doses=euds(:,k);
                % regression using exact EUD
                [b,dev,s]=glmfit(doses,[ptcomp pttotal],'binomial','link','logit');
                EPGobj.mLogitMat(k).b=b;
                EPGobj.mLogitMat(k).dev=dev;
                EPGobj.mLogitMat(k).stats=s;
            end
            
        end
        function EPGobj = fLogitAnalysisBin_EUD(EPGobj)
            EPGobj.mLogitMatBin = repmat( struct('b',[],'dev',[],'stats',[]), [size(EPGobj.mLgN,1),1] );
            
            % using bins
            if isempty(EPGobj.mEUD_DoseBins) || isempty(EPGobj.mEUD_AtlasTotal)
                error('mEUD_DoseBins or mEUD_AtlasTotal not determined before method "LogitAnalysisBin_EUD" is called (classComplicationGroup)');
            end
            
            dosebins = (EPGobj.mEUD_DoseBins(1:end-1)+EPGobj.mEUD_DoseBins(2:end))/2; % dose bins are at the middle of the intervals
            pttotal = ones( EPGobj.mEUD_AtlasTotal(1,1), 1 ); % each patient has his own row
            ptcomp = true( EPGobj.mEUD_AtlasTotal(1,1), 1 ); % allocate space for complication of each patient
            doseval = -inf(EPGobj.mNumInGrp,1); % allocate space for dose of each patient.

            % patient complication and censor info
            pta = EPGobj.mEUD_AtlasTotal; % patient total from atlas
            pca = EPGobj.mEUD_AtlasComp; % patient complication from atlas
            pta = pta - pca; % censored info of patients
            
            pta(1:end-1,:) = abs(diff(pta,1,1)); pta(end,:) = 0; % locations of patients with censor info
            pca(1:end-1,:) = abs(diff(pca,1,1)); pca(end,:) = 0; % locations of patient complications
            
            warning('off','stats:glmfit:IterationLimit');
            for k=1:length(EPGobj.mLgN) % Logistic Regression for each ln(n)
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
                EPGobj.mLogitMatBin(k).b=b;
                EPGobj.mLogitMatBin(k).dev=dev;
                EPGobj.mLogitMatBin(k).stats=s;
            end
            warning('on','stats:glmfit:IterationLimit');



%             EPGobj.mLogitMatBin = repmat( struct('b',[],'dev',[],'stats',[]), [size(EPGobj.mLgN,1),1] );
%             
%             % using bins
%             if isempty(EPGobj.mEUD_DoseBins)
%                 error('mEUD_DoseBins not determined before method "LogitAnalysis_EUD" is called (classComplicationGroup)');
%             end
%             dosebins = (EPGobj.mEUD_DoseBins(1:end-1)+EPGobj.mEUD_DoseBins(2:end))/2; % dose bins are at the middle of the intervals
%             pttotal = ones( EPGobj.mNumInGrp, 1 ); % each patient has his own row
%             ptcomp = ones( EPGobj.mNumInGrp, 1 ); ptcomp(EPGobj.mGrp.mFlgCensor)=0; % allocate space for complication of each patient
%             doseval = -inf(EPGobj.mNumInGrp,1); % allocate space for dose of each patient.
%             euds = [EPGobj.mGrp.EUD]';
%             warning('off','stats:glmfit:IterationLimit');
%             for k=1:length(EPGobj.mLgN)
%                 % determine the dose for each patient
%                 for m = 1:length(EPGobj.mEUD_DoseBins)-1
%                     f = euds(:,k)>=EPGobj.mEUD_DoseBins(m);
%                     doseval(f) = dosebins(m);
%                 end
%                 f = euds(:,k)>=EPGobj.mEUD_DoseBins(end); % patient with dose outside the required dose bins should be removed
%                 doseval(f) = min(EPGobj.mEUD_DoseBins)-1;
%                 
%                 f = doseval>=min(EPGobj.mEUD_DoseBins); % pick up the patients whose dose fall into the dose bins
%                 % logistic regression
%                 [b,dev,s]=glmfit(doseval(f),[ptcomp(f) pttotal(f)],'binomial','link','logit');
%                 EPGobj.mLogitMatBin(k).b=b;
%                 EPGobj.mLogitMatBin(k).dev=dev;
%                 EPGobj.mLogitMatBin(k).stats=s;
%             end
%             warning('on','stats:glmfit:IterationLimit');
        end
        function EPGobj = fLogitHosmerLemeshowTestExact_EUD(EPGobj,loga)
            % loga determination
            st = [EPGobj.mLogitMat];
            dpf = [st.dev]; % deviations
            [~,loc] = min(dpf); % the minimum deviation corresponds to the maximum likelihood
            disp(['The best loga of Logistic Regression is: ',num2str(-EPGobj.mLgN(loc))]);
            if ~exist('loga','var') || ischar(loga)
                loga = -EPGobj.mLgN(loc);
            end
            [~,loc] = min(abs(EPGobj.mLgN - (-loga))); % find the loga whose goodness of fit is computed
            disp(['The loga for the goodness of fit computation is: ',num2str(-EPGobj.mLgN(loc))]);

            % group number
            if isequal(EPGobj.mLogitHosmerLemeshow.n, 0)
                if EPGobj.mNumInGrp<100
                    EPGobj.mLogitHosmerLemeshow.n = 5;
                else
                    EPGobj.mLogitHosmerLemeshow.n = 10;
                end
            end
%             EPGobj.mLogitHosmerLemeshow.n = EPGobj.mNumInGrp;

            % prepare goodness of fit computation
            flg=[EPGobj.mGrp.mFlgCensor]; % censor flags of patients
            euds = [EPGobj.mGrp.EUD]; euds = euds(loc,:); % the gEUDs of that n
            st = EPGobj.mLogitMat(loc); % the fitting result of that n
            disp(['The b0 and b1 are: ',num2str(st.b')]);
            % goodness of fit
            [medianVal,numComp,numTotal,BetaInv84,BetaInv16] = EventObserved(flg,euds,EPGobj.mLogitHosmerLemeshow.n); % observations
            [rpb] = glmval(st.b, medianVal,'logit',st.stats); % predictions
            numE= rpb.*numTotal; % expectations
            EPGobj.mLogitHosmerLemeshow.Chi2 = sum( (numComp-numE).^2 ./ (numTotal.*rpb.*(1-rpb)) );
            EPGobj.mLogitHosmerLemeshow.df = EPGobj.mLogitHosmerLemeshow.n - (length(st.b)+1); % (length(st.b)+1) because the loga is not counted in st.b
            EPGobj.mLogitHosmerLemeshow.p_value = 1 - chi2cdf( EPGobj.mLogitHosmerLemeshow.Chi2, EPGobj.mLogitHosmerLemeshow.df );
            disp(EPGobj.mLogitHosmerLemeshow);
        end
        function EPGobj = fLogitHosmerLemeshowTestBin_EUD(EPGobj,loga)
            % loga determination
            st = [EPGobj.mLogitMatBin];
            dpf = [st.dev]; % deviations
            [~,loc] = min(dpf); % the minimum deviation corresponds to the maximum likelihood
            disp(['The best loga of Logistic Regression is: ',num2str(-EPGobj.mLgN(loc))]);
            if ~exist('loga','var') || ischar(loga)
                loga = -EPGobj.mLgN(loc);
            end
            [~,loc] = min(abs(EPGobj.mLgN - (-loga))); % find the loga whose goodness of fit is computed
            disp(['The loga for the goodness of fit computation is: ',num2str(-EPGobj.mLgN(loc))]);

            % group number
            if isequal(EPGobj.mLogitHosmerLemeshowBin.n, 0)
                if EPGobj.mNumInGrp<100
                    EPGobj.mLogitHosmerLemeshowBin.n = 5;
                else
                    EPGobj.mLogitHosmerLemeshowBin.n = 10;
                end
            end
%             EPGobj.mLogitHosmerLemeshowBin.n = EPGobj.mNumInGrp;

            % prepare goodness of fit computation
            flg=[EPGobj.mGrp.mFlgCensor]; % censor flags of patients
            euds = [EPGobj.mGrp.EUD]; euds = euds(loc,:); % the gEUDs of that n
            st = EPGobj.mLogitMatBin(loc); % the fitting result of that n
            disp(['The b0 and b1 are: ',num2str(st.b')]);
            % goodness of fit
            [medianVal,numComp,numTotal,BetaInv84,BetaInv16] = EventObserved(flg,euds,EPGobj.mLogitHosmerLemeshowBin.n); % observations
            [rpb] = glmval(st.b, medianVal,'logit',st.stats); % predictions
            numE= rpb.*numTotal; % expectations
            EPGobj.mLogitHosmerLemeshowBin.Chi2 = sum( (numComp-numE).^2 ./ (numTotal.*rpb.*(1-rpb)) );
            EPGobj.mLogitHosmerLemeshowBin.df = EPGobj.mLogitHosmerLemeshowBin.n - (length(st.b)+1); % (length(st.b)+1) because the loga is not counted in st.b
            EPGobj.mLogitHosmerLemeshowBin.p_value = 1 - chi2cdf( EPGobj.mLogitHosmerLemeshowBin.Chi2, EPGobj.mLogitHosmerLemeshowBin.df );
            disp(EPGobj.mLogitHosmerLemeshowBin);
        end
        function EPGobj = fLogitGTestExact_EUD(EPGobj,loga)
            % loga determination
            st = [EPGobj.mLogitMat];
            dpf = [st.dev]; % deviations
            [~,loc] = min(dpf); % the minimum deviation corresponds to the maximum likelihood
            disp(['The best loga of Logistic Regression is: ',num2str(-EPGobj.mLgN(loc))]);
            if ~exist('loga','var') || ischar(loga)
                loga = -EPGobj.mLgN(loc);
            end
            [~,loc] = min(abs(EPGobj.mLgN - (-loga))); % find the loga whose goodness of fit is computed
            disp(['The loga for the goodness of fit computation is: ',num2str(-EPGobj.mLgN(loc))]);

            % group number
            if isequal(EPGobj.mLogitGTest.n, 0)
                if EPGobj.mNumInGrp<100
                    EPGobj.mLogitGTest.n = 5;
                else
                    EPGobj.mLogitGTest.n = 10;
                end
            end
%             EPGobj.mLogitGTest.n = EPGobj.mNumInGrp;

            % prepare goodness of fit computation
            flg=[EPGobj.mGrp.mFlgCensor]; % censor flags of patients
            euds = [EPGobj.mGrp.EUD]; euds = euds(loc,:); % the gEUDs of that n
            st = EPGobj.mLogitMat(loc); % the fitting result of that n
            disp(['The b0 and b1 are: ',num2str(st.b')]);
            % goodness of fit
            [medianVal,numComp,numTotal,BetaInv84,BetaInv16] = EventObserved(flg,euds,EPGobj.mLogitGTest.n); % observations
            [rpb] = glmval(st.b, medianVal,'logit',st.stats); % predictions
            numE= rpb.*numTotal; % expectations
            f = find(numComp);
            EPGobj.mLogitGTest.Chi2 = 2 * sum( numComp(f) .* log(numComp(f)./numE(f)) );
            EPGobj.mLogitGTest.df = EPGobj.mLogitGTest.n - (length(st.b)+1); % (length(st.b)+1) because the loga is not counted in st.b
            EPGobj.mLogitGTest.p_value = 1 - chi2cdf( EPGobj.mLogitGTest.Chi2, EPGobj.mLogitGTest.df );
            disp(EPGobj.mLogitGTest);
        end
        function EPGobj = fLogitGTestBin_EUD(EPGobj,loga)
            % loga determination
            st = [EPGobj.mLogitMatBin];
            dpf = [st.dev]; % deviations
            [~,loc] = min(dpf); % the minimum deviation corresponds to the maximum likelihood
            disp(['The best loga of Logistic Regression is: ',num2str(-EPGobj.mLgN(loc))]);
            if ~exist('loga','var') || ischar(loga)
                loga = -EPGobj.mLgN(loc);
            end
            [~,loc] = min(abs(EPGobj.mLgN - (-loga))); % find the loga whose goodness of fit is computed
            disp(['The loga for the goodness of fit computation is: ',num2str(-EPGobj.mLgN(loc))]);

            % group number
            if isequal(EPGobj.mLogitGTestBin.n, 0)
                if EPGobj.mNumInGrp<100
                    EPGobj.mLogitGTestBin.n = 5;
                else
                    EPGobj.mLogitGTestBin.n = 10;
                end
            end
%             EPGobj.mLogitGTestBin.n = EPGobj.mNumInGrp;

            % prepare goodness of fit computation
            flg=[EPGobj.mGrp.mFlgCensor]; % censor flags of patients
            euds = [EPGobj.mGrp.EUD]; euds = euds(loc,:); % the gEUDs of that n
            st = EPGobj.mLogitMatBin(loc); % the fitting result of that n
            disp(['The b0 and b1 are: ',num2str(st.b')]);
            % goodness of fit
            [medianVal,numComp,numTotal,BetaInv84,BetaInv16] = EventObserved(flg,euds,EPGobj.mLogitGTestBin.n); % observations
            [rpb] = glmval(st.b, medianVal,'logit',st.stats); % predictions
            numE= rpb.*numTotal; % expectations
            f = find(numComp);
            EPGobj.mLogitGTestBin.Chi2 = 2 * sum( numComp(f) .* log(numComp(f)./numE(f)) );
            EPGobj.mLogitGTestBin.df = EPGobj.mLogitGTestBin.n - (length(st.b)+1); % (length(st.b)+1) because the loga is not counted in st.b
            EPGobj.mLogitGTestBin.p_value = 1 - chi2cdf( EPGobj.mLogitGTestBin.Chi2, EPGobj.mLogitGTestBin.df );
            disp(EPGobj.mLogitGTestBin);
        end
        function EPGobj = fLogitPearsonTestExact_EUD(EPGobj,loga)
            % loga determination
            st = [EPGobj.mLogitMat];
            dpf = [st.dev]; % deviations
            [~,loc] = min(dpf); % the minimum deviation corresponds to the maximum likelihood
            disp(['The best loga of Logistic Regression is: ',num2str(-EPGobj.mLgN(loc))]);
            if ~exist('loga','var') || ischar(loga)
                loga = -EPGobj.mLgN(loc);
            end
            [~,loc] = min(abs(EPGobj.mLgN - (-loga))); % find the loga whose goodness of fit is computed
            disp(['The loga for the goodness of fit computation is: ',num2str(-EPGobj.mLgN(loc))]);

            % group number
            if isequal(EPGobj.mLogitPearson.n, 0)
                if EPGobj.mNumInGrp<100
                    EPGobj.mLogitPearson.n = 5;
                else
                    EPGobj.mLogitPearson.n = 10;
                end
            end

            % prepare goodness of fit computation
            flg=[EPGobj.mGrp.mFlgCensor]; % censor flags of patients
            euds = [EPGobj.mGrp.EUD]; euds = euds(loc,:); % the gEUDs of that n
            st = EPGobj.mLogitMat(loc); % the fitting result of that n
            disp(['The b0 and b1 are: ',num2str(st.b')]);
            % goodness of fit
            [medianVal,numComp,numTotal,BetaInv84,BetaInv16] = EventObserved(flg,euds,EPGobj.mLogitPearson.n); % observations
            [rpb] = glmval(st.b, medianVal,'logit',st.stats); % predictions
            numE= rpb.*numTotal; % expectations
            EPGobj.mLogitPearson.Chi2 = sum( (numComp-numE).^2 ./ numE );
            EPGobj.mLogitPearson.df = EPGobj.mLogitPearson.n - (length(st.b)+1); % (length(st.b)+1) because the loga is not counted in st.b
            EPGobj.mLogitPearson.p_value = 1 - chi2cdf( EPGobj.mLogitPearson.Chi2, EPGobj.mLogitPearson.df );
            disp(EPGobj.mLogitPearson);
        end
        function EPGobj = fLogitPearsonTestBin_EUD(EPGobj,loga)
            % loga determination
            st = [EPGobj.mLogitMatBin];
            dpf = [st.dev]; % deviations
            [~,loc] = min(dpf); % the minimum deviation corresponds to the maximum likelihood
            disp(['The best loga of Logistic Regression is: ',num2str(-EPGobj.mLgN(loc))]);
            if ~exist('loga','var') || ischar(loga)
                loga = -EPGobj.mLgN(loc);
            end
            [~,loc] = min(abs(EPGobj.mLgN - (-loga))); % find the loga whose goodness of fit is computed
            disp(['The loga for the goodness of fit computation is: ',num2str(-EPGobj.mLgN(loc))]);

            % group number
            if isequal(EPGobj.mLogitPearsonBin.n, 0)
                if EPGobj.mNumInGrp<100
                    EPGobj.mLogitPearsonBin.n = 5;
                else
                    EPGobj.mLogitPearsonBin.n = 10;
                end
            end

            % prepare goodness of fit computation
            flg=[EPGobj.mGrp.mFlgCensor]; % censor flags of patients
            euds = [EPGobj.mGrp.EUD]; euds = euds(loc,:); % the gEUDs of that n
            st = EPGobj.mLogitMatBin(loc); % the fitting result of that n
            disp(['The b0 and b1 are: ',num2str(st.b')]);
            % goodness of fit
            [medianVal,numComp,numTotal,BetaInv84,BetaInv16] = EventObserved(flg,euds,EPGobj.mLogitPearsonBin.n); % observations
            [rpb] = glmval(st.b, medianVal,'logit',st.stats); % predictions
            numE= rpb.*numTotal; % expectations
            EPGobj.mLogitPearsonBin.Chi2 = sum( (numComp-numE).^2 ./ numE );
            EPGobj.mLogitPearsonBin.df = EPGobj.mLogitPearsonBin.n - (length(st.b)+1); % (length(st.b)+1) because the loga is not counted in st.b
            EPGobj.mLogitPearsonBin.p_value = 1 - chi2cdf( EPGobj.mLogitPearsonBin.Chi2, EPGobj.mLogitPearsonBin.df );
            disp(EPGobj.mLogitPearsonBin);
        end
        function EPGobj = fLogitAnalysisGridExact_EUD(EPGobj)
            % preparation
            euds = [EPGobj.mGrp.EUD]';
            flg = [EPGobj.mGrp.mFlgCensor]';
            
            b0 = EPGobj.mLogitGridBetaRange{1};
            b1 = EPGobj.mLogitGridBetaRange{2};

            % for each mLgN,b0, and b1, compute the log likelihood
            loglikelihood = -inf(length(b0),length(b1),length(EPGobj.mLgN));
            for kk = 1:length(EPGobj.mLgN)
                for jj = 1:length(b1)
                    for ii = 1:length(b0)
                        pr = exp(b0(ii)+b1(jj)*euds(:,kk));
                        pr = pr./(1+pr); % logistic probability
                        pr(flg) = 1-pr(flg); % non-complication patients
                        pr = log(pr); % log likelihood of each patients
                        loglikelihood(ii,jj,kk) = sum(pr); % loglikelihood of all
                    end
                end
            end
            EPGobj.mLogitGrid = struct('b0',b0,'b1',b1,'loglikelihood',loglikelihood);
        end
        function EPGobj = fLogitAnalysisGridBin_EUD(EPGobj)
            % parse atlas (part)
            dosebins = (EPGobj.mEUD_DoseBins(1:end-1)+EPGobj.mEUD_DoseBins(2:end))/2; % dose bins are at the middle of the intervals
            ptcomp = true( EPGobj.mEUD_AtlasTotal(1,1), 1 ); % allocate space for complication of each patient
            doseval = -inf(EPGobj.mNumInGrp,1); % allocate space for dose of each patient.

            pta = EPGobj.mEUD_AtlasTotal; % patient total from atlas
            pca = EPGobj.mEUD_AtlasComp; % patient complication from atlas
            pta = pta - pca; % censored info of patients
            pta(1:end-1,:) = abs(diff(pta,1,1)); pta(end,:) = 0; % locations of patients with censor info
            pca(1:end-1,:) = abs(diff(pca,1,1)); pca(end,:) = 0; % locations of patient complications

            % preparation
            b0 = EPGobj.mLogitGridBetaRange{1};
            b1 = EPGobj.mLogitGridBetaRange{2};

            % for each mLgN, TD50, and m, compute the log likelihood
            loglikelihood = -inf(length(b0),length(b1),length(EPGobj.mLgN));
            for kk = 1:length(EPGobj.mLgN)
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
                for jj = 1:length(b1)
                    for ii = 1:length(b0)
                        pr = exp(b0(ii)+b1(jj)*doseval); % logistic probability
                        pr = pr./(1+pr);
                        pr(flg) = 1-pr(flg); % non-complication patients
                        pr = log(pr); % log likelihood of each patients
                        loglikelihood(ii,jj,kk) = sum(pr); % loglikelihood of all
                    end
                end
            end
            EPGobj.mLogitGridBin = struct('b0',b0,'b1',b1,'loglikelihood',loglikelihood);
        end
        function EPGobj = fLogitGOFSimulationExact_EUD(EPGobj)
            % loga determination
            st = [EPGobj.mLogitMat];
            dpf = [st.dev]; % deviations
            [~,loc] = min(dpf); % the minimum deviation corresponds to the maximum likelihood
            disp(['The best loga of Logistic Regression is: ',num2str(-EPGobj.mLgN(loc))]);
            % fitting result
            st = st(loc);
            % euds
            euds = [EPGobj.mGrp.EUD]; euds = euds(loc,:)'; % the gEUDs of that n

            % expectation for each patient
            numE = glmval(st.b,euds,'logit',st.stats);

            % SSR at simulated obserations
            numSim = 100000; % number of simulations
            ssr = zeros(numSim,1); % SSR
            for k = 1:numSim
                % generate observations
                numO = rand([EPGobj.mNumInGrp,1]); % random number with uniform distribution on (0,1)
                f = numO<=numE; numO(f) = 1; numO(~f) = 0; % the lower the expectation, the lower the probability it has a complication
                % SSR
                ssr(k) = sum( (numO-numE).^2 );
            end

            % p-value
            flg=[EPGobj.mGrp.mFlgCensor]'; % censor flags of patients
            ssrObsv = sum( (~flg-numE).^2 ); % the actual ssr
            f = ssrObsv <= ssr; % the simulatations that observed SSR is less than simulation
            p = min(sum(f),sum(~f)); % smaller one shows the extreme situations at the two ends
            p = p/numSim;

            % save
            EPGobj.mLogitGOFSim.SSRSim = ssr;
            EPGobj.mLogitGOFSim.SSRObserve = ssrObsv;
            EPGobj.mLogitGOFSim.p_value = p;
        end
        function EPGobj = fLogitGOFSimulationBin_EUD(EPGobj)
            % loga determination
            st = [EPGobj.mLogitMatBin];
            dpf = [st.dev]; % deviations
            [~,loc] = min(dpf); % the minimum deviation corresponds to the maximum likelihood
            disp(['The best loga of Logistic Regression is: ',num2str(-EPGobj.mLgN(loc))]);
            % fitting result
            st = st(loc);
            % euds
            euds = [EPGobj.mGrp.EUD]; euds = euds(loc,:)'; % the gEUDs of that n

            % expectation for each patient
            numE = glmval(st.b,euds,'logit',st.stats);

            % SSR at simulated obserations
            numSim = 100000; % number of simulations
            ssr = zeros(numSim,1); % SSR
            for k = 1:numSim
                % generate observations
                numO = rand([EPGobj.mNumInGrp,1]); % random number with uniform distribution on (0,1)
                f = numO<=numE; numO(f) = 1; numO(~f) = 0; % the lower the expectation, the lower the probability it has a complication
                % SSR
                ssr(k) = sum( (numO-numE).^2 );
            end

            % p-value
            flg=[EPGobj.mGrp.mFlgCensor]'; % censor flags of patients
            ssrObsv = sum( (~flg-numE).^2 ); % the actual ssr
            f = ssrObsv <= ssr; % the simulatations that observed SSR is less than simulation
            p = min(sum(f),sum(~f)); % smaller one shows the extreme situations at the two ends
            p = p/numSim;

            % save
            EPGobj.mLogitGOFSimBin.SSRSim = ssr;
            EPGobj.mLogitGOFSimBin.SSRObserve = ssrObsv;
            EPGobj.mLogitGOFSimBin.p_value = p;
        end

        function EPGobj = fLymanAnalysisGridExact_EUD(EPGobj)
            % preparation
            euds = [EPGobj.mGrp.EUD]';
            flg = [EPGobj.mGrp.mFlgCensor]';
            
            TD50 = EPGobj.mLymanGridTD50Range;
            m = EPGobj.mLymanGridMRange;
%             m = [0:0.01:1, 1.1:0.1:2, 3:1:10]; m(1) = 0.001;
%             lgm = -1:0.01:1;
%             m = 10.^lgm';

            % for each mLgN,TD50, and m, compute the log likelihood
            loglikelihood = -inf(length(TD50),length(m),length(EPGobj.mLgN));
            for kk = 1:length(EPGobj.mLgN)
                for jj = 1:length(m)
                    for ii = 1:length(TD50)
                        pr = normcdf((euds(:,kk)-TD50(ii))/(m(jj)*TD50(ii)),0,1); % Lyman probability
                        pr(flg) = 1-pr(flg); % non-complication patients
                        pr = log(pr); % log likelihood of each patients
                        loglikelihood(ii,jj,kk) = sum(pr); % loglikelihood of all
                    end
                end
            end
            EPGobj.mLymanGrid = struct('TD50',TD50,'m',m,'loglikelihood',loglikelihood);
        end
        function EPGobj = fLymanAnalysisGridBin_EUD(EPGobj)
            % preparation
            dosebins = (EPGobj.mEUD_DoseBins(1:end-1)+EPGobj.mEUD_DoseBins(2:end))/2; % dose bins are at the middle of the intervals
            ptcomp = true( EPGobj.mEUD_AtlasTotal(1,1), 1 ); % allocate space for complication of each patient
            doseval = -inf(EPGobj.mNumInGrp,1); % allocate space for dose of each patient.

            pta = EPGobj.mEUD_AtlasTotal; % patient total from atlas
            pca = EPGobj.mEUD_AtlasComp; % patient complication from atlas
            pta = pta - pca; % censored info of patients
            pta(1:end-1,:) = abs(diff(pta,1,1)); pta(end,:) = 0; % locations of patients with censor info
            pca(1:end-1,:) = abs(diff(pca,1,1)); pca(end,:) = 0; % locations of patient complications

            TD50 = EPGobj.mLymanGridTD50Range;
            m = EPGobj.mLymanGridMRange;

            % for each mLgN, TD50, and m, compute the log likelihood
            loglikelihood = -inf(length(TD50),length(m),length(EPGobj.mLgN));
            for kk = 1:length(EPGobj.mLgN)
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
                        pr = log(pr); % log likelihood of each patients
                        loglikelihood(ii,jj,kk) = sum(pr); % loglikelihood of all
                    end
                end
            end
            EPGobj.mLymanGridBin = struct('TD50',TD50,'m',m,'loglikelihood',loglikelihood);
        end
        function EPGobj = fLymanHosmerLemeshowTestAnalysisExact_EUD(EPGobj,loga)
            % loga determination
            [~,loc] = max(EPGobj.mLymanGrid.loglikelihood(:));
            [~,~,loc] = ind2sub(size(EPGobj.mLymanGrid.loglikelihood),loc);
            disp(['the "log10(a)" for the maximum likelihood in Lyman model is: ',num2str(-EPGobj.mLgN(loc))]);
            if ~exist('loga','var') || ischar(loga)
                loga = -EPGobj.mLgN(loc);
            end
            [~,loc] = min(abs(EPGobj.mLgN+loga)); 
            disp(['The loga for the goodness of fit computation is: ',num2str(-EPGobj.mLgN(loc))]);
            ll = EPGobj.mLymanGrid.loglikelihood(:,:,loc);
            mx = max(ll(:));
            [dd,mm] = find(ll == mx,1);
            TD50 = EPGobj.mLymanGrid.TD50(dd);
            m = EPGobj.mLymanGrid.m(mm);
            disp(['TD50 & m are: ',num2str([TD50,m])]);

            % group number
            if isequal(EPGobj.mLymanHosmerLemeshow.n, 0)
                if EPGobj.mNumInGrp<100
                    EPGobj.mLymanHosmerLemeshow.n = 5;
                else
                    EPGobj.mLymanHosmerLemeshow.n = 10;
                end
            end
%             EPGobj.mLymanHosmerLemeshow.n = EPGobj.mNumInGrp;

            % prepare goodness of fit computation
            flg=[EPGobj.mGrp.mFlgCensor]; % censor flags of patients
            euds = [EPGobj.mGrp.EUD]; euds = euds(loc,:); % the gEUDs of that n
            % goodness of fit
            [medianVal,numComp,numTotal,~,~] = EventObserved(flg,euds,EPGobj.mLymanHosmerLemeshow.n); % observations
            rpb = normcdf((medianVal-TD50)/(m*TD50),0,1); % Lyman probability
            numE= rpb.*numTotal; % expectations
            EPGobj.mLymanHosmerLemeshow.Chi2 = sum( (numComp-numE).^2 ./ (numTotal.*rpb.*(1-rpb)) );
            EPGobj.mLymanHosmerLemeshow.df = EPGobj.mLymanHosmerLemeshow.n - 2 - 1;
            EPGobj.mLymanHosmerLemeshow.p_value = 1 - chi2cdf( EPGobj.mLymanHosmerLemeshow.Chi2, EPGobj.mLymanHosmerLemeshow.df );
            disp(EPGobj.mLymanHosmerLemeshow);
        end
        function EPGobj = fLymanHosmerLemeshowTestAnalysisBin_EUD(EPGobj,loga)
            % loga determination
            [~,loc] = max(EPGobj.mLymanGridBin.loglikelihood(:));
            [~,~,loc] = ind2sub(size(EPGobj.mLymanGridBin.loglikelihood),loc);
            disp(['the "log10(a)" for the maximum likelihood in Lyman model is: ',num2str(-EPGobj.mLgN(loc))]);
            if ~exist('loga','var') || ischar(loga)
                loga = -EPGobj.mLgN(loc);
            end
            [~,loc] = min(abs(EPGobj.mLgN+loga)); % the n whose corresponding responding function will be ploted
            disp(['The loga for the goodness of fit computation is: ',num2str(-EPGobj.mLgN(loc))]);
            ll = EPGobj.mLymanGridBin.loglikelihood(:,:,loc);
            mx = max(ll(:));
            [dd,mm] = find(ll == mx,1);
            TD50 = EPGobj.mLymanGridBin.TD50(dd);
            m = EPGobj.mLymanGridBin.m(mm);
            disp(['TD50 & m are: ',num2str([TD50,m])]);

            % group number
            if isequal(EPGobj.mLymanHosmerLemeshowBin.n, 0)
                if EPGobj.mNumInGrp<100
                    EPGobj.mLymanHosmerLemeshowBin.n = 5;
                else
                    EPGobj.mLymanHosmerLemeshowBin.n = 10;
                end
            end
%             EPGobj.mLymanHosmerLemeshowBin.n = EPGobj.mNumInGrp;

            % prepare goodness of fit computation
            flg=[EPGobj.mGrp.mFlgCensor]; % censor flags of patients
            euds = [EPGobj.mGrp.EUD]; euds = euds(loc,:); % the gEUDs of that n
            % goodness of fit
            [medianVal,numComp,numTotal,BetaInv84,BetaInv16] = EventObserved(flg,euds,EPGobj.mLymanHosmerLemeshowBin.n); % observations
            rpb = normcdf((medianVal-TD50)/(m*TD50),0,1); % Lyman probability
            numE= rpb.*numTotal; % expectations
            EPGobj.mLymanHosmerLemeshowBin.Chi2 = sum( (numComp-numE).^2 ./ (numTotal.*rpb.*(1-rpb)) );
            EPGobj.mLymanHosmerLemeshowBin.df = EPGobj.mLymanHosmerLemeshowBin.n - 2 - 1;
            EPGobj.mLymanHosmerLemeshowBin.p_value = 1 - chi2cdf( EPGobj.mLymanHosmerLemeshowBin.Chi2, EPGobj.mLymanHosmerLemeshowBin.df );
            disp(EPGobj.mLymanHosmerLemeshowBin);
        end
        function EPGobj = fLymanGTestAnalysisExact_EUD(EPGobj,loga)
            % loga determination
            [~,loc] = max(EPGobj.mLymanGrid.loglikelihood(:));
            [~,~,loc] = ind2sub(size(EPGobj.mLymanGrid.loglikelihood),loc);
            disp(['the "log10(a)" for the maximum likelihood in Lyman model is: ',num2str(-EPGobj.mLgN(loc))]);
            if ~exist('loga','var') || ischar(loga)
                loga = -EPGobj.mLgN(loc);
            end
            [~,loc] = min(abs(EPGobj.mLgN+loga)); % the n whose corresponding responding function will be ploted
            disp(['The loga for the goodness of fit computation is: ',num2str(-EPGobj.mLgN(loc))]);
            ll = EPGobj.mLymanGrid.loglikelihood(:,:,loc);
            mx = max(ll(:));
            [dd,mm] = find(ll == mx,1);
            TD50 = EPGobj.mLymanGrid.TD50(dd);
            m = EPGobj.mLymanGrid.m(mm);
            disp(['TD50 & m are: ',num2str([TD50,m])]);

            % group number
            if isequal(EPGobj.mLymanGTest.n, 0)
                if EPGobj.mNumInGrp<100
                    EPGobj.mLymanGTest.n = 5;
                else
                    EPGobj.mLymanGTest.n = 10;
                end
            end
%             EPGobj.mLymanGTest.n = EPGobj.mNumInGrp;

            % prepare goodness of fit computation
            flg=[EPGobj.mGrp.mFlgCensor]; % censor flags of patients
            euds = [EPGobj.mGrp.EUD]; euds = euds(loc,:); % the gEUDs of that n
            % goodness of fit
            [medianVal,numComp,numTotal,BetaInv84,BetaInv16] = EventObserved(flg,euds,EPGobj.mLymanGTest.n); % observations
            rpb = normcdf((medianVal-TD50)/(m*TD50),0,1); % Lyman probability
            numE= rpb.*numTotal; % expectations
            f = find(numComp);
            EPGobj.mLymanGTest.Chi2 = 2 * sum( numComp(f) .* log(numComp(f)./numE(f)) );
            EPGobj.mLymanGTest.df = EPGobj.mLymanGTest.n - 2 - 1;
            EPGobj.mLymanGTest.p_value = 1 - chi2cdf( EPGobj.mLymanGTest.Chi2, EPGobj.mLymanGTest.df );
            disp(EPGobj.mLymanGTest);
        end
        function EPGobj = fLymanGTestAnalysisBin_EUD(EPGobj,loga)
            % loga determination
            [~,loc] = max(EPGobj.mLymanGridBin.loglikelihood(:));
            [~,~,loc] = ind2sub(size(EPGobj.mLymanGridBin.loglikelihood),loc);
            disp(['the "log10(a)" for the maximum likelihood in Lyman model is: ',num2str(-EPGobj.mLgN(loc))]);
            if ~exist('loga','var') || ischar(loga)
                loga = -EPGobj.mLgN(loc);
            end
            [~,loc] = min(abs(EPGobj.mLgN+loga)); % the n whose corresponding responding function will be ploted
            disp(['The loga for the goodness of fit computation is: ',num2str(-EPGobj.mLgN(loc))]);
            ll = EPGobj.mLymanGridBin.loglikelihood(:,:,loc);
            mx = max(ll(:));
            [dd,mm] = find(ll == mx,1);
            TD50 = EPGobj.mLymanGridBin.TD50(dd);
            m = EPGobj.mLymanGridBin.m(mm);
            disp(['TD50 & m are: ',num2str([TD50,m])]);

            % group number
            if isequal(EPGobj.mLymanGTestBin.n, 0)
                if EPGobj.mNumInGrp<100
                    EPGobj.mLymanGTestBin.n = 5;
                else
                    EPGobj.mLymanGTestBin.n = 10;
                end
            end
%             EPGobj.mLymanGTestBin.n = EPGobj.mNumInGrp;

            % prepare goodness of fit computation
            flg=[EPGobj.mGrp.mFlgCensor]; % censor flags of patients
            euds = [EPGobj.mGrp.EUD]; euds = euds(loc,:); % the gEUDs of that n
            % goodness of fit
            [medianVal,numComp,numTotal,BetaInv84,BetaInv16] = EventObserved(flg,euds,EPGobj.mLymanGTestBin.n); % observations
            rpb = normcdf((medianVal-TD50)/(m*TD50),0,1); % Lyman probability
            numE= rpb.*numTotal; % expectations
            f = find(numComp);
            EPGobj.mLymanGTestBin.Chi2 = 2 * sum( numComp(f) .* log(numComp(f)./numE(f)) );
            EPGobj.mLymanGTestBin.df = EPGobj.mLymanGTestBin.n - 2 - 1;
            EPGobj.mLymanGTestBin.p_value = 1 - chi2cdf( EPGobj.mLymanGTestBin.Chi2, EPGobj.mLymanGTestBin.df );
            disp(EPGobj.mLymanGTestBin);
        end
        function EPGobj = fLymanPearsonTestAnalysisExact_EUD(EPGobj,loga)
            % loga determination
            [~,loc] = max(EPGobj.mLymanGrid.loglikelihood(:));
            [~,~,loc] = ind2sub(size(EPGobj.mLymanGrid.loglikelihood),loc);
            disp(['the "log10(a)" for the maximum likelihood in Lyman model is: ',num2str(-EPGobj.mLgN(loc))]);
            if ~exist('loga','var') || ischar(loga)
                loga = -EPGobj.mLgN(loc);
            end
            [~,loc] = min(abs(EPGobj.mLgN+loga)); % the n whose corresponding responding function will be ploted
            disp(['The loga for the goodness of fit computation is: ',num2str(-EPGobj.mLgN(loc))]);
            ll = EPGobj.mLymanGrid.loglikelihood(:,:,loc);
            mx = max(ll(:));
            [dd,mm] = find(ll == mx,1);
            TD50 = EPGobj.mLymanGrid.TD50(dd);
            m = EPGobj.mLymanGrid.m(mm);
            disp(['TD50 & m are: ',num2str([TD50,m])]);

            % group number
            if isequal(EPGobj.mLymanPearson.n, 0)
                if EPGobj.mNumInGrp<100
                    EPGobj.mLymanPearson.n = 5;
                else
                    EPGobj.mLymanPearson.n = 10;
                end
            end

            % prepare goodness of fit computation
            flg=[EPGobj.mGrp.mFlgCensor]; % censor flags of patients
            euds = [EPGobj.mGrp.EUD]; euds = euds(loc,:); % the gEUDs of that n
            % goodness of fit
            [medianVal,numComp,numTotal,BetaInv84,BetaInv16] = EventObserved(flg,euds,EPGobj.mLymanPearson.n); % observations
            rpb = normcdf((medianVal-TD50)/(m*TD50),0,1); % Lyman probability
            numE= rpb.*numTotal; % expectations
            EPGobj.mLymanPearson.Chi2 = sum( (numComp-numE).^2 ./ numE );
            EPGobj.mLymanPearson.df = EPGobj.mLymanPearson.n - 2 - 1;
            EPGobj.mLymanPearson.p_value = 1 - chi2cdf( EPGobj.mLymanPearson.Chi2, EPGobj.mLymanPearson.df );
            disp(EPGobj.mLymanPearson);
        end
        function EPGobj = fLymanPearsonTestAnalysisBin_EUD(EPGobj,loga)
            % loga determination
            [~,loc] = max(EPGobj.mLymanGridBin.loglikelihood(:));
            [~,~,loc] = ind2sub(size(EPGobj.mLymanGridBin.loglikelihood),loc);
            disp(['the "log10(a)" for the maximum likelihood in Lyman model is: ',num2str(-EPGobj.mLgN(loc))]);
            if ~exist('loga','var') || ischar(loga)
                loga = -EPGobj.mLgN(loc);
            end
            [~,loc] = min(abs(EPGobj.mLgN+loga)); % the n whose corresponding responding function will be ploted
            disp(['The loga for the goodness of fit computation is: ',num2str(-EPGobj.mLgN(loc))]);
            ll = EPGobj.mLymanGridBin.loglikelihood(:,:,loc);
            mx = max(ll(:));
            [dd,mm] = find(ll == mx,1);
            TD50 = EPGobj.mLymanGridBin.TD50(dd);
            m = EPGobj.mLymanGridBin.m(mm);
            disp(['TD50 & m are: ',num2str([TD50,m])]);

            % group number
            if isequal(EPGobj.mLymanPearsonBin.n, 0)
                if EPGobj.mNumInGrp<100
                    EPGobj.mLymanPearsonBin.n = 5;
                else
                    EPGobj.mLymanPearsonBin.n = 10;
                end
            end

            % prepare goodness of fit computation
            flg=[EPGobj.mGrp.mFlgCensor]; % censor flags of patients
            euds = [EPGobj.mGrp.EUD]; euds = euds(loc,:); % the gEUDs of that n
            % goodness of fit
            [medianVal,numComp,numTotal,BetaInv84,BetaInv16] = EventObserved(flg,euds,EPGobj.mLymanPearsonBin.n); % observations
            rpb = normcdf((medianVal-TD50)/(m*TD50),0,1); % Lyman probability
            numE= rpb.*numTotal; % expectations
            EPGobj.mLymanPearsonBin.Chi2 = sum( (numComp-numE).^2 ./ numE );
            EPGobj.mLymanPearsonBin.df = EPGobj.mLymanPearsonBin.n - 2 - 1;
            EPGobj.mLymanPearsonBin.p_value = 1 - chi2cdf( EPGobj.mLymanPearsonBin.Chi2, EPGobj.mLymanPearsonBin.df );
            disp(EPGobj.mLymanPearsonBin);
        end
        function EPGobj = fLymanGOFAnalysisSimulationExact_EUD(EPGobj)
            % loga determination
            [~,loc] = max(EPGobj.mLymanGrid.loglikelihood(:));
            [dd,mm,loc] = ind2sub(size(EPGobj.mLymanGrid.loglikelihood),loc);
            disp(['the "log10(a)" for the maximum likelihood in Lyman model is: ',num2str(-EPGobj.mLgN(loc))]);
            % fitting result
            TD50 = EPGobj.mLymanGrid.TD50(dd);
            m = EPGobj.mLymanGrid.m(mm);
            % euds
            euds = [EPGobj.mGrp.EUD]; euds = euds(loc,:)'; % the gEUDs of that n

            % expectation for each patient
            numE = normcdf((euds-TD50)/(m*TD50),0,1); % Lyman probability

            % SSR at simulated obserations
            numSim = 100000; % number of simulations
            ssr = zeros(numSim,1); % SSR
            for k = 1:numSim
                % generate observations
                numO = rand([EPGobj.mNumInGrp,1]); % random number with uniform distribution on (0,1)
                f = numO<=numE; numO(f) = 1; numO(~f) = 0; % the lower the expectation, the lower the probability it has a complication
                % SSR
                ssr(k) = sum( (numO-numE).^2 );
            end

            % p-value
            flg=[EPGobj.mGrp.mFlgCensor]'; % censor flags of patients
            ssrObsv = sum( (~flg-numE).^2 ); % the actual ssr
            f = ssrObsv <= ssr; % the simulatations that observed SSR is less than simulation
            p = min(sum(f),sum(~f)); % smaller one shows the extreme situations at the two ends
            p = p/numSim;

            % save
            EPGobj.mLymanGOFSim.SSRSim = ssr;
            EPGobj.mLymanGOFSim.SSRObserve = ssrObsv;
            EPGobj.mLymanGOFSim.p_value = p;
        end
        function EPGobj = fLymanGOFAnalysisSimulationBin_EUD(EPGobj)
            % loga determination
            [~,loc] = max(EPGobj.mLymanGridBin.loglikelihood(:));
            [dd,mm,loc] = ind2sub(size(EPGobj.mLymanGridBin.loglikelihood),loc);
            disp(['the "log10(a)" for the maximum likelihood in Lyman model is: ',num2str(-EPGobj.mLgN(loc))]);
            % fitting result
            TD50 = EPGobj.mLymanGridBin.TD50(dd);
            m = EPGobj.mLymanGridBin.m(mm);
            % euds
            euds = [EPGobj.mGrp.EUD]; euds = euds(loc,:)'; % the gEUDs of that n

            % expectation for each patient
            numE = normcdf((euds-TD50)/(m*TD50),0,1); % Lyman probability

            % SSR at simulated obserations
            numSim = 100000; % number of simulations
            ssr = zeros(numSim,1); % SSR
            for k = 1:numSim
                % generate observations
                numO = rand([EPGobj.mNumInGrp,1]); % random number with uniform distribution on (0,1)
                f = numO<=numE; numO(f) = 1; numO(~f) = 0; % the lower the expectation, the lower the probability it has a complication
                % SSR
                ssr(k) = sum( (numO-numE).^2 );
            end

            % p-value
            flg=[EPGobj.mGrp.mFlgCensor]'; % censor flags of patients
            ssrObsv = sum( (~flg-numE).^2 ); % the actual ssr
            f = ssrObsv <= ssr; % the simulatations that observed SSR is less than simulation
            p = min(sum(f),sum(~f)); % smaller one shows the extreme situations at the two ends
            p = p/numSim;

            % save
            EPGobj.mLymanGOFSimBin.SSRSim = ssr;
            EPGobj.mLymanGOFSimBin.SSRObserve = ssrObsv;
            EPGobj.mLymanGOFSimBin.p_value = p;
        end
    end
    
    methods % EUD plot and atlas writing
        function fEUDCurvesFig_a_EUD(EPGobj)
            f=[EPGobj.mGrp.mFlgCensor]; g=find(f);
            a1 = gca;
            a2 = copyobj(a1,gcf);
            set(a2,'Color','none');
            set(a2,'Xtick',[]);
            hold on;
            for m = 1:length(g)
                plot(a1,EPGobj.mGrp(g(m)).EUD,EPGobj.mGrp(g(m)).mLgN,'b');
            end
            g=find(~f);
            for m = 1:length(g)
                plot(a1,EPGobj.mGrp(g(m)).EUD,EPGobj.mGrp(g(m)).mLgN,'r');
            end
            set(a2,'YAxisLocation','right');
            set(a2,'YTickLabel',num2str(EPGobj.mLgN(1:2:end)));
            set(a1,'YTickLabel',num2str(EPGobj.mLgN(end:-2:1)));
            set(a1,'xminortick','on','yminortick','on');
            set(a1,'box','on');
            xlabel(a1,'EUD'); ylabel(a1,'log_1_0(a)');
        end
        function fEUDCurvesSummary_a_EUD(EPGobj)
            % prepare
            f = [EPGobj.mGrp.mFlgCensor]; % censor info
            numcensor = sum(f); % number of patients censored
            numcensor95 = round((1-0.95)*numcensor/2);
            numcensor68 = round((1-0.68)*numcensor/2);
            numcomp = sum(~f); % number of patients with complication
            numcomp95 = round((1-0.95)*numcomp/2);
            numcomp68 = round((1-0.68)*numcomp/2);
            
            eud = [EPGobj.mGrp.EUD]'; % all doses
            eud_censor  = zeros(length(EPGobj.mLgN),5); % 5 columns for -95%, -68%, median, 68%, 95% lines
            eud_comp = eud_censor;

            % eud computation
            for kk = 1:length(EPGobj.mLgN)
                % eud of censored patients
                dose = sort(eud(f,kk)); % dose of censored patients
                if numcensor95>0
                    eud_censor(kk,1) = dose(numcensor95);
                    eud_censor(kk,5) = dose(numcensor-numcensor95);
                else
                    eud_censor(kk,1) = -inf;
                    eud_censor(kk,5) = -inf;
                end
                if numcensor68>0
                    eud_censor(kk,2) = dose(numcensor68);
                    eud_censor(kk,4) = dose(numcensor-numcensor68);
                else
                    eud_censor(kk,2) = -inf;
                    eud_censor(kk,4) = -inf;
                end
                eud_censor(kk,3) = median(dose);

                % eud of complication patients
                dose = sort(eud(~f,kk)); % dose of complication patients
                if numcomp95>0
                    eud_comp(kk,1) = dose(numcomp95);
                    eud_comp(kk,5) = dose(numcomp-numcomp95);
                else
                    eud_comp(kk,1) = -inf;
                    eud_comp(kk,5) = -inf;
                end
                if numcomp68>0
                    eud_comp(kk,2) = dose(numcomp68);
                    eud_comp(kk,4) = dose(numcomp-numcomp68);
                else
                    eud_comp(kk,2) = -inf;
                    eud_comp(kk,4) = -inf;
                end
                eud_comp(kk,3) = median(dose);
            end
            
            % plot the curves
            a1 = gca;
            a2 = copyobj(a1,gcf);
            set(a2,'Color','none');
            set(a2,'Xtick',[]);

            hold on;
            plot(eud_censor(:,3),EPGobj.mLgN,'b*-');
            plot(eud_censor(:,[2,4]),EPGobj.mLgN,'b');
            plot(eud_censor(:,[1,5]),EPGobj.mLgN,'b--');
            plot(eud_comp(:,3),EPGobj.mLgN,'r*-');
            plot(eud_comp(:,[2,4]),EPGobj.mLgN,'r');
            plot(eud_comp(:,[1,5]),EPGobj.mLgN,'r--');
            hold off;

            set(a2,'YAxisLocation','right');
            set(a2,'YTickLabel',num2str(EPGobj.mLgN(1:2:end)));
            set(a1,'YTickLabel',num2str(EPGobj.mLgN(end:-2:1)));
            set(a1,'xminortick','on','yminortick','on');
            set(a1,'box','on');
            xlabel(a1,'Dose (Gy)'); ylabel(a1,'log_1_0(a)');
        end
        function fAtlasFig_EUD(EPGobj)
            if isempty(EPGobj.mEUD_AtlasTotal)
                disp('Empty member "mEUD_AtlasTotal", cannot display its figure.'); return;
            end
            dosestep = 5;
            doses=EPGobj.mEUD_DoseBins; x=mod(doses,dosestep)==0;
            [xx,yy]=ndgrid(doses(x),1:length(EPGobj.mLgN)); xx=num2cell(xx); yy=num2cell(yy);
            strComp=(EPGobj.mEUD_AtlasComp(x,:)); strTotal=(EPGobj.mEUD_AtlasTotal(x,:)); strAtlas=arrayfun(@(a,b) strcat(num2str(a),'/',num2str(b)),strComp,strTotal,'UniformOutpu',false);
            cellfun(@(a,b,c) text(a,b,c,'fontsize',16),xx,yy,strAtlas);
            set(gca,'XLim',[xx{1,1}-1,xx{end,1}+2]); set(gca,'YLim',[1-1,length(EPGobj.mLgN)+1]);
            set(gca,'YTick',1:2:length(EPGobj.mLgN)); set(gca,'YTickLabel',EPGobj.mLgN(1:2:end));
            pos=get(gcf,'Position'); set(gcf,'Position',[pos(1:2), 900/3*4, 900]);
            set(gca,'xminortick','on','yminortick','on');
            set(gca,'box','on');
%             title([EPGobj.xlsSheet,', number of patients and numbers with severe pneumonitis treated to at least a given EUD'],'fontsize',16);
            xlabel('gEUD doses (Gy)'); ylabel('log_1_0a'); % set(gca,'fontsize',16);
        end
        function fAtlasRotatedFig_EUD(EPGobj,fonts,ticks)
            if ~exist('fonts','var')
                fonts = 10;
            end
            if ~exist('ticks','var')
                ticks = 10;
            end
            
            if isempty(EPGobj.mEUD_AtlasTotal)
                disp('Empty member "mEUD_AtlasTotal", cannot display its figure.'); return;
            end
            % check columns with informaiton
            % complication part
            f = diff(EPGobj.mEUD_AtlasComp);
            for k = 1:size(EPGobj.mEUD_AtlasComp,1) % upper side
                if any(f(k,:))
                    colcompu = k;
                    break;
                end
            end
            for k = size(EPGobj.mEUD_AtlasComp,1)-1 : -1 : 1 % lower side
                if any(f(k,:))
                    colcompl = k+1;
                    break;
                end
            end
            % total part
            f = diff(EPGobj.mEUD_AtlasTotal);
            for k = 1:size(EPGobj.mEUD_AtlasTotal,1) % upper side
                if any(f(k,:))
                    coltotalu = k;
                    break;
                end
            end
            for k = size(EPGobj.mEUD_AtlasTotal,1)-1 : -1 : 1 % right side
                if any(f(k,:))
                    coltotall = k+1;
                    break;
                end
            end
            % combination of the comp and total
            colu = min( colcompu, coltotalu );
            coll = max( colcompl, coltotall );
            % prepare to write the table in a figure
%             doses = zeros(size(EPGobj.mEUD_DoseBins));
            doses=EPGobj.mEUD_DoseBins(colu:coll);
            [xx,yy]=ndgrid(doses,1:length(EPGobj.mLgN)); xx=num2cell(xx); yy=num2cell(yy);
            strComp=EPGobj.mEUD_AtlasComp(colu:coll,:); strTotal=EPGobj.mEUD_AtlasTotal(colu:coll,:);
            strAtlas=arrayfun(@(a,b) strcat(num2str(a),'/',num2str(b)),strComp,strTotal,'UniformOutpu',false);
%             figure(1);
            clf reset;
            cellfun(@(a,b,c) text(b,a,c,'fontsize',fonts),xx,yy,strAtlas);
            set(gca,'YLim',[xx{1,1}-1,xx{end,1}+2]);
            set(gca,'XLim',[1-1,length(EPGobj.mLgN)+1]);
            set(gca,'XTick',1:2:length(EPGobj.mLgN)); set(gca,'XTickLabel',EPGobj.mLgN(1:2:end));
            set(gca,'fontsize',ticks);
            set(gca,'xminortick','on','yminortick','on');
            set(gca,'Position',[0.04,0.04,0.955,0.955]);
            set(gca,'box','on');
%             pos=get(gcf,'Position'); set(gcf,'Position',[pos(1:2), 900/3*4, 900]);
%             title([EPGobj.xlsSheet,', number of patients and numbers with severe pneumonitis treated to at least a given EUD'],'fontsize',16);
            xlabel('gEUD doses (Gy)'); ylabel('log_1_0a'); % set(gca,'fontsize',16);
        end
        function fAtlasCompactFig_EUD(EPGobj,fonts,ticks)
            if ~exist('fonts','var')
                fonts = 10;
            end
            if ~exist('ticks','var')
                ticks = 10;
            end
            
            if isempty(EPGobj.mEUD_AtlasTotal)
                disp('Empty member "mEUD_AtlasTotal", cannot display its figure.'); return;
            end
            
            % check duplicated columns for each n to determine the talbe size
            dupcol = logical(diff(EPGobj.mEUD_AtlasComp));
            dupcol = dupcol | logical(diff(EPGobj.mEUD_AtlasTotal));
            shiftpos = zeros(2,length(EPGobj.mLgN));
            for n = 1:length(EPGobj.mLgN)
                f = find(dupcol(:,n));
                shiftpos(1,n) = f(1);
                shiftpos(2,n) = f(end);
            end
            f = diff(shiftpos); fy = max(f)+1;

            % prepare data
            strComp = zeros(fy+1,length(EPGobj.mLgN));
            strTotal = zeros(fy+1,length(EPGobj.mLgN));
            for n = 1:length(EPGobj.mLgN)
                strComp(:,n) = EPGobj.mEUD_AtlasComp(shiftpos(1,n):shiftpos(1,n)+fy,n);
                strTotal(:,n) = EPGobj.mEUD_AtlasTotal(shiftpos(1,n):shiftpos(1,n)+fy,n);
            end
            strAtlas = arrayfun(@(a,b) strcat(num2str(a),'/',num2str(b)),strComp,strTotal,'UniformOutput',false);
%             strshift = num2cell(shiftpos(1,:))';
            strshift = arrayfun(@(a) num2str(EPGobj.mEUD_DoseBins(a)),shiftpos(1,:),'UniformOutput',false);

            % plot the table
            clf reset;
            [x,y] = ndgrid(0:fy,1:length(EPGobj.mLgN));
            xx = num2cell(x); yy = num2cell(y);
            cellfun(@(a,b,c) text(b,a,c,'fontsize',fonts),xx,yy,strAtlas);
            set(gca,'YLim',[xx{1,1}-1,xx{end,1}+2]);
            set(gca,'XLim',[1-1,length(EPGobj.mLgN)+1]);
            set(gca,'XTick',1:1:length(EPGobj.mLgN)); set(gca,'XTickLabel',EPGobj.mLgN(end:-1:1));
            set(gca,'fontsize',ticks);

%             axis manual;
            cellfun(@(b,x) text(b,-4,x,'fontsize',fonts),yy(1,:),strshift);
            set(gca,'Position',[0.04,0.1,0.956,0.89]);
            set(gca,'xminortick','off','yminortick','on');
            set(gca,'box','on');
            xlabel('log_1_0a'); ylabel('gEUD doses (Gy)'); % set(gca,'fontsize',16);
        end
        
        function loga = fLogitLikelyhoodExactFig_a_EUD(EPGobj,loga,strMarker,lw)
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
            st = [EPGobj.mLogitMat];
            dpf = [st.dev]; % deviations
            st =[st.stats];
            df = [st.dfe]; % degree of freedom
            dpf = dpf./df; % deviations per degree of freedom
            [~,loc] = min(dpf); % the minimum deviation corresponds to the maximum likelihood
            disp(['best log10(a) of Logistic Regression of exact gEUD is: ',num2str(-EPGobj.mLgN(loc))]);
            if ~exist('loga','var') || ischar(loga)
                loga = -EPGobj.mLgN(loc);
            end
            disp(['the log10(a) in coefficient searching is: ',num2str(loga)]);
            [~,loc] = min(abs(EPGobj.mLgN+loga));
            disp('the corresponding coefficients, sd, and 95% CI are:');
            disp(num2str([st(loc).beta, st(loc).se, st(loc).beta-1.96*st(loc).se, st(loc).beta+1.96*st(loc).se]));

            loglikelyhood = -0.5*dpf;
            [mx,loc] = max(loglikelyhood); % the maximum loglikelyhood
            loglikelyhood68 = repmat(mx-0.5* 1 /df(loc),size(EPGobj.mLgN));
            loglikelyhood95 = repmat(mx-0.5* (1.96*2) /df(loc),size(EPGobj.mLgN));
            hold on;
            plot(-EPGobj.mLgN, loglikelyhood,strMarker,'LineWidth',lw);
            plot(-EPGobj.mLgN, loglikelyhood68,strcat(strMarker(1),'--'),'LineWidth',1);
            plot(-EPGobj.mLgN, loglikelyhood95,strcat(strMarker(1),'-'),'LineWidth',1);
            plot(-EPGobj.mLgN(loc), loglikelyhood(loc),strMarker,'LineWidth',lw+2);
            set(gca,'xminortick','on','yminortick','on');
            set(gca,'box','on');
            xlabel('log_1_0a'); ylabel('loglikely per degree of freedom');

%             % disp the result from atlas
%             if ~isempty(EPGobj.mLogitMatBin)
%                 st = [EPGobj.mLogitMatBin];
%                 disp(['the coefficients for the atlas are: ',num2str(st(loc).stats.beta')]);
%             end
        end
        function fLogitLikelyhoodBinFig_a_EUD(EPGobj,loga,strMarker,lw)
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
            st = [EPGobj.mLogitMatBin];
            dpf = [st.dev]; % deviations
            st =[st.stats];
            df = [st.dfe]; % degree of freedom
            dpf = dpf./df; % deviations per degree of freedom
            [~,loc] = min(dpf); % the minimum deviation corresponds to the maximum likelihood
            disp(['best log10(a) of Logistic Regression of gEUD Atlas is: ',num2str(-EPGobj.mLgN(loc))]);
            if ~exist('loga','var') || ischar(loga)
                loga = -EPGobj.mLgN(loc);
            end
            disp(['the log10(a) in coefficient searching is: ',num2str(loga)]);
            [~,loc] = min(abs(EPGobj.mLgN+loga));
            disp('the corresponding coefficients, sd, and 95% CI are:');
            disp(num2str([st(loc).beta, st(loc).se, st(loc).beta-1.96*st(loc).se, st(loc).beta+1.96*st(loc).se]));

            loglikelyhood = -0.5*dpf;
            [mx,loc] = max(loglikelyhood); % the maximum loglikelyhood
            loglikelyhood68 = repmat(mx-0.5* 1/df(loc),size(EPGobj.mLgN));
            loglikelyhood95 = repmat(mx-0.5* (1.96*2) /df(loc),size(EPGobj.mLgN));
            hold on;
            plot(-EPGobj.mLgN, loglikelyhood,strMarker,'LineWidth',lw);
%             plot(-EPGobj.mLgN, loglikelyhood68,strcat(strMarker(1),'--'),'LineWidth',1);
%             plot(-EPGobj.mLgN, loglikelyhood95,strcat(strMarker(1),'-'),'LineWidth',1);
%             plot(-EPGobj.mLgN(loc), loglikelyhood(loc),strMarker,'LineWidth',lw+2);
            set(gca,'xminortick','on','yminortick','on');
            set(gca,'box','on');
            xlabel('log_1_0a'); ylabel('loglikely per degree of freedom');
        end
        function fLogitPvalueExactFig_a_EUD(EPGobj,strMarker,lw)
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
            st = [EPGobj.mLogitMat];
            st =[st.stats];
            pvalue = [st.p];
            pvalue = pvalue(2,:); % the p-value corresponding to gEUD
            [~,loc] = min(pvalue); % the location of mininum p-value
            semilogy(-EPGobj.mLgN, pvalue,strMarker,'LineWidth',lw);
            hold on;
            semilogy(-EPGobj.mLgN(loc), pvalue(loc),strMarker,'LineWidth',lw+2);
            semilogy(-EPGobj.mLgN, repmat(0.05,size(EPGobj.mLgN)), strcat(strMarker(1),'--'),'LineWidth',1);
            set(gca,'xminortick','on','yminortick','on');
            set(gca,'box','on');
            xlabel('log_1_0a'); ylabel('p-value');
        end
        function fLogitPvalueBinFig_a_EUD(EPGobj,strMarker,lw)
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
            st = [EPGobj.mLogitMatBin];
            st =[st.stats];
            pvalue = [st.p];
            pvalue = pvalue(2,:); % the p-value corresponding to gEUD
            semilogy(-EPGobj.mLgN, pvalue,strMarker,'LineWidth',lw);
            hold on;
            semilogy(-EPGobj.mLgN, repmat(0.05,size(EPGobj.mLgN)), strcat(strMarker(1),'--'),'LineWidth',1);
            set(gca,'xminortick','on','yminortick','on');
            set(gca,'box','on');
            xlabel('log_1_0a'); ylabel('p-value');
        end
        function loga = fLogitRespondingCurveExactFig_a_EUD(EPGobj,loga,strMarker,lw)
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
            st = [EPGobj.mLogitMat];
            dpf = [st.dev]; % deviations
            [~,loc] = min(dpf); % the minimum deviation corresponds to the maximum likelihood
            disp(['best log10(a) of Logistic Regression of exact gEUD is: ',num2str(-EPGobj.mLgN(loc))]);
            if ~exist('loga','var') || ischar(loga)
                loga = -EPGobj.mLgN(loc);
            end
            [~,loc] = min(abs(EPGobj.mLgN+loga)); % the n whose corresponding responding function will be ploted
            disp(['the log10(a) in responding curve is: ',num2str(loga)]);

            % responding curve
            st = EPGobj.mLogitMat(loc); % the fitting result of that n
            euds = [EPGobj.mGrp.EUD]; euds = euds(loc,:); % the gEUDs of that n
            doses = (0:max(euds))'; % doses (gEUD) whose RP probability will be computed
            [rpb,rplo,rphi] = glmval(st.b, doses,'logit',st.stats); % the responding function values at doses
            disp(['the beta values are: ',num2str(st.b')]);

            % plot
            hold on;
            plot(doses,rpb,strMarker,'LineWidth',lw*2); % responding function
            plot(doses,rpb-rplo,strMarker,'LineWidth',lw); % low CI curve
            plot(doses,rpb+rphi,strMarker,'LineWidth',lw); % high CI curve
            set(gca,'xminortick','on','yminortick','on');
            set(gca,'box','on');
            xlabel('gEUD'); ylabel('RP probability');
        end
        function loga = fLogitRespondingCurveBinFig_a_EUD(EPGobj,loga,strMarker,lw)
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
            st = [EPGobj.mLogitMatBin];
            dpf = [st.dev]; % deviations
            [~,loc] = min(dpf); % the minimum deviation corresponds to the maximum likelihood
            disp(['best log10(a) of Logistic Regression of gEUD Atlas is: ',num2str(-EPGobj.mLgN(loc))]);
            if ~exist('loga','var') || ischar(loga)
                loga = -EPGobj.mLgN(loc);
            end
            disp(['the log10(a) in responding curve is: ',num2str(loga)]);

            % responding curve
            [~,loc] = min(abs(EPGobj.mLgN+loga)); % the n whose corresponding responding function will be ploted
            st = EPGobj.mLogitMatBin(loc); % the fitting result of that n
            euds = [EPGobj.mGrp.EUD]; euds = euds(loc,:); % the gEUDs of that n
            doses = (0:max(euds))'; % doses (gEUD) whose RP probability will be computed
            [rpb,rplo,rphi] = glmval(st.b, doses,'logit',st.stats); % the responding function values at doses
            disp(['the beta values are: ',num2str(st.b')]);

            % plot
            hold on;
            plot(doses,rpb,strMarker,'LineWidth',lw*2); % responding function
            plot(doses,rpb-rplo,strMarker,'LineWidth',lw); % low CI curve
            plot(doses,rpb+rphi,strMarker,'LineWidth',lw); % high CI curve
            set(gca,'xminortick','on','yminortick','on');
            set(gca,'box','on');
            xlabel('gEUD'); ylabel('RP probability');
        end

        function fLogitGOFFig(EPGobj,loga)
            % compute the p-value curves
            nmn = 5; nmx = 30;
            nGrp = [nmn:nmx, EPGobj.mNumInGrp]';
            nGrp = [nmn:nmx]';
            pValExact = zeros(length(nGrp),3);
            pValBin = zeros(length(nGrp),3);
            for n = 1:length(nGrp)
                % assign number of groups
                EPGobj.mLogitHosmerLemeshow.n = nGrp(n);
                EPGobj.mLogitHosmerLemeshowBin.n = nGrp(n);
                EPGobj.mLogitGTest.n = nGrp(n);
                EPGobj.mLogitGTestBin.n = nGrp(n);
                EPGobj.mLogitPearson.n = nGrp(n);
                EPGobj.mLogitPearsonBin.n = nGrp(n);

                % compute the goodness of fit
                EPGobj = EPGobj.LogitHosmerLemeshowTestAnalysisExact_EUD(loga);
                EPGobj = EPGobj.LogitHosmerLemeshowTestAnalysisBin_EUD(loga);
                EPGobj = EPGobj.LogitGTestAnalysisExact_EUD(loga);
                EPGobj = EPGobj.LogitGTestAnalysisBin_EUD(loga);
                EPGobj = EPGobj.LogitPearsonTestAnalysisExact_EUD(loga);
                EPGobj = EPGobj.LogitPearsonTestAnalysisBin_EUD(loga);

                % save
                pValExact(n,1) = EPGobj.mLogitHosmerLemeshow.p_value;
                pValExact(n,2) = EPGobj.mLogitGTest.p_value;
                pValExact(n,3) = EPGobj.mLogitPearson.p_value;
                pValBin(n,1) = EPGobj.mLogitHosmerLemeshowBin.p_value;
                pValBin(n,2) = EPGobj.mLogitGTestBin.p_value;
                pValBin(n,3) = EPGobj.mLogitPearsonBin.p_value;
            end

            % plot p-value curves
            hold on;
            plot(nGrp,pValExact(:,1),'b--',nGrp,pValBin(:,1),'b-');
            plot(nGrp,pValExact(:,2),'r--',nGrp,pValBin(:,2),'r-');
            plot(nGrp,pValExact(:,3),'k--',nGrp,pValBin(:,3),'k-');
            set(gca,'xminortick','on','yminortick','on');
            set(gca,'box','on');
            xlabel('number of groups'); ylabel('p-value of goodness of fit');
        end

        function loga = fLogitGridRespondingCurveExactFig_a_EUD(EPGobj,loga,strMarker,lw)
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
            [mx,loc] = max(EPGobj.mLogitGrid.loglikelihood(:));
            [~,~,loc] = ind2sub(size(EPGobj.mLogitGrid.loglikelihood),loc);
            disp(['the maximum log likelihood and its "log10(a)" in Logit of exact gEUD are: ',num2str([mx, -EPGobj.mLgN(loc)])]);
            if ~exist('loga','var') || ischar(loga)
                loga = -EPGobj.mLgN(loc);
            end
            disp(['the plotted response curve is at log10(a) = ',num2str(loga)]);

            % coefficients for the Lyman model
            [~,n] = min(abs(EPGobj.mLgN+loga)); % the n whose corresponding responding function will be ploted
            ll = EPGobj.mLogitGrid.loglikelihood(:,:,n); % log likelihood of log10(a) = loga
            mx = max(ll(:));
            [xx,yy] = find(ll == mx); % the coefficients
            b0 = EPGobj.mLogitGrid.b0(xx);
            b1 = EPGobj.mLogitGrid.b1(yy);
            disp(['the b0 and b1 for the response function are: ',num2str([b0,b1])]);

            % curves
            euds = [EPGobj.mGrp.EUD]; euds = euds(n,:); % the gEUDs of that n
            doses = (0:max(euds))'; % doses (gEUD) whose RP probability will be computed
            rpb = exp(b0+b1*doses);
            rpb = rpb./(1+rpb); % logistic probability
            % 95% CI
            low95 = mx-0.5* (1.96*2); % CI level
            cxy = contourc(EPGobj.mLogitGrid.b0,EPGobj.mLogitGrid.b1,ll',[low95,low95]); % coordinates of CI contours at CI level
            CImx = zeros(length(doses),1); % upper end of CI for response at each gEUD
            CImn = zeros(length(doses),1); % lower end of CI for response at each gEUD
            b0 = cxy(1,2:end);
            b1 = cxy(2,2:end); 
            for d = 1:length(doses)
                CI = exp(b0+b1*doses(d));
                CI = CI./(1+CI); % probability at each iso point
                CImx(d) = max(CI);
                CImn(d) = min(CI);
            end
            % plot
            hold on;
            plot(doses,rpb,strMarker,'LineWidth',lw*2); % responding function
            plot(doses,CImx,strMarker,'LineWidth',lw); % responding function
            plot(doses,CImn,strMarker,'LineWidth',lw); % responding function
            set(gca,'xminortick','on','yminortick','on');
            set(gca,'box','on');
            xlabel('gEUD'); ylabel('RP probability');
        end
        function loga = fLogitGridRespondingCurveBinFig_a_EUD(EPGobj,loga,strMarker,lw)
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
            [mx,loc] = max(EPGobj.mLogitGridBin.loglikelihood(:));
            [~,~,loc] = ind2sub(size(EPGobj.mLogitGridBin.loglikelihood),loc);
            disp(['the maximum log likelihood and its "log10(a)" in Logit of Bin gEUD are: ',num2str([mx, -EPGobj.mLgN(loc)])]);
            if ~exist('loga','var') || ischar(loga)
                loga = -EPGobj.mLgN(loc);
            end
            disp(['the plotted response curve is at log10(a) = ',num2str(loga)]);

            % coefficients for the Lyman model
            [~,n] = min(abs(EPGobj.mLgN+loga)); % the n whose corresponding responding function will be ploted
            ll = EPGobj.mLogitGridBin.loglikelihood(:,:,n); % log likelihood of log10(a) = loga
            mx = max(ll(:));
            [xx,yy] = find(ll == mx); % the coefficients
            b0 = EPGobj.mLogitGridBin.b0(xx);
            b1 = EPGobj.mLogitGridBin.b1(yy);
            disp(['the b0 and b1 for the response function are: ',num2str([b0,b1])]);

            % curves
            euds = [EPGobj.mGrp.EUD]; euds = euds(n,:); % the gEUDs of that n
            doses = (0:max(euds))'; % doses (gEUD) whose RP probability will be computed
            rpb = exp(b0+b1*doses);
            rpb = rpb./(1+rpb); % logistic probability
            % 95% CI
            low95 = mx- 0.5* (1.96*2); % CI level
            cxy = contourc(EPGobj.mLogitGridBin.b0,EPGobj.mLogitGridBin.b1,ll',[low95,low95]); % coordinates of CI contours at CI level
            CImx = zeros(length(doses),1); % upper end of CI for response at each gEUD
            CImn = zeros(length(doses),1); % lower end of CI for response at each gEUD
            b0 = cxy(1,2:end);
            b1 = cxy(2,2:end); 
            for d = 1:length(doses)
                CI = exp(b0+b1*doses(d));
                CI = CI./(1+CI); % probability at each iso point
                CImx(d) = max(CI);
                CImn(d) = min(CI);
            end
            % plot
            hold on;
            plot(doses,rpb,strMarker,'LineWidth',lw*2); % responding function
%             plot(doses,CImx,strMarker,'LineWidth',lw); % responding function
%             plot(doses,CImn,strMarker,'LineWidth',lw); % responding function
            set(gca,'xminortick','on','yminortick','on');
            set(gca,'box','on');
            xlabel('gEUD'); ylabel('RP probability');
        end

        function fLogitGridExactFig_a_loglikelhood_EUD(EPGobj,strMarker,lw)
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
            [mx,loc] = max(EPGobj.mLogitGrid.loglikelihood(:));
            [~,~,loc] = ind2sub(size(EPGobj.mLogitGrid.loglikelihood),loc);
            disp(['the maximum log likelihood and its "log10(a)" in Logistic Regression of exact gEUD are: ',num2str([mx, -EPGobj.mLgN(loc)])]);

            % plot the maximum log likelihood curve w.r.t. log10(a)
            llmx = zeros(size(EPGobj.mLgN)); %log likelihood's max for each log10(a)
            for kk = 1:length(llmx)
                ll = EPGobj.mLogitGrid.loglikelihood(:,:,kk);
                llmx(kk) = max(ll(:));
            end
            llmx = llmx / (EPGobj.mNumInGrp-2);
            mx = mx / (EPGobj.mNumInGrp-2);
            loglikelyhood68 = repmat(mx-0.5* 1 /(EPGobj.mNumInGrp-2),size(llmx));
            loglikelyhood95 = repmat(mx-0.5* (1.96*2) /(EPGobj.mNumInGrp-2),size(llmx));
            hold on;
            plot(-EPGobj.mLgN,llmx,strMarker,'LineWidth',lw);
            plot(-EPGobj.mLgN, loglikelyhood68,strcat(strMarker(1),'--'),'LineWidth',1);
            plot(-EPGobj.mLgN, loglikelyhood95,strcat(strMarker(1),'-'),'LineWidth',1);
            if length(strMarker)==1
                strMarker = [strMarker,'*'];
            elseif strMarker(2) == '-' || strMarker(2) == '.'
                strMarker(2) = '*';
            end
            plot(-EPGobj.mLgN(loc),mx,strMarker,'LineWidth',lw+2);
            set(gca,'xminortick','on','yminortick','on');
            set(gca,'box','on');
            xlabel('log_1_0(a)'); ylabel('log likelihood per degree of freedom');
        end
        function fLogitGridBinFig_a_loglikelihood_EUD(EPGobj,strMarker,lw)
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
            [mx,loc] = max(EPGobj.mLogitGridBin.loglikelihood(:));
            [~,~,loc] = ind2sub(size(EPGobj.mLogitGridBin.loglikelihood),loc);
            disp(['the maximum log likelihood and its "log10(a)" in Logistic Regression of gEUD Atlas are: ',num2str([mx, -EPGobj.mLgN(loc)])]);

            % plot the maximum log likelihood curve w.r.t. log10(a)
            llmx = zeros(size(EPGobj.mLgN)); %log likelihood's max for each log10(a)
            for kk = 1:length(llmx)
                ll = EPGobj.mLogitGridBin.loglikelihood(:,:,kk);
                llmx(kk) = max(ll(:));
            end
            llmx = llmx / (EPGobj.mNumInGrp-2);
            mx = mx / (EPGobj.mNumInGrp-2);
            loglikelyhood68 = repmat(mx-0.5* 1 /(EPGobj.mNumInGrp-2),size(llmx));
            loglikelyhood95 = repmat(mx-0.5* (1.96*2) /(EPGobj.mNumInGrp-2),size(llmx));
            hold on;
            plot(-EPGobj.mLgN,llmx,strMarker,'LineWidth',lw);
            plot(-EPGobj.mLgN(loc),mx,strMarker,'LineWidth',lw+2);
%             plot(-EPGobj.mLgN, loglikelyhood68,strcat(strMarker(1),'--'),'LineWidth',1);
%             plot(-EPGobj.mLgN, loglikelyhood95,strcat(strMarker(1),'-'),'LineWidth',1);
            set(gca,'xminortick','on','yminortick','on');
            set(gca,'box','on');
            xlabel('log_1_0(a)'); ylabel('log likelihood per degree of freedom');
        end
        function fLogitGridExactFig_b0_loglikelihood_EUD(EPGobj,strMarker,lw)
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

            % find the best "b0"
            [mx,loc] = max(EPGobj.mLogitGrid.loglikelihood(:));
            [loc,~,~] = ind2sub(size(EPGobj.mLogitGrid.loglikelihood),loc);
            disp(['the maximum log likelihood and its b0 in Logistic Regression of exact gEUD are: ',num2str([mx, EPGobj.mLogitGrid.b0(loc)])]);

            % plot the maximum log likelihood curve w.r.t. TD50
            llmx = zeros(size(EPGobj.mLogitGrid.b0)); %log likelihood's max for each log10(a)
            for kk = 1:length(llmx)
                ll = EPGobj.mLogitGrid.loglikelihood(kk,:,:);
                llmx(kk) = max(ll(:));
            end
            llmx = llmx / (EPGobj.mNumInGrp-2);
            mx = mx / (EPGobj.mNumInGrp-2);
            loglikelyhood68 = repmat(mx-0.5* 1 /(EPGobj.mNumInGrp-2),size(llmx));
            loglikelyhood95 = repmat(mx-0.5* (1.96*2) /(EPGobj.mNumInGrp-2),size(llmx));
            hold on;
            plot(EPGobj.mLogitGrid.b0,llmx,strMarker,'LineWidth',lw);
            plot(EPGobj.mLogitGrid.b0, loglikelyhood68,strcat(strMarker(1),'--'),'LineWidth',1);
            plot(EPGobj.mLogitGrid.b0, loglikelyhood95,strcat(strMarker(1),'-'),'LineWidth',1);
            if length(strMarker)==1
                strMarker = [strMarker,'*'];
            elseif strMarker(2) == '-' || strMarker(2) == '.'
                strMarker(2) = '*';
            end
            plot(EPGobj.mLogitGrid.b0(loc),mx,strMarker,'LineWidth',lw+2);
            hold off;
            set(gca,'xminortick','on','yminortick','on');
            set(gca,'box','on');
            xlabel('b0'); ylabel('log likelihood per degree of freedom');
        end
        function fLogitGridBinFig_b0_loglikelihood_EUD(EPGobj,strMarker,lw)
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

            % find the best "b0"
            [mx,loc] = max(EPGobj.mLogitGridBin.loglikelihood(:));
            [loc,~,~] = ind2sub(size(EPGobj.mLogitGridBin.loglikelihood),loc);
            disp(['the maximum log likelihood and its b0 in Logistic Regression of gEUD Atlas are: ',num2str([mx, EPGobj.mLogitGridBin.b0(loc)])]);

            % plot the maximum log likelihood curve w.r.t. TD50
            llmx = zeros(size(EPGobj.mLogitGridBin.b0)); %log likelihood's max for each log10(a)
            for kk = 1:length(llmx)
                ll = EPGobj.mLogitGridBin.loglikelihood(kk,:,:);
                llmx(kk) = max(ll(:));
            end
            llmx = llmx / (EPGobj.mNumInGrp-2);
            mx = mx / (EPGobj.mNumInGrp-2);
            loglikelyhood68 = repmat(mx-0.5* 1 /(EPGobj.mNumInGrp-2),size(llmx));
            loglikelyhood95 = repmat(mx-0.5* (1.96*2) /(EPGobj.mNumInGrp-2),size(llmx));
            hold on;
            plot(EPGobj.mLogitGridBin.b0,llmx,strMarker,'LineWidth',lw);
%             plot(EPGobj.mLogitGridBin.b0, loglikelyhood68,strcat(strMarker(1),'--'),'LineWidth',1);
%             plot(EPGobj.mLogitGridBin.b0, loglikelyhood95,strcat(strMarker(1),'-'),'LineWidth',1);
            if length(strMarker)==1
                strMarker = [strMarker,'*'];
            elseif strMarker(2) == '-' || strMarker(2) == '.'
                strMarker(2) = '*';
            end
            plot(EPGobj.mLogitGridBin.b0(loc),mx,strMarker,'LineWidth',lw+2);
            hold off;
            set(gca,'xminortick','on','yminortick','on');
            set(gca,'box','on');
            xlabel('b0'); ylabel('log likelihood per degree of freedom');
        end
        function fLogitGridExactFig_b1_loglikelihood_EUD(EPGobj,strMarker,lw)
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
            [mx,loc] = max(EPGobj.mLogitGrid.loglikelihood(:));
            [~,loc,~] = ind2sub(size(EPGobj.mLogitGrid.loglikelihood),loc);
            disp(['the maximum log likelihood and its "b1" in Logistic Regression of exact gEUD are: ',num2str([mx, EPGobj.mLogitGrid.b1(loc)])]);

            % plot the maximum log likelihood curve w.r.t. TD50
            llmx = zeros(size(EPGobj.mLogitGrid.b1)); %log likelihood's max for each log10(a)
            for kk = 1:length(llmx)
                ll = EPGobj.mLogitGrid.loglikelihood(:,kk,:);
                llmx(kk) = max(ll(:));
            end
            llmx = llmx / (EPGobj.mNumInGrp-2);
            mx = mx / (EPGobj.mNumInGrp-2);
            loglikelyhood68 = repmat(mx-0.5* 1 /(EPGobj.mNumInGrp-2),size(llmx));
            loglikelyhood95 = repmat(mx-0.5* (1.96*2) /(EPGobj.mNumInGrp-2),size(llmx));
            hold on;
            plot(EPGobj.mLogitGrid.b1,llmx,strMarker,'LineWidth',lw);
            plot(EPGobj.mLogitGrid.b1, loglikelyhood68,strcat(strMarker(1),'--'),'LineWidth',1);
            plot(EPGobj.mLogitGrid.b1, loglikelyhood95,strcat(strMarker(1),'-'),'LineWidth',1);
            if length(strMarker)==1
                strMarker = [strMarker,'*'];
            elseif strMarker(2) == '-' || strMarker(2) == '.'
                strMarker(2) = '*';
            end
            plot(EPGobj.mLogitGrid.b1(loc),mx,strMarker,'LineWidth',lw+2);
            hold off;
            set(gca,'xminortick','on','yminortick','on');
            set(gca,'box','on');
            xlabel('b1'); ylabel('log likelihood per degree of freedom');
        end
        function fLogitGridBinFig_b1_loglikelihood_EUD(EPGobj,strMarker,lw)
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
            [mx,loc] = max(EPGobj.mLogitGridBin.loglikelihood(:));
            [~,loc,~] = ind2sub(size(EPGobj.mLogitGridBin.loglikelihood),loc);
            disp(['the maximum log likelihood and its "b1" in Logistic Regression of gEUD Atlas are: ',num2str([mx, EPGobj.mLogitGridBin.b1(loc)])]);

            % plot the maximum log likelihood curve w.r.t. TD50
            llmx = zeros(size(EPGobj.mLogitGridBin.b1)); %log likelihood's max for each log10(a)
            for kk = 1:length(llmx)
                ll = EPGobj.mLogitGridBin.loglikelihood(:,kk,:);
                llmx(kk) = max(ll(:));
            end
            llmx = llmx / (EPGobj.mNumInGrp-2);
            mx = mx / (EPGobj.mNumInGrp-2);
            loglikelyhood68 = repmat(mx-0.5* 1 /(EPGobj.mNumInGrp-2),size(llmx));
            loglikelyhood95 = repmat(mx-0.5* (1.96*2) /(EPGobj.mNumInGrp-2),size(llmx));
            hold on;
            plot(EPGobj.mLogitGridBin.b1,llmx,strMarker,'LineWidth',lw);
%             plot(EPGobj.mLogitGridBin.b1, loglikelyhood68,strcat(strMarker(1),'--'),'LineWidth',1);
%             plot(EPGobj.mLogitGridBin.b1, loglikelyhood95,strcat(strMarker(1),'-'),'LineWidth',1);
            if length(strMarker)==1
                strMarker = [strMarker,'*'];
            elseif strMarker(2) == '-' || strMarker(2) == '.'
                strMarker(2) = '*';
            end
            plot(EPGobj.mLogitGridBin.b1(loc),mx,strMarker,'LineWidth',lw+2);
            hold off;
            set(gca,'xminortick','on','yminortick','on');
            set(gca,'box','on');
            xlabel('b1'); ylabel('log likelihood per degree of freedom');
        end

        function fLogitGridExactFig_b0_b1_EUD(EPGobj,loga,zdepth)
            if ~exist('zdepth','var')
                zdepth = 10;
            end
            
            % find the best "a"
            [mx,loc] = max(EPGobj.mLogitGrid.loglikelihood(:));
            [~,~,loc] = ind2sub(size(EPGobj.mLogitGrid.loglikelihood),loc);
            disp(['the "log10(a)" for the maximum likelihood in Logistic Regression of exact gEUD is: ',num2str(-EPGobj.mLgN(loc))]);
            if ~exist('loga','var') || ischar(loga)
                loga = -EPGobj.mLgN(loc);
            end
            [~,loc] = min(abs(EPGobj.mLgN+loga)); % the n whose corresponding responding function will be ploted
            disp(['the "log10(a)" for the map of likelihood in Logistic Regression is: ',num2str(-EPGobj.mLgN(loc))]);

            
            % b0 and b1 of the best "log10(a)"
            ll = EPGobj.mLogitGrid.loglikelihood(:,:,loc);
            mx = max(ll(:));
            [dd,mm] = find(ll == mx,1);
            disp(['b0 & b1 are: ',num2str([EPGobj.mLogitGrid.b0(dd),EPGobj.mLogitGrid.b1(mm)])]);
            
            % map of the grid of the best "a"
            if 0 % show the mesh surface
                ll( ll < mx-zdepth ) = -inf;
                surf(EPGobj.mLogitGrid.b1,EPGobj.mLogitGrid.b0,ll,'EdgeColor','none');
                ylabel('b0'); xlabel('b1 (Gy^-^1)'); zlabel('Log likelihood');
            else % show the contours
                % compute the 68% and 95% CIs
                low68 = mx-0.5* (1);
                low95 = mx-0.5* (1.96^2);
                low99 = mx-0.5* (3^2);
%                 low68 = mx-0.5* (2.30);
%                 low95 = mx-0.5* (6.17);
%                 low99 = mx-0.5* (9.21);
                % display the contours
                [c,h] = contour(EPGobj.mLogitGrid.b0,EPGobj.mLogitGrid.b1,ll',[low99,low95,low68]);
                hold on; plot(EPGobj.mLogitGrid.b0(dd),EPGobj.mLogitGrid.b1(mm),'k*'); hold off;
%                 text_handle = clabel(c,h,'LabelSpacing',100000);%'String',{'low 68%';'low 95%';'low 99.7%'},
%                 set(text_handle,'BackgroundColor',[1 1 .6], 'Edgecolor',[.7 .7 .7], 'String',{'low 68%';'low 95%';'low 99.7%'});
                set(gca,'xminortick','on','yminortick','on');
                set(gca,'box','on');
                xlabel('b0'); ylabel('b1 (Gy^-^1)');
            end

%             % the coresponding parameters computed from atlas
%             if ~isempty(EPGobj.mLogitGridBin)
%                 ll = EPGobj.mLogitGridBin.loglikelihood(:,:,loc); % log likelihood of log10(a) = loga
%                 mx = max(ll(:)); % the best point for this log10(a)
%                 [dd,mm] = find(ll == mx); % the coefficients
%                 disp(['b0 & b1 for the same log10(a) in the atlas are: ', num2str([EPGobj.mLogitGridBin.b0(dd),EPGobj.mLogitGridBin.b1(mm)])]);
%             end
        end
        function fLogitGridBinFig_b0_b1_EUD(EPGobj,loga,zdepth)
            if ~exist('zdepth','var')
                zdepth = 10;
            end
            
            % find the best "a"
            [mx,loc] = max(EPGobj.mLogitGridBin.loglikelihood(:));
            [~,~,loc] = ind2sub(size(EPGobj.mLogitGridBin.loglikelihood),loc);
            disp(['the "log10(a)" for the maximum likelihood in Logistic Regression of gEUD Atlas is: ',num2str(-EPGobj.mLgN(loc))]);
            if ~exist('loga','var') || ischar(loga)
                loga = -EPGobj.mLgN(loc);
            end
            [~,loc] = min(abs(EPGobj.mLgN+loga)); % the n whose corresponding responding function will be ploted
            disp(['the "log10(a)" for the map of likelihood in Logistic Regression is: ',num2str(-EPGobj.mLgN(loc))]);

            
            % b0 and b1 of the best "log10(a)"
            ll = EPGobj.mLogitGridBin.loglikelihood(:,:,loc);
            mx = max(ll(:));
            [dd,mm] = find(ll == mx,1);
            disp(['b0 & b1 are: ',num2str([EPGobj.mLogitGridBin.b0(dd),EPGobj.mLogitGridBin.b1(mm)])]);
            
            % map of the grid of the best "a"
            if 0 % show the mesh surface
                ll( ll < mx-zdepth ) = -inf;
                surf(EPGobj.mLogitGridBin.b1,EPGobj.mLogitGridBin.b0,ll,'EdgeColor','none');
                ylabel('b0'); xlabel('b1 (Gy^-^1)'); zlabel('Log likelihood');
            else % show the contours
                % compute the 68% and 95% CIs
                low68 = mx-0.5* (1);
                low95 = mx-0.5* (1.96^2);
                low99 = mx-0.5* (3^2);
%                 low68 = mx-0.5* (2.30);
%                 low95 = mx-0.5* (6.17);
%                 low99 = mx-0.5* (9.21);
                % display the contours
                [c,h] = contour(EPGobj.mLogitGridBin.b0,EPGobj.mLogitGridBin.b1,ll',[low99,low95,low68]);
                hold on; plot(EPGobj.mLogitGridBin.b0(dd),EPGobj.mLogitGridBin.b1(mm),'k*'); hold off;
%                 text_handle = clabel(c,h,'LabelSpacing',100000);%'String',{'low 68%';'low 95%';'low 99.7%'},
%                 set(text_handle,'BackgroundColor',[1 1 .6], 'Edgecolor',[.7 .7 .7], 'String',{'low 68%';'low 95%';'low 99.7%'});
                set(gca,'xminortick','on','yminortick','on');
                set(gca,'box','on');
                xlabel('b0'); ylabel('b1 (Gy^-^1)');
            end
        end
        function fLogitGridExactFig_b0_a_EUD(EPGobj,zdepth)
            if ~exist('zdepth','var')
                zdepth = 10;
            end
            
            % find the best "m"
            [mx,loc] = max(EPGobj.mLogitGrid.loglikelihood(:));
            [~,loc,~] = ind2sub(size(EPGobj.mLogitGrid.loglikelihood),loc);
            disp(['the "b1" for the maximum likelihood in Logistic Regression of exact gEUD is: ',num2str(EPGobj.mLogitGrid.b1(loc))]);
            
            % TD50 and log10(a) corresponding to the best m
            ll = squeeze(EPGobj.mLogitGrid.loglikelihood(:,loc,:));
            [dd,aa] = find(ll == mx,1);
            disp(['b0 & log10(a) are: ',num2str([EPGobj.mLogitGrid.b0(dd),-EPGobj.mLgN(aa)])]);
            
            % map of the grid of the best "m"
            if 0 % show the mesh surface
                ll( ll < mx-zdepth ) = -inf;
                surf(-EPGobj.mLgN,EPGobj.mLogitGrid.b0,ll,'EdgeColor','none');
                ylabel('b0'); xlabel('log_1_0(a)'); zlabel('Log likelihood');
            else % show the contours
                % compute the 68% and 95% CIs
                low68 = mx-0.5* (1);
                low95 = mx-0.5* (1.96^2);
                low99 = mx-0.5* (3^2);
%                 low68 = mx-0.5* (2.30);
%                 low95 = mx-0.5* (6.17);
%                 low99 = mx-0.5* (9.21);
                % display the contours
                contour(-EPGobj.mLgN,EPGobj.mLogitGrid.b0,ll,[low99,low95,low68]);
                hold on; plot(-EPGobj.mLgN(aa),EPGobj.mLogitGrid.b0(dd),'k*'); hold off;
                set(gca,'xminortick','on','yminortick','on');
                set(gca,'box','on');
                xlabel('log_1_0(a)'); ylabel('b0');
            end
        end
        function fLogitGridBinFig_b0_a_EUD(EPGobj,zdepth)
            if ~exist('zdepth','var')
                zdepth = 10;
            end
            
            % find the best "m"
            [mx,loc] = max(EPGobj.mLogitGridBin.loglikelihood(:));
            [~,loc,~] = ind2sub(size(EPGobj.mLogitGridBin.loglikelihood),loc);
            disp(['the "b1" for the maximum likelihood in Logistic Regression of gEUD Atlas is: ',num2str(EPGobj.mLogitGridBin.b1(loc))]);
            
%             % reconstruct the log likelihood matrix using image reconstruction
%             ll = -inf(size(EPGobj.mLogitGridBin.loglikelihood));
%             [mx,loc1] = max(EPGobj.mLogitGridBin.loglikelihood(:));
%             ll(loc1) = mx;
%             EPGobj.mLogitGridBin.loglikelihood = imreconstruct(ll,EPGobj.mLogitGridBin.loglikelihood,26);
            
            % TD50 and log10(a) corresponding to the best m
            ll = squeeze(EPGobj.mLogitGridBin.loglikelihood(:,loc,:));
            [dd,aa] = find(ll == mx,1);
            disp(['b0 & log10(a) are: ',num2str([EPGobj.mLogitGridBin.b0(dd),-EPGobj.mLgN(aa)])]);
            
            % map of the grid of the best "m"
            if 0 % show the mesh surface
                ll( ll < mx-zdepth ) = -inf;
                surf(-EPGobj.mLgN,EPGobj.mLogitGridBin.b0,ll,'EdgeColor','none');
                ylabel('b0'); xlabel('log_1_0(a)'); zlabel('Log likelihood');
            else % show the contours
                % compute the 68% and 95% CIs
                low68 = mx-0.5* (1);
                low95 = mx-0.5* (1.96^2);
                low99 = mx-0.5* (3^2);
%                 low68 = mx-0.5* (2.30);
%                 low95 = mx-0.5* (6.17);
%                 low99 = mx-0.5* (9.21);
                % display the contours
                contour(-EPGobj.mLgN,EPGobj.mLogitGridBin.b0,ll,[low99,low95,low68]);
                hold on; plot(-EPGobj.mLgN(aa),EPGobj.mLogitGridBin.b0(dd),'k*'); hold off;
                set(gca,'xminortick','on','yminortick','on');
                set(gca,'box','on');
                xlabel('log_1_0(a)'); ylabel('b0');
            end
        end
        function fLogitGridExactFig_b1_a_EUD(EPGobj,zdepth)
            if ~exist('zdepth','var')
                zdepth = 10;
            end
            
            % find the best "m"
            [mx,loc] = max(EPGobj.mLogitGrid.loglikelihood(:));
            [loc,~,~] = ind2sub(size(EPGobj.mLogitGrid.loglikelihood),loc);
            disp(['the b0 for the maximum likelihood in Logistic Regression of exact gEUD is: ',num2str(EPGobj.mLogitGrid.b0(loc))]);
            
            % TD50 and log10(a) corresponding to the best m
            ll = squeeze(EPGobj.mLogitGrid.loglikelihood(loc,:,:));
            [mm,aa] = find(ll == mx,1);
            disp(['b1 & log10(a) are: ',num2str([EPGobj.mLogitGrid.b1(mm),-EPGobj.mLgN(aa)])]);
            
            % map of the grid of the best "m"
            if 0 % show the mesh surface
                ll( ll < mx-zdepth ) = -inf;
                surf(-EPGobj.mLgN,EPGobj.mLogitGrid.b1,ll,'EdgeColor','none');
                ylabel('b1 (Gy^-^1)'); xlabel('log_1_0(a)'); zlabel('Log likelihood');
            else % show the contours
                % compute the 68% and 95% CIs
                low68 = mx-0.5* (1);
                low95 = mx-0.5* (1.96^2);
                low99 = mx-0.5* (3^2);
%                 low68 = mx-0.5* (2.30);
%                 low95 = mx-0.5* (6.17);
%                 low99 = mx-0.5* (9.21);
                % display the contours
                contour(-EPGobj.mLgN,EPGobj.mLogitGrid.b1,ll,[low99,low95,low68]);
                hold on; plot(-EPGobj.mLgN(aa),EPGobj.mLogitGrid.b1(mm),'k*'); hold off;
                set(gca,'xminortick','on','yminortick','on');
                set(gca,'box','on');
                xlabel('log_1_0(a)'); ylabel('b1 (Gy^-^1)');
            end
        end
        function fLogitGridBinFig_b1_a_EUD(EPGobj,zdepth)
            if ~exist('zdepth','var')
                zdepth = 10;
            end
            
            % find the best "m"
            [mx,loc] = max(EPGobj.mLogitGridBin.loglikelihood(:));
            [loc,~,~] = ind2sub(size(EPGobj.mLogitGridBin.loglikelihood),loc);
            disp(['the b0 for the maximum likelihood in Logistic Regression of gEUD Atlas is: ',num2str(EPGobj.mLogitGridBin.b0(loc))]);
            
            % TD50 and log10(a) corresponding to the best m
            ll = squeeze(EPGobj.mLogitGridBin.loglikelihood(loc,:,:));
            [mm,aa] = find(ll == mx,1);
            disp(['b1 & log10(a) are: ',num2str([EPGobj.mLogitGridBin.b1(mm),-EPGobj.mLgN(aa)])]);
            
            % map of the grid of the best "m"
            if 0 % show the mesh surface
                ll( ll < mx-zdepth ) = -inf;
                surf(-EPGobj.mLgN,EPGobj.mLogitGridBin.b1,ll,'EdgeColor','none');
                ylabel('b1 (Gy^-^1)'); xlabel('log_1_0(a)'); zlabel('Log likelihood');
            else % show the contours
                % compute the 68% and 95% CIs
                low68 = mx-0.5* (1);
                low95 = mx-0.5* (1.96^2);
                low99 = mx-0.5* (3^2);
%                 low68 = mx-0.5* (2.30);
%                 low95 = mx-0.5* (6.17);
%                 low99 = mx-0.5* (9.21);
                % display the contours
                contour(-EPGobj.mLgN,EPGobj.mLogitGridBin.b1,ll,[low99,low95,low68]);
                hold on; plot(-EPGobj.mLgN(aa),EPGobj.mLogitGridBin.b1(mm),'k*'); hold off;
                set(gca,'xminortick','on','yminortick','on');
                set(gca,'box','on');
                xlabel('log_1_0(a)'); ylabel('b1 (Gy^-^1)');
            end
        end

        function fLymanGridExactFig_TD50_m_EUD(EPGobj,loga,zdepth)
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
            [mx,loc] = max(EPGobj.mLymanGrid.loglikelihood(:));
            [~,~,loc] = ind2sub(size(EPGobj.mLymanGrid.loglikelihood),loc);
%             st = [EPGobj.mLogitMat];
%             dpf = [st.dev]; % deviations
%             st =[st.stats];
%             df = [st.dfe]; % degree of freedom
%             dpf = dpf./df; % deviations per degree of freedom
%             loglikelyhood = -0.5*dpf;
%             [~,loc] = max(loglikelyhood); % the maximum loglikelyhood
            disp(['the "log10(a)" for the maximum likelihood in Lyman model of exact gEUD is: ',num2str(-EPGobj.mLgN(loc))]);
            if ~exist('loga','var') || ischar(loga)
                loga = -EPGobj.mLgN(loc);
            end
            [~,loc] = min(abs(EPGobj.mLgN+loga)); % the n whose corresponding responding function will be ploted
            disp(['the "log10(a)" for the map of likelihood in Lyman model is: ',num2str(-EPGobj.mLgN(loc))]);

            
            % TD50 and gamma of the best "log10(a)"
%             lymangrid = EPGobj.mLymanGrid;
            ll = EPGobj.mLymanGrid.loglikelihood(:,:,loc);
            mx = max(ll(:));
            [dd,mm] = find(ll == mx,1);
            disp(['TD50 & m are: ',num2str([EPGobj.mLymanGrid.TD50(dd),EPGobj.mLymanGrid.m(mm)])]);
            
            % map of the grid of the best "a"
            if 0 % show the mesh surface
                ll( ll < mx-zdepth ) = -inf;
                surf(EPGobj.mLymanGrid.m,EPGobj.mLymanGrid.TD50,ll,'EdgeColor','none');
                ylabel('TD50 (Gy)'); xlabel('m'); zlabel('Log likelihood');
            else % show the contours
                % compute the 68% and 95% CIs
                low68 = mx-0.5* (1);
                low95 = mx-0.5* (1.96^2);
                low99 = mx-0.5* (3^2);
%                 low68 = mx-0.5* (2.30);
%                 low95 = mx-0.5* (6.17);
%                 low99 = mx-0.5* (9.21);
                % display the contours
                [c,h] = contour(EPGobj.mLymanGrid.TD50,EPGobj.mLymanGrid.m,ll',[low99,low95,low68]);
                hold on; plot(EPGobj.mLymanGrid.TD50(dd),EPGobj.mLymanGrid.m(mm),'k*'); hold off;
%                 text_handle = clabel(c,h,'LabelSpacing',100000);%'String',{'low 68%';'low 95%';'low 99.7%'},
%                 set(text_handle,'BackgroundColor',[1 1 .6], 'Edgecolor',[.7 .7 .7], 'String',{'low 68%';'low 95%';'low 99.7%'});
                set(gca,'xminortick','on','yminortick','on');
                set(gca,'box','on');
                xlabel('TD50 (Gy)'); ylabel('m');
            end

%             % the coresponding parameters computed from atlas
%             if ~isempty(EPGobj.mLymanGridBin)
%                 ll = EPGobj.mLymanGridBin.loglikelihood(:,:,loc); % log likelihood of log10(a) = loga
%                 mx = max(ll(:)); % the best point for this log10(a)
%                 [dd,mm] = find(ll == mx); % the coefficients
%                 disp(['TD50 & m for the same log10(a) in the atlas are: ', num2str([EPGobj.mLymanGridBin.TD50(dd),EPGobj.mLymanGridBin.m(mm)])]);
%             end
        end
        function fLymanGridBinFig_TD50_m_EUD(EPGobj,loga,zdepth)
            if ~exist('zdepth','var')
                zdepth = 10;
            end
            
            % find the best "a"
            [mx,loc] = max(EPGobj.mLymanGridBin.loglikelihood(:));
            [~,~,loc] = ind2sub(size(EPGobj.mLymanGridBin.loglikelihood),loc);
            disp(['the "log10(a)" for the maximum likelihood in Lyman model of gEUD Atlas is: ',num2str(-EPGobj.mLgN(loc))]);
            if ~exist('loga','var') || ischar(loga)
                loga = -EPGobj.mLgN(loc);
            end
            [~,loc] = min(abs(EPGobj.mLgN+loga)); % the n whose corresponding responding function will be ploted
            disp(['the "log10(a)" for the map of likelihood in Lyman model is: ',num2str(-EPGobj.mLgN(loc))]);
            
            % TD50 and gamma of the best "log10(a)"
            ll = EPGobj.mLymanGridBin.loglikelihood(:,:,loc);
            mx = max(ll(:));
            [dd,mm] = find(ll == mx,1);
            disp(['TD50 & m are: ',num2str([EPGobj.mLymanGridBin.TD50(dd),EPGobj.mLymanGridBin.m(mm)])]);
            
            % map of the grid of the best "a"
            if 0 % show the mesh surface
                ll( ll < mx-zdepth ) = -inf;
                surf(EPGobj.mLymanGridBin.m,EPGobj.mLymanGridBin.TD50,ll,'EdgeColor','none');
                ylabel('TD50 (Gy)'); xlabel('m'); zlabel('Log likelihood');
            else % show the contours
                % compute the 68% and 95% CIs
                low68 = mx-0.5* (1);
                low95 = mx-0.5* (1.96^2);
                low99 = mx-0.5* (3^2);
%                 low68 = mx-0.5* (2.30);
%                 low95 = mx-0.5* (6.17);
%                 low99 = mx-0.5* (9.21);
                % display the contours
                contour(EPGobj.mLymanGridBin.TD50,EPGobj.mLymanGridBin.m,ll',[low99,low95,low68]);
                hold on; plot(EPGobj.mLymanGridBin.TD50(dd),EPGobj.mLymanGridBin.m(mm),'k*'); hold off;
                set(gca,'xminortick','on','yminortick','on');
                set(gca,'box','on');
                xlabel('TD50 (Gy)'); ylabel('m');
            end
        end
        function fLymanGridExactFig_TD50_a_EUD(EPGobj,zdepth)
            if ~exist('zdepth','var')
                zdepth = 10;
            end
            
            % find the best "m"
            [mx,loc] = max(EPGobj.mLymanGrid.loglikelihood(:));
            [~,loc,~] = ind2sub(size(EPGobj.mLymanGrid.loglikelihood),loc);
            disp(['the "m" for the maximum likelihood in Lyman model of exact gEUD is: ',num2str(EPGobj.mLymanGrid.m(loc))]);
            
            % TD50 and log10(a) corresponding to the best m
            ll = squeeze(EPGobj.mLymanGrid.loglikelihood(:,loc,:));
            [dd,aa] = find(ll == mx,1);
            disp(['TD50 & a are: ',num2str([EPGobj.mLymanGrid.TD50(dd),-EPGobj.mLgN(aa)])]);
            
            % map of the grid of the best "m"
            if 0 % show the mesh surface
                ll( ll < mx-zdepth ) = -inf;
                surf(-EPGobj.mLgN,EPGobj.mLymanGrid.TD50,ll,'EdgeColor','none');
                ylabel('TD50 (Gy)'); xlabel('log_1_0(a)'); zlabel('Log likelihood');
            else % show the contours
                % compute the 68% and 95% CIs
                low68 = mx-0.5* (1);
                low95 = mx-0.5* (1.96^2);
                low99 = mx-0.5* (3^2);
%                 low68 = mx-0.5* (2.30);
%                 low95 = mx-0.5* (6.17);
%                 low99 = mx-0.5* (9.21);
                % display the contours
                contour(-EPGobj.mLgN,EPGobj.mLymanGrid.TD50,ll,[low99,low95,low68]);
                hold on; plot(-EPGobj.mLgN(aa),EPGobj.mLymanGrid.TD50(dd),'k*'); hold off;
                set(gca,'xminortick','on','yminortick','on');
                set(gca,'box','on');
                xlabel('log_1_0(a)'); ylabel('TD50 (Gy)');
            end
        end
        function fLymanGridBinFig_TD50_a_EUD(EPGobj,zdepth)
            if ~exist('zdepth','var')
                zdepth = 10;
            end
            
            % find the best "m"
            [mx,loc] = max(EPGobj.mLymanGridBin.loglikelihood(:));
            [~,loc,~] = ind2sub(size(EPGobj.mLymanGridBin.loglikelihood),loc);
            disp(['the "m" for the maximum likelihood in Lyman model of gEUD Atlas is: ',num2str(EPGobj.mLymanGridBin.m(loc))]);
            
            % TD50 and log10(a) corresponding to the best m
            ll = squeeze(EPGobj.mLymanGridBin.loglikelihood(:,loc,:));
            [dd,aa] = find(ll == mx,1);
            disp(['TD50 & a are: ',num2str([EPGobj.mLymanGridBin.TD50(dd),-EPGobj.mLgN(aa)])]);
            
            % map of the grid of the best "m"
            if 0 % show the mesh surface
                ll( ll < mx-zdepth ) = -inf;
                surf(-EPGobj.mLgN,EPGobj.mLymanGridBin.TD50,ll,'EdgeColor','none');
                ylabel('TD50 (Gy)'); xlabel('log_1_0(a)'); zlabel('Log likelihood');
            else % show the contours
                % compute the 68% and 95% CIs
                low68 = mx-0.5* (1);
                low95 = mx-0.5* (1.96^2);
                low99 = mx-0.5* (3^2);
%                 low68 = mx-0.5* (2.30);
%                 low95 = mx-0.5* (6.17);
%                 low99 = mx-0.5* (9.21);
                % display the contours
                contour(-EPGobj.mLgN,EPGobj.mLymanGridBin.TD50,ll,[low99,low95,low68]);
                hold on; plot(-EPGobj.mLgN(aa),EPGobj.mLymanGridBin.TD50(dd),'k*'); hold off;
                set(gca,'xminortick','on','yminortick','on');
                set(gca,'box','on');
                xlabel('log_1_0(a)'); ylabel('TD50 (Gy)');
            end
        end
        function fLymanGridExactFig_m_a_EUD(EPGobj,zdepth)
            if ~exist('zdepth','var')
                zdepth = 10;
            end
            
            % find the best "m"
            [mx,loc] = max(EPGobj.mLymanGrid.loglikelihood(:));
            [loc,~,~] = ind2sub(size(EPGobj.mLymanGrid.loglikelihood),loc);
            disp(['the "TD50" for the maximum likelihood in Lyman model of exact gEUD is: ',num2str(EPGobj.mLymanGrid.TD50(loc))]);
            
            % TD50 and log10(a) corresponding to the best m
            ll = squeeze(EPGobj.mLymanGrid.loglikelihood(loc,:,:));
            [mm,aa] = find(ll == mx,1);
            disp(['m & a are: ',num2str([EPGobj.mLymanGrid.m(mm),-EPGobj.mLgN(aa)])]);
            
            % map of the grid of the best "m"
            if 0 % show the mesh surface
                ll( ll < mx-zdepth ) = -inf;
                surf(-EPGobj.mLgN,EPGobj.mLymanGrid.m,ll,'EdgeColor','none');
                ylabel('m'); xlabel('log_1_0(a)'); zlabel('Log likelihood');
            else % show the contours
                % compute the 68% and 95% CIs
                low68 = mx-0.5* (1);
                low95 = mx-0.5* (1.96^2);
                low99 = mx-0.5* (3^2);
%                 low68 = mx-0.5* (2.30);
%                 low95 = mx-0.5* (6.17);
%                 low99 = mx-0.5* (9.21);
                % display the contours
                contour(-EPGobj.mLgN,EPGobj.mLymanGrid.m,ll,[low99,low95,low68]);
                hold on; plot(-EPGobj.mLgN(aa),EPGobj.mLymanGrid.m(mm),'k*'); hold off;
                set(gca,'xminortick','on','yminortick','on');
                set(gca,'box','on');
                xlabel('log_1_0(a)'); ylabel('m');
            end
        end
        function fLymanGridBinFig_m_a_EUD(EPGobj,zdepth)
            if ~exist('zdepth','var')
                zdepth = 10;
            end
            
            % find the best "m"
            [mx,loc] = max(EPGobj.mLymanGridBin.loglikelihood(:));
            [loc,~,~] = ind2sub(size(EPGobj.mLymanGridBin.loglikelihood),loc);
            disp(['the "TD50" for the maximum likelihood in Lyman model of gEUD Atlas is: ',num2str(EPGobj.mLymanGridBin.TD50(loc))]);
            
            % TD50 and log10(a) corresponding to the best m
            ll = squeeze(EPGobj.mLymanGridBin.loglikelihood(loc,:,:));
            [mm,aa] = find(ll == mx,1);
            disp(['m & a are: ',num2str([EPGobj.mLymanGridBin.m(mm),-EPGobj.mLgN(aa)])]);
            
            % map of the grid of the best "m"
            if 0 % show the mesh surface
                ll( ll < mx-zdepth ) = -inf;
                surf(-EPGobj.mLgN,EPGobj.mLymanGridBin.m,ll,'EdgeColor','none');
                ylabel('m'); xlabel('log_1_0(a)'); zlabel('Log likelihood');
            else % show the contours
                % compute the 68% and 95% CIs
                low68 = mx-0.5* (1);
                low95 = mx-0.5* (1.96^2);
                low99 = mx-0.5* (3^2);
%                 low68 = mx-0.5* (2.30);
%                 low95 = mx-0.5* (6.17);
%                 low99 = mx-0.5* (9.21);
                % display the contours
                contour(-EPGobj.mLgN,EPGobj.mLymanGridBin.m,ll,[low99,low95,low68]);
                hold on; plot(-EPGobj.mLgN(aa),EPGobj.mLymanGridBin.m(mm),'k*'); hold off;
                set(gca,'xminortick','on','yminortick','on');
                set(gca,'box','on');
                xlabel('log_1_0(a)'); ylabel('m');
            end
        end
        
        function fLymanGridExactFig_a_loglikelihood(EPGobj,strMarker,lw)
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
            [mx,loc] = max(EPGobj.mLymanGrid.loglikelihood(:));
            [~,~,loc] = ind2sub(size(EPGobj.mLymanGrid.loglikelihood),loc);
            disp(['the maximum log likelihood and its corresponding "log10(a)" in Lyman model of exact gEUD are: ',num2str([mx, -EPGobj.mLgN(loc)])]);

            % plot the maximum log likelihood curve w.r.t. log10(a)
            llmx = zeros(size(EPGobj.mLgN)); %log likelihood's max for each log10(a)
            for kk = 1:length(llmx)
                ll = EPGobj.mLymanGrid.loglikelihood(:,:,kk);
                llmx(kk) = max(ll(:));
            end
            llmx = llmx / (EPGobj.mNumInGrp-2);
            mx = mx / (EPGobj.mNumInGrp-2);
            loglikelyhood68 = repmat(mx-0.5* 1 /(EPGobj.mNumInGrp-2),size(llmx));
            loglikelyhood95 = repmat(mx-0.5* (1.96*2) /(EPGobj.mNumInGrp-2),size(llmx));
            hold on;
            plot(-EPGobj.mLgN,llmx,strMarker,'LineWidth',lw);
            plot(-EPGobj.mLgN, loglikelyhood68,strcat(strMarker(1),'--'),'LineWidth',1);
            plot(-EPGobj.mLgN, loglikelyhood95,strcat(strMarker(1),'-'),'LineWidth',1);
            if length(strMarker)==1
                strMarker = [strMarker,'*'];
            elseif strMarker(2) == '-' || strMarker(2) == '.'
                strMarker(2) = '*';
            end
            plot(-EPGobj.mLgN(loc),mx,strMarker,'LineWidth',lw+2);
            set(gca,'xminortick','on','yminortick','on');
            set(gca,'box','on');
            xlabel('log_1_0(a)'); ylabel('log likelihood per degree of freedom');
        end
        function fLymanGridBinFig_a_loglikelihood(EPGobj,strMarker,lw)
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
            [mx,loc] = max(EPGobj.mLymanGridBin.loglikelihood(:));
            [~,~,loc] = ind2sub(size(EPGobj.mLymanGridBin.loglikelihood),loc);
            disp(['the maximum log likelihood and its corresponding "log10(a)" in Lyman model of gEUD Atlas are: ',num2str([mx, -EPGobj.mLgN(loc)])]);

            % plot the maximum log likelihood curve w.r.t. log10(a)
            llmx = zeros(size(EPGobj.mLgN)); %log likelihood's max for each log10(a)
            for kk = 1:length(llmx)
                ll = EPGobj.mLymanGridBin.loglikelihood(:,:,kk);
                llmx(kk) = max(ll(:));
            end
            llmx = llmx / (EPGobj.mNumInGrp-2);
            mx = mx / (EPGobj.mNumInGrp-2);
            loglikelyhood68 = repmat(mx-0.5* 1 /(EPGobj.mNumInGrp-2),size(llmx));
            loglikelyhood95 = repmat(mx-0.5* (1.96*2) /(EPGobj.mNumInGrp-2),size(llmx));
            hold on;
            plot(-EPGobj.mLgN,llmx,strMarker,'LineWidth',lw);
            plot(-EPGobj.mLgN(loc),mx,strMarker,'LineWidth',lw+2);
%             plot(-EPGobj.mLgN, loglikelyhood68,strcat(strMarker(1),'--'),'LineWidth',1);
%             plot(-EPGobj.mLgN, loglikelyhood95,strcat(strMarker(1),'-'),'LineWidth',1);
            set(gca,'xminortick','on','yminortick','on');
            set(gca,'box','on');
            xlabel('log_1_0(a)'); ylabel('log likelihood per degree of freedom');
        end
        function fLymanGridExactFig_TD50_loglikelihood(EPGobj,strMarker,lw)
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
            [mx,loc] = max(EPGobj.mLymanGrid.loglikelihood(:));
            [loc,~,~] = ind2sub(size(EPGobj.mLymanGrid.loglikelihood),loc);
            disp(['the maximum log likelihood and its corresponding "TD50" in Lyman model of exact gEUD are: ',num2str([mx, EPGobj.mLymanGrid.TD50(loc)])]);

            % plot the maximum log likelihood curve w.r.t. TD50
            llmx = zeros(size(EPGobj.mLymanGrid.TD50)); %log likelihood's max for each log10(a)
            for kk = 1:length(llmx)
                ll = EPGobj.mLymanGrid.loglikelihood(kk,:,:);
                llmx(kk) = max(ll(:));
            end
            llmx = llmx / (EPGobj.mNumInGrp-2);
            mx = mx / (EPGobj.mNumInGrp-2);
            loglikelyhood68 = repmat(mx-0.5* 1 /(EPGobj.mNumInGrp-2),size(llmx));
            loglikelyhood95 = repmat(mx-0.5* (1.96*2) /(EPGobj.mNumInGrp-2),size(llmx));
            hold on;
            plot(EPGobj.mLymanGrid.TD50,llmx,strMarker,'LineWidth',lw);
            plot(EPGobj.mLymanGrid.TD50, loglikelyhood68,strcat(strMarker(1),'--'),'LineWidth',1);
            plot(EPGobj.mLymanGrid.TD50, loglikelyhood95,strcat(strMarker(1),'-'),'LineWidth',1);
            if length(strMarker)==1
                strMarker = [strMarker,'*'];
            elseif strMarker(2) == '-' || strMarker(2) == '.'
                strMarker(2) = '*';
            end
            plot(EPGobj.mLymanGrid.TD50(loc),mx,strMarker,'LineWidth',lw+2);
            set(gca,'xminortick','on','yminortick','on');
            set(gca,'box','on');
            xlabel('TD50'); ylabel('log likelihood per degree of freedom');
        end
        function fLymanGridBinFig_TD50_loglikelihood(EPGobj,strMarker,lw)
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
            [mx,loc] = max(EPGobj.mLymanGridBin.loglikelihood(:));
            [loc,~,~] = ind2sub(size(EPGobj.mLymanGridBin.loglikelihood),loc);
            disp(['the maximum log likelihood and its corresponding "TD50" in Lyman model of gEUD Atlas are: ',num2str([mx, EPGobj.mLymanGridBin.TD50(loc)])]);

            % plot the maximum log likelihood curve w.r.t. TD50
            llmx = zeros(size(EPGobj.mLymanGridBin.TD50)); %log likelihood's max for each log10(a)
            for kk = 1:length(llmx)
                ll = EPGobj.mLymanGridBin.loglikelihood(kk,:,:);
                llmx(kk) = max(ll(:));
            end
            llmx = llmx / (EPGobj.mNumInGrp-2);
            mx = mx / (EPGobj.mNumInGrp-2);
            loglikelyhood68 = repmat(mx-0.5* 1 /(EPGobj.mNumInGrp-2),size(llmx));
            loglikelyhood95 = repmat(mx-0.5* (1.96*2) /(EPGobj.mNumInGrp-2),size(llmx));
            hold on;
            plot(EPGobj.mLymanGridBin.TD50,llmx,strMarker,'LineWidth',lw);
            plot(EPGobj.mLymanGridBin.TD50(loc),mx,strMarker,'LineWidth',lw+2);
%             plot(EPGobj.mLymanGrid.TD50, loglikelyhood68,strcat(strMarker(1),'--'),'LineWidth',1);
%             plot(EPGobj.mLymanGrid.TD50, loglikelyhood95,strcat(strMarker(1),'-'),'LineWidth',1);
            set(gca,'xminortick','on','yminortick','on');
            set(gca,'box','on');
            xlabel('TD50'); ylabel('log likelihood per degree of freedom');
        end
        function fLymanGridExactFig_m_loglikelihood(EPGobj,strMarker,lw)
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
            [mx,loc] = max(EPGobj.mLymanGrid.loglikelihood(:));
            [~,loc,~] = ind2sub(size(EPGobj.mLymanGrid.loglikelihood),loc);
            disp(['the maximum log likelihood and its corresponding "m" in Lyman model of exact gEUD are: ',num2str([mx, EPGobj.mLymanGrid.m(loc)])]);

            % plot the maximum log likelihood curve w.r.t. TD50
            llmx = zeros(size(EPGobj.mLymanGrid.m)); %log likelihood's max for each log10(a)
            for kk = 1:length(llmx)
                ll = EPGobj.mLymanGrid.loglikelihood(:,kk,:);
                llmx(kk) = max(ll(:));
            end
            llmx = llmx / (EPGobj.mNumInGrp-2);
            mx = mx / (EPGobj.mNumInGrp-2);
            loglikelyhood68 = repmat(mx-0.5* 1 /(EPGobj.mNumInGrp-2),size(llmx));
            loglikelyhood95 = repmat(mx-0.5* (1.96*2) /(EPGobj.mNumInGrp-2),size(llmx));
            hold on;
            plot(EPGobj.mLymanGrid.m,llmx,strMarker,'LineWidth',lw);
            plot(EPGobj.mLymanGrid.m, loglikelyhood68,strcat(strMarker(1),'--'),'LineWidth',1);
            plot(EPGobj.mLymanGrid.m, loglikelyhood95,strcat(strMarker(1),'-'),'LineWidth',1);
            if length(strMarker)==1
                strMarker = [strMarker,'*'];
            elseif strMarker(2) == '-' || strMarker(2) == '.'
                strMarker(2) = '*';
            end
            plot(EPGobj.mLymanGrid.m(loc),mx,strMarker,'LineWidth',lw+2);
            set(gca,'xminortick','on','yminortick','on');
            set(gca,'box','on');
            xlabel('m'); ylabel('log likelihood per degree of freedom');
        end
        function fLymanGridBinFig_m_loglikelihood(EPGobj,strMarker,lw)
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
            [mx,loc] = max(EPGobj.mLymanGridBin.loglikelihood(:));
            [~,loc,~] = ind2sub(size(EPGobj.mLymanGridBin.loglikelihood),loc);
            disp(['the maximum log likelihood and its corresponding "m" in Lyman model of gEUD Atlas are: ',num2str([mx, EPGobj.mLymanGridBin.m(loc)])]);

            % plot the maximum log likelihood curve w.r.t. TD50
            llmx = zeros(size(EPGobj.mLymanGridBin.m)); %log likelihood's max for each log10(a)
            for kk = 1:length(llmx)
                ll = EPGobj.mLymanGridBin.loglikelihood(:,kk,:);
                llmx(kk) = max(ll(:));
            end
            llmx = llmx / (EPGobj.mNumInGrp-2);
            mx = mx / (EPGobj.mNumInGrp-2);
            loglikelyhood68 = repmat(mx-0.5* 1 /(EPGobj.mNumInGrp-2),size(llmx));
            loglikelyhood95 = repmat(mx-0.5* (1.96*2) /(EPGobj.mNumInGrp-2),size(llmx));
            hold on;
            plot(EPGobj.mLymanGridBin.m,llmx,strMarker,'LineWidth',lw);
            plot(EPGobj.mLymanGridBin.m(loc),mx,strMarker,'LineWidth',lw+2);
%             plot(EPGobj.mLymanGrid.m, loglikelyhood68,strcat(strMarker(1),'--'),'LineWidth',1);
%             plot(EPGobj.mLymanGrid.m, loglikelyhood95,strcat(strMarker(1),'-'),'LineWidth',1);
            set(gca,'xminortick','on','yminortick','on');
            set(gca,'box','on');
            xlabel('m'); ylabel('log likelihood per degree of freedom');
        end
        
        function fLymanGridExactFig_TD50_loglikelihoodAtLoga_EUD(EPGobj,loga,strMarker,lw)
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
            [~,loc] = max(EPGobj.mLymanGrid.loglikelihood(:));
            [~,~,loc] = ind2sub(size(EPGobj.mLymanGrid.loglikelihood),loc);
            disp(['the "log10(a)" for the maximum likelihood in Lyman model of exact gEUD is: ',num2str(-EPGobj.mLgN(loc))]);
            if ~exist('loga','var') || ischar(loga)
                loga = -EPGobj.mLgN(loc);
            end
            [~,loc] = min(abs(EPGobj.mLgN+loga)); % the n whose corresponding responding function will be ploted
            disp(['the "log10(a)" for the coefficients is: ',num2str(-EPGobj.mLgN(loc))]);

            % plot the maximum log likelihood curve w.r.t. TD50 under a spicific loga
            llmx = max(EPGobj.mLymanGrid.loglikelihood(:,:,loc),[],2);
            llmx = llmx / (EPGobj.mNumInGrp-2);
            [mx,loc] = max(llmx);
            disp(['the best TD50 is: ',num2str(EPGobj.mLymanGrid.TD50(loc))]);
            loglikelyhood68 = repmat(mx-0.5* 1 /(EPGobj.mNumInGrp-2),size(llmx));
            loglikelyhood95 = repmat(mx-0.5* (1.96*2) /(EPGobj.mNumInGrp-2),size(llmx));
            hold on;
            plot(EPGobj.mLymanGrid.TD50,llmx,strMarker,'LineWidth',lw);
            plot(EPGobj.mLymanGrid.TD50, loglikelyhood68,strcat(strMarker(1),'--'),'LineWidth',1);
            plot(EPGobj.mLymanGrid.TD50, loglikelyhood95,strcat(strMarker(1),'-'),'LineWidth',1);
            if length(strMarker)==1
                strMarker = [strMarker,'*'];
            elseif strMarker(2) == '-' || strMarker(2) == '.'
                strMarker(2) = '*';
            end
            plot(EPGobj.mLymanGrid.TD50(loc),mx,strMarker,'LineWidth',lw+2);
            set(gca,'xminortick','on','yminortick','on');
            set(gca,'box','on');
            xlabel('TD50'); ylabel('log likelihood per degree of freedom');
        end
        function fLymanGridExactFig_m_loglikelihoodAtLoga_EUD(EPGobj,loga,strMarker,lw)
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
            [~,loc] = max(EPGobj.mLymanGrid.loglikelihood(:));
            [~,~,loc] = ind2sub(size(EPGobj.mLymanGrid.loglikelihood),loc);
            disp(['the "log10(a)" for the maximum likelihood in Lyman model of exact gEUD is: ',num2str(-EPGobj.mLgN(loc))]);
            if ~exist('loga','var') || ischar(loga)
                loga = -EPGobj.mLgN(loc);
            end
            [~,loc] = min(abs(EPGobj.mLgN+loga)); % the n whose corresponding responding function will be ploted
            disp(['the "log10(a)" for the coefficients is: ',num2str(-EPGobj.mLgN(loc))]);

            % plot the maximum log likelihood curve w.r.t. TD50 under a spicific loga
            llmx = max(EPGobj.mLymanGrid.loglikelihood(:,:,loc),[],1);
            llmx = llmx / (EPGobj.mNumInGrp-2);
            [mx,loc] = max(llmx);
            disp(['the best m is: ',num2str(EPGobj.mLymanGrid.m(loc))]);
            loglikelyhood68 = repmat(mx-0.5* 1 /(EPGobj.mNumInGrp-2),size(llmx));
            loglikelyhood95 = repmat(mx-0.5* (1.96*2) /(EPGobj.mNumInGrp-2),size(llmx));
            hold on;
            plot(EPGobj.mLymanGrid.m,llmx,strMarker,'LineWidth',lw);
            plot(EPGobj.mLymanGrid.m, loglikelyhood68,strcat(strMarker(1),'--'),'LineWidth',1);
            plot(EPGobj.mLymanGrid.m, loglikelyhood95,strcat(strMarker(1),'-'),'LineWidth',1);
            if length(strMarker)==1
                strMarker = [strMarker,'*'];
            elseif strMarker(2) == '-' || strMarker(2) == '.'
                strMarker(2) = '*';
            end
            plot(EPGobj.mLymanGrid.m(loc),mx,strMarker,'LineWidth',lw+2);
            set(gca,'xminortick','on','yminortick','on');
            set(gca,'box','on');
            xlabel('m'); ylabel('log likelihood per degree of freedom');
        end

        function fLymanGridGOFFig(EPGobj,loga)
            % compute the p-value curves
            nmn = 5; nmx = 30;
            nGrp = [nmn:nmx, EPGobj.mNumInGrp]';
            nGrp = [nmn:nmx]';
            pValExact = zeros(length(nGrp),3);
            pValBin = zeros(length(nGrp),3);
            for n = 1:length(nGrp)
                % assign number of groups
                EPGobj.mLymanHosmerLemeshow.n = nGrp(n);
                EPGobj.mLymanHosmerLemeshowBin.n = nGrp(n);
                EPGobj.mLymanGTest.n = nGrp(n);
                EPGobj.mLymanGTestBin.n = nGrp(n);
                EPGobj.mLymanPearson.n = nGrp(n);
                EPGobj.mLymanPearsonBin.n = nGrp(n);

                % compute the goodness of fit
                EPGobj = EPGobj.LymanHosmerLemeshowTestAnalysisExact_EUD(loga);
                EPGobj = EPGobj.LymanHosmerLemeshowTestAnalysisBin_EUD(loga);
                EPGobj = EPGobj.LymanGTestAnalysisExact_EUD(loga);
                EPGobj = EPGobj.LymanGTestAnalysisBin_EUD(loga);
                EPGobj = EPGobj.LymanPearsonTestAnalysisExact_EUD(loga);
                EPGobj = EPGobj.LymanPearsonTestAnalysisBin_EUD(loga);

                % save
                pValExact(n,1) = EPGobj.mLymanHosmerLemeshow.p_value;
                pValExact(n,2) = EPGobj.mLymanGTest.p_value;
                pValExact(n,3) = EPGobj.mLymanPearson.p_value;
                pValBin(n,1) = EPGobj.mLymanHosmerLemeshowBin.p_value;
                pValBin(n,2) = EPGobj.mLymanGTestBin.p_value;
                pValBin(n,3) = EPGobj.mLymanPearsonBin.p_value;
            end

            % plot p-value curves
            hold on;
            plot(nGrp,pValExact(:,1),'b--',nGrp,pValBin(:,1),'b-');
            plot(nGrp,pValExact(:,2),'r--',nGrp,pValBin(:,2),'r-');
            plot(nGrp,pValExact(:,3),'k--',nGrp,pValBin(:,3),'k-');
            set(gca,'xminortick','on','yminortick','on');
            set(gca,'box','on');
            xlabel('number of groups'); ylabel('p-value of goodness of fit');
        end

        function loga = fLymanGridResponseExactFig_a_EUD(EPGobj,loga,strMarker,lw)
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
            [mx,loc] = max(EPGobj.mLymanGrid.loglikelihood(:));
            [~,~,loc] = ind2sub(size(EPGobj.mLymanGrid.loglikelihood),loc);
            disp(['the maximum log likelihood and its "log10(a)" in Lyman model of exact gEUD are: ',num2str([mx, -EPGobj.mLgN(loc)])]);
            if ~exist('loga','var') || ischar(loga)
                loga = -EPGobj.mLgN(loc);
            end
            disp(['the plotted response curve is at log10(a) = ',num2str(loga)]);

            % coefficients for the Lyman model
            [~,n] = min(abs(EPGobj.mLgN+loga)); % the n whose corresponding responding function will be ploted
            ll = EPGobj.mLymanGrid.loglikelihood(:,:,n); % log likelihood of log10(a) = loga
            mx = max(ll(:));
            [xx,yy] = find(ll == mx); % the coefficients
            TD50 = EPGobj.mLymanGrid.TD50(xx);
            m = EPGobj.mLymanGrid.m(yy);
            disp(['the TD50 and m for the response function are: ',num2str([TD50,m])]);

            % curves
            euds = [EPGobj.mGrp.EUD]; euds = euds(n,:); % the gEUDs of that n
            doses = (0:max(euds))'; % doses (gEUD) whose RP probability will be computed
            rpb = normcdf((doses-TD50)/(m*TD50),0,1); % Lyman probability
            % 95% CI
            low95 = mx-0.5* (1.96*2); % CI level
            cxy = contourc(EPGobj.mLymanGrid.TD50,EPGobj.mLymanGrid.m,ll',[low95,low95]); % coordinates of CI contours at CI level
            CImx = zeros(length(doses),1); % upper end of CI for response at each gEUD
            CImn = zeros(length(doses),1); % lower end of CI for response at each gEUD
            TD50 = cxy(1,2:end);
            m = cxy(2,2:end); 
            for d = 1:length(doses)
                CI = normcdf((doses(d)-TD50)./(m.*TD50),0,1); % probability at each iso point
                CImx(d) = max(CI);
                CImn(d) = min(CI);
            end
            % plot
            hold on;
            plot(doses,rpb,strMarker,'LineWidth',lw*2); % responding function
            plot(doses,CImx,strMarker,'LineWidth',lw); % responding function
            plot(doses,CImn,strMarker,'LineWidth',lw); % responding function
            set(gca,'xminortick','on','yminortick','on');
            set(gca,'box','on');
            xlabel('gEUD'); ylabel('RP probability');
        end
        function loga = fLymanGridResponseBinFig_a_EUD(EPGobj,loga,strMarker,lw)
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
            [mx,loc] = max(EPGobj.mLymanGridBin.loglikelihood(:));
            [~,~,loc] = ind2sub(size(EPGobj.mLymanGridBin.loglikelihood),loc);
            disp(['the maximum log likelihood and its "log10(a)" in Lyman model of Bin gEUD are: ',num2str([mx, -EPGobj.mLgN(loc)])]);
            if ~exist('loga','var') || ischar(loga)
                loga = -EPGobj.mLgN(loc);
            end
            disp(['the plotted response curve is at log10(a) = ',num2str(loga)]);

            % coefficients for the Lyman model
            [~,n] = min(abs(EPGobj.mLgN+loga)); % the n whose corresponding responding function will be ploted
            ll = EPGobj.mLymanGridBin.loglikelihood(:,:,n); % log likelihood of log10(a) = loga
            mx = max(ll(:));
            [xx,yy] = find(ll == mx); % the coefficients
            TD50 = EPGobj.mLymanGridBin.TD50(xx);
            m = EPGobj.mLymanGridBin.m(yy);
            disp(['the TD50 and m for the response function are: ',num2str([TD50,m])]);

            % curves
            euds = [EPGobj.mGrp.EUD]; euds = euds(n,:); % the gEUDs of that n
            doses = (0:max(euds))'; % doses (gEUD) whose RP probability will be computed
            rpb = normcdf((doses-TD50)/(m*TD50),0,1); % Lyman probability
            % 95% CI
            low95 = mx-0.5* (1.96*2); % CI level
            cxy = contourc(EPGobj.mLymanGridBin.TD50,EPGobj.mLymanGridBin.m,ll',[low95,low95]); % coordinates of CI contours at CI level
            CImx = zeros(length(doses),1); % upper end of CI for response at each gEUD
            CImn = zeros(length(doses),1); % lower end of CI for response at each gEUD
            TD50 = cxy(1,2:end);
            m = cxy(2,2:end); 
            for d = 1:length(doses)
                CI = normcdf((doses(d)-TD50)./(m.*TD50),0,1); % probability at each iso point
                CImx(d) = max(CI);
                CImn(d) = min(CI);
            end
            % plot
            hold on;
            plot(doses,rpb,strMarker,'LineWidth',lw*2); % responding function
            plot(doses,CImx,strMarker,'LineWidth',lw); % responding function
            plot(doses,CImn,strMarker,'LineWidth',lw); % responding function
            set(gca,'xminortick','on','yminortick','on');
            set(gca,'box','on');
            xlabel('gEUD'); ylabel('RP probability');
        end
        
        function [medianeud,betainv84,betainv16] = fComplicationObservedFig_EUD(EPGobj,loga,numIntervals,strMarker,lw)
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
                st = [EPGobj.mLogitMatBin];
                dpf = [st.dev]; % deviations
                [~,loc] = min(dpf); % the minimum deviation corresponds to the maximum likelihood
                loga = -EPGobj.mLgN(loc);
            end
            if ~exist('numIntervals','var')
                numIntervals = 4;
            end
            
            % parameters
            if ~exist('loga','var') || ischar(loga)
                st = [EPGobj.mLogitMat];
                dpf = [st.dev]; % deviations
                [~,loc] = min(dpf); % the minimum deviation corresponds to the maximum likelihood
                loga = -EPGobj.mLgN(loc);
            end
            [~,n] = min(abs(EPGobj.mLgN + loga)); % the n whose gEUDs will be ploted
            euds = [EPGobj.mGrp.EUD]; euds = euds(n,:); % the gEUDs of that n
            if ~exist('numIntervals','var')
                numIntervals = 4;
            end

            % grouping patients
            flg=[EPGobj.mGrp.mFlgCensor]; % censor flags of patients
            [medianeud,numcomp,numtotal,betainv84,betainv16] = EventObserved(flg,euds,numIntervals);
            prob = numcomp./numtotal;
            % plot
            errorbar(medianeud,prob,max(0,prob-betainv16),max(0,betainv84-prob),strMarker,'LineWidth',lw);
            set(gca,'xminortick','on','yminortick','on');
            set(gca,'box','on');
            xlabel('gEUD'); ylabel('RP rate observed');
        end
        function fProbabilityFig_EUD(EPGobj)
            % map of the probability of have complication rate of at least value (e.g. 20%)
            if isempty(EPGobj.mBetaCumMat)
                disp('Empty member (BetaCumulativeMat), can not display its figure.'); return;
            end

            % prepare
            cm = colormap(jet(300)); cm=cm(1:256,:); %cm(end,:) = 0.5;
            imgmsk = EPGobj.mEUD_AtlasTotal > 0; imgmsk=imgmsk';
            eudmx = EPGobj.mEUD_DoseBins(end); eudmn = EPGobj.mEUD_DoseBins(1);

            % image data 
            img=1-EPGobj.mBetaCumMat; img=img';
            mx = ceil(max(img(imgmsk))*256);
            img(~imgmsk) = NaN;
            colormap(cm(1:mx,:));
            contourf(img); axis xy;
            colorbar;

%             set(gca,'Position',[0.05,0.05,0.8,0.9]);
            
            set(gca,'YTick',1:length(EPGobj.mLgN)); set(gca,'YTickLabel',-EPGobj.mLgN);
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
%         set(gca,'YTick',1:length(CGobjs(k).mLgN)); set(gca,'YTickLabel',CGobjs(k).mLgN);

%             img=1-EPGobj.BetaCumulativeMat'; img(1,end)=1; img(2,end)=0;
%             contourf(img); %contourf(flipud(rot90(EPGobj.BetaCumulativeMat,1)));
%             xsteps=0:EPGobj.GyStep:max(EPGobj.EUD(:)); xtickvec=10:10:xsteps(end); xtickstep=(size(EPGobj.BetaCumulativeMat,1)/xsteps(end)*10);
%             set(gca,'XTick',xtickstep:xtickstep:size(EPGobj.BetaCumulativeMat,1)); set(gca,'XTickLabel',xtickvec);
%             ytickvec=EPGobj.log10n; ytickstep=size(EPGobj.BetaCumulativeMat,2)/length(ytickvec);
%             set(gca,'YTick',ytickstep:ytickstep:size(EPGobj.BetaCumulativeMat,2)); set(gca,'YTickLabel',ytickvec);
%             %             pmin=min(EPGobj.BetaCumulativeMat(:)); pmax=max(EPGobj.BetaCumulativeMat(:)); %bstep=(pmax-pmin)/10*256;
%             colorbar; %colorbar('YTickLabel',round([pmin:(pmax-pmin)/9:pmax]*10)/10);
%             title([EPGobj.xlsSheet,', the probability that observed complications arise from true rate > 20%']);
        end
        function fLow68pctConfidenceFig_EUD(EPGobj)
            % map of lower 68% confidence limit on complication probability
            if isempty(EPGobj.mBetaCumMat)
                disp('Empty member (BetaCumulativeMat), can not display its figure.'); return;
            end

            % prepare
            cm = colormap(jet(300)); cm=cm(1:256,:); %cm(end,:) = 0.5;
            imgmsk = EPGobj.mEUD_AtlasTotal > 0; imgmsk=imgmsk';
            eudmx = EPGobj.mEUD_DoseBins(end); eudmn = EPGobj.mEUD_DoseBins(1);

            img=EPGobj.mBetaInvMat';
            mx = ceil(max(img(imgmsk))*256);
            img(~imgmsk) = NaN;
            colormap(cm(1:mx,:));
            colormap(cm);
            contourf(img); axis xy;
            colorbar;

%             set(gca,'Position',[0.05,0.05,0.8,0.9]);
           
            set(gca,'YTick',1:length(EPGobj.mLgN)); set(gca,'YTickLabel',-EPGobj.mLgN);
            xlim = get(gca,'XLim'); % x limit
            xticklabel = 0:5:eudmx; xtick = diff(xlim)/(eudmx-eudmn)*xticklabel+xlim(1);
            set(gca,'XTick',xtick);
            set(gca,'XTickLabel',xticklabel);
            
            set(gca,'xminortick','on','yminortick','on');
            set(gca,'box','on');
            xlabel('EUD doses (Gy)'); ylabel('log a');
            
%             if isempty(EPGobj.BetaInverseMat)
%                 disp('Empty member (BetaInverseMat), can not display its figure.'); return;
%             end
%             img=EPGobj.BetaInverseMat'; img(1,end)=1; img(2,end)=0;
%             contourf(img); %contourf(flipud(rot90(EPGobj.BetaInverseMat,1)));
%             xsteps=0:EPGobj.GyStep:max(EPGobj.EUD(:)); xtickvec=10:10:xsteps(end); xtickstep=(size(EPGobj.BetaInverseMat,1)/xsteps(end)*10);
%             set(gca,'XTick',xtickstep:xtickstep:size(EPGobj.BetaInverseMat,1)); set(gca,'XTickLabel',xtickvec);
%             ytickvec=EPGobj.log10n; ytickstep=size(EPGobj.BetaInverseMat,2)/length(ytickvec);
%             set(gca,'YTick',ytickstep:ytickstep:size(EPGobj.BetaInverseMat,2)); set(gca,'YTickLabel',ytickvec);
%             colorbar;
%             title([EPGobj.xlsSheet,', lower 68% confidence limit on complication probability']);
        end
    end
    
    methods(Static)
        function EPGobj = fCombineAtlas_EUD(CGobj1,CGobj2)
            % integrity check
            % mLgN
            if any(abs(CGobj1.mLgN-CGobj2.mLgN)>1e-6)
                error('The log10(n) is not equal in the two CGobjs so they can not be combined');
            end
            % dose bins in atlas computation
            if max(CGobj1.mEUD_DoseBins) >= max(CGobj2.mEUD_DoseBins)
                dosebins1 = CGobj1.mEUD_DoseBins;
                dosebins2 = CGobj2.mEUD_DoseBins;
            else
                dosebins1 = CGobj2.mEUD_DoseBins;
                dosebins2 = CGobj1.mEUD_DoseBins;
            end
            d = abs( dosebins1(1:size(dosebins2,1)) - dosebins2 );
            if any(d>1e-6)
                error('The dose bins in atlas computation were not consistant so they are not combined');
            end
            
            % combination
            EPGobj = CGobj1;
            EPGobj = EPGobj.AddPatient(CGobj2.mGrp);
            
            EPGobj.mEUD_DoseBins = dosebins1;
            size1 = size(CGobj1.mEUD_AtlasTotal); size2 = size(CGobj2.mEUD_AtlasTotal);
            EPGobj.mEUD_AtlasTotal = zeros(max(size1,size2));
            EPGobj.mEUD_AtlasTotal(1:size1(1),1:size1(2)) = CGobj1.mEUD_AtlasTotal;
            EPGobj.mEUD_AtlasTotal(1:size2(1),1:size2(2)) = EPGobj.mEUD_AtlasTotal(1:size2(1),1:size2(2))+CGobj2.mEUD_AtlasTotal;
            EPGobj.mEUD_AtlasComp = zeros(max(size1,size2));
            EPGobj.mEUD_AtlasComp(1:size1(1),1:size1(2)) = CGobj1.mEUD_AtlasComp;
            EPGobj.mEUD_AtlasComp(1:size2(1),1:size2(2)) = EPGobj.mEUD_AtlasComp(1:size2(1),1:size2(2))+CGobj2.mEUD_AtlasComp;
            
%             % compute probabilities
%             EPGobj=EPGobj.BetaCumulativeProbability_EUD();
%             EPGobj=EPGobj.BetaInverseProbability_EUD();
%             EPGobj=EPGobj.LogitAnalysis_EUD();
        end
    end
    
    methods % unused
        function EPGobj = fLymanAnalysisFittingExact_EUD(EPGobj)
            % preparation
            fttp = fittype('normcdf((eud-TD50)/(m*TD50),0,1)', 'independent','eud', 'coefficients',{'TD50','m'});
            fttp = fittype('normcdf(eud,TD50,m)', 'independent','eud', 'coefficients',{'TD50','m'});
            ftopt = fitoptions('method','LinearLeastSquares','algorith','Levenberg-Marquardt','Display','on');
            euds = [EPGobj.mGrp.EUD]';
            flg = double([EPGobj.mGrp.mFlgCensor]');

%             warning('off','MATLAB:singularMatrix');
            % for each mLgN, compute the parameters of Lyman model, TD50 and m
            kk = 1; % the first fit
            [dose,indx] = sort(euds(:,kk),'ascend');
            comp = ~flg(indx); f = find(comp); % rearrange the complication info so it corresponds to the dose
            comp = cumsum(comp)/EPGobj.mNumInGrp; % cumulative complication for dose
            % found the best grid as start point
                llh = EPGobj.mLymanGrid.loglikelihood(:,:,kk); %log likeli-hood from grid computation
                mn = max(llh(:)); [dd,mm]=find(llh==mn);
            [fittedobj,goodness,output] = fit(dose,comp,fttp,fitoptions(ftopt,'StartPoint',[EPGobj.mLymanGrid.TD50(dd(1)),EPGobj.mLymanGrid.m(mm(1))]+10));
            [fittedobj,goodness,output] = fit(euds(:,:,kk),flg,fttp,fitoptions(ftopt,'StartPoint',[EPGobj.mLymanGrid.TD50(dd(1)),EPGobj.mLymanGrid.m(mm(1))]+10));
            [fittedobj,goodness,output] = fit(euds(:,:,kk),flg,fttp,'StartPoint',[EPGobj.mLymanGrid.TD50(dd(1)),EPGobj.mLymanGrid.m(mm(1))]+10);
            fitobj = cell(size(EPGobj.mLgN));
            fitobj{kk} = fittedobj;
            goodness = repmat(goodness,size(EPGobj.mLgN));
            output = repmat(output,size(EPGobj.mLgN));
            for kk = 2:size(EPGobj.mLgN,1)
                llh = EPGobj.mLymanGrid.loglikelihood(:,:,kk); %log likeli-hood from grid computation
                mn = max(llh(:)); [dd,mm]=find(llh==mn);

                [fittedobj,goodness(kk),output(kk)] = fit(euds(:,kk),flg,fttp,fitoptions(ftopt,'StartPoint',[EPGobj.mLymanGrid.TD50(dd(1)),EPGobj.mLymanGrid.m(mm(1))]));
                fittedobj
                fitobj{kk} = fittedobj;
            end
%             warning('on','MATLAB:singularMatrix');
            EPGobj.mLyman = struct('FitObj',fitobj,'Goodness',goodness,'Output',output);
        end
        function fWriteXls_EUD(EPGobj) % write results from EUD into a spread sheet
            warning('off','MATLAB:xlswrite:AddSheet');
            if EPGobj.flgColOutput
                % EUD
                xlswrite(EPGobj.xlsFile_output,EPGobj.log10n',strcat(EPGobj.xlsSheet,'_EUD'),'B1');
                if ~isempty(EPGobj.PatientCode)
                    xlswrite(EPGobj.xlsFile_output,EPGobj.PatientCode(EPGobj.PatientRows),strcat(EPGobj.xlsSheet,'_EUD'),'A2');
                end
                if ~isempty(EPGobj.EUD)
                    xlswrite(EPGobj.xlsFile_output,EPGobj.EUD(EPGobj.PatientRows,:),strcat(EPGobj.xlsSheet,'_EUD'),'B2');
                end
                
                % total patients at log10n and dose
                xlswrite(EPGobj.xlsFile_output,EPGobj.log10n', strcat(EPGobj.xlsSheet,'_Total'),'B1');
                if ~isempty(EPGobj.mEUD_AtlasTotal);
%                     doses = (0:size(EPGobj.mEUD_AtlasTotal,1)-1)*EPGobj.GyStep;
                    xlswrite(EPGobj.xlsFile_output,EPGobj.AtlasBins,strcat(EPGobj.xlsSheet,'_Total'),'A2');
                    xlswrite(EPGobj.xlsFile_output,EPGobj.mEUD_AtlasTotal,strcat(EPGobj.xlsSheet,'_Total'),'B2');
                end
                
                % complication patients at log10n and dose
                xlswrite(EPGobj.xlsFile_output,EPGobj.log10n', strcat(EPGobj.xlsSheet,'_Comp'),'B1');
                if ~isempty(EPGobj.mEUD_AtlasComp)
%                     doses = (0:size(EPGobj.mEUD_AtlasComp,1)-1)*EPGobj.GyStep;
                    xlswrite(EPGobj.xlsFile_output,EPGobj.AtlasBins,strcat(EPGobj.xlsSheet,'_Comp'),'A2');
                    xlswrite(EPGobj.xlsFile_output,EPGobj.mEUD_AtlasComp,strcat(EPGobj.xlsSheet,'_Comp'),'B2');
                end
                
                % Beta probability
                if ~isempty(EPGobj.BetaCumulativeMat)
%                     doses = (0:size(EPGobj.mEUD_AtlasComp,1)-1)*EPGobj.GyStep;
                    for k=1:length(EPGobj.mBetaCumTh)
                        xlswrite(EPGobj.xlsFile_output, EPGobj.log10n', strcat(EPGobj.xlsSheet,'_prob_',num2str(EPGobj.mBetaCumTh(k))),'B1');
                        xlswrite(EPGobj.xlsFile_output,EPGobj.AtlasBins,strcat(EPGobj.xlsSheet,'_prob_',num2str(EPGobj.mBetaCumTh(k))),'A2');
                        xlswrite(EPGobj.xlsFile_output,EPGobj.BetaCumulativeMat(:,:,k),strcat(EPGobj.xlsSheet,'_prob_',num2str(EPGobj.mBetaCumTh(k))),'B2');
                    end
                end
                if ~isempty(EPGobj.BetaInverseMat)
%                     doses = (0:size(EPGobj.mEUD_AtlasComp,1)-1)*EPGobj.GyStep;
                    for k=1:length(EPGobj.mBetaInvTh)
                        xlswrite(EPGobj.xlsFile_output,EPGobj.log10n', strcat(EPGobj.xlsSheet,'_Low_',num2str(EPGobj.mBetaInvTh(k))),'B1');
                        xlswrite(EPGobj.xlsFile_output,EPGobj.AtlasBins,strcat(EPGobj.xlsSheet,'_Low_',num2str(EPGobj.mBetaInvTh(k))),'A2');
                        xlswrite(EPGobj.xlsFile_output,EPGobj.BetaCumulativeMat(:,:,k),strcat(EPGobj.xlsSheet,'_Low_',num2str(EPGobj.mBetaInvTh(k))),'B2');
                    end
                end
            else
                % EUD
                xlswrite(EPGobj.xlsFile_output, EPGobj.log10n, strcat(EPGobj.xlsSheet,'_EUD'),'A2');
                if ~isempty(EPGobj.PatientCode)
                    xlswrite(EPGobj.xlsFile_output,EPGobj.PatientCode(EPGobj.PatientRows)',strcat(EPGobj.xlsSheet,'_EUD'),'B1');
                end
                if ~isempty(EPGobj.EUD)
                    xlswrite(EPGobj.xlsFile_output,EPGobj.EUD(EPGobj.PatientRows,:)',strcat(EPGobj.xlsSheet,'_EUD'),'B2');
                end
                
                % total patients at log10n and dose
                xlswrite(EPGobj.xlsFile_output, EPGobj.log10n, strcat(EPGobj.xlsSheet,'_Total'),'A2');
                if ~isempty(EPGobj.mEUD_AtlasTotal);
%                     doses = (0:size(EPGobj.mEUD_AtlasTotal,1)-1)*EPGobj.GyStep;
                    xlswrite(EPGobj.xlsFile_output,EPGobj.AtlasBins',strcat(EPGobj.xlsSheet,'_Total'),'B1');
                    xlswrite(EPGobj.xlsFile_output,EPGobj.mEUD_AtlasTotal',strcat(EPGobj.xlsSheet,'_Total'),'B2');
                end
                
                % complication patients at log10n and dose
                xlswrite(EPGobj.xlsFile_output, EPGobj.log10n, strcat(EPGobj.xlsSheet,'_Comp'),'A2');
                if ~isempty(EPGobj.mEUD_AtlasComp)
%                     doses = (0:size(EPGobj.mEUD_AtlasComp,1)-1)*EPGobj.GyStep;
                    xlswrite(EPGobj.xlsFile_output,EPGobj.AtlasBins',strcat(EPGobj.xlsSheet,'_Comp'),'B1');
                    xlswrite(EPGobj.xlsFile_output,EPGobj.mEUD_AtlasComp',strcat(EPGobj.xlsSheet,'_Comp'),'B2');
                end
                
                % Beta probability
                if ~isempty(EPGobj.BetaCumulativeMat)
%                     doses = (0:size(EPGobj.mEUD_AtlasComp,1)-1)*EPGobj.GyStep;
                    for k=1:length(EPGobj.mBetaCumTh)
                        xlswrite(EPGobj.xlsFile_output, EPGobj.log10n, strcat(EPGobj.xlsSheet,'_prob_',num2str(EPGobj.mBetaCumTh(k))),'A2');
                        xlswrite(EPGobj.xlsFile_output,EPGobj.AtlasBins',strcat(EPGobj.xlsSheet,'_prob_',num2str(EPGobj.mBetaCumTh(k))),'B1');
                        xlswrite(EPGobj.xlsFile_output,EPGobj.BetaCumulativeMat(:,:,k)',strcat(EPGobj.xlsSheet,'_prob_',num2str(EPGobj.mBetaCumTh(k))),'B2');
                    end
                end
                if ~isempty(EPGobj.BetaInverseMat)
%                     doses = (0:size(EPGobj.mEUD_AtlasComp,1)-1)*EPGobj.GyStep;
                    for k=1:length(EPGobj.mBetaInvTh)
                        xlswrite(EPGobj.xlsFile_output, EPGobj.log10n, strcat(EPGobj.xlsSheet,'_Low_',num2str(EPGobj.mBetaInvTh(k))),'A2');
                        xlswrite(EPGobj.xlsFile_output,EPGobj.AtlasBins',strcat(EPGobj.xlsSheet,'_Low_',num2str(EPGobj.mBetaInvTh(k))),'B1');
                        xlswrite(EPGobj.xlsFile_output,EPGobj.BetaInverseMat(:,:,k)',strcat(EPGobj.xlsSheet,'_Low_',num2str(EPGobj.mBetaInvTh(k))),'B2');
                    end
                end
            end
            warning('on','MATLAB:xlswrite:AddSheet');
        end
    end
end