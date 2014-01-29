classdef classOutcomeAnalysis
    properties
        % general
        mGrp = classOutcomeIndividual.empty(1,0); % patient group of classOutcomeIndividual objects
        mNumInGrp = 0; % number of patients in the grou

        mBeta2Alpha = 0; % beta to alpha ratio
        mLymanN

        mStepDose % dose step size for computation based on DVH (in Gy)
        mStepVol % volume step size for computation based on DVH (in cc)
        mStepTime % time step size in month for dynamic analysis
        mBinsDose % the dose bins can be generated from mStepDose, or overwritten in case irregular bin steps were used
        mBinsVol % volume bins
        mBinsTime % the time bins can be overwritten in case irregular bins steps were used


        % Universal Survival Curve
        mUscAlpha
        mUscRangeD0
        mUscRangeDq
        mUscLogLikelihoods
        mUscDt
        
        mAtlasTotal_DVH % matrix for total patients computed from DVH
        mAtlasComp_DVH % matrix for patients with complications computed from DVH
        
        mAtlasTotal % matrix for total patients computed from DVH
        mAtlasComp % matrix for patients with complications computed from DVH
        mBetaCumulativeMat % based on EUD, the probability that the complication rate is larger than a specific number mBetaCumulativeThreshold
        mBetaCumulativeThreshold = 0.2;
        mBetaInverseMat % based on EUD, The chi-square where the probability is mBetaInverseThreshold
        mBetaInverseThreshold = 0.16;

        mKaplanMeierCompMat % patient complication Kaplan Meier curve at each atlas grid point
        mKaplanMeierCompSample % sample of patient complication curve at specific time point at each atlas grid point
        mKaplanMeierCompFromAtlas % patient complication info at binned time point in atlas form.

        mKaplanMeierCompOverall = classKaplanMeierCurve.empty(1,0); % overall complication curve
        mKaplanMeierSurvivalOverall = classKaplanMeierCurve.empty(1,0); % overall survival curve
        mKaplanMeierRelapseFree = classKaplanMeierCurve.empty(1,0); % relapse free survival curve
        
        mCoxMinSize = 2; % minimum sample size for Cox model
        mLogRankMinSize =2; % minimum sample size for log rank test
        mCoxParameter % cox model pairs (cell collumn 1 -- label of the Cox model fitting, cell collumn 2 -- fitting result)
        mLogRank % log rank test for (Di,Vj) (cell 1 -- label of the test (Dx, Vx, DVx), cell 2 -- n1, censor1, n2, censor2, p-value for each grid point

        mLogisticRegressionMat % logistic regression model using exact EUD
        mLogisticRegressionMatBin % logistic regression model using binned EUD

        mLogisticRegressionHosmerLemeshow = struct('n',0,'Chi2',[],'df',[],'p_value',[]); % goodness of fit for logistic regression (binary only) using Hosmer Lemshow test, ...
        % structure: n (number of risk groups), Chi2 (test result which is chi-square), df (degree of freedom, usually n-2), p-value (p-value of the chi-square)
        mLogisticRegressionHosmerLemeshowBin = struct('n',0,'Chi2',[],'df',[],'p_value',[]);
        mLogisticRegressionGTest = struct('n',0,'Chi2',[],'df',[],'p_value',[]); % goodness of fit for logistic regression (binary only) using G-test, ...
        % structure: n (number of risk groups), Chi2 (test result which is chi-square), df (degree of freedom, usually n-2), p-value (p-value of the chi-square)
        mLogisticRegressionGTestBin = struct('n',0,'Chi2',[],'df',[],'p_value',[]); % goodness of fit for logistic regression (binary only), ...
        mLogisticRegressionPearson = struct('n',0,'Chi2',[],'df',[],'p_value',[]); % goodness of fit for logistic regression (binary only) using Pearson test, ...
        % structure: n (number of risk groups), Chi2 (test result which is chi-square), df (degree of freedom, usually n-2), p-value (p-value of the chi-square)
        mLogisticRegressionPearsonBin = struct('n',0,'Chi2',[],'df',[],'p_value',[]); % goodness of fit for logistic regression (binary only), ...

        mLogisticRegressionGridBetaRange % beta ranges for grid searching
        mLogisticRegressionGrid % logistic regression model using exact EUD for grid
        mLogisticRegressionGridBin % logistic regression model using binned EUD for grid
        mLogisticRegressionGridMixtureModel % logistic regression model using exact EUD for grid
        mLogisticRegressionGoodnessOfFitSim = struct('SSRSim',[],'SSRObserve',[],'p_value',[]);  % goodness of fit using simulations
        mLogisticRegressionGoodnessOfFitSimBin = struct('SSRSim',[],'SSRObserve',[],'p_value',[]); % goodness of fit using simulations

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
        mLymanGoodnessOfFitSim = struct('SSRSim',[],'SSRObserve',[],'p_value',[]);  % goodness of fit using simulations
        mLymanGoodnessOfFitSimBin = struct('SSRSim',[],'SSRObserve',[],'p_value',[]); % goodness of fit using simulations
    end

    methods % set member vaulues
        function OCobj = classOutcomeAnalysis()
            %             OCobj.mGrp = classOutcomeIndividual.empty(1,0);
        end
        function OCobj = set.mBeta2Alpha(OCobj,mBeta2Alpha)
            if any (mBeta2Alpha < 0)
                error('Beta to Alpha ratio is negative, in set.mBeta2Alpha of classOutcomeAnalysis');
            end
            if length(mBeta2Alpha)>1
                warning('classOutcomeAnalysis:mBeta2Alpha','Parameter for mBeta2Alpha not a scalar, choose the first element');
            end
            for m = 1:OCobj.mNumInGrp
                OCobj.mGrp(m).mBeta2Alpha = mBeta2Alpha;
            end
            OCobj.mBeta2Alpha = mBeta2Alpha;
        end
        function OCobj = set.mBinsDose(OCobj,dosebins)
            if isempty(dosebins)
                OCobj.mBinsDose = dosebins;
                return;
            end

            euds = [OCobj.mGrp.mEUD];
            if ~isempty(euds)
                if max(euds(:)) > dosebins(end)
                    disp('assigned dose bins can not cover the maximum EUD');
                end
            end
            
            if any(diff(dosebins)<=0)
                disp('dose bins is not monotonic increasing when assigned to mBinsDose in classOutcomeAnalysis');
                disp('thus they are sorted');
                dosebins = sort( dosebins );
                dosebins(diff(dosebins)==0)=[];
            end
            OCobj.mBinsDose = dosebins(:);
        end
        function OCobj = set.mBinsVol(OCobj,volbins)
            if any(diff(volbins)<=0)
                disp('volume bins is not monotonic increasing when assigned to mBinsVol in classOutcomeAnalysis');
                disp('thus they are sorted');
                volbins = sort( volbins );
                volbins(diff(bolbins)==0)=[];
            end
            OCobj.mBinsVol = volbins(:);
        end
        function OCobj = set.mBinsTime(OCobj,timebins)
            if any( diff( timebins ) <=0 )
                disp('time bins is not strictly monotonic increasing when assign to mBinsTime in classOutcomeAnalysis');
                disp('thus they are sorted');
                timebins = sort(timebins);
                timebins(diff(timebins)==0)=[];
            end
            OCobj.mBinsTime = timebins(:);
        end
        function OCobj = set.mBetaCumulativeThreshold(OCobj, betacumthreshold)
            OCobj.mBetaCumulativeThreshold = betacumthreshold(:);
        end
        function OCobj = set.mBetaInverseThreshold(OCobj,betainvthreshold)
            OCobj.mBetaInverseThreshold = betainvthreshold(:);
        end
    end

    methods % patient operations like adding and removing
        function OCobj = fSetBeta2AlphaStatic(OCobj,beta2alpha)
            OCobj.mBeta2Alpha = beta2alpha;
        end
        function OCobj = fAddPatient(OCobj,ptobjs)
            % check if the new patient info matches the class definition
            if ~all(isa(ptobjs,'classOutcomeIndividual'))
                disp('Not an instance of patient individual class when adding a new patient');
                return;
            end
            
            % make sure the log10(n) and the mBeta2Alpha are the same for all patients
            for m = 1:length(ptobjs)
                ptobjs(m).mLymanN = OCobj.mLymanN;
                %ptobjs(m).mBeta2Alpha = OCobj.mBeta2Alpha;
            end

            % for empty patient object, add the patient directly
            if isempty([OCobj.mGrp.mID])
                OCobj.mGrp = ptobjs(:);
                OCobj.mNumInGrp = length(OCobj.mGrp);
                return;
            end
            
            % update existing patients in the group
            [fg,g] = ismember( {ptobjs.mID}, {OCobj.mGrp.mID} );
            f=find(g);
            for m = 1:length(f)
                OCobj.mGrp(g(f(m))) = ptobjs(f(m));
            end
            
            % add new patient
            OCobj.mNumInGrp = OCobj.mNumInGrp + sum(~fg);
            OCobj.mGrp(end+1:OCobj.mNumInGrp) = ptobjs(~fg);
        end
        function OCobj = fRemovePatient(OCobj,idx)
            if ~exist('idx','var') % if idx is not passed on, remove all patients (default)
                idx = true(OCobj.mNumInGrp,1);
            end
            OCobj.mGrp(idx) = [];
            OCobj.mNumInGrp = length(OCobj.mGrp);
        end
        function OCobj = fLinearQuartraticCorrection(OCobj)
            for k = 1:OCobj.mNumInGrp
                OCobj.mGrp(k) = OCobj.mGrp(k).fLinearQuartraticCorrection();
            end
        end
        
        function [OCobj, frac_lq,frac_full_lq] = fUSCCorrection(OCobj,dq,alfad0,a2b)
            tot_num_lq=0;
            tot_num_bins=0;
            tot_num_full_lq=0;
            for k = 1:OCobj.mNumInGrp
                [OCobj.mGrp(k), num_lq,num_bins,full_lq] = OCobj.mGrp(k).fUSCCorrection(dq,alfad0,a2b);
                tot_num_lq=num_lq+tot_num_lq;
                tot_num_bins=num_bins+tot_num_bins;
                tot_num_full_lq=full_lq+tot_num_full_lq;
            end
            frac_lq =tot_num_lq/tot_num_bins;
            frac_full_lq = tot_num_full_lq/(OCobj.mNumInGrp);
        end
        function flg = fPatientsWithComplicationData(OCobj)
            pt = OCobj.mGrp;
            f1 = cellfun(@(x) isempty(x)|isinf(x)|isnan(x),{pt.mDateBaseline}); % patients with no baseline date
            f2 = cellfun(@(x) isempty(x)|isinf(x)|isnan(x),{pt.mDateComp}); % patients with no complication date
            flg = find(~f2);  % patients with complciation date were not censored, so this information is restated
            for k = 1:length(flg)
                pt(flg(k)).mFlgCensor = 0;
            end
            
            f3 = cellfun(@(x) isempty(x)|isinf(x)|isnan(x),{pt.mDateLastFollowup}); % patients with no last follow up date
            f4 = cellfun(@(x) isempty(x)|isinf(x)|isnan(x),{pt.mFlgCensor}); % patients with no censor info
            flg = f1 | (f2&f3) | f4; % patients missing at least one data
            flg = ~flg; % patients with all data
        end
    end

    methods % general operations
        function OCobj = fOverallCompCurve(OCobj)
            % select patients with data
            f = OCobj.fPatientsWithComplicationData();
            pt = OCobj.mGrp(f);
            num = sum(f); % number of patients

            % complication time
            f2 = ~cellfun(@(x) isempty(x)|isinf(x)|isnan(x),{pt.mDateComp}); % patients with complication date
            f3 = ~cellfun(@(x) isempty(x)|isinf(x)|isnan(x),{pt.mDateLastFollowup}); % patients with last follow up date
            complicationdays = inf(num,1);
            lastfollowup = inf(num,1);
            complicationdays(f2) = ([pt(f2).mDateComp] - [pt(f2).mDateBaseline])' / 30;
            lastfollowup(f3) = ([pt(f3).mDateLastFollowup] - [pt(f3).mDateBaseline])' / 30;
            lastfollowup = min( lastfollowup, complicationdays );
            flgcensor = [pt.mFlgCensor]';

            % overall complication curve
            sa = classKaplanMeierCurve();
            sa.mSurvivalTime = {lastfollowup};
            sa.mFlgCensor = {flgcensor};
            sa = sa.fCalculateSurvivalCurve();

            OCobj.mKaplanMeierCompOverall = sa;
        end
        function OCobj = fOverallSurviveCurve(OCobj)
            % select patients with data
            pt = OCobj.mGrp;
            f1 = cellfun(@(x) isempty(x)|isinf(x)|isnan(x),{pt.mDateBaseline}); % patients with no baseline date
            f2 = cellfun(@(x) isempty(x)|isinf(x)|isnan(x),{pt.mDateDeath}); % patients with no death date
            f3 = cellfun(@(x) isempty(x)|isinf(x)|isnan(x),{pt.mDateLastFollowup}); % patients with no last follow-up date
            f4 = cellfun(@(x) isempty(x)|isinf(x)|isnan(x),{pt.mDateComp}); % patients with no complication date
            flg = f1 | (f2 & f3 & f4); % patient with out death date or complication data, so shall not be included in survival curve
            pt(flg) = [];

            num = sum(~flg); % number of patients

            % survival time
            f2 = ~cellfun(@(x) isempty(x)|isinf(x)|isnan(x),{pt.mDateDeath}); % patients with death date
            f3 = ~cellfun(@(x) isempty(x)|isinf(x)|isnan(x),{pt.mDateLastFollowup}); % patients with last follow up date
            f4 = ~cellfun(@(x) isempty(x)|isinf(x)|isnan(x),{pt.mDateComp}); % patients with complication date
            survivaltime = inf(num,1);
            survivaltime(f4) = ([pt(f4).mDateComp] - [pt(f4).mDateBaseline])' / 30;
            survivaltime(f3) = ([pt(f3).mDateLastFollowup] - [pt(f3).mDateBaseline])' / 30;
            survivaltime(f2) = ([pt(f2).mDateDeath] - [pt(f2).mDateBaseline])' / 30;
            
            flgcensor = ~f2'; %overall survival events is death
            

            % overall survival
            sa = classKaplanMeierCurve();
            sa.mSurvivalTime = {survivaltime};
            sa.mFlgCensor = {flgcensor};
            sa = sa.fCalculateSurvivalCurve();

            OCobj.mKaplanMeierSurvivalOverall = sa;
            
            % 2013_01_10 NOT calculating relapse for CW study
            % to implement, fix 'CW Pain recurrence date' in master list
            % to have 'date of survivial' or 'last follow-up date' when no
            % recurrence  (btw, only 2 recurrances)
            
            % free of relapse survival
%             f2 = ~cellfun(@(x) isempty(x)|isinf(x)|isnan(x),{pt.mDateRelapse}); % patients with relapse date
%             flgcensor(:) = 1;
%             flgcensor(f2) = 0;
%             survivaltime(f2) = ([pt(f2).mDateRelapse] - [pt(f2).mDateBaseline])'/30;
%             sa.mSurvivalTime = {survivaltime};
%             sa.mFlgCensor = {flgcensor};
%             sa = sa.fCalculateSurvivalCurve();

            OCobj.mKaplanMeierRelapseFree = sa;
        end
        
        function OCobj = fCalculateDoseBins(OCobj)
            doses = cellfun(@(x) x(end), {OCobj.mGrp.mDoseBins_LQ}); dmax = max(doses);
            step = dmax/200;
            OCobj.mBinsDose = ( 0 : step : dmax+step )';
            %OCobj.mBinsDose = ( 0 : OCobj.mStepDose : dmax+OCobj.mStepDose )';
        end
        function OCobj = fCalculateDoseBinsLog(OCobj)
            doses = cellfun(@(x) x(end), {OCobj.mGrp.mDoseBins_LQ}); dmax = log10(max(doses));
            doses = cellfun(@(x) x(1), {OCobj.mGrp.mDoseBins_LQ}); dmin = min(doses); dmin = log10(max(0.1,dmin));
            dosebinslog = dmin : OCobj.mStepDose : (dmax+OCobj.mStepDose);
            OCobj.mBinsDose = [0, 10.^dosebinslog]';
        end
        function OCobj = fCalculateVolBins(OCobj)
            vols = cellfun(@(x) x(1), {OCobj.mGrp.mVolCum}); vmax = max(vols);
            step = vmax/200;
            OCobj.mBinsVol = ( 0 : step : vmax+step )';
            
            %OCobj.mBinsVol = ( 0 : OCobj.mStepVol : vmax+OCobj.mStepVol )';
        end
        function OCobj = fCalculateVolBinsLog(OCobj)
            vols = cellfun(@(x) x(1), {OCobj.mGrp.mVolCum}); vmax = log10(max(vols));
            vols = cellfun(@(x) min(x(x>0)), {OCobj.mGrp.mVolCum}); vmin = min(vols); vmin = log10(max(0.1,vmin));
            vmin = fix(vmin/OCobj.mStepVol)*OCobj.mStepVol;
            volbinslog = vmin : OCobj.mStepVol : (vmax + OCobj.mStepVol);
            OCobj.mBinsVol = [0, 10.^volbinslog]';
        end
        function OCobj = fCalulateTimeBins(OCobj)
            timebins = ( [OCobj.mGrp.mDateLastFollowup] - [OCobj.mGrp.mDateBaseline] ) / 30; tmax = max(timebins);
            OCobj.mBinsTime = ( 0 : OCobj.mStepTime : tmax+OCobj.mStepTime )';
        end
    end

    methods % DVH operations
        function OCobj = fKaplanMeierCompMat_DVH(OCobj)
            % select patients with data
            f = OCobj.fPatientsWithComplicationData();
            CG = OCobj.fRemovePatient(~f);

            % complication time
            f2 = ~cellfun('isempty',{CG.mGrp.mDateComp}); % patients with no complication date
            f3 = ~cellfun('isempty',{CG.mGrp.mDateLastFollowup}); % patients with no last follow up date
            complicationdays = inf(CG.mNumInGrp,1);
            lastfollowup = inf(CG.mNumInGrp,1);
            complicationdays(f2) = ([CG.mGrp(f2).mDateComp] - [CG.mGrp(f2).mDateBaseline])' / 30;
            lastfollowup(f3) = ([CG.mGrp(f3).mDateLastFollowup] - [CG.mGrp(f3).mDateBaseline])' / 30;
            lastfollowup = min( lastfollowup, complicationdays );
            flgcensor = [CG.mGrp.mFlgCensor]';

            % populate to all grid points
            sa = classKaplanMeierCurve(); % a new clean object
            pse = repmat(sa, [length(CG.mBinsDose), length(CG.mBinsVol)]); % patient survival exact

            % Vx
            vols = zeros( CG.mNumInGrp, 1 );
            for ii = 1:length(CG.mBinsDose)
                % volume of each patient at current dose (x)
                vols(:) = 0;
                for jj = 1:CG.mNumInGrp
                    vols(jj) = CG.mGrp(jj).fVolAtDose( CG.mBinsDose(ii) );
                end
                vols(vols==0) = -1; % exclude zero volume patients
                
                % matrix at each (Di, Vj)
                for jj = 1:length(CG.mBinsVol) % for each volume point under dose x
                    f = vols >= CG.mBinsVol(jj); % patient at the grid point
                    if ~any(f) % no patient for the grid, clean its contents
                        continue;
                    end

                    % use clasSurvivalAnalysis to compute the Kelplan Meier curve
                    pse(ii,jj).mSurvivalTime = {lastfollowup(f)};
                    pse(ii,jj).mFlgCensor = {flgcensor(f)};
                    pse(ii,jj) = pse(ii,jj).fCalculateSurvivalCurve();
                end
            end
            OCobj.mKaplanMeierCompMat = pse;
        end
        function OCobj = fCoxModel_DVH(OCobj)
            % prepare
            if isempty(OCobj.mBinsVol)
                OCobj = OCobj.fCalculateVolBins();
            end
            if isempty(OCobj.mBinsDose)
                OCobj = OCobj.fCalculateDoseBins();
            end

            % select patients with data
            f = OCobj.fPatientsWithComplicationData();
            CG = OCobj.fRemovePatient(~f);
            
            % survival/complication time
            f2 = ~cellfun('isempty',{CG.mGrp.mDateComp}); % patients with no complication date
            f3 = ~cellfun('isempty',{CG.mGrp.mDateLastFollowup}); % patients with no last follow up date
            compdate = inf(CG.mNumInGrp,1);
            lastfollowup = inf(CG.mNumInGrp,1);
            compdate(f2) = ([CG.mGrp(f2).mDateComp] - [CG.mGrp(f2).mDateBaseline])' / 30;
            lastfollowup(f3) = ([CG.mGrp(f3).mDateLastFollowup] - [CG.mGrp(f3).mDateBaseline])' / 30;
            compdate = min( lastfollowup, compdate );
            flgcensor = [CG.mGrp.mFlgCensor]';

            warning('off');
            % dmax
                dmx = cellfun(@(x) x(end), {CG.mGrp.mDoseBins_LQ}');

                [~,logl,h,stats]=coxphfit(dmx,compdate,'baseline',0,'censoring',flgcensor);
                stats.logl=logl; stats.h=h;
                stats.data_exposure = dmx; stats.data_hazard = compdate;
                
                    if isempty(OCobj.mCoxParameter)
                            OCobj.mCoxParameter{end+1,1}='Dmax';
                            OCobj.mCoxParameter{end,2}=stats;
                    else
                        f = cellfun(@(x) strcmpi('Dmax',x),OCobj.mCoxParameter(:,1));
                        if any(f)
                            OCobj.mCoxParameter{f,2}=stats;
                        else
                            OCobj.mCoxParameter{end+1,1}='Dmax';
                            OCobj.mCoxParameter{end,2}=stats;
                        end
                    end

            % cm2cw
                cm2cw = [CG.mGrp.mDistanceToChestWall]';
                if ~isempty(cm2cw)
                [~,logl,h,stats]=coxphfit(cm2cw,compdate,'baseline',0,'censoring',flgcensor);
                stats.logl=logl; stats.h=h;
                stats.data_exposure = cm2cw; stats.data_hazard = compdate;
                
                    if isempty(OCobj.mCoxParameter)
                            OCobj.mCoxParameter{end+1,1}='CM2CW';
                            OCobj.mCoxParameter{end,2}=stats;
                    else
                        f = cellfun(@(x) strcmpi('CM2CW',x),OCobj.mCoxParameter(:,1));
                        if any(f)
                            OCobj.mCoxParameter{f,2}=stats;
                        else
                            OCobj.mCoxParameter{end+1,1}='CM2CW';
                            OCobj.mCoxParameter{end,2}=stats;
                        end
                    end
                end
            % fx
                fxnum=[CG.mGrp.mFxNum]'; % all fractions
                flgfx = length(unique(fxnum))>1; % if there is more than one fraction numbers, do the cox model

                %flgfx=false;
                if flgfx
                    [~,logl,h,stats]=coxphfit([CG.mGrp.mFxNum]',compdate,'baseline',0,'censoring',flgcensor);
                    stats.logl=logl; stats.h=h;
                    stats.data_exposure = [CG.mGrp.mFxNum]'; stats.data_hazard = compdate;
                    
                    if isempty(OCobj.mCoxParameter)
                            OCobj.mCoxParameter{end+1,1}='Fx';
                            OCobj.mCoxParameter{end,2}=stats;
                    else
                        f = cellfun(@(x) strcmpi('Fx',x),OCobj.mCoxParameter(:,1));
                        if any(f)
                            OCobj.mCoxParameter{f,2}=stats;
                        else
                            OCobj.mCoxParameter{end+1,1}='Fx';
                            OCobj.mCoxParameter{end,2}=stats;
                        end
                    end
                end
                
               
                
                % age
                age=[CG.mGrp.mAgeAtTx]';
                [~,logl,h,stats]=coxphfit(age,compdate,'baseline',0,'censoring',flgcensor);
                stats.logl=logl; stats.h=h;
                stats.data_exposure = age; stats.data_hazard = compdate;
                    
                if isempty(OCobj.mCoxParameter)
                        OCobj.mCoxParameter{end+1,1}='Age';
                        OCobj.mCoxParameter{end,2}=stats;
                else
                    f = cellfun(@(x) strcmpi('Age',x),OCobj.mCoxParameter(:,1));
                    if any(f)
                        OCobj.mCoxParameter{f,2}=stats;
                    else
                        OCobj.mCoxParameter{end+1,1}='Age';
                        OCobj.mCoxParameter{end,2}=stats;
                    end
                end
            
                 % gender
                gender=[CG.mGrp.mGender]';
                [~,logl,h,stats]=coxphfit(gender,compdate,'baseline',0,'censoring',flgcensor);
                stats.logl=logl; stats.h=h;
                stats.data_exposure = gender; stats.data_hazard = compdate;
                    
                if isempty(OCobj.mCoxParameter)
                        OCobj.mCoxParameter{end+1,1}='Gender';
                        OCobj.mCoxParameter{end,2}=stats;
                else
                    f = cellfun(@(x) strcmpi('Gender',x),OCobj.mCoxParameter(:,1));
                    if any(f)
                        OCobj.mCoxParameter{f,2}=stats;
                    else
                        OCobj.mCoxParameter{end+1,1}='Gender';
                        OCobj.mCoxParameter{end,2}=stats;
                    end
                end
                
                 % kps
                kps=[CG.mGrp.mKPS]';
                if ~isempty(kps)
                [~,logl,h,stats]=coxphfit(kps,compdate,'baseline',0,'censoring',flgcensor);
                stats.logl=logl; stats.h=h;
                stats.data_exposure = kps; stats.data_hazard = compdate;
                    
                if isempty(OCobj.mCoxParameter)
                        OCobj.mCoxParameter{end+1,1}='KPS';
                        OCobj.mCoxParameter{end,2}=stats;
                else
                    f = cellfun(@(x) strcmpi('KPS',x),OCobj.mCoxParameter(:,1));
                    if any(f)
                        OCobj.mCoxParameter{f,2}=stats;
                    else
                        OCobj.mCoxParameter{end+1,1}='KPS';
                        OCobj.mCoxParameter{end,2}=stats;
                    end
                end
                end
        
                  % dperfx
                DperFx=[CG.mGrp.mDosePerFx]';
                [~,logl,h,stats]=coxphfit(DperFx,compdate,'baseline',0,'censoring',flgcensor);
                stats.logl=logl; stats.h=h;
                stats.data_exposure = DperFx; stats.data_hazard = compdate;
                    
                if isempty(OCobj.mCoxParameter)
                        OCobj.mCoxParameter{end+1,1}='DperFx';
                        OCobj.mCoxParameter{end,2}=stats;
                else
                    f = cellfun(@(x) strcmpi('DperFx',x),OCobj.mCoxParameter(:,1));
                    if any(f)
                        OCobj.mCoxParameter{f,2}=stats;
                    else
                        OCobj.mCoxParameter{end+1,1}='DperFx';
                        OCobj.mCoxParameter{end,2}=stats;
                    end
                end
        
        % bmi
                bmi=[CG.mGrp.mBMI]';
                if ~isempty(bmi)
                % remove patients with no bmi info
                bmi_idx = bmi>0;
                [~,logl,h,stats]=coxphfit(bmi(bmi_idx),compdate(bmi_idx),'baseline',0,'censoring',flgcensor(bmi_idx));
                stats.logl=logl; stats.h=h;
                stats.data_exposure = bmi; stats.data_hazard = compdate;
                    
                if isempty(OCobj.mCoxParameter)
                        OCobj.mCoxParameter{end+1,1}='BMI';
                        OCobj.mCoxParameter{end,2}=stats;
                else
                    f = cellfun(@(x) strcmpi('BMI',x),OCobj.mCoxParameter(:,1));
                    if any(f)
                        OCobj.mCoxParameter{f,2}=stats;
                    else
                        OCobj.mCoxParameter{end+1,1}='BMI';
                        OCobj.mCoxParameter{end,2}=stats;
                    end
                end
                end
        
            % delivered dose
                tx = [CG.mGrp.mDoseTx]';
                [~,logl,h,stats]=coxphfit(tx,compdate,'baseline',0,'censoring',flgcensor);
                stats.logl=logl; stats.h=h;
                stats.data_exposure = tx; stats.data_hazard = compdate;

                    if isempty(OCobj.mCoxParameter)
                            OCobj.mCoxParameter{end+1,1}='Tx'; % dose treated
                            OCobj.mCoxParameter{end,2}=stats;
                    else
                        f = cellfun(@(x) strcmpi('Tx',x),OCobj.mCoxParameter(:,1));
                        if any(f)
                            OCobj.mCoxParameter{f,2}=stats;
                        else
                            OCobj.mCoxParameter{end+1,1}='Tx';
                            OCobj.mCoxParameter{end,2}=stats;
                        end
                    end

                % fx + prescription dose
                if flgfx
                    [~,logl,h,stats]=coxphfit([tx,fxnum],compdate,'baseline',0,'censoring',flgcensor);
                    stats.logl=logl; stats.h=h;
                    stats.data_exposure = [tx,fxnum]; stats.data_hazard = compdate;
                    
                        f = cellfun(@(x) strcmpi('TxFx',x),OCobj.mCoxParameter(:,1));
                        if any(f)
                            OCobj.mCoxParameter{f,2}=stats;
                        else
                            OCobj.mCoxParameter{end+1,1}='TxFx';
                            OCobj.mCoxParameter{end,2}=stats;
                        end
                end
               %euds
%                 euds=[CG.mGrp.mEUD]';
%                 eudnum = size(euds,2);
%                 CoxEUDs=repmat(stats,[eudnum,1]);
%                 for i=1:eudnum %loop over each eud value
%                     cur_euds =euds(:,i);
%                     
%                      try
%                         if length(unique(cur_euds))>1
%                             [~,logl,h,stats]=coxphfit(cur_euds,compdate,'baseline',0,'censoring',flgcensor);
%                             stats.logl=logl; stats.h=h;
%                         else
%                             stats.beta=-inf; stats.logl=-inf; stats.h=-inf;
%                         end
%                     catch
%                         stats.beta=-inf; stats.logl=-inf; stats.h=-inf;
%                     end
%                     stats.data_exposure = cur_euds; stats.data_hazard = compdate;
%                     CoxEUDs(i)=stats;
%                 end  
            % Dx
                %tmp
                
                %volnum=200;
                %vbins = [round(linspace(1,max(CG.mBinsVol),min(length(CG.mBinsVol),200)))];
                
                volnum=length(CG.mBinsVol);
                vbins = [round(linspace(1,max(CG.mBinsVol),volnum))];
                %vbins = [1:volnum];
                
                
                
                %volnum=min(length(CG.mBinsVol),200);
                %vbins = [round(linspace(1,max(CG.mBinsVol),volnum))];
        
                stats.bin=0;
                CoxDx=repmat(stats,[volnum,1]);
                CoxDVx=repmat(stats,[volnum,1]);
                CoxDxFx=repmat(stats,[volnum,1]);
                CoxDVxFx=repmat(stats,[volnum,1]);
                CoxDVxTx=repmat(stats,[volnum,1]);
                CoxDVxDperFx=repmat(stats,[volnum,1]);
                CoxDVxBMI=repmat(stats,[volnum,1]);
                CoxDVxCM2CW=repmat(stats,[volnum,1]);
                
                % check for Dx one by one
                dv = zeros(CG.mNumInGrp,1);
                for v=1:volnum
                   
                    % doses under v
                    dv(:) = 0;
                    for k = 1:CG.mNumInGrp
                        %dv(k) = CG.mGrp(k).fDoseAtVol( CG.mBinsVol(vbins(v)) );
                        dv(k) = CG.mGrp(k).fDoseAtVol( CG.mBinsVol(v) );
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
                    stats.bin=vbins(v);
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
                        stats.bin=vbins(v);
                        CoxDVxFx(v)=stats;
                    end
                    % Cox model for patients with non-zero data
%                     g=find(dv); % patients with non-zero doses
%                     
% %                     
% %                     if v>650 && v<700,
% %                             disp([num2str(v),...
% %                             ', length(g): ',num2str(length(g)),...
% %                             ' min size: ',num2str(CG.mCoxMinSize)]);
% %                     end
% %                    
%                     % check if sample size is too small
%                     if length(g)<CG.mCoxMinSize
%                         stats.beta=-inf; stats.logl=-inf; stats.h=-inf;
%                         stats.bin=-inf;
%                         CoxDx(v)=stats; CoxDxFx(v)=stats;
%                         continue;
%                     end
                     g=find(dv);
                      if length(g)<=(2*(CG.mCoxMinSize)+2)
                        stats.beta=-inf; stats.logl=-inf; stats.h=-inf;
                        stats.bin = -inf;
                        CoxDx(v)=stats; CoxDxFx(v)=stats;
                        continue;
                    end
                    
                       g2=find(dv.*~flgcensor); % non-zero complications dose cases
                    % check if sample size is too small
                    if length(g2)<CG.mCoxMinSize
                        stats.beta=-inf; stats.logl=-inf; stats.h=-inf;
                        stats.bin = -inf;
                        CoxDx(v)=stats; CoxDxFx(v)=stats;
                        continue;
                    end
                    
                    
                    % otherwise calcualte the cox model
                    % Dx
                    [~,logl,h,stats]=coxphfit(dv(g),compdate(g),'baseline',0,'censoring',flgcensor(g));
                    stats.logl=logl; stats.h=h;
                    stats.data_exposure = dv(g); stats.data_hazard = compdate(g);
                    stats.bin=vbins(v);
                    CoxDx(v)=stats;
                    
                    % Dx + fx
                    if flgfx
                        [~,logl,h,stats]=coxphfit([dv(g),fxnum(g)],compdate(g),'baseline',0,'censoring',flgcensor(g));
                        stats.logl=logl; stats.h=h;
                        stats.data_exposure = [dv(g),fxnum(g)]; stats.data_hazard = compdate(g);
                        stats.bin=vbins(v);
                        CoxDxFx(v)=stats;
                    end
                    
                    
                    %here
                     % Dx + Tx
                        try
                            [~,logl,h,stats]=coxphfit([dv,tx],compdate,'baseline',0,'censoring',flgcensor);
                            stats.logl=logl; stats.h=h;
                        catch
                            stats.beta=-inf; stats.logl=-inf; stats.h=-inf;
                        end
                        stats.data_exposure = [dv,tx]; stats.data_hazard = compdate;
                        stats.bin=vbins(v);
                        CoxDVxTx(v)=stats;
                   
                        % Dx + Dose/Frac
                        try
                            [~,logl,h,stats]=coxphfit([dv,DperFx],compdate,'baseline',0,'censoring',flgcensor);
                            stats.logl=logl; stats.h=h;
                        catch
                            stats.beta=-inf; stats.logl=-inf; stats.h=-inf;
                        end
                        stats.data_exposure = [dv,DperFx]; stats.data_hazard = compdate;
                        stats.bin=vbins(v);
                        CoxDVxDperFx(v)=stats;

                    
                      % Dx + BMI
                      if ~isempty(bmi)
                        try
                            [~,logl,h,stats]=coxphfit([dv(bmi_idx),bmi(bmi_idx)],compdate(bmi_idx),'baseline',0,'censoring',flgcensor(bmi_idx));
                            stats.logl=logl; stats.h=h;
                        catch
                            stats.beta=-inf; stats.logl=-inf; stats.h=-inf;
                        end
                        stats.data_exposure = [dv(bmi_idx),bmi(bmi_idx)]; stats.data_hazard = compdate(bmi_idx);
                        stats.bin=vbins(v);
                        CoxDVxBMI(v)=stats;
                      end
                    
                    % Vx + cm2cw
                    if ~isempty(cm2cw)
                    try
                            [~,logl,h,stats]=coxphfit([dv,cm2cw],compdate,'baseline',0,'censoring',flgcensor);
                            stats.logl=logl; stats.h=h;
                        catch
                            stats.beta=-inf; stats.logl=-inf; stats.h=-inf;
                        end
                        stats.data_exposure = [dv,cm2cw]; stats.data_hazard = compdate;
                        stats.bin=vbins(v);
                        CoxDVxCM2CW(v)=stats;
                    end                        
                end
                
                % save results
                if isempty(OCobj.mCoxParameter)
                    OCobj.mCoxParameter{end+1,1}='Dx';
                    OCobj.mCoxParameter{end,2}=CoxDx;
                else
                    f = cellfun(@(x) strcmpi('Dx',x),OCobj.mCoxParameter(:,1));
                    if any(f)
                        OCobj.mCoxParameter{f,2}=CoxDx;
                    else
                        OCobj.mCoxParameter{end+1,1}='Dx';
                        OCobj.mCoxParameter{end,2}=CoxDx;
                    end
                end
                
%                 f = cellfun(@(x) strcmpi('EUDs',x),OCobj.mCoxParameter(:,1));
%                 if any(f)
%                     OCobj.mCoxParameter{f,2}=CoxEUDs;
%                 else
%                     OCobj.mCoxParameter{end+1,1}='EUDs';
%                     OCobj.mCoxParameter{end,2}=CoxEUDs;
%                 end
                
                f = cellfun(@(x) strcmpi('DVx',x),OCobj.mCoxParameter(:,1));
                if any(f)
                    OCobj.mCoxParameter{f,2}=CoxDVx;
                else
                    OCobj.mCoxParameter{end+1,1}='DVx';
                    OCobj.mCoxParameter{end,2}=CoxDVx;
                end
                if flgfx
                    f = cellfun(@(x) strcmpi('DVxFx',x),OCobj.mCoxParameter(:,1));
                    if any(f)
                        OCobj.mCoxParameter{f,2}=CoxDVxFx;
                    else
                        OCobj.mCoxParameter{end+1,1}='DVxFx';
                        OCobj.mCoxParameter{end,2}=CoxDVxFx;
                    end
                    f = cellfun(@(x) strcmpi('DxFx',x),OCobj.mCoxParameter(:,1));
                    if any(f)
                        OCobj.mCoxParameter{f,2}=CoxDxFx;
                    else
                        OCobj.mCoxParameter{end+1,1}='DxFx';
                        OCobj.mCoxParameter{end,2}=CoxDxFx;
                    end
                
                    % VDxTx
                     f = cellfun(@(x) strcmpi('DVxTx',x),OCobj.mCoxParameter(:,1));
                    if any(f)
                        OCobj.mCoxParameter{f,2}=CoxDVxTx;
                    else
                        OCobj.mCoxParameter{end+1,1}='DVxTx';
                        OCobj.mCoxParameter{end,2}=CoxDVxTx;
                    end
                    
                      % DVxTx
                     f = cellfun(@(x) strcmpi('DVxDperFx',x),OCobj.mCoxParameter(:,1));
                    if any(f)
                        OCobj.mCoxParameter{f,2}=CoxDVxDperFx;
                    else
                        OCobj.mCoxParameter{end+1,1}='DVxDperFx';
                        OCobj.mCoxParameter{end,2}=CoxDVxDperFx;
                    end
                     % DVxBMI
                     f = cellfun(@(x) strcmpi('DVxBMI',x),OCobj.mCoxParameter(:,1));
                    if any(f)
                        OCobj.mCoxParameter{f,2}=CoxDVxBMI;
                    else
                        OCobj.mCoxParameter{end+1,1}='DVxBMI';
                        OCobj.mCoxParameter{end,2}=CoxDVxBMI;
                    end
                    % DVxCM2CW
                     f = cellfun(@(x) strcmpi('DVxCM2CW',x),OCobj.mCoxParameter(:,1));
                    if any(f)
                        OCobj.mCoxParameter{f,2}=CoxDVxCM2CW;
                    else
                        OCobj.mCoxParameter{end+1,1}='DVxCM2CW';
                        OCobj.mCoxParameter{end,2}=CoxDVxCM2CW;
                    end
                    end
            % Vx
                % tmp 2013_01_31 instead of 1 Gy bins
                % divide up  equally (say 200 bins)
                
                %dosenum=min(length(CG.mBinsDose),200);
                %
                dosenum=length(CG.mBinsDose);
                %dbins = [1:dosenum];
                dbins = [round(linspace(1,max(CG.mBinsDose),dosenum))];
                
                stats.bin=0;
                CoxVx=repmat(stats,[dosenum,1]);
                CoxVDx=repmat(stats,[dosenum,1]);
                CoxVxFx=repmat(stats,[dosenum,1]);
                CoxVDxFx=repmat(stats,[dosenum,1]);
                CoxVDxDMax=repmat(stats,[dosenum,1]);
                CoxVDxTx=repmat(stats,[dosenum,1]);
                CoxVDxDperFx=repmat(stats,[dosenum,1]);
                CoxVDxBMI=repmat(stats,[dosenum,1]);
                CoxVDxCM2CW=repmat(stats,[dosenum,1]);
                
                vd=zeros(CG.mNumInGrp,1);
                % check for Vx one by one
                for d=1:dosenum
                    
                   
                    
                    
                    % volumes under d
                    vd(:)=0;
                    for k=1:CG.mNumInGrp
                        vd(k) = CG.mGrp(k).fVolAtDose( CG.mBinsDose(d) );
                        %vd(k) = CG.mGrp(k).fVolAtDose( CG.mBinsDose(dbins(d)) );
                    end
                    
                    % Cox model for patients with non-zero data
                    g=find(vd); % non-zeros dose cases
%                     vd(vd==0) = -1; % exclude zero volume patients
                    % check if sample size is too small
                        
                        % length(g) <= (2*min+2) -> 6 for ILung(pultox) only
                        % length(g) < 2*min+1 works for everything else
                    if length(g)<=(2*(CG.mCoxMinSize)+2)
                        stats.beta=-inf; stats.logl=-inf; stats.h=-inf;
                        stats.bin = -inf;
                        CoxVx(d)=stats; CoxVxFx(d)=stats;
                        continue;
                    end
                    
                       g2=find(vd.*~flgcensor); % non-zero complications dose cases
                    % check if sample size is too small
                    if length(g2)<CG.mCoxMinSize
                        stats.beta=-inf; stats.logl=-inf; stats.h=-inf;
                        stats.bin = -inf;
                        CoxVx(d)=stats; CoxVxFx(d)=stats;
                        continue;
                    end
                  
                    
                    
                    % Cox model for all patients
                    % Vx
                    try
                        if length(unique(vd))>1
                             [~,logl,h,stats]=coxphfit(vd,compdate,'baseline',0,'censoring',flgcensor);
                            %[~,logl,h,stats]=coxphfit(vd,compdate,'censoring',flgcensor);
                            stats.logl=logl; stats.h=h;
                        else
                            stats.beta=-inf; stats.logl=-inf; stats.h=-inf;
                        end
                    catch
                        stats.beta=-inf; stats.logl=-inf; stats.h=-inf;
                    end
                    stats.data_exposure = vd; stats.data_hazard = compdate;
                    stats.bin = dbins(d);
                    %stats.bin = d;
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
                        stats.bin = dbins(d);
                        CoxVDxFx(d)=stats;
                    end
                    
  
                    % otherwise calcualte the cox model
                    % Vx
                    [~,logl,h,stats]=coxphfit(vd(g),compdate(g),'baseline',0,'censoring',flgcensor(g));
                    stats.logl=logl; stats.h=h;
                    stats.data_exposure = vd(g); stats.data_hazard = compdate(g);
                    stats.bin = dbins(d);
                    CoxVx(d)=stats;
                    
                   
                    
                    % Vx + fx
                    if flgfx
                        [~,logl,h,stats]=coxphfit([vd(g),fxnum(g)],compdate(g),'baseline',0,'censoring',flgcensor(g));
                        stats.logl=logl; stats.h=h;
                        stats.data_exposure = [vd(g),fxnum(g)]; stats.data_hazard = compdate(g);
                        stats.bin = dbins(d);
                        CoxVxFx(d)=stats;
                    end
                    
                     % Vx + dmax
                        try
                            [~,logl,h,stats]=coxphfit([vd,dmx],compdate,'baseline',0,'censoring',flgcensor);
                            stats.logl=logl; stats.h=h;
                        catch
                            stats.beta=-inf; stats.logl=-inf; stats.h=-inf;
                        end
                        stats.data_exposure = [vd,tx]; stats.data_hazard = compdate;
                        stats.bin = dbins(d);
                        CoxVDxDMax(d)=stats;

                    
                     % Vx + Tx
                        try
                            [~,logl,h,stats]=coxphfit([vd,tx],compdate,'baseline',0,'censoring',flgcensor);
                            stats.logl=logl; stats.h=h;
                        catch
                            stats.beta=-inf; stats.logl=-inf; stats.h=-inf;
                        end
                        stats.data_exposure = [vd,tx]; stats.data_hazard = compdate;
                        stats.bin = dbins(d);
                        CoxVDxTx(d)=stats;
                   
                        % Vx + Dose/Frac
                    if flgfx
                        try
                            [~,logl,h,stats]=coxphfit([vd,DperFx],compdate,'baseline',0,'censoring',flgcensor);
                            stats.logl=logl; stats.h=h;
                        catch
                            stats.beta=-inf; stats.logl=-inf; stats.h=-inf;
                        end
                        stats.data_exposure = [vd,DperFx]; stats.data_hazard = compdate;
                        stats.bin = dbins(d);
                        CoxVDxDperFx(d)=stats;
                    end
                
                      % Vx + BMI
                      if ~isempty(bmi)
                        try
                            [~,logl,h,stats]=coxphfit([vd(bmi_idx),bmi(bmi_idx)],compdate(bmi_idx),'baseline',0,'censoring',flgcensor(bmi_idx));
                            stats.logl=logl; stats.h=h;
                        catch
                            stats.beta=-inf; stats.logl=-inf; stats.h=-inf;
                        end
                        stats.data_exposure = [vd(bmi_idx),bmi(bmi_idx)]; stats.data_hazard = compdate(bmi_idx);
                        stats.bin = dbins(d);
                        CoxVDxBMI(d)=stats;
                
                      end
                      
                      if ~isempty(cm2cw)
                    % Vx + cm2cw
                        try
                            [~,logl,h,stats]=coxphfit([vd,cm2cw],compdate,'baseline',0,'censoring',flgcensor);
                            stats.logl=logl; stats.h=h;
                        catch
                            stats.beta=-inf; stats.logl=-inf; stats.h=-inf;
                        end
                        stats.data_exposure = [vd,cm2cw]; stats.data_hazard = compdate;
                        stats.bin = dbins(d);
                        CoxVDxCM2CW(d)=stats;
                    
                      end    
                end
                
                                
                
                % save results
                f = cellfun(@(x) strcmpi('VDx',x),OCobj.mCoxParameter(:,1));
                if any(f)
                    OCobj.mCoxParameter{f,2}=CoxVDx;
                else
                    OCobj.mCoxParameter{end+1,1}='VDx';
                    OCobj.mCoxParameter{end,2}=CoxVDx;
                end
                f = cellfun(@(x) strcmpi('Vx',x),OCobj.mCoxParameter(:,1));
                if any(f)
                    OCobj.mCoxParameter{f,2}=CoxVx;
                else
                    OCobj.mCoxParameter{end+1,1}='Vx';
                    OCobj.mCoxParameter{end,2}=CoxVx;
                end
                if flgfx
                    f = cellfun(@(x) strcmpi('VDxFx',x),OCobj.mCoxParameter(:,1));
                    if any(f)
                        OCobj.mCoxParameter{f,2}=CoxVDxFx;
                    else
                        OCobj.mCoxParameter{end+1,1}='VDxFx';
                        OCobj.mCoxParameter{end,2}=CoxVDxFx;
                    end
                    f = cellfun(@(x) strcmpi('VxFx',x),OCobj.mCoxParameter(:,1));
                    if any(f)
                        OCobj.mCoxParameter{f,2}=CoxVxFx;
                    else
                        OCobj.mCoxParameter{end+1,1}='VxFx';
                        OCobj.mCoxParameter{end,2}=CoxVxFx;
                    end
                        f = cellfun(@(x) strcmpi('VDxDperFx',x),OCobj.mCoxParameter(:,1));
                    if any(f)
                        OCobj.mCoxParameter{f,2}=CoxVDxDperFx;
                    else
                        OCobj.mCoxParameter{end+1,1}='VDxDperFx';
                        OCobj.mCoxParameter{end,2}=CoxVDxDperFx;
                    end
                end
                    % VDxDMax
                     f = cellfun(@(x) strcmpi('VDxDmax',x),OCobj.mCoxParameter(:,1));
                    if any(f)
                        OCobj.mCoxParameter{f,2}=CoxVDxDMax;
                    else
                        OCobj.mCoxParameter{end+1,1}='VDxDMax';
                        OCobj.mCoxParameter{end,2}=CoxVDxDMax;
                    end
                    % VDxTx
                     f = cellfun(@(x) strcmpi('VDxTx',x),OCobj.mCoxParameter(:,1));
                    if any(f)
                        OCobj.mCoxParameter{f,2}=CoxVDxTx;
                    else
                        OCobj.mCoxParameter{end+1,1}='VDxTx';
                        OCobj.mCoxParameter{end,2}=CoxVDxTx;
                    end
                    
                      % VDxTx
                 
                     % VDxBMI
                     f = cellfun(@(x) strcmpi('VDxBMI',x),OCobj.mCoxParameter(:,1));
                    if any(f)
                        OCobj.mCoxParameter{f,2}=CoxVDxBMI;
                    else
                        OCobj.mCoxParameter{end+1,1}='VDxBMI';
                        OCobj.mCoxParameter{end,2}=CoxVDxBMI;
                    end
                
                    % VDxCM2CW
                     f = cellfun(@(x) strcmpi('VDxCM2CW',x),OCobj.mCoxParameter(:,1));
                    if any(f)
                        OCobj.mCoxParameter{f,2}=CoxVDxCM2CW;
                    else
                        OCobj.mCoxParameter{end+1,1}='VDxCM2CW';
                        OCobj.mCoxParameter{end,2}=CoxVDxCM2CW;
                    end
                
                
                warning('on');
        end
        
        function OCobj = fUniversalSurvivalCurve(OCobj)
            
           %set mUniversalSurvivalCurveLogLikelihood
           %assume mUniversalSurvivalCurveAlpha has been set
           if isempty(OCobj.mUscAlpha)
              error('classOutcomeAnalysis:fUniversalSurvivalCurve, alpha not set!');
           end

           if isempty(OCobj.mUscRangeD0)
               OCobj.mUscRangeD0 = 0.2:0.2:4;
           end
           if isempty(OCobj.mUscRangeDq)
               OCobj.mUscRangeDq = 0.2:0.2:4;
           end
            
           %% Loop over D0 Dq values, calculate DT, compute USC BED, fit outcome data           
           OCobj.mUscLogLikelihoods = inf(length(OCobj.mUscRangeD0),length(OCobj.mUscRangeDq));
           OCobj.mUscDt = inf(length(OCobj.mUscRangeD0),length(OCobj.mUscRangeDq));
           
           alpha = OCobj.mUscAlpha;
           for i=1:length(OCobj.mUscRangeD0)
               cur_d0 = OCobj.mUscRangeD0(i);
               for j=1:length(OCobj.mUscRangeDq)
                   cur_dq = OCobj.mUscRangeDq(j);
                   
                   % Transition Dose
                   cur_dt = (2*cur_dq)/(1-(alpha*cur_d0));
                   
                   % USC BED
                    lq_flag = (EPIobj.mDoseBins_org./100)<cur_dt;
                   
                   
               end
           end
                   
            
        end
        function OCobj = fLogRankDx_DVH(OCobj)
            % prepare
            if isempty(OCobj.mBinsVol)
                OCobj = OCobj.fCalculateVolBins();
            end
            if isempty(OCobj.mBinsDose)
                OCobj = OCobj.fCalculateDoseBins();
            end

            % select patients with data
            f = OCobj.fPatientsWithComplicationData();
            CG = OCobj.fRemovePatient(~f);
            
            % survival/complication time
            f2 = ~cellfun('isempty',{CG.mGrp.mDateComp}); % patients with no complication date
            f3 = ~cellfun('isempty',{CG.mGrp.mDateLastFollowup}); % patients with no last follow up date
            compdate = inf(CG.mNumInGrp,1);
            lastfollowup = inf(CG.mNumInGrp,1);
            compdate(f2) = ([CG.mGrp(f2).mDateComp] - [CG.mGrp(f2).mDateBaseline])' / 30;
            lastfollowup(f3) = ([CG.mGrp(f3).mDateLastFollowup] - [CG.mGrp(f3).mDateBaseline])' / 30;
            compdate = min( lastfollowup, compdate );
            flgcensor = [CG.mGrp.mFlgCensor]';

                volnum=length(CG.mBinsVol);
                dosenum=length(CG.mBinsDose);
                Dxmat = ones(dosenum,volnum,6); % n1,c1,n2,c2,p,flg (0 -- positive corelation, 1 -- negative corelation, 2 -- not available)
                Dxmat(:,:,6) = 2; % default is not available
                
                sa=classKaplanMeierCurve(); % initialize a survivalanalysis obj
                
            % Dx
                dv=zeros(CG.mNumInGrp,1);
                numstart=CG.mLogRankMinSize;
                for v=1:volnum
                    % doses under v
                    dv(:) = 0;
                    for k = 1:CG.mNumInGrp
                        dv(k) = CG.mGrp(k).fDoseAtVol( CG.mBinsVol(v) );
                    end
                    
                    % Dx at (di,vj)
                    g=find(dv); % patients with non-zero doses, used to remove patients (as opposed to DVx)
                    flg_dosebelow2=-1;
                    numend=length(g)-numstart;
                    for d=1:dosenum
                        % check the sample size
                        flg_dosebelow1=dv(g)<=CG.mBinsDose(d); f=length(find(flg_dosebelow1));
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
                if isempty(OCobj.mLogRank)
                        OCobj.mLogRank{end+1,1}='Dx';
                        OCobj.mLogRank{end,2}=Dxmat;
                else
                    f = cellfun(@(x) strcmpi('Dx',x),OCobj.mLogRank(:,1));
                    if any(f)
                        OCobj.mLogRank{f,2}=Dxmat;
                    else
                        OCobj.mLogRank{end+1,1}='Dx';
                        OCobj.mLogRank{end,2}=Dxmat;
                    end
                end
        end
        function OCobj = fLogRankVx_DVH(OCobj)
            % prepare
            if isempty(OCobj.mBinsVol)
                OCobj = OCobj.fCalculateVolBins();
            end
            if isempty(OCobj.mBinsDose)
                OCobj = OCobj.fCalculateDoseBins();
            end

            % select patients with data
            f = OCobj.fPatientsWithComplicationData();
            CG = OCobj.fRemovePatient(~f);
            
            % survival/complication time
            f2 = ~cellfun('isempty',{CG.mGrp.mDateComp}); % patients with no complication date
            f3 = ~cellfun('isempty',{CG.mGrp.mDateLastFollowup}); % patients with no last follow up date
            compdate = inf(CG.mNumInGrp,1);
            lastfollowup = inf(CG.mNumInGrp,1);
            compdate(f2) = ([CG.mGrp(f2).mDateComp] - [CG.mGrp(f2).mDateBaseline])' / 30;
            lastfollowup(f3) = ([CG.mGrp(f3).mDateLastFollowup] - [CG.mGrp(f3).mDateBaseline])' / 30;
            compdate = min( lastfollowup, compdate );
            flgcensor = [CG.mGrp.mFlgCensor]';

                volnum=length(CG.mBinsVol);
                dosenum=length(CG.mBinsDose);
                Vxmat = ones(dosenum,volnum,6); % (n1,c1,n2,c2,p,flg (0 -- positive corelation, 1 -- negative corelation, 2 -- not available))
                Vxmat(:,:,6) = 2;
                
                sa=classKaplanMeierCurve(); % initialize a survivalanalysis obj
                
            % Vx
                vd=zeros(CG.mNumInGrp,1); % volume v at dose d
                numstart=CG.mLogRankMinSize;
                for d=1:dosenum
                    % volume under d
                    vd(:)=0;
                    for k=1:CG.mNumInGrp
                        vd(k) = CG.mGrp(k).fVolAtDose( CG.mBinsDose(d) );
                    end
                    g=find(vd); % non-zeros volume cases
                    
                    % (di,vj)
                    flg_volbelow2=-1;
                    numend=length(g)-numstart;
                    for v=1:volnum
                        % check smaple size
                        flg_volbelow1=vd(g)<=CG.mBinsVol(v); f=length(find(flg_volbelow1)); % group DVHs by (d,v)
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
                if isempty(OCobj.mLogRank)
                        OCobj.mLogRank{end+1,1}='Vx';
                        OCobj.mLogRank{end,2}=Vxmat;
                else
                    f = cellfun(@(x) strcmpi('Vx',x),OCobj.mLogRank(:,1));
                    if any(f)
                        OCobj.mLogRank{f,2}=Vxmat;
                    else
                        OCobj.mLogRank{end+1,1}='Vx';
                        OCobj.mLogRank{end,2}=Vxmat;
                    end
                end
        end
       function OCobj = fLogRankVDx_DVH(OCobj)
            % prepare
            if isempty(OCobj.mBinsVol)
                OCobj = OCobj.fCalculateVolBins();
            end
            if isempty(OCobj.mBinsDose)
                OCobj = OCobj.fCalculateDoseBins();
            end

            % select patients with data
            f = OCobj.fPatientsWithComplicationData();
            CG = OCobj.fRemovePatient(~f);
            
            % survival/complication time
            f2 = ~cellfun('isempty',{CG.mGrp.mDateComp}); % patients with no complication date
            f3 = ~cellfun('isempty',{CG.mGrp.mDateLastFollowup}); % patients with no last follow up date
            compdate = inf(CG.mNumInGrp,1);
            lastfollowup = inf(CG.mNumInGrp,1);
            compdate(f2) = ([CG.mGrp(f2).mDateComp] - [CG.mGrp(f2).mDateBaseline])' / 30;
            lastfollowup(f3) = ([CG.mGrp(f3).mDateLastFollowup] - [CG.mGrp(f3).mDateBaseline])' / 30;
            compdate = min( lastfollowup, compdate );
            flgcensor = [CG.mGrp.mFlgCensor]';
            
            %tmp for a2b
            volnum=length(CG.mBinsVol);
            %vbins = [1:volnum];
            dosenum=length(CG.mBinsDose);
            %dbins = [1:dosenum];
            
            %dosenum=min(length(CG.mBinsDose),200);
            %dbins = [round(linspace(1,max(CG.mBinsDose),dosenum))];
            %dbins = [round(linspace(1,max(CG.mBinsDose),200))];
            %volnum=min(length(CG.mBinsVol),200);
            %vbins = [round(linspace(1,max(CG.mBinsVol),volnum))];
            %vbins = [round(linspace(1,max(CG.mBinsVol),200))];
            
            % (n1,c1,n2,c2,p,flg,cox_hr (flg=0 -- positive corelation, 1 -- negative corelation, 2 -- not available))
            VDxmat = ones(dosenum,volnum,7); 
            VDxmat(:,:,6) = 2;
                
            sa=classKaplanMeierCurve(); % initialize a survivalanalysis obj
                
            % Vx
            vd=zeros(CG.mNumInGrp,1); % volume v at dose d
            numstart=CG.mLogRankMinSize;
            for d=1:dosenum
                    % volume under d
                    vd(:)=0;
                    for k=1:CG.mNumInGrp
                        %tmp
                        %vd(k) = CG.mGrp(k).fVolAtDose( CG.mBinsDose(d) );
                        vd(k) = CG.mGrp(k).fVolAtDose( CG.mBinsDose(d) );
                    end
                    
                    % EDT: 2013-05-03
                    g=find(vd>-1); % non-zeros volume cases, if -1 <- inlclude all
                    %g = intersect(g,find(flgfx')); %from
                    %fPlotKaplanMeierCurve, NOT same flgfx as above!
                    
                    % (di,vj)
                    flg_volbelow2=-1;
                    
                    % EDT: 2013-05-03
                    numend=length(g)-numstart;
                    %numend=CG.mNumInGrp-numstart;
                    for v=1:volnum
                        % check smaple size
                        flg_volbelow1=vd(g)<=CG.mBinsVol(v); f=length(find(flg_volbelow1)); % group DVHs by (d,v)
                        if f<numstart || f>numend % one group has too less patients, or the volume at the dose is too small, skip it
                            continue;
                        end
                        % check the change of grouping
                        if isequal(flg_volbelow1,flg_volbelow2) % if it is the same grouping, skip the computation and save the result directly
                            VDxmat(d,v,:)=VDxmat(d,v-1,:);
                            continue;
                        end
                        flg_volbelow2=flg_volbelow1; % keep the current grouping

                        
                        % assign properties of object sa
                        survivedate={compdate(g(flg_volbelow1)); compdate(g(~flg_volbelow1))}; % survive time of each group
                        fcensor={flgcensor(g(flg_volbelow1)); flgcensor(g(~flg_volbelow1))}; % censor flag for each group
                        
                        cox_beta=coxphfit(~flg_volbelow1,compdate(g),'baseline',0,'censoring',flgcensor(g));
                        cox_hr = exp(cox_beta);
%     
%                         if d==56 && v==14
%                             disp(['']);
%                         end
                        
                        sa.mSurvivalTime=survivedate;
                        sa.mFlgCensor=fcensor;
                        % compute survival curves and compare them
                        sa=sa.fCalculateSurvivalCurve();
                        sa=sa.fCombineSurvivalTime();
                        sa=sa.fCompareSurvivalByLogrank();
                        VDxmat(d,v,1:5)=[length(survivedate{1,1}),sum(fcensor{1,1}), length(survivedate{2,1}),sum(fcensor{2,1}),sa.mpValue];
                        
                        %[~, twoyr_below] =min(abs(sa.mSurvivalTimeSorted{1} - 24));
                        %asympt_below=1-sa.mSurvivalCurve{1}(length(sa.mSurvivalTimeSorted{1}));
                        %[~, twoyr_above] =min(abs(sa.mSurvivalTimeSorted{2} - 24));
                        %asympt_above=1-sa.mSurvivalCurve{2}(length(sa.mSurvivalTimeSorted{2}));
                        
                        % =1 for positive corr at end, =0 for negative
                        %VDxmat(d,v,6)= (asympt_below<asympt_above);
                        
                        % below: 1 = anti-corr, 0 = pos-corr (area of
                        % survival (!) curves
                        VDxmat(d,v,6)=sa.mCurveArea(1)<sa.mCurveArea(2); % the group with lower volume had worse survival curve, record it
                        VDxmat(d,v,7)=cox_hr;
                    end
                end
       
            % save reslt
                if isempty(OCobj.mLogRank)
                        OCobj.mLogRank{end+1,1}='VDx';
                        OCobj.mLogRank{end,2}=VDxmat;
                else
                    f = cellfun(@(x) strcmpi('VDx',x),OCobj.mLogRank(:,1));
                    if any(f)
                        OCobj.mLogRank{f,2}=VDxmat;
                    else
                        OCobj.mLogRank{end+1,1}='VDx';
                        OCobj.mLogRank{end,2}=VDxmat;
                    end
                end
       end 
 function OCobj = fLogRankDVx_DVH(OCobj)
            % prepare
            if isempty(OCobj.mBinsVol)
                OCobj = OCobj.fCalculateVolBins();
            end
            if isempty(OCobj.mBinsDose)
                OCobj = OCobj.fCalculateDoseBins();
            end

            % select patients with data
            f = OCobj.fPatientsWithComplicationData();
            CG = OCobj.fRemovePatient(~f);
            
            % survival/complication time
            f2 = ~cellfun('isempty',{CG.mGrp.mDateComp}); % patients with no complication date
            f3 = ~cellfun('isempty',{CG.mGrp.mDateLastFollowup}); % patients with no last follow up date
            compdate = inf(CG.mNumInGrp,1);
            lastfollowup = inf(CG.mNumInGrp,1);
            compdate(f2) = ([CG.mGrp(f2).mDateComp] - [CG.mGrp(f2).mDateBaseline])' / 30;
            lastfollowup(f3) = ([CG.mGrp(f3).mDateLastFollowup] - [CG.mGrp(f3).mDateBaseline])' / 30;
            compdate = min( lastfollowup, compdate );
            flgcensor = [CG.mGrp.mFlgCensor]';
            
            volnum=length(CG.mBinsVol);
            dosenum=length(CG.mBinsDose);
                        
            % (n1,c1,n2,c2,p,flg,cox_hr (flg=0 -- positive corelation, 1 -- negative corelation, 2 -- not available))
            DVxmat = ones(dosenum,volnum,7); 
            DVxmat(:,:,6) = 2;
                        
            sa=classKaplanMeierCurve(); % initialize a survivalanalysis obj
                
            % Vx
            dv=zeros(CG.mNumInGrp,1); % volume v at dose d
            numstart=CG.mLogRankMinSize;
            for v=1:volnum
                    % volume under d
                    dv(:)=0;
                    for k=1:CG.mNumInGrp
                        %vd(k) = CG.mGrp(k).fVolAtDose( CG.mBinsDose(d) );
                        dv(k) = CG.mGrp(k).fDoseAtVol( CG.mBinsVol(v) );
                    end
                    
                    % EDT: 2013-05-03
                    g=find(dv>-1); % non-zeros volume cases, if -1 <- inlclude all
                    %g = intersect(g,find(flgfx')); %from
                    %fPlotKaplanMeierCurve, NOT same flgfx as above!
                    
                    % (di,vj)
                    flg_volbelow2=-1;
                    
                    % EDT: 2013-05-03
                    numend=length(g)-numstart;
                    %numend=CG.mNumInGrp-numstart;
                    for d=1:dosenum
                        % check smaple size
                        flg_volbelow1=dv(g)<=CG.mBinsDose(d); f=length(find(flg_volbelow1)); % group DVHs by (d,v)
                        if f<numstart || f>numend % one group has too less patients, or the volume at the dose is too small, skip it
                            continue;
                        end
                        % check the change of grouping
                        if isequal(flg_volbelow1,flg_volbelow2) % if it is the same grouping, skip the computation and save the result directly
                            DVxmat(d,v,:)=DVxmat(d-1,v,:);
                            continue;
                        end
                        flg_volbelow2=flg_volbelow1; % keep the current grouping

                        
                        % assign properties of object sa
                        survivedate={compdate(g(flg_volbelow1)); compdate(g(~flg_volbelow1))}; % survive time of each group
                        fcensor={flgcensor(g(flg_volbelow1)); flgcensor(g(~flg_volbelow1))}; % censor flag for each group
                        
                        cox_beta=coxphfit(~flg_volbelow1,compdate(g),'baseline',0,'censoring',flgcensor(g));
                        cox_hr = exp(cox_beta);
%     
%                         if d==56 && v==14
%                             disp(['']);
%                         end
                        
                        sa.mSurvivalTime=survivedate;
                        sa.mFlgCensor=fcensor;
                        % compute survival curves and compare them
                        sa=sa.fCalculateSurvivalCurve();
                        sa=sa.fCombineSurvivalTime();
                        sa=sa.fCompareSurvivalByLogrank();
                        DVxmat(d,v,1:5)=[length(survivedate{1,1}),sum(fcensor{1,1}),length(survivedate{2,1}),sum(fcensor{2,1}),sa.mpValue];
                        DVxmat(d,v,6)=sa.mCurveArea(1)<sa.mCurveArea(2); % the group with lower volume had worse survival curve, record it
                        DVxmat(d,v,7)=cox_hr;
                    end
                end
       
            % save reslt
                if isempty(OCobj.mLogRank)
                        OCobj.mLogRank{end+1,1}='DVx';
                        OCobj.mLogRank{end,2}=DVxmat;
                else
                    f = cellfun(@(x) strcmpi('DVx',x),OCobj.mLogRank(:,1));
                    if any(f)
                        OCobj.mLogRank{f,2}=DVxmat;
                    else
                        OCobj.mLogRank{end+1,1}='DVx';
                        OCobj.mLogRank{end,2}=DVxmat;
                    end
                end
        end 
        function OCobj = fLogRankCM2CW_DVH(OCobj)
            % prepare
            if isempty(OCobj.mBinsVol)
                OCobj = OCobj.fCalculateVolBins();
            end
            if isempty(OCobj.mBinsDose)
                OCobj = OCobj.fCalculateDoseBins();
            end

            % select patients with data
            f = OCobj.fPatientsWithComplicationData();
            CG = OCobj.fRemovePatient(~f);
            
            % survival/complication time
            f2 = ~cellfun('isempty',{CG.mGrp.mDateComp}); % patients with no complication date
            f3 = ~cellfun('isempty',{CG.mGrp.mDateLastFollowup}); % patients with no last follow up date
            compdate = inf(CG.mNumInGrp,1);
            lastfollowup = inf(CG.mNumInGrp,1);
            compdate(f2) = ([CG.mGrp(f2).mDateComp] - [CG.mGrp(f2).mDateBaseline])' / 30;
            lastfollowup(f3) = ([CG.mGrp(f3).mDateLastFollowup] - [CG.mGrp(f3).mDateBaseline])' / 30;
            compdate = min( lastfollowup, compdate );
            flgcensor = [CG.mGrp.mFlgCensor]';

            cm2cw = [CG.mGrp.mDistanceToChestWall]';
            unique_cm2cw = unique(cm2cw);
            distnum=length(unique_cm2cw);
            CM2CWmat = ones(distnum,6); % (n1,c1,n2,c2,p,flg,sa (flg: 0 -- positive corelation, 1 -- negative corelation, 2 -- not available))
            CM2CWmat(:,6) = 2;
                
            sa=classKaplanMeierCurve(); % initialize a survivalanalysis obj
            CM2CWsa =repmat({sa},distnum,1);
            
            numstart=CG.mLogRankMinSize;
            % g=find(vd); % non-zeros volume cases

            flg_volbelow2=-1;
            numend=length(cm2cw)-numstart;
            
            for v=1:length(unique_cm2cw)
                flg_volbelow1=cm2cw<=unique_cm2cw(v); f=length(find(flg_volbelow1)); % group DVHs by (d,v)
                        if f<numstart || f>numend % one group has too less patients, or the volume at the dose is too small, skip it
                            continue;
                        end
                        % check the change of grouping
                        if isequal(flg_volbelow1,flg_volbelow2) % if it is the same grouping, skip the computation and save the result directly
                            CM2CWmat(v,:)=CM2CWmat(v-1,:);
                            continue;
                        end
                        flg_volbelow2=flg_volbelow1; % keep the current grouping
                        % HR
                        cox_beta=coxphfit(flg_volbelow1,compdate,'baseline',0,'censoring',flgcensor);
                        cox_hr = exp(cox_beta);
    
                        % assign properties of object sa
                        survivedate={compdate(flg_volbelow1); compdate(~flg_volbelow1)}; % survive time of each group
                        fcensor={flgcensor(flg_volbelow1); flgcensor(~flg_volbelow1)}; % censor flag for each group
                        sa.mSurvivalTime=survivedate;
                        sa.mFlgCensor=fcensor;
                        sa.mHR = cox_hr;
                        % compute survival curves and compare them
                        sa=sa.fCalculateSurvivalCurve();
                        sa=sa.fCombineSurvivalTime();
                        sa=sa.fCompareSurvivalByLogrank();
                        CM2CWmat(v,1:5)=[length(survivedate{1,1}),sum(fcensor{1,1}), length(survivedate{2,1}),sum(fcensor{2,1}),sa.mpValue];
                        CM2CWmat(v,6)=sa.mCurveArea(1)<sa.mCurveArea(2); % the group with lower volume had worse survival curve, record it
                        CM2CWsa{v}=sa;
            end
        
            % save reslt
                if isempty(OCobj.mLogRank)
                        OCobj.mLogRank{end+1,1}='CM2CW';
                        OCobj.mLogRank{end,2}={CM2CWmat CM2CWsa};
                else
                    f = cellfun(@(x) strcmpi('CM2CW',x),OCobj.mLogRank(:,1));
                    if any(f)
                        OCobj.mLogRank{f,2}={CM2CWmat CM2CWsa};
                    else
                        OCobj.mLogRank{end+1,1}='CM2CW';
                        OCobj.mLogRank{end,2}={CM2CWmat CM2CWsa};
                    end
                end
        end
       function OCobj = fLogRankBMI_DVH(OCobj)
            % prepare
            if isempty(OCobj.mBinsVol)
                OCobj = OCobj.fCalculateVolBins();
            end
            if isempty(OCobj.mBinsDose)
                OCobj = OCobj.fCalculateDoseBins();
            end

            % select patients with data
            f = OCobj.fPatientsWithComplicationData();
            CG = OCobj.fRemovePatient(~f);
            
            % survival/complication time
            f2 = ~cellfun('isempty',{CG.mGrp.mDateComp}); % patients with no complication date
            f3 = ~cellfun('isempty',{CG.mGrp.mDateLastFollowup}); % patients with no last follow up date
            compdate = inf(CG.mNumInGrp,1);
            lastfollowup = inf(CG.mNumInGrp,1);
            compdate(f2) = ([CG.mGrp(f2).mDateComp] - [CG.mGrp(f2).mDateBaseline])' / 30;
            lastfollowup(f3) = ([CG.mGrp(f3).mDateLastFollowup] - [CG.mGrp(f3).mDateBaseline])' / 30;
            compdate = min( lastfollowup, compdate );
            flgcensor = [CG.mGrp.mFlgCensor]';

            bmi = [CG.mGrp.mBMI]';
            bmi_idx = bmi>0;
            
            compdate = compdate(bmi_idx);
            flgcensor = flgcensor(bmi_idx);
            
            bmi = bmi(bmi_idx);
            unique_bmi = unique(bmi);
            distnum=length(unique_bmi);
            BMImat = ones(distnum,6); % (n1,c1,n2,c2,p,flg,sa (flg: 0 -- positive corelation, 1 -- negative corelation, 2 -- not available))
            BMImat(:,6) = 2;
                
            sa=classKaplanMeierCurve(); % initialize a survivalanalysis obj
            BMIsa =repmat({sa},distnum,1);
            
            numstart=CG.mLogRankMinSize;
            % g=find(vd); % non-zeros volume cases

            flg_volbelow2=-1;
            numend=length(bmi)-numstart;
            
            for v=1:length(unique_bmi)
                flg_volbelow1=bmi<=unique_bmi(v); f=length(find(flg_volbelow1)); % group DVHs by (d,v)
                        if f<numstart || f>numend % one group has too less patients, or the volume at the dose is too small, skip it
                            continue;
                        end
                        % check the change of grouping
                        if isequal(flg_volbelow1,flg_volbelow2) % if it is the same grouping, skip the computation and save the result directly
                            BMImat(v,:)=BMImat(v-1,:);
                            continue;
                        end
                        flg_volbelow2=flg_volbelow1; % keep the current grouping

                        % assign properties of object sa
                        survivedate={compdate(flg_volbelow1); compdate(~flg_volbelow1)}; % survive time of each group
                        fcensor={flgcensor(flg_volbelow1); flgcensor(~flg_volbelow1)}; % censor flag for each group
                        sa.mSurvivalTime=survivedate;
                        sa.mFlgCensor=fcensor;
                        % compute survival curves and compare them
                        sa=sa.fCalculateSurvivalCurve();
                        sa=sa.fCombineSurvivalTime();
                        sa=sa.fCompareSurvivalByLogrank();
                        BMImat(v,1:5)=[length(survivedate{1,1}),sum(fcensor{1,1}), length(survivedate{2,1}),sum(fcensor{2,1}),sa.mpValue];
                        BMImat(v,6)=sa.mCurveArea(1)<sa.mCurveArea(2); % the group with lower volume had worse survival curve, record it
                        BMIsa{v}=sa;
                    end
                                
            % save reslt
                if isempty(OCobj.mLogRank)
                        OCobj.mLogRank{end+1,1}='BMI';
                        OCobj.mLogRank{end,2}={BMImat BMIsa};
                else
                    f = cellfun(@(x) strcmpi('BMI',x),OCobj.mLogRank(:,1));
                    if any(f)
                        OCobj.mLogRank{f,2}={BMImat BMIsa};
                    else
                        OCobj.mLogRank{end+1,1}='BMI';
                        OCobj.mLogRank{end,2}={BMImat BMIsa};
                    end
                end
        end
          
        function OCobj = fLogRankV39BMI_DVH(OCobj)
            % prepare
            if isempty(OCobj.mBinsVol)
                OCobj = OCobj.fCalculateVolBins();
            end
            if isempty(OCobj.mBinsDose)
                OCobj = OCobj.fCalculateDoseBins();
            end

            % select patients with data
            f = OCobj.fPatientsWithComplicationData();
            CG = OCobj.fRemovePatient(~f);
            
            % survival/complication time
            f2 = ~cellfun('isempty',{CG.mGrp.mDateComp}); % patients with no complication date
            f3 = ~cellfun('isempty',{CG.mGrp.mDateLastFollowup}); % patients with no last follow up date
            compdate = inf(CG.mNumInGrp,1);
            lastfollowup = inf(CG.mNumInGrp,1);
            compdate(f2) = ([CG.mGrp(f2).mDateComp] - [CG.mGrp(f2).mDateBaseline])' / 30;
            lastfollowup(f3) = ([CG.mGrp(f3).mDateLastFollowup] - [CG.mGrp(f3).mDateBaseline])' / 30;
            compdate = min( lastfollowup, compdate );
            flgcensor = [CG.mGrp.mFlgCensor]';

            dose=39;         
            [~,fdose] = min(abs(CG.mBinsDose - dose));
            vd=zeros(CG.mNumInGrp,1); % volume v at dose d
            vd(:)=0;
            for k=1:CG.mNumInGrp
                vd(k) = CG.mGrp(k).fVolAtDose( CG.mBinsDose(fdose) );
            end
                        
            bmi = [CG.mGrp.mBMI]';
            bmi_idx = bmi>0;
  
            compdate = compdate(bmi_idx);
            flgcensor = flgcensor(bmi_idx);
            
            vd_bmi_split = 0.0207*vd(bmi_idx)+0.0408*bmi(bmi_idx);
            
            
            unique_split = unique(vd_bmi_split);
            distnum=length(unique_split);
            VDxBMImat = ones(distnum,6); % (n1,c1,n2,c2,p,flg,sa (flg: 0 -- positive corelation, 1 -- negative corelation, 2 -- not available))
            VDxBMImat(:,6) = 2;
                
            sa=classKaplanMeierCurve(); % initialize a survivalanalysis obj
            VDxBMIsa =repmat({sa},distnum,1);
            
            numstart=CG.mLogRankMinSize;
            % g=find(vd); % non-zeros volume cases

            flg_volbelow2=-1;
            numend=length(vd_bmi_split)-numstart;
            
            for v=1:length(unique_split)
                flg_volbelow1=vd_bmi_split<=unique_split(v); f=length(find(flg_volbelow1)); % group DVHs by (d,v)
                        if f<numstart || f>numend % one group has too less patients, or the volume at the dose is too small, skip it
                            continue;
                        end
                        % check the change of grouping
                        if isequal(flg_volbelow1,flg_volbelow2) % if it is the same grouping, skip the computation and save the result directly
                            VDxBMImat(v,:)=VDxBMImat(v-1,:);
                            continue;
                        end
                        flg_volbelow2=flg_volbelow1; % keep the current grouping

                        % assign properties of object sa
                      
                        survivedate={compdate(flg_volbelow1); compdate(~flg_volbelow1)}; % survive time of each group
                        fcensor={flgcensor(flg_volbelow1); flgcensor(~flg_volbelow1)}; % censor flag for each group
                        sa.mSurvivalTime=survivedate;
                        sa.mFlgCensor=fcensor;
                        % compute survival curves and compare them
                        sa=sa.fCalculateSurvivalCurve();
                        sa=sa.fCombineSurvivalTime();
                        sa=sa.fCompareSurvivalByLogrank();
                        VDxBMImat(v,1:5)=[length(survivedate{1,1}),sum(fcensor{1,1}), length(survivedate{2,1}),sum(fcensor{2,1}),sa.mpValue];
                        VDxBMImat(v,6)=sa.mCurveArea(1)<sa.mCurveArea(2); % the group with lower volume had worse survival curve, record it
                        VDxBMIsa{v}=sa;
                    end
                                
            % save reslt
                if isempty(OCobj.mLogRank)
                        OCobj.mLogRank{end+1,1}='V39BMI';
                        OCobj.mLogRank{end,2}={VDxBMImat VDxBMIsa};
                else
                    f = cellfun(@(x) strcmpi('V39BMI',x),OCobj.mLogRank(:,1));
                    if any(f)
                        OCobj.mLogRank{f,2}={VDxBMImat VDxBMIsa};
                    else
                        OCobj.mLogRank{end+1,1}='V39BMI';
                        OCobj.mLogRank{end,2}={VDxBMImat VDxBMIsa};
                    end
                end
        end 
        
        function OCobj = fLogRankV30BMI_DVH(OCobj)
            % prepare
            if isempty(OCobj.mBinsVol)
                OCobj = OCobj.fCalculateVolBins();
            end
            if isempty(OCobj.mBinsDose)
                OCobj = OCobj.fCalculateDoseBins();
            end

            % select patients with data
            f = OCobj.fPatientsWithComplicationData();
            CG = OCobj.fRemovePatient(~f);
            
            % survival/complication time
            f2 = ~cellfun('isempty',{CG.mGrp.mDateComp}); % patients with no complication date
            f3 = ~cellfun('isempty',{CG.mGrp.mDateLastFollowup}); % patients with no last follow up date
            compdate = inf(CG.mNumInGrp,1);
            lastfollowup = inf(CG.mNumInGrp,1);
            compdate(f2) = ([CG.mGrp(f2).mDateComp] - [CG.mGrp(f2).mDateBaseline])' / 30;
            lastfollowup(f3) = ([CG.mGrp(f3).mDateLastFollowup] - [CG.mGrp(f3).mDateBaseline])' / 30;
            compdate = min( lastfollowup, compdate );
            flgcensor = [CG.mGrp.mFlgCensor]';

            dose=30;         
            [~,fdose] = min(abs(CG.mBinsDose - dose));
            vd=zeros(CG.mNumInGrp,1); % volume v at dose d
            vd(:)=0;
            for k=1:CG.mNumInGrp
                vd(k) = CG.mGrp(k).fVolAtDose( CG.mBinsDose(fdose) );
            end
                        
            bmi = [CG.mGrp.mBMI]';
            bmi_idx = bmi>0;
  
            compdate = compdate(bmi_idx);
            flgcensor = flgcensor(bmi_idx);
            
            vd_bmi_split = 0.0131*vd(bmi_idx)+0.0426*bmi(bmi_idx);
            
            
            unique_split = unique(vd_bmi_split);
            distnum=length(unique_split);
            VDxBMImat = ones(distnum,6); % (n1,c1,n2,c2,p,flg,sa (flg: 0 -- positive corelation, 1 -- negative corelation, 2 -- not available))
            VDxBMImat(:,6) = 2;
                
            sa=classKaplanMeierCurve(); % initialize a survivalanalysis obj
            VDxBMIsa =repmat({sa},distnum,1);
            
            numstart=CG.mLogRankMinSize;
            % g=find(vd); % non-zeros volume cases

            flg_volbelow2=-1;
            numend=length(vd_bmi_split)-numstart;
            
            for v=1:length(unique_split)
                flg_volbelow1=vd_bmi_split<=unique_split(v); f=length(find(flg_volbelow1)); % group DVHs by (d,v)
                        if f<numstart || f>numend % one group has too less patients, or the volume at the dose is too small, skip it
                            continue;
                        end
                        % check the change of grouping
                        if isequal(flg_volbelow1,flg_volbelow2) % if it is the same grouping, skip the computation and save the result directly
                            VDxBMImat(v,:)=VDxBMImat(v-1,:);
                            continue;
                        end
                        flg_volbelow2=flg_volbelow1; % keep the current grouping

                        % assign properties of object sa
                      
                        survivedate={compdate(flg_volbelow1); compdate(~flg_volbelow1)}; % survive time of each group
                        fcensor={flgcensor(flg_volbelow1); flgcensor(~flg_volbelow1)}; % censor flag for each group
                        sa.mSurvivalTime=survivedate;
                        sa.mFlgCensor=fcensor;
                        % compute survival curves and compare them
                        sa=sa.fCalculateSurvivalCurve();
                        sa=sa.fCombineSurvivalTime();
                        sa=sa.fCompareSurvivalByLogrank();
                        VDxBMImat(v,1:5)=[length(survivedate{1,1}),sum(fcensor{1,1}), length(survivedate{2,1}),sum(fcensor{2,1}),sa.mpValue];
                        VDxBMImat(v,6)=sa.mCurveArea(1)<sa.mCurveArea(2); % the group with lower volume had worse survival curve, record it
                        VDxBMIsa{v}=sa;
                    end
                                
            % save reslt
                if isempty(OCobj.mLogRank)
                        OCobj.mLogRank{end+1,1}='V30BMI';
                        OCobj.mLogRank{end,2}={VDxBMImat VDxBMIsa};
                else
                    f = cellfun(@(x) strcmpi('V30BMI',x),OCobj.mLogRank(:,1));
                    if any(f)
                        OCobj.mLogRank{f,2}={VDxBMImat VDxBMIsa};
                    else
                        OCobj.mLogRank{end+1,1}='V30BMI';
                        OCobj.mLogRank{end,2}={VDxBMImat VDxBMIsa};
                    end
                end
        end
  
        function OCobj = fLogRankV30BMITx_DVH(OCobj)
            % prepare
            if isempty(OCobj.mBinsVol)
                OCobj = OCobj.fCalculateVolBins();
            end
            if isempty(OCobj.mBinsDose)
                OCobj = OCobj.fCalculateDoseBins();
            end

            % select patients with data
            f = OCobj.fPatientsWithComplicationData();
            CG = OCobj.fRemovePatient(~f);
            
            % survival/complication time
            f2 = ~cellfun('isempty',{CG.mGrp.mDateComp}); % patients with no complication date
            f3 = ~cellfun('isempty',{CG.mGrp.mDateLastFollowup}); % patients with no last follow up date
            compdate = inf(CG.mNumInGrp,1);
            lastfollowup = inf(CG.mNumInGrp,1);
            compdate(f2) = ([CG.mGrp(f2).mDateComp] - [CG.mGrp(f2).mDateBaseline])' / 30;
            lastfollowup(f3) = ([CG.mGrp(f3).mDateLastFollowup] - [CG.mGrp(f3).mDateBaseline])' / 30;
            compdate = min( lastfollowup, compdate );
            flgcensor = [CG.mGrp.mFlgCensor]';

            dose=30;
            [~,fdose] = min(abs(CG.mBinsDose - dose));
            vd=zeros(CG.mNumInGrp,1); % volume v at dose d
            vd(:)=0;
            for k=1:CG.mNumInGrp
                vd(k) = CG.mGrp(k).fVolAtDose( CG.mBinsDose(fdose) );
            end
                        
            bmi = [CG.mGrp.mBMI]';
            bmi_idx = bmi>0;
            
            compdate = compdate(bmi_idx);
            flgcensor = flgcensor(bmi_idx);
            
            tx = [CG.mGrp.mDoseTx]';

            vd_bmi_tx_split = 0.011902*vd(bmi_idx)+...
                0.042235*bmi(bmi_idx)+...
                0.000527*tx(bmi_idx);
            
            
            unique_split = unique(vd_bmi_tx_split);
            distnum=length(unique_split);
            V30BMITxmat = ones(distnum,6); % (n1,c1,n2,c2,p,flg,sa (flg: 0 -- positive corelation, 1 -- negative corelation, 2 -- not available))
            V30BMITxmat(:,6) = 2;
                
            sa=classKaplanMeierCurve(); % initialize a survivalanalysis obj
            V30BMITxsa =repmat({sa},distnum,1);
            
            numstart=CG.mLogRankMinSize;
            % g=find(vd); % non-zeros volume cases

            flg_volbelow2=-1;
            numend=length(vd_bmi_tx_split)-numstart;
            
            for v=1:length(unique_split)
                flg_volbelow1=vd_bmi_tx_split<=unique_split(v); f=length(find(flg_volbelow1)); % group DVHs by (d,v)
                        if f<numstart || f>numend % one group has too less patients, or the volume at the dose is too small, skip it
                            continue;
                        end
                        % check the change of grouping
                        if isequal(flg_volbelow1,flg_volbelow2) % if it is the same grouping, skip the computation and save the result directly
                            V30BMITxmat(v,:)=V30BMITxmat(v-1,:);
                            continue;
                        end
                        flg_volbelow2=flg_volbelow1; % keep the current grouping

                        % assign properties of object sa
                        survivedate={compdate(flg_volbelow1); compdate(~flg_volbelow1)}; % survive time of each group
                        fcensor={flgcensor(flg_volbelow1); flgcensor(~flg_volbelow1)}; % censor flag for each group
                        sa.mSurvivalTime=survivedate;
                        sa.mFlgCensor=fcensor;
                        % compute survival curves and compare them
                        sa=sa.fCalculateSurvivalCurve();
                        sa=sa.fCombineSurvivalTime();
                        sa=sa.fCompareSurvivalByLogrank();
                        V30BMITxmat(v,1:5)=[length(survivedate{1,1}),sum(fcensor{1,1}), length(survivedate{2,1}),sum(fcensor{2,1}),sa.mpValue];
                        V30BMITxmat(v,6)=sa.mCurveArea(1)<sa.mCurveArea(2); % the group with lower volume had worse survival curve, record it
                        V30BMITxsa{v}=sa;
                    end
                                
            % save reslt
                if isempty(OCobj.mLogRank)
                        OCobj.mLogRank{end+1,1}='V30BMITX';
                        OCobj.mLogRank{end,2}={V30BMITxmat V30BMITxsa};
                else
                    f = cellfun(@(x) strcmpi('V30BMITX',x),OCobj.mLogRank(:,1));
                    if any(f)
                        OCobj.mLogRank{f,2}={V30BMITxmat V30BMITxsa};
                    else
                        OCobj.mLogRank{end+1,1}='V30BMITX';
                        OCobj.mLogRank{end,2}={V30BMITxmat V30BMITxsa};
                    end
                end
        end 
        function OCobj = fLogRankDVx_DVH_old(OCobj)
            % prepare
            if isempty(OCobj.mBinsVol)
                OCobj = OCobj.fCalculateVolBins();
            end
            if isempty(OCobj.mBinsDose)
                OCobj = OCobj.fCalculateDoseBins();
            end

            % select patients with data
            f = OCobj.fPatientsWithComplicationData();
            CG = OCobj.fRemovePatient(~f);
            
            % survival/complication time
            f2 = ~cellfun('isempty',{CG.mGrp.mDateComp}); % patients with no complication date
            f3 = ~cellfun('isempty',{CG.mGrp.mDateLastFollowup}); % patients with no last follow up date
            compdate = inf(CG.mNumInGrp,1);
            lastfollowup = inf(CG.mNumInGrp,1);
            compdate(f2) = ([CG.mGrp(f2).mDateComp] - [CG.mGrp(f2).mDateBaseline])' / 30;
            lastfollowup(f3) = ([CG.mGrp(f3).mDateLastFollowup] - [CG.mGrp(f3).mDateBaseline])' / 30;
            compdate = min( lastfollowup, compdate );
            flgcensor = [CG.mGrp.mFlgCensor]';

                volnum=length(CG.mBinsVol);
                dosenum=length(CG.mBinsDose);
                DVxmat = ones(dosenum,volnum,6); % n1,c1,n2,c2,p,flg (0 -- positive corelation, 1 -- negative corelation, 2 -- not available)
                DVxmat(:,:,6) = 2; % default is not available
                
                sa=classKaplanMeierCurve(); % initialize a survivalanalysis obj
                
            % Dx
                dv=zeros(CG.mNumInGrp,1);
                numstart=CG.mLogRankMinSize;
                numend=CG.mNumInGrp-numstart;
                for v=1:volnum
                    % doses under v
                    dv(:) = 0;
                    for k = 1:CG.mNumInGrp
                        dv(k) = CG.mGrp(k).fDoseAtVol( CG.mBinsVol(v) );
                    end
                    
                    % DVx at (dj,vj)
                    flg_dosebelow2=-1;
                    for d=1:dosenum
                        % check the sample size
                        flg_dosebelow1=dv<=CG.mBinsDose(d); f=length(find(flg_dosebelow1));
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
                if isempty(OCobj.mLogRank)
                        OCobj.mLogRank{end+1,1}='DVx';
                        OCobj.mLogRank{end,2}=DVxmat;
                else
                    f = cellfun(@(x) strcmpi('DVx',x),OCobj.mLogRank(:,1));
                    if any(f)
                        OCobj.mLogRank{f,2}=DVxmat;
                    else
                        OCobj.mLogRank{end+1,1}='DVx';
                        OCobj.mLogRank{end,2}=DVxmat;
                    end
                end
        end
        function [allCox,flgCox,flgAnti] = fCoxParameter_DVH(OCobj,strCoxVx)
            f = cellfun(@(x) strcmpi(strCoxVx,x),OCobj.mCoxParameter(:,1)); % search the label
            if any(f)
                allCox = OCobj.mCoxParameter{f,2}; % extract Cox model result
                flgCox = ~arrayfun( @(y) any(structfun(@(x) any(isempty(x(:)))|any(isinf(x(:))), y)), allCox); % some fields are empty or infinite, indicating no data for those values
                % correlation
                f = [allCox.beta]';
                flgAnti = f<0;
            else
                allCox = []; flgCox = []; flgAnti = [];
            end
        end
    end

    methods % DVH atlas
        function OCobj = fKaplanMeierCompMatAtSampleTime_DVH(OCobj)
            % prepare
            psem = zeros( length(OCobj.mBinsDose), length(OCobj.mBinsVol), length(OCobj.mBinsTime) );
%             pse = OCobj.mKaplanMeierCompMat(:); 
%             pse = reshape(pse,size(OCobj.mKaplanMeierCompMat));
            pse = OCobj.mKaplanMeierCompMat;
            tb = OCobj.mBinsTime(:);
            db = OCobj.mBinsDose(:);
            vb = OCobj.mBinsVol(:);

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
            OCobj.mKaplanMeierCompSample = psem;
        end
        function OCobj = fAtlasAlongTime_DVH(OCobj)
            % prepare
            OCobj.mAtlasTotal = zeros( length(OCobj.mBinsDose), length(OCobj.mBinsVol), length(OCobj.mBinsTime) );
            OCobj.mAtlasComp = zeros( length(OCobj.mBinsDose), length(OCobj.mBinsVol), length(OCobj.mBinsTime) );
            
            % select patients with data
            f = OCobj.fPatientsWithComplicationData();
            pt = OCobj.mGrp(f);
            num = sum(f); % number of patients

            % complication time
            f2 = ~cellfun(@(x) isempty(x)|isinf(x)|isnan(x),{pt.mDateComp}); % patients with complication date
            f3 = ~cellfun(@(x) isempty(x)|isinf(x)|isnan(x),{pt.mDateLastFollowup}); % patients with last follow up date
            complicationdays = inf(num,1);
            lastfollowup = inf(num,1);
            complicationdays(f2) = ([pt(f2).mDateComp] - [pt(f2).mDateBaseline])' / 30; % convert to months
            lastfollowup(f3) = ([pt(f3).mDateLastFollowup] - [pt(f3).mDateBaseline])' / 30;
            lastfollowup = min( lastfollowup, complicationdays );

            % Vx
            vols = zeros( num, 1 );
            for ii = 1:length(OCobj.mBinsDose)
                % volume of each patient at current dose (x)
                vols(:) = 0;
                for jj = 1:num
                    vols(jj) = pt(jj).fVolAtDose( OCobj.mBinsDose(ii) );
                end
                vols(vols==0) = -1; % exclude zero volume patients
                
                % matrix at each (Di, Vj, Tk)
                for jj = 1:length(OCobj.mBinsVol) % for each volume point under dose x
                    f = find( vols >= OCobj.mBinsVol(jj) ); % patient at the grid point
                    % for each time point
                    for kk = 1:length(OCobj.mBinsTime)
                        % patients with last followup or complications before the time point
                        f1 = find( lastfollowup(f) <= OCobj.mBinsTime(kk) );
                        % total patients
                        OCobj.mAtlasTotal(ii,jj,kk) = length(f1);
                        % patients with complications
                        f1 = find( complicationdays(f) <= OCobj.mBinsTime(kk) );
                        OCobj.mAtlasComp(ii,jj,kk) = length(f1);
                    end
                end
            end
        end
        function OCobj = fKaplanMeierCompFromAtlas_DVH(OCobj)
            % prepare
            OCobj.mKaplanMeierCompFromAtlas = ones( length(OCobj.mBinsDose), length(OCobj.mBinsVol), length(OCobj.mBinsTime) );
            
            % survival along the time bins
            for kk = 2:length(OCobj.mBinsTime)
                % patients at risk
                ptrisk = OCobj.mAtlasTotal(:,:,end) - OCobj.mAtlasTotal(:,:,kk-1);
                % patient with complication
                ptcomp = OCobj.mAtlasComp(:,:,end) - OCobj.mAtlasComp(:,:,kk-1);
                % rate of patients without complication
                OCobj.mKaplanMeierCompFromAtlas(:,:,kk) = OCobj.mKaplanMeierCompFromAtlas(:,:,kk-1).*(1-ptcomp./ptrisk);
            end
        end

        function OCobj = fCoxModelAtlas_DVH(OCobj)
            if isempty(OCobj.mAtlasTotal)
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
            flgcensor = false(OCobj.mAtlasTotal(1,1,end),1);
            vx = zeros(size(flgcensor));
            dx = zeros(size(flgcensor));
            compdate = zeros(size(flgcensor));
            
            [dimd,dimv,dimt] = size( OCobj.mAtlasTotal ); % dimensions of dose, volume, and time
        
            dosebins = OCobj.mBinsDose; dosebins(1:end-1) = (dosebins(1:end-1) + dosebins(2:end)) / 2; % dose bins
            volbins = OCobj.mBinsVol; volbins(1:end-1) = (volbins(1:end-1) + volbins(2:end)) / 2; % volume bins
            timebins = OCobj.mBinsTime; timebins(2:end) = (timebins(1:end-1) + timebins(2:end)) / 2; % time bins
            
            % Cox model
            [~,logl,h,stats]=coxphfit([0; 1],[0; 1],'baseline',0,'censoring',[0; 1]);
            stats.logl=logl; stats.h=h;
            stats.data_exposure = [0; 1]; stats.data_hazard = [0; 1];
            % Vx
            CoxVxAtlas = repmat(stats,[dimd,1]);
            for d = 1:dimd % Cox model for each dose (x)
                % extract the matrix corresponding to the dose
                vc = squeeze( OCobj.mAtlasComp(d,:,:) ); % complication at different volume and time points
                vt = squeeze( OCobj.mAtlasTotal(d,:,:) ); % complication and censored at different time points
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
                dc = squeeze( OCobj.mAtlasComp(:,v,:) ); % complication at different volume and time points
                dt = squeeze( OCobj.mAtlasTotal(:,v,:) ); % complication and censored at different time points
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
            f = cellfun(@(x) strcmpi('VxAtlas',x),OCobj.mCoxParameter(:,1));
            if any(f)
                OCobj.mCoxParameter{f,2}=CoxVxAtlas;
            else
                OCobj.mCoxParameter{end+1,1}='VxAtlas';
                OCobj.mCoxParameter{end,2}=CoxVxAtlas;
            end
            f = cellfun(@(x) strcmpi('DxAtlas',x),OCobj.mCoxParameter(:,1));
            if any(f)
                OCobj.mCoxParameter{f,2}=CoxDxAtlas;
            else
                OCobj.mCoxParameter{end+1,1}='DxAtlas';
                OCobj.mCoxParameter{end,2}=CoxDxAtlas;
            end
            
            warning('on');
%             warning('on','MATLAB:singularMatrix');
%             warning('on','stats:coxphfit:FitWarning');
%             warning('on','stats:coxphfit:RankDeficient');
%             warning('on','stats:coxphfit:IterOrEvalLimit');
        end
    end

    methods % DVH plot and atlas writing
        function fCoxRiskVDxFig_DVH(OCobj,d,t)
            numintv = 1; % group number, i.e., interval numbers
            % Vx
            % search the best Cox model
            [allCox,flgCox,flganti] = fCoxParameter_DVH(OCobj,'VDx'); % find availabe Cox models
            flgCox(flganti)=false; % anti-correlations were not be considered
            logl = [allCox.logl]'; logl(~flgCox) = -inf; % log likelihood of Cox model, anti-correlation points not counted
            [~,doseloc]=max(logl); % the best fitting of Cox model
            disp(' ');
            disp(['Best Cox Model is at dose :',num2str(OCobj.mBinsDose(doseloc))]);
            disp(allCox(doseloc));
            if exist('d','var')
                if d~=-1
                    doseloc = OCobj.mBinsDose == d;
                end
            end
            allCox = allCox(doseloc);
%             disp(' ');
            disp(['Cox model at dose: ',num2str(OCobj.mBinsDose(doseloc))]);
            disp(allCox);

            flgcensor = [OCobj.mGrp.mFlgCensor]';% the overall patient complication info
            % specify event time
            comptime = ([OCobj.mGrp.mDateComp]' - [OCobj.mGrp.mDateBaseline]')/30;
%             disp(' ');
            disp(['median complication time: ',num2str(median(comptime(~flgcensor)))]);
            if ~exist('t','var')
                t = median(comptime(~flgcensor));
            end
            disp(['Cox model analysis at time: ',num2str(t)]);

            % volumes of patients at best Vx
            Vx=zeros(OCobj.mNumInGrp,1);
            d = OCobj.mBinsDose(doseloc);
            for k=1:OCobj.mNumInGrp
                Vx(k) = OCobj.mGrp(k).fVolAtDose( d );
            end
            % observed volume at best Vx in groups
%             flgcensor(comptime>t) = 1; % by the time t some patients might not develop complications so they shall be excluded.
%             for k = 1:OCobj.mNumInGrp
%                 OCobj.mGrp(k).mFlgCensor = flgcensor(k);
%             end
            [sortQ,indxQ,indxorg] = EqualIntervals(Vx,numintv);
            meanvol = zeros(numintv,1);
            prob = zeros(numintv,1);
            stdvol = zeros(numintv,1);
            stdprob = zeros(numintv,1);
%             betainv84 = zeros(numintv,1);
%             betainv16 = zeros(numintv,1);
%             f = ~flgcensor(indxorg);
            for m = 1 : numintv
                %vv=sortQ(indxQ==m); vv(vv==0)=[];
                %meanvol(m) = median(vv);
                %meanvol(m) = mean(vv);
                
                meanvol(m) = mean(sortQ(indxQ==m));
                stdvol(m) = std(sortQ(indxQ==m));
%                 numcomp = sum(f(indxQ==m)); numtotal = sum(indxQ==m);
%                 prob(m) = numcomp/numtotal;
%                 betainv84(m) = betainv( .84, numcomp+1, numtotal - numcomp + 1 );
%                 betainv16(m) = betainv( .16, numcomp+1, numtotal - numcomp + 1 );
            end
            
            
            stdvolL=meanvol-stdvol; stdvolL=max(0,min(inf,stdvolL));
            stdvolU=meanvol+stdvol; stdvolU=max(0,min(inf,stdvolU));
            
            % servivial curves for each group
%             str = 'rbkc';
            CoxComplicationTime = cell(numintv,1);
            CoxComplicationCurve = cell(numintv,1);
            CG = repmat(OCobj,[numintv,1]);
            for m = 1:numintv
                % compute the complication of the group
                CG(m) = OCobj.fRemovePatient(indxorg(indxQ ~= m));
                
%                 for k = 1:CG(m).mNumInGrp
%                     Vx(k) = CG(m).mGrp(k).fVolAtDose(d);
%                 end
%                 hold on; plot(m,Vx(1:CG(m).mNumInGrp),[str(m),'*']); hold off;
                
                CG(m) = CG(m).fOverallCompCurve();
                [CoxComplicationTime{m},CoxComplicationCurve{m}] = fCoxRiskVDxComplicationFig_DVH(CG(m),allCox,meanvol(m));

                f = find(CG(m).mKaplanMeierCompOverall.mSurvivalTimeSorted{1} <= t);
                prob(m) = CG(m).mKaplanMeierCompOverall.mSurvivalCurve{1}(f(end));
                stdprob(m) = CG(m).mKaplanMeierCompOverall.mSurvivalCurveVariance{1}(f(end));
            end
            prob = 1-prob;
            stdprobL = prob-stdprob; stdprobL=max(0,min(1,stdprobL));
            stdprobU = prob+stdprob; stdprobU=max(0,min(1,stdprobU));
            
            % observed risk ratio
            disp(' ');
            disp(['observed risk and the ratio of last and first quartile: ',num2str([prob', prob(end)/prob(1)])]);

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
%             hold on;
%             plot(meanvol,prob,'r*','linewidth',1,'markersize',12); hold off;
            hold on; ploterr(meanvol,prob,{stdvolL,stdvolU},{stdprobL,stdprobU},'r*'); hold off;
%             hold on; ploterr(meanvol,prob,stdvol,stdprob,'r*','LineWidth',1,'markersize',12); hold off;



            % plot Cox survival curves at different volume resolutions
            % survival curve from Cox model at volume vol
            CoxCompTime_HighResolution = allCox.h(1,1):0.1:allCox.h(end,1);
            h = interp1(allCox.h(:,1),allCox.h(:,2),CoxCompTime_HighResolution,'linear','extrap');
            expbetax = exp(allCox.beta*meanvol);
            CoxCompCurve_HighResolution = zeros(length(CoxCompTime_HighResolution),numintv);
            for m = 1:numintv
               CoxCompCurve_HighResolution(:,m) = exp( -h * expbetax(m) );
               %CoxCompCurve_HighResolution(:,m) = exp( -h );
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

                plot(CoxCompTime_HighResolution, 1-CoxCompCurve_HighResolution(:,m),[str(mod(m,length(str))+1),'--'],'LineWidth',1); % high resolution Cox curve

                % agrees with
                %plot([CoxComplicationTime{m}], [1-CoxComplicationCurve{m}],'r--','LineWidth',1); % high resolution Cox curve 
                
                sa = CG(m).mKaplanMeierCompOverall; 
                stairs(sa.mSurvivalTimeSorted{1}, 1-sa.mSurvivalCurve{1}, [str(mod(m,length(str))+1),'-']);
                plot(sa.mSurvivalTimeSorted{1}(sa.mCensorStatistics{1}(:,1)),...
                    1-sa.mSurvivalCurve{1}(sa.mCensorStatistics{1}(:,1)),'+');
            end
            hold off;



            set(gca,'xminortick','on','yminortick','on');
            set(gca,'box','on');
            xlabel('Month'); ylabel('Probability of Complication');
        end
        function [CoxComplicationTime,CoxComplicationCurve] = fCoxRiskVDxComplicationFig_DVH(OCobj,CoxPar,vol) % plot K-M survival curves for patients in the OCobj and that predicted by Cox model "CoxPar" at specified volume "vol"
%             % K-M survival curve
%             OCobj = OCobj.ComplicationCurves_DVH();

            % h(t) clean up
            f = find(diff(CoxPar.h(:,1))==0); % find duplicate time values of h(t)
            while ~isempty(f)
                CoxPar.h(f,1) = CoxPar.h(f,1)-eps*10; % adjust it a bit to avoid ambiguius
                f = find(diff(CoxPar.h(:,1))==0); % find duplicate time values of h(t)
            end

            % survival curve from Cox model at volume vol
            CoxComplicationTime = OCobj.mKaplanMeierCompOverall.mSurvivalTimeSorted{1};
            CoxComplicationTime(OCobj.mKaplanMeierCompOverall.mCensorStatistics{1}(:,1)) = [];
%             CoxComplicationTime(CoxComplicationTime<CoxPar.h(1,1)) = [];
            h = interp1(CoxPar.h(:,1),CoxPar.h(:,2),CoxComplicationTime,'linear','extrap');
            expbetax = exp(CoxPar.beta*vol);
            CoxComplicationCurve = exp( -h * expbetax );
        end
        function fCoxFig_DVH(OCobj,parm)
            % this is only a demo of how to use the cox model to plot response function

            % Vx
            % search the best Cox model
            [allCox,flgCox,flganti] = fCoxParameter_DVH(OCobj,parm); % find availabe Cox models
            
            flgCox(flganti)=false; % anti-correlations were not be considered
            
            logl = [allCox.logl]'; logl(~flgCox) = -inf; % log likelihood of Cox model, anti-correlation points not counted
            [mx,doseloc]=max(logl); % the best fitting of Cox model
            lowCI68 = mx - 0.5; % 68% confidence
            lowCI95 = mx - 1.96; % 95% confidence

%             num = cellfun(@(x) size(x,1),{allCox.data_exposure});
            
            figure(1); clf reset; plot(OCobj.mBinsDose(flgCox), [allCox(flgCox).logl],'.-');
            hold on; plot(OCobj.mBinsDose(flgCox),repmat(lowCI68,size(OCobj.mBinsDose(flgCox))),'r--'); hold off;
            hold on; plot(OCobj.mBinsDose(flgCox),repmat(lowCI95,size(OCobj.mBinsDose(flgCox))),'c--'); hold off;
            set(gca,'fontsize',14);
            set(gca,'xminortick','on','yminortick','on');
            set(gca,'box','on');
            xlabel('Dose (Gy)','fontsize',18); ylabel('log likelihood','fontsize',18);
            figure(2); clf reset; semilogy(OCobj.mBinsDose(flgCox),[allCox(flgCox).p],'.-');
            hold on; semilogy(OCobj.mBinsDose(flgCox),repmat(0.05,size(OCobj.mBinsDose(flgCox))),'r--'); hold off;
            set(gca,'fontsize',14);
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
%             disp([OCobj.mBinsDose(doseloc),m,allCox.p]); % display the selected Cox model
            
%             % median patient complication time
%             f = OCobj.fPatientsWithComplicationData(); % select patients with data
%             CG = OCobj.fRemovePatient(~f);
%             f2 = ~cellfun('isempty',{CG.mGrp.mDateComp}); % patients with no complication date
%             compdate = inf(CG.mNumInGrp,1);
%             compdate(f2) = ([CG.mGrp(f2).mDateComp] - [CG.mGrp(f2).mDateBaseline])' / 30;
%             t = median(compdate(isfinite(compdate)));
% %             t = 50;
%             flg = compdate > t;
% 
% %             f3 = ~cellfun('isempty',{CG.mGrp.mDateLastFollowup}); % patients with no last follow up date
% %             lastfollowup = inf(CG.mNumInGrp,1);
% %             lastfollowup(f3) = ([CG.mGrp(f3).mDateLastFollowup] - [CG.mGrp(f3).mDateBaseline])' / 30;
% %             compdate = min( lastfollowup, compdate );
% 
%             % observed data at time t
%             % volumes of patients
%             vd=zeros(CG.mNumInGrp,1);
%             d = CG.mBinsDose(doseloc);
%             for k=1:CG.mNumInGrp
%                 vd(k) = CG.mGrp(k).fVolAtDose( d );
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
        function fDVHCurvesSummary_DVH(OCobj)
            % prepare
            f = [OCobj.mGrp.mFlgCensor]; % censor info
            dosebins = cat(1,OCobj.mGrp.mDoseBins_LQ); % all doses
            dosebins = (0:max(dosebins))'; % dose bins for complication patients
            vol_center = zeros(length(dosebins),5); % 5 columns for -95%, -68%, median, 68%, 95% lines
            vol_comp = vol_center;

            % volume computation
            vol = -inf(OCobj.mNumInGrp,1);
            for kk = 1:length(dosebins)
                % volumes of each patient at dose kk
                vol(:) = -inf;
                for mm = 1:OCobj.mNumInGrp
                    vol(mm) = OCobj.mGrp(mm).fVolAtDose( dosebins(kk) );
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
            h_cens=plot(dosebins,vol_center(:,3),'b.-');
            h_cens_68=plot(dosebins,vol_center(:,[2,4]),'b');
            h_cens_95=plot(dosebins,vol_center(:,[1,5]),'b--');
            h_comp=plot(dosebins,vol_comp(:,3),'r.-');
            h_comp_68=plot(dosebins,vol_comp(:,[2,4]),'r');
            h_comp_95=plot(dosebins,vol_comp(:,[1,5]),'r--');
            
            lgnd=legend([h_cens h_cens_68(1) h_cens_95(1)...
                h_comp h_comp_68(1) h_comp_95(1)],...
                'W/out comp. median',...
                '  68% CL',...
                '  95% CL',...
                'With comp. median',...
                '  68% CL',...
                '  95% CL');
            set(lgnd,'FontSize',16);
            set(gca,'xminortick','on','yminortick','on');
            set(gca,'box','on');
            xlabel('Dose (Gy)','FontSize',18);
            ylabel('Volume (cc)','FontSize',18);
            hold off;
        end
    end

    methods % EUD operations
        function OCobj = fCalculateEUDBins(OCobj)
            euds = [OCobj.mGrp.mEUD]'; dmax = max(euds(:));
            OCobj.mBinsDose = (0:OCobj.mStepDose:dmax+OCobj.mStepDose)';
        end
        function OCobj = fCalculateEUDBinsLog(OCobj)
            euds = [OCobj.mGrp.mEUD]'; dmax = log10(max(euds(:)));
%             f = euds(:)>0; dmin = log10(min(euds(f)));
%             dosebinslog = ((dmin-OCobj.mStepDose): OCobj.mStepDose:(dmax+OCobj.mStepDose))';
            dosebinslog = 0 : OCobj.mStepDose : (dmax+OCobj.mStepDose)';
            OCobj.mBinsDose = 10.^dosebinslog;
        end
        function OCobj = fCalculateEUD(OCobj)
            for k = 1:OCobj.mNumInGrp
                OCobj.mGrp(k) = OCobj.mGrp(k).fCalculateEUD();
            end
        end
        
        function OCobj = fCrudeAtlas_EUD(OCobj)
            % extract all EUDs from all patients
                euds = [OCobj.mGrp.mEUD]'; % matrix (pt,euds)

            % generate dose bins if mBinsDose is not specifically assigned
                if isempty(OCobj.mBinsDose)
                    error('mBinsDose not determined before method "CrudeAtlas_EUD" is called (classOutcomeAnalysis)');
                end
                numdosebins=size(OCobj.mBinsDose,1);
                
            % censor and complication info
                flgcensor = [OCobj.mGrp.mFlgCensor]; flgcomp = ~flgcensor; % a patient either was censored or had complication
                
            % for each log10(n) and each dose step, compute the total patients and their complications
                OCobj.mAtlasTotal = zeros( numdosebins, size(OCobj.mLymanN,1) );
                OCobj.mAtlasComp = OCobj.mAtlasTotal; % (numdosebins X mLymanN)
                for n = 1:size(OCobj.mLymanN,1)
                    for m = 1:numdosebins
                        % f indices of euds >= current bin dose
                        f = find( euds(:,n) >= OCobj.mBinsDose(m) );
                        % g indices of pts with complication and eud >
                        % current bin dose
                        g = find( flgcomp(f) );
                        OCobj.mAtlasTotal(m,n) = length(f);
                        OCobj.mAtlasComp(m,n) = length(g);
                    end
                end
        end
        
        function OCobj = fCrudeAtlas_DVH(OCobj,nfx) % !! functionality moved to scripts/GenerateDvhAtlas.m
            
            % generate dose bins if mBinsDose is not specifically assigned
            if isempty(OCobj.mBinsDose)
                error('mBinsDose not determined before method "CrudeAtlas_EUD" is called (classOutcomeAnalysis)');
            end
            numdosebins=size(OCobj.mBinsDose,1);
            
            OCgrp = [OCobj.mGrp];
            fx_flg = ones(length(OCgrp),1);
            
            if nfx>0 %nfx= -1, all fx schemes
                
                grp_flgs = [OCgrp.mFxNum];
                fx_flg = grp_flgs==nfx;
                
                if sum(fx_flg)==0
                    disp(['Invalid nfx arg']);
                    return;
                end
                
                OCgrp = OCgrp(fx_flg);
                
            end
            
            
            % censor and complication info
                flgcensor = [OCgrp.mFlgCensor]; flgcomp = ~flgcensor; % a patient either was censored or had complication
                
            vol_bins = [0:5:100];
            % prepare
            %OCobj.mAtlasTotal_DVH = zeros( length(OCobj.mBinsDose), length(OCobj.mBinsVol));
            %OCobj.mAtlasComp_DVH = zeros( length(OCobj.mBinsDose), length(OCobj.mBinsVol));
            
            OCobj.mAtlasTotal_DVH = zeros( length(OCobj.mBinsDose), length(vol_bins));
            OCobj.mAtlasComp_DVH = zeros( length(OCobj.mBinsDose), length(vol_bins));
            
            
            % assuming all have comp data
            %f = OCobj.fPatientsWithComplicationData();
            pt = OCgrp;
            num = length(pt); % number of patients

            % Vx
            vols = zeros( num, 1 );
            for ii = 1:length(OCobj.mBinsDose)
                % volume of each patient at current dose (x)
                vols(:) = 0;
                for jj = 1:num
                    cur_vol = pt(jj).fVolAtDose( OCobj.mBinsDose(ii) );
                    vols(jj) = cur_vol/max(pt(jj).mVolCum);
                    vols(jj) = vols(jj).*100;
                end
                vols(vols==0) = -1; % exclude zero volume patients
                
                % matrix at each (Di, Vj, Tk)
                for kk = 1:length(vol_bins) % for each volume point under dose x
                    %f = find( vols >= OCobj.mBinsVol(jj) ); % patient at the grid point
                %f = find( vols >= (OCobj.mBinsVol(jj)/max(OCobj.mBinsVol)) ); % patient at the grid point
                    
                f = find( vols >= vol_bins(kk)); % patient at the grid point
                    g = find( flgcomp(f) );

                    % for each time point
                        % total patients
                        OCobj.mAtlasTotal_DVH(ii,kk) = length(f);
                        % patients with complications
                        OCobj.mAtlasComp_DVH(ii,kk) = length(g);
                end
            end
        end
                
              
         function OCobj = fBetaCumulativeProbability_DVH(OCobj)
            OCobj.mBetaCumulativeMat = zeros( [size(OCobj.mAtlasTotal_DVH), size(OCobj.mBetaCumulativeThreshold,1)] );
            for k = 1:size(OCobj.mBetaCumulativeThreshold,1)
                OCobj.mBetaCumulativeMat(:,:,k) = betacdf( OCobj.mBetaCumulativeThreshold(k), OCobj.mAtlasComp_DVH+1, OCobj.mAtlasTotal_DVH - OCobj.mAtlasComp_DVH + 1 );
            end
        end
        
        
        
        
        function OCobj = fBetaCumulativeProbability_EUD(OCobj)
            OCobj.mBetaCumulativeMat = zeros( [size(OCobj.mAtlasTotal), size(OCobj.mBetaCumulativeThreshold,1)] );
            for k = 1:size(OCobj.mBetaCumulativeThreshold,1)
                OCobj.mBetaCumulativeMat(:,:,k) = betacdf( OCobj.mBetaCumulativeThreshold(k), OCobj.mAtlasComp+1, OCobj.mAtlasTotal - OCobj.mAtlasComp + 1 );
            end
        end
        function OCobj = fBetaInverseProbability_EUD(OCobj)
            OCobj.mBetaInverseMat = zeros( [size(OCobj.mAtlasTotal), size(OCobj.mBetaInverseThreshold,1)] );
            for k=1:length(OCobj.mBetaInverseThreshold)
                OCobj.mBetaInverseMat(:,:,k) = betainv( OCobj.mBetaInverseThreshold(k), OCobj.mAtlasComp+1, OCobj.mAtlasTotal - OCobj.mAtlasComp + 1 );
            end
        end
        
         function OCobj = fBetaInverseProbability_DVH(OCobj)
            OCobj.mBetaInverseMat = zeros( [size(OCobj.mAtlasTotal_DVH), size(OCobj.mBetaInverseThreshold,1)] );
            for k=1:length(OCobj.mBetaInverseThreshold)
                OCobj.mBetaInverseMat(:,:,k) = betainv( OCobj.mBetaInverseThreshold(k), OCobj.mAtlasComp_DVH+1, OCobj.mAtlasTotal_DVH - OCobj.mAtlasComp_DVH + 1 );
            end
        end
        
        
        function OCobj = fLogisticRegressionExact_EUD(OCobj)
            OCobj.mLogisticRegressionMat = repmat( struct('b',[],'dev',[],'stats',[]), [size(OCobj.mLymanN,1),1] );
            % mLogisticRegressionMat: (nLymanN X 1) matrix that holds
            % results of lyman fit for each logn value
            
            % using exact EUD
            euds = [OCobj.mGrp.mEUD]';
            pttotal = ones(OCobj.mNumInGrp,1);
            ptcomp = ones(OCobj.mNumInGrp,1); ptcomp([OCobj.mGrp.mFlgCensor])=0;
            for k=1:size(OCobj.mLymanN,1) % loop over each logn value, calculate fit of gEUDs
                doses=euds(:,k);
                % regression using exact EUD
                [b,dev,s]=glmfit(doses,[ptcomp pttotal],'binomial','link','logit');
                OCobj.mLogisticRegressionMat(k).b=b;
                OCobj.mLogisticRegressionMat(k).dev=dev;
                OCobj.mLogisticRegressionMat(k).stats=s;
            end
        end
        function OCobj = fLogisticRegressionGridExact_EUD(OCobj)
            % preparation
            euds = [OCobj.mGrp.mEUD]';
            flg = [OCobj.mGrp.mFlgCensor]';
            
            b0 = OCobj.mLogisticRegressionGridBetaRange{1};
            b1 = OCobj.mLogisticRegressionGridBetaRange{2};
            
            pttotal = ones(OCobj.mNumInGrp,1);
            ptcomp = ones(OCobj.mNumInGrp,1); 
            ptcomp([OCobj.mGrp.mFlgCensor])=0; 
            
            % for each mLymanN,b0, and b1, compute the log likelihood
            loglikelihood = -inf(length(b0),length(b1),length(OCobj.mLymanN));
            
            
            pvals = -inf(length(OCobj.mLymanN),1);    

            for kk = 1:length(OCobj.mLymanN)
                %inefficient...
                [~,~,s]=glmfit(euds(:,kk),[ptcomp pttotal],'binomial','link','logit');
                pval = s.p;
                pvals(kk) = pval(2);
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
            OCobj.mLogisticRegressionGrid = struct('b0',b0,'b1',b1,'loglikelihood',loglikelihood,'pvals',pvals);
        end
        function OCobj = fLogisticRegressionGridExactMixtureModel_EUD(OCobj,do_print,fig_loc)
            % sa is kaplan meier class
            
            % preparation
            euds = [OCobj.mGrp.mEUD]';
            flg = [OCobj.mGrp.mFlgCensor]';
            
                        pttotal = ones(OCobj.mNumInGrp,1);

            ptcomp = ones(OCobj.mNumInGrp,1); 
            ptcomp([OCobj.mGrp.mFlgCensor])=0;

            f2 = ~cellfun('isempty',{OCobj.mGrp.mDateComp}); % patients with no complication date
            f3 = ~cellfun('isempty',{OCobj.mGrp.mDateLastFollowup}); % patients with no last follow up date
            compdate = inf(OCobj.mNumInGrp,1);
            lastfollowup = inf(OCobj.mNumInGrp,1);
            compdate(f2) = ([OCobj.mGrp(f2).mDateComp] - [OCobj.mGrp(f2).mDateBaseline])' / 30;
            lastfollowup(f3) = ([OCobj.mGrp(f3).mDateLastFollowup] - [OCobj.mGrp(f3).mDateBaseline])' / 30;
            compdate = min( lastfollowup, compdate );
            
            % fit with just complication dates only
            compdates = compdate(~flg);
            compdates(compdates<=0)=0.00001;
            mu_sig = lognfit(sort(compdates));
            
            if do_print %deb
                
                screen_size=get(0,'ScreenSize');
                cur_fig=figure(1000);clf reset;hold on;
                set(cur_fig,'Position',[0 0 screen_size(3)/2 screen_size(4)/2]);
                x_axis = linspace(0,max(compdates)+1,500);
                
                [ax,h1,h2]=plotyy(x_axis,lognpdf(x_axis,mu_sig(1),mu_sig(2)),...
                    x_axis,logncdf(x_axis,mu_sig(1),mu_sig(2)));hold on;
                set(h1,'LineWidth',2);
                set(h2,'LineWidth',2);
                %set(ax(1),'YLim',[0.001 1]);
                set(get(ax(1),'Ylabel'),'String','P','FontSize',15);
                
                %set(ax(2),'YLim',[1 10]);
                set(get(ax(2),'Ylabel'),'String','P','FontSize',15);
                %semilogy(x_axis,repmat(0.05,length(x_axis)),'r--','LineWidth',1);
                set(ax(1),'XLim',[min(x_axis) max(x_axis)]);
                set(ax(2),'XLim',[min(x_axis) max(x_axis)]);
                set(ax,'FontSize',14);
                set(get(ax(2),'XLabel'),'String','Complication Date','FontSize',15);
                title('Log-Normal fits to complication data','FontSize',15);
                lgnd=legend([h1 h2],'PDF','CDF','Location','Best');
                set(lgnd,'FontSize',15);
                
                set(cur_fig,'Color','w');
                export_fig(cur_fig,[fig_loc,'log_normal_fits'],'-pdf');
            end;
            [sort_dates,dates_ind] = sort(compdate);
            f_t = lognpdf(sort_dates,mu_sig(1),mu_sig(2));
            F_t = logncdf(sort_dates,mu_sig(1),mu_sig(2));
            
%             f_t = lognpdf(compdate,mu_sig(1),mu_sig(2));
%             F_t = logncdf(compdate,mu_sig(1),mu_sig(2));
%             
            f_t(f_t==0)=0.0001;
            F_t(F_t==0)=0.0001;
            
            b0 = OCobj.mLogisticRegressionGridBetaRange{1};
            b1 = OCobj.mLogisticRegressionGridBetaRange{2};
            
            pvals = -inf(length(OCobj.mLymanN),1);
%             tstats = -inf(length(OCobj.mLymanN),1);
%             
       
            loglikelihood = -inf(length(b0),length(b1),length(OCobj.mLymanN));
            for ii = 1:length(OCobj.mLymanN)
%                  %inefficient...
%                 [~,~,s]=glmfit(euds(:,ii),[ptcomp pttotal],'binomial','link','logit');
%                 cur_se = s.se(2);
%                 
% %                 pval = s.p;
% %                 pvals(ii) = pval(2);
%                 
%                 cur_tstats = zeros(length(b0),length(b1));
%                 cur_pvals = zeros(length(b0),length(b1));
                for jj=1:length(b0)
                    for kk=1:length(b1)
                        cur_euds = euds(:,ii);
                        sorted_euds = cur_euds(dates_ind);
                        tmp_p = exp(b0(jj)+b1(kk)*sorted_euds)./(1+exp(b0(jj)+b1(kk)*sorted_euds));
                        ll = sum(ptcomp(dates_ind).*log(tmp_p.*f_t)+(1-ptcomp(dates_ind)).*log((1-(tmp_p.*F_t))));
                        loglikelihood(jj,kk,ii) = ll;
                        
%                         cur_tstats(jj,kk) = b1(kk)/cur_se;
%                         cur_pvals(jj,kk) = 2 * normcdf(-abs(cur_tstats(jj,kk)));
                        
                    end;
                end;
                % only store t/pvalue for max llhd for given n
%                 tmp_llhds = loglikelihood(:,:,ii);
%                 [~,loc] = max(tmp_llhds(:));
%                 [b0_ind,b1_ind] = ind2sub(size(tmp_llhds),loc);
%                 pvals(ii) = cur_pvals(b0_ind,b1_ind);
%                 tstats(ii) = cur_tstats(b0_ind,b1_ind);
%                 
            end;
            OCobj.mLogisticRegressionGridMixtureModel = struct('b0',b0,'b1',b1,'loglikelihood',loglikelihood,'pvals',pvals);
        end
        
        function OCobj = fLogisticRegressionGridExactMixtureModel_2_EUD(OCobj,do_print,fig_loc)
            % sa is kaplan meier class
            
            % preparation
            euds = [OCobj.mGrp.mEUD]';
            flg = [OCobj.mGrp.mFlgCensor]';
            
                        pttotal = ones(OCobj.mNumInGrp,1);

            ptcomp = ones(OCobj.mNumInGrp,1); 
            ptcomp([OCobj.mGrp.mFlgCensor])=0;

            f2 = ~cellfun('isempty',{OCobj.mGrp.mDateComp}); % patients with no complication date
            f3 = ~cellfun('isempty',{OCobj.mGrp.mDateLastFollowup}); % patients with no last follow up date
            compdate = inf(OCobj.mNumInGrp,1);
            lastfollowup = inf(OCobj.mNumInGrp,1);
            compdate(f2) = ([OCobj.mGrp(f2).mDateComp] - [OCobj.mGrp(f2).mDateBaseline])' / 30;
            lastfollowup(f3) = ([OCobj.mGrp(f3).mDateLastFollowup] - [OCobj.mGrp(f3).mDateBaseline])' / 30;
            compdate = min( lastfollowup, compdate );
            
            % fit with just complication dates only
            compdates = compdate(~flg);
            compdates(compdates<=0)=0.00001;
           % mu_sig = lognfit(sort(compdates));
            
            if 0 %deb
                
                screen_size=get(0,'ScreenSize');
                cur_fig=figure(1000);clf reset;hold on;
                set(cur_fig,'Position',[0 0 screen_size(3)/2 screen_size(4)/2]);
                x_axis = linspace(0,max(compdates)+1,500);
                
                [ax,h1,h2]=plotyy(x_axis,lognpdf(x_axis,mu_sig(1),mu_sig(2)),...
                    x_axis,logncdf(x_axis,mu_sig(1),mu_sig(2)));hold on;
                set(h1,'LineWidth',2);
                set(h2,'LineWidth',2);
                %set(ax(1),'YLim',[0.001 1]);
                set(get(ax(1),'Ylabel'),'String','P','FontSize',15);
                
                %set(ax(2),'YLim',[1 10]);
                set(get(ax(2),'Ylabel'),'String','P','FontSize',15);
                %semilogy(x_axis,repmat(0.05,length(x_axis)),'r--','LineWidth',1);
                set(ax(1),'XLim',[min(x_axis) max(x_axis)]);
                set(ax(2),'XLim',[min(x_axis) max(x_axis)]);
                set(ax,'FontSize',14);
                set(get(ax(2),'XLabel'),'String','Complication Date','FontSize',15);
                title('Log-Normal fits to complication data','FontSize',15);
                lgnd=legend([h1 h2],'PDF','CDF','Location','Best');
                set(lgnd,'FontSize',15);
                
                set(cur_fig,'Color','w');
                export_fig(cur_fig,[fig_loc,'log_normal_fits'],'-pdf');
            end;

             [sort_dates,dates_ind] = sort(compdate);
%             f_t = lognpdf(sort_dates,mu_sig(1),mu_sig(2));
%             F_t = logncdf(sort_dates,mu_sig(1),mu_sig(2));
%             f_t(f_t==0)=0.0001;
%             F_t(F_t==0)=0.0001;
% 
        %% Get cumulative incidenc
        
            function f=loggauss(x,xdata)
                mu=x(1);
                sigma=x(2);
                f=normcdf(log(xdata),mu,sigma);
            end
            
        test_sigma=1;
        test_mu=1;
        x0=[test_mu test_sigma];
        
        
         sa=OCobj.mKaplanMeierCompOverall;
        comp_time = sa.mSurvivalTimeSorted{1};
        comp_time(comp_time<=0)=0.001;
        comp_cinc = 1-sa.mSurvivalCurve{1};
 %       comp_inc = [0;diff(comp_cinc)];
        
        x=lsqcurvefit(@loggauss,x0,comp_time,comp_cinc);
   
        F_t = logncdf(sort_dates,x(1),x(2));     
     
        F_t = max(F_t)-F_t;

        
            b0 = OCobj.mLogisticRegressionGridBetaRange{1};
            b1 = OCobj.mLogisticRegressionGridBetaRange{2};
            pvals = -inf(length(OCobj.mLymanN),1);
            loglikelihood = -inf(length(b0),length(b1),length(OCobj.mLymanN));
            for ii = 1:length(OCobj.mLymanN)
%                  %inefficient...
%                 [~,~,s]=glmfit(euds(:,ii),[ptcomp pttotal],'binomial','link','logit');
%                 cur_se = s.se(2);
%                 4
% %                 pval = s.p;
% %                 pvals(ii) = pval(2);
%                 
%                 cur_tstats = zeros(length(b0),length(b1));
%                 cur_pvals = zeros(length(b0),length(b1));
                for jj=1:length(b0)
                    for kk=1:length(b1)
                        cur_euds = euds(:,ii);
                        sorted_euds = cur_euds(dates_ind);
                        tmp_p = exp(b0(jj)+b1(kk)*sorted_euds)./(1+exp(b0(jj)+b1(kk)*sorted_euds));
                        %ll = sum(ptcomp(dates_ind).*log(tmp_p.*f_t)+(1-ptcomp(dates_ind)).*log((1-(tmp_p.*F_t))));
                        %ll = sum(ptcomp(dates_ind).*log(tmp_p)+(1-ptcomp(dates_ind)).*log((1-(tmp_p.*F_t))));
                        ll = sum(ptcomp(dates_ind).*log(tmp_p)+(1-ptcomp(dates_ind)).*log((1-tmp_p)+tmp_p.*F_t));
                        loglikelihood(jj,kk,ii) = ll;
                        
%                         cur_tstats(jj,kk) = b1(kk)/cur_se;
%                         cur_pvals(jj,kk) = 2 * normcdf(-abs(cur_tstats(jj,kk)));
                        
                    end;
                end;
                % only store t/pvalue for max llhd for given n
%                 tmp_llhds = loglikelihood(:,:,ii);
%                 [~,loc] = max(tmp_llhds(:));
%                 [b0_ind,b1_ind] = ind2sub(size(tmp_llhds),loc);
%                 pvals(ii) = cur_pvals(b0_ind,b1_ind);
%                 tstats(ii) = cur_tstats(b0_ind,b1_ind);
%                 
            end;
            OCobj.mLogisticRegressionGridMixtureModel = struct('b0',b0,'b1',b1,'loglikelihood',loglikelihood,'pvals',pvals);
        end
        
        function OCobj = fLogisticRegressionGoodnessOfFitSimulationExact_EUD(OCobj)
            % loga determination
            st = [OCobj.mLogisticRegressionMat];
            dpf = [st.dev]; % deviations
            [~,loc] = min(dpf); % the minimum deviation corresponds to the maximum likelihood
            % fitting result
            st = st(loc);
            % euds
            euds = [OCobj.mGrp.mEUD]; euds = euds(loc,:)'; % the gEUDs of that n

            % expectation for each patient
            numE = glmval(st.b,euds,'logit',st.stats);

            % SSR at simulated obserations
            numSim = 100000; % number of simulations
            ssr = zeros(numSim,1); % SSR
            for k = 1:numSim
                % generate observations
                numO = rand([OCobj.mNumInGrp,1]); % random number with uniform distribution on (0,1)
                f = numO<=numE; numO(f) = 1; numO(~f) = 0; % the lower the expectation, the lower the probability it has a complication
                % SSR
                ssr(k) = sum( (numO-numE).^2 );
            end

            % G-value
            flg=[OCobj.mGrp.mFlgCensor]'; % censor flags of patients
            ssrObsv = sum( (~flg-numE).^2 ); % the actual ssr
            f = ssr > ssrObsv; % the proportion that simulated ssr is greater than the observed SSR
            p = sum(f)/numSim; % p close to 1 means overfitting, p close to 0 means the model dose not fit the data

            disp(['The best loga of Logistic Regression is: ',num2str(-log10(OCobj.mLymanN(loc))),...
                'with G=',num2str(p)]);

            % save
            OCobj.mLogisticRegressionGoodnessOfFitSim.SSRSim = ssr;
            OCobj.mLogisticRegressionGoodnessOfFitSim.SSRObserve = ssrObsv;
            OCobj.mLogisticRegressionGoodnessOfFitSim.p_value = p;
        end

        function OCobj = fLymanAnalysisGridExact_EUD(OCobj,binning)
            % preparation
            euds = [OCobj.mGrp.mEUD]';
            flg = [OCobj.mGrp.mFlgCensor]';
            
            if isempty(OCobj.mLymanGridTD50Range)
                % fine
                if isequal(binning,'fine')
                    OCobj.mLymanGridTD50Range = (0:0.1:max(euds(:)))';
                elseif isequal(binning,'med')
                %med
                    OCobj.mLymanGridTD50Range = (0:0.5:max(euds(:)))';
                elseif isequal(binning,'coarse')
                %coarse
                    OCobj.mLymanGridTD50Range = (0:1:max(euds(:)))';
                else
                    OCobj.mLymanGridTD50Range = (0:0.1:max(euds(:)))';
                end
            end
            if isempty(OCobj.mLymanGridMRange)
                if isequal(binning,'fine')
                    OCobj.mLymanGridMRange = (0:0.01:1.5)';
                    OCobj.mLymanGridMRange(1) = 0.0001;
                elseif isequal(binning,'med')
                    OCobj.mLymanGridMRange = (0:0.05:1.5)';
                    OCobj.mLymanGridMRange(1) = 0.001;
                elseif isequal(binning,'coarse')
                    OCobj.mLymanGridMRange = (0:0.1:1.5)';
                    OCobj.mLymanGridMRange(1) = 0.001;
                else
                    OCobj.mLymanGridMRange = (0:0.01:1.5)';
                    OCobj.mLymanGridMRange(1) = 0.0001;
                end
               end
            TD50 = OCobj.mLymanGridTD50Range;
            m = OCobj.mLymanGridMRange;
            pttotal = ones(OCobj.mNumInGrp,1);
            ptcomp = ones(OCobj.mNumInGrp,1); 
            ptcomp([OCobj.mGrp.mFlgCensor])=0;    
            % for each mLymanN,TD50, and m, compute the log likelihood
            loglikelihood = -inf(length(TD50),length(m),length(OCobj.mLymanN));
            
            % test
            pvals = -inf(length(OCobj.mLymanN),1);
            %pvals = -inf(length(TD50),length(m),length(OCobj.mLymanN));
            
            for kk = 1:length(OCobj.mLymanN)
                % pvals from probit fit
                %[~,~,s2]=glmfit(euds(:,kk),[ptcomp
                %pttotal],'binomial','link','probit');
                %cur_pvals = s2.p;    
                %pvals(kk)=cur_pvals(2);
                
                for jj = 1:length(m)
                    for ii = 1:length(TD50)
                
                        pr = normcdf((euds(:,kk)-TD50(ii))/(m(jj)*TD50(ii)),0,1); % Lyman probability
                        pr(flg) = 1-pr(flg); % non-complication patients
                        pr = log(pr); % log likelihood of each patients
                        loglikelihood(ii,jj,kk) = sum(pr); % loglikelihood of all
                    
                    end
                end
                    % find best TD50 and m
              % loga determination
            [~,loc] = max(loglikelihood(:));
            [dd,mm,loc] = ind2sub(size(loglikelihood),loc);
            % fitting result
            cur_TD50 = TD50(dd);
            cur_m = m(mm);
            [~,~,s2]=glmfit(((euds(:,kk)-cur_TD50)/(cur_m*cur_TD50)),...
                            [ptcomp pttotal],'binomial','link','probit');        

            cur_pvals = s2.p;
            pvals(kk)=cur_pvals(2);   
            
                        
            end
            
             
%mLymanGrid.TD50
                        
            OCobj.mLymanGrid = struct('TD50',TD50,'m',m,'loglikelihood',loglikelihood,'pvals',pvals);
        end
        function OCobj = fLymanGoodnessOfFitSimulationExact_EUD(OCobj)
           
            
            % loga determination
            [~,loc] = max(OCobj.mLymanGrid.loglikelihood(:));
            [dd,mm,loc] = ind2sub(size(OCobj.mLymanGrid.loglikelihood),loc);
           
            
            
              % find the best "a"
            [mx,loc] = max(OCobj.mLymanGrid.loglikelihood(:));
            [~,~,loc] = ind2sub(size(OCobj.mLymanGrid.loglikelihood),loc);
            disp(['the maximum log likelihood and its corresponding "log10(a)" in Lyman model of exact gEUD are: ',num2str([mx, -OCobj.mLymanN(loc)])]);
            
            % fitting result
            TD50 = OCobj.mLymanGrid.TD50(dd);
            m = OCobj.mLymanGrid.m(mm);
            % euds
            euds = [OCobj.mGrp.mEUD]; euds = euds(loc,:)'; % the gEUDs of that n

            % calculate p-value from probit
            %pttotal = ones(OCobj.mNumInGrp,1);
            %ptcomp = ones(OCobj.mNumInGrp,1); 
            %ptcomp([OCobj.mGrp.mFlgCensor])=0;
            %[~,~,s2]=glmfit(euds,[ptcomp pttotal],'binomial','link','probit');
            %p2_vals=s2.p;
            
            % expectation for each patient
            numE = normcdf((euds-TD50)/(m*TD50),0,1); % Lyman probability

            % SSR at simulated obserations
            numSim = 1000000; % number of simulations
            ssr = zeros(numSim,1); % SSR
            for k = 1:numSim
                % generate observations
                numO = rand([OCobj.mNumInGrp,1]); % random number with uniform distribution on (0,1)
                f = numO<=numE; numO(f) = 1; numO(~f) = 0; % the lower the expectation, the lower the probability it has a complication
                % SSR
                ssr(k) = sum( (numO-numE).^2 );
            end

            % G-value
            flg=[OCobj.mGrp.mFlgCensor]'; % censor flags of patients
            ssrObsv = sum( (~flg-numE).^2 ); % the actual ssr
            f = ssr > ssrObsv; % the proportion that simulated ssr is greater than the observed SSR
            p = sum(f)/numSim; % p close to 1 means overfitting, p close to 0 means the model dose not fit the data

            disp(['the "log10(a)" for the maximum likelihood in Lyman model is: ',num2str(-OCobj.mLymanN(loc)),...
                ' with G = ',num2str(p)]);

            % save
            OCobj.mLymanGoodnessOfFitSim.SSRSim = ssr;
            OCobj.mLymanGoodnessOfFitSim.SSRObserve = ssrObsv;
            OCobj.mLymanGoodnessOfFitSim.p_value = p;
            %OCobj.mLymanGoodnessOfFitSim.p2_value = p2_vals(2);
        end
    
    
    function OCobj = fLymanGoodnessOfFitSimulationExact_plot_EUD(OCobj)
           
    pvals = inf(21,1);
       
    for i=1:21,
        loc = i;
        
         loga = -OCobj.mLymanN(loc);
           
        % coefficients for the Lyman model
        [~,loc] = min(abs(OCobj.mLymanN+loga)); % the n whose corresponding responding function will be ploted
        ll = OCobj.mLymanGrid.loglikelihood(:,:,loc); % log likelihood of log10(a) = loga
        mx = max(ll(:));
        [xx,yy] = find(ll == mx); % the coefficients
        TD50 = OCobj.mLymanGrid.TD50(xx);
        m = OCobj.mLymanGrid.m(yy);
        
%         [dd,mm,~] = ind2sub(size(OCobj.mLymanGrid.loglikelihood),loc);
%         % fitting result
%         TD50 = OCobj.mLymanGrid.TD50(dd);
%         m = OCobj.mLymanGrid.m(mm);
            % euds
            euds = [OCobj.mGrp.mEUD]; euds = euds(loc,:)'; % the gEUDs of that n

            % calculate p-value from probit
            %pttotal = ones(OCobj.mNumInGrp,1);
            %ptcomp = ones(OCobj.mNumInGrp,1); 
            %ptcomp([OCobj.mGrp.mFlgCensor])=0;
            %[~,~,s2]=glmfit(euds,[ptcomp pttotal],'binomial','link','probit');
            %p2_vals=s2.p;
            
            % expectation for each patient
            numE = normcdf((euds-TD50)/(m*TD50),0,1); % Lyman probability

            % SSR at simulated obserations
            numSim = 100000; % number of simulations
            ssr = zeros(numSim,1); % SSR
            for k = 1:numSim
                % generate observations
                numO = rand([OCobj.mNumInGrp,1]); % random number with uniform distribution on (0,1)
                f = numO<=numE; numO(f) = 1; numO(~f) = 0; % the lower the expectation, the lower the probability it has a complication
                % SSR
                ssr(k) = sum( (numO-numE).^2 );
            end

            % G-value
            flg=[OCobj.mGrp.mFlgCensor]'; % censor flags of patients
            ssrObsv = sum( (~flg-numE).^2 ); % the actual ssr
            f = ssr > ssrObsv; % the proportion that simulated ssr is greater than the observed SSR
            p = sum(f)/numSim; % p close to 1 means overfitting, p close to 0 means the model dose not fit the data

            % save
            pvals(i)=p;
                        %OCobj.mLymanGoodnessOfFitSim.p2_value = p2_vals(2);
    
    end
    
    
    cur_fig=figure(1111);clf reset;
    screen_size=get(0,'ScreenSize');
    set(cur_fig,'Position',[0 0 screen_size(3) screen_size(4)]);
    plot(-OCobj.mLymanN,pvals,'ko--','MarkerSize',18,'LineWidth',2);
    grid on;
    ylabel('Goodness of Fit','FontSize',24);
    xlabel('log_1_0(a)','FontSize',24);
    set(gca,'FontSize',14);
    set(cur_fig,'Color','w');
    fig_loc = 'Z:/elw/MATLAB/meta_analy/slides/figures/latest/';
    export_fig(cur_fig,[fig_loc,'gof'],'-pdf');
    
    
    end
    
    end

    methods % EUD atlas
        function OCobj = fLogisticRegressionBin_EUD(OCobj)
            OCobj.mLogisticRegressionMatBin = repmat( struct('b',[],'dev',[],'stats',[]), [size(OCobj.mLymanN,1),1] );
            
            % using bins
            if isempty(OCobj.mBinsDose) || isempty(OCobj.mAtlasTotal)
                error('mBinsDose or mAtlasTotal not determined before method "LogitAnalysisBin_EUD" is called (classOutcomeAnalysis)');
            end
            
            dosebins = (OCobj.mBinsDose(1:end-1)+OCobj.mBinsDose(2:end))/2; % dose bins are at the middle of the intervals
            pttotal = ones( OCobj.mAtlasTotal(1,1), 1 ); % each patient has his own row
            ptcomp = true( OCobj.mAtlasTotal(1,1), 1 ); % allocate space for complication of each patient
            doseval = -inf(OCobj.mNumInGrp,1); % allocate space for dose of each patient.

            % patient complication and censor info
            pta = OCobj.mAtlasTotal; % patient total from atlas
            pca = OCobj.mAtlasComp; % patient complication from atlas
            pta = pta - pca; % censored info of patients
            
            pta(1:end-1,:) = abs(diff(pta,1,1)); pta(end,:) = 0; % locations of patients with censor info
            pca(1:end-1,:) = abs(diff(pca,1,1)); pca(end,:) = 0; % locations of patient complications
            
            warning('off','stats:glmfit:IterationLimit');
            for k=1:length(OCobj.mLymanN) % Logistic Regression for each ln(n)
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
                OCobj.mLogisticRegressionMatBin(k).b=b;
                OCobj.mLogisticRegressionMatBin(k).dev=dev;
                OCobj.mLogisticRegressionMatBin(k).stats=s;
            end
            warning('on','stats:glmfit:IterationLimit');



%             OCobj.mLogisticRegressionMatBin = repmat( struct('b',[],'dev',[],'stats',[]), [size(OCobj.mLymanN,1),1] );
%             
%             % using bins
%             if isempty(OCobj.mBinsDose)
%                 error('mBinsDose not determined before method "LogitAnalysis_EUD" is called (classOutcomeAnalysis)');
%             end
%             dosebins = (OCobj.mBinsDose(1:end-1)+OCobj.mBinsDose(2:end))/2; % dose bins are at the middle of the intervals
%             pttotal = ones( OCobj.mNumInGrp, 1 ); % each patient has his own row
%             ptcomp = ones( OCobj.mNumInGrp, 1 ); ptcomp(OCobj.mGrp.mFlgCensor)=0; % allocate space for complication of each patient
%             doseval = -inf(OCobj.mNumInGrp,1); % allocate space for dose of each patient.
%             euds = [OCobj.mGrp.mEUD]';
%             warning('off','stats:glmfit:IterationLimit');
%             for k=1:length(OCobj.mLymanN)
%                 % determine the dose for each patient
%                 for m = 1:length(OCobj.mBinsDose)-1
%                     f = euds(:,k)>=OCobj.mBinsDose(m);
%                     doseval(f) = dosebins(m);
%                 end
%                 f = euds(:,k)>=OCobj.mBinsDose(end); % patient with dose outside the required dose bins should be removed
%                 doseval(f) = min(OCobj.mBinsDose)-1;
%                 
%                 f = doseval>=min(OCobj.mBinsDose); % pick up the patients whose dose fall into the dose bins
%                 % logistic regression
%                 [b,dev,s]=glmfit(doseval(f),[ptcomp(f) pttotal(f)],'binomial','link','logit');
%                 OCobj.mLogisticRegressionMatBin(k).b=b;
%                 OCobj.mLogisticRegressionMatBin(k).dev=dev;
%                 OCobj.mLogisticRegressionMatBin(k).stats=s;
%             end
%             warning('on','stats:glmfit:IterationLimit');
        end
        
        function OCobj = fLogisticRegressionGoodnessOfFitSimulationBin_EUD(OCobj)
            % loga determination
            st = [OCobj.mLogisticRegressionMatBin];
            dpf = [st.dev]; % deviations
            [~,loc] = min(dpf); % the minimum deviation corresponds to the maximum likelihood
            % fitting result
            st = st(loc);
            % euds
            euds = [OCobj.mGrp.mEUD]; euds = euds(loc,:)'; % the gEUDs of that n

            % expectation for each patient
            numE = glmval(st.b,euds,'logit',st.stats);

            % SSR at simulated obserations
            numSim = 1000000; % number of simulations
            ssr = zeros(numSim,1); % SSR
            for k = 1:numSim
                % generate observations
                numO = rand([OCobj.mNumInGrp,1]); % random number with uniform distribution on (0,1)
                f = numO<=numE; numO(f) = 1; numO(~f) = 0; % the lower the expectation, the lower the probability it has a complication
                % SSR
                ssr(k) = sum( (numO-numE).^2 );
            end

            % p-value
            flg=[OCobj.mGrp.mFlgCensor]'; % censor flags of patients
            ssrObsv = sum( (~flg-numE).^2 ); % the actual ssr
            f = ssr > ssrObsv; % the proportion that simulated ssr is greater than the observed SSR
            p = sum(f)/numSim; % p close to 1 means overfitting, p close to 0 means the model dose not fit the data

            disp(['The best loga of Logistic Regression based on atlas is: ',num2str(-OCobj.mLymanN(loc)),...
                'with p=',num2str(p)]);

            % save
            OCobj.mLogisticRegressionGoodnessOfFitSimBin.SSRSim = ssr;
            OCobj.mLogisticRegressionGoodnessOfFitSimBin.SSRObserve = ssrObsv;
            OCobj.mLogisticRegressionGoodnessOfFitSimBin.p_value = p;
        end
        function OCobj = fLogisticRegressionGridBin_EUD(OCobj)
            % parse atlas (part)
            dosebins = (OCobj.mBinsDose(1:end-1)+OCobj.mBinsDose(2:end))/2; % dose bins are at the middle of the intervals
            ptcomp = true( OCobj.mAtlasTotal(1,1), 1 ); % allocate space for complication of each patient
            doseval = -inf(OCobj.mNumInGrp,1); % allocate space for dose of each patient.

            pta = OCobj.mAtlasTotal; % patient total from atlas
            pca = OCobj.mAtlasComp; % patient complication from atlas
            pta = pta - pca; % censored info of patients
            pta(1:end-1,:) = abs(diff(pta,1,1)); pta(end,:) = 0; % locations of patients with censor info
            pca(1:end-1,:) = abs(diff(pca,1,1)); pca(end,:) = 0; % locations of patient complications

            % preparation
            b0 = OCobj.mLogisticRegressionGridBetaRange{1};
            b1 = OCobj.mLogisticRegressionGridBetaRange{2};

            % for each mLymanN, TD50, and m, compute the log likelihood
            loglikelihood = -inf(length(b0),length(b1),length(OCobj.mLymanN));
            for kk = 1:length(OCobj.mLymanN)
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
            OCobj.mLogisticRegressionGridBin = struct('b0',b0,'b1',b1,'loglikelihood',loglikelihood);
        end

        function OCobj = fLymanAnalysisGridBin_EUD(OCobj)
            % preparation
            dosebins = (OCobj.mBinsDose(1:end-1)+OCobj.mBinsDose(2:end))/2; % dose bins are at the middle of the intervals
            ptcomp = true( OCobj.mAtlasTotal(1,1), 1 ); % allocate space for complication of each patient
            doseval = -inf(OCobj.mNumInGrp,1); % allocate space for dose of each patient.
            pttotal = ones(OCobj.mNumInGrp,1);
            
            pta = OCobj.mAtlasTotal; % patient total from atlas
            pca = OCobj.mAtlasComp; % patient complication from atlas
            pta = pta - pca; % censored info of patients
            pta(1:end-1,:) = abs(diff(pta,1,1)); pta(end,:) = 0; % locations of patients with censor info
            pca(1:end-1,:) = abs(diff(pca,1,1)); pca(end,:) = 0; % locations of patient complications

            TD50 = OCobj.mLymanGridTD50Range;
            m = OCobj.mLymanGridMRange;
            
            pvals = -inf(length(OCobj.mLymanN),1);    
            % for each mLymanN, TD50, and m, compute the log likelihood
            loglikelihood = -inf(length(TD50),length(m),length(OCobj.mLymanN));
            for kk = 1:length(OCobj.mLymanN)
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

                % calc p-values from probit
                [~,~,s2]=glmfit(doseval,[ptcomp pttotal],'binomial','link','probit');
                cur_pvals = s2.p;    
                pvals(kk)=cur_pvals(2);
                
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
            OCobj.mLymanGridBin = struct('TD50',TD50,'m',m,'loglikelihood',loglikelihood,'pvals',pvals);
        end
        function OCobj = fLymanGoodnessOfFitSimulationBin_EUD(OCobj)
            % loga determination
            [~,loc] = max(OCobj.mLymanGridBin.loglikelihood(:));
            [dd,mm,loc] = ind2sub(size(OCobj.mLymanGridBin.loglikelihood),loc);
            % fitting result
            TD50 = OCobj.mLymanGridBin.TD50(dd);
            m = OCobj.mLymanGridBin.m(mm);
            % euds
            euds = [OCobj.mGrp.mEUD]; euds = euds(loc,:)'; % the gEUDs of that n

            % expectation for each patient
            numE = normcdf((euds-TD50)/(m*TD50),0,1); % Lyman probability

            % SSR at simulated obserations
            numSim = 1000000; % number of simulations
            ssr = zeros(numSim,1); % SSR
            for k = 1:numSim
                % generate observations
                numO = rand([OCobj.mNumInGrp,1]); % random number with uniform distribution on (0,1)
                f = numO<=numE; numO(f) = 1; numO(~f) = 0; % the lower the expectation, the lower the probability it has a complication
                % SSR
                ssr(k) = sum( (numO-numE).^2 );
            end

            % p-value
            flg=[OCobj.mGrp.mFlgCensor]'; % censor flags of patients
            ssrObsv = sum( (~flg-numE).^2 ); % the actual ssr
            f = ssr > ssrObsv; % the proportion that simulated ssr is greater than the observed SSR
            p = sum(f)/numSim; % p close to 1 means overfitting, p close to 0 means the model dose not fit the data

            disp(['the "log10(a)" for the maximum likelihood in Lyman model is: ',num2str(-OCobj.mLymanN(loc)),...
                'with p=',num2str(p)]);

            % save
            OCobj.mLymanGoodnessOfFitSimBin.SSRSim = ssr;
            OCobj.mLymanGoodnessOfFitSimBin.SSRObserve = ssrObsv;
            OCobj.mLymanGoodnessOfFitSimBin.p_value = p;
        end
    end

    methods % EUD plot and atlas writing
        function fEUDCurvesFig_a_EUD(OCobj)
            f=[OCobj.mGrp.mFlgCensor]; g=find(f);
            a1 = gca;
            a2 = copyobj(a1,gcf);
            set(a2,'Color','none');
            set(a2,'Xtick',[]);
            hold on;
            for m = 1:length(g)
                plot(a1,OCobj.mGrp(g(m)).mEUD,OCobj.mGrp(g(m)).mLymanN,'b');
            end
            g=find(~f);
            for m = 1:length(g)
                plot(a1,OCobj.mGrp(g(m)).mEUD,OCobj.mGrp(g(m)).mLymanN,'r');
            end
            set(a2,'YAxisLocation','right');
            set(a2,'YTickLabel',num2str(OCobj.mLymanN(1:2:end)));
            set(a1,'YTickLabel',num2str(OCobj.mLymanN(end:-2:1)));
            set(a1,'xminortick','on','yminortick','on');
            set(a1,'box','on');
            xlabel(a1,'EUD'); ylabel(a1,'log_1_0(a)');
        end
        function fEUDCurvesSummary_a_EUD(OCobj)
            % prepare
            f = [OCobj.mGrp.mFlgCensor]; % censor info
            numcensor = sum(f); % number of patients censored
            numcensor95 = round((1-0.95)*numcensor/2);
            numcensor68 = round((1-0.68)*numcensor/2);
            numcomp = sum(~f); % number of patients with complication
            numcomp95 = round((1-0.95)*numcomp/2);
            numcomp68 = round((1-0.68)*numcomp/2);
            
            eud = [OCobj.mGrp.mEUD]'; % all doses
            eud_censor  = zeros(length(OCobj.mLymanN),5); % 5 columns for -95%, -68%, median, 68%, 95% lines
            eud_comp = eud_censor;

            % eud computation
            for kk = 1:length(OCobj.mLymanN)
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
            plot(eud_censor(:,3),OCobj.mLymanN,'b*-');
            plot(eud_censor(:,[2,4]),OCobj.mLymanN,'b');
            plot(eud_censor(:,[1,5]),OCobj.mLymanN,'b--');
            plot(eud_comp(:,3),OCobj.mLymanN,'r*-');
            plot(eud_comp(:,[2,4]),OCobj.mLymanN,'r');
            plot(eud_comp(:,[1,5]),OCobj.mLymanN,'r--');
            hold off;

            set(a2,'YAxisLocation','right');
            set(a2,'YTickLabel',num2str(OCobj.mLymanN(1:2:end)));
            set(a1,'YTickLabel',num2str(OCobj.mLymanN(end:-2:1)));
            set(a1,'xminortick','on','yminortick','on');
            set(a1,'box','on');
            xlabel(a1,'Dose (Gy)'); ylabel(a1,'log_1_0(a)');
        end
        function fAtlasFig_EUD(OCobj)
            if isempty(OCobj.mAtlasTotal)
                disp('Empty member "mAtlasTotal", cannot display its figure.'); return;
            end
            dosestep = 5;
            doses=OCobj.mBinsDose; x=mod(doses,dosestep)==0;
            [xx,yy]=ndgrid(doses(x),1:length(OCobj.mLymanN)); xx=num2cell(xx); yy=num2cell(yy);
            strComp=(OCobj.mAtlasComp(x,:)); strTotal=(OCobj.mAtlasTotal(x,:)); strAtlas=arrayfun(@(a,b) strcat(num2str(a),'/',num2str(b)),strComp,strTotal,'UniformOutpu',false);
            cellfun(@(a,b,c) text(a,b,c,'fontsize',16),xx,yy,strAtlas);
            set(gca,'XLim',[xx{1,1}-1,xx{end,1}+2]); set(gca,'YLim',[1-1,length(OCobj.mLymanN)+1]);
            set(gca,'YTick',1:2:length(OCobj.mLymanN)); set(gca,'YTickLabel',OCobj.mLymanN(1:2:end));
            pos=get(gcf,'Position'); set(gcf,'Position',[pos(1:2), 900/3*4, 900]);
            set(gca,'xminortick','on','yminortick','on');
            set(gca,'box','on');
%             title([OCobj.xlsSheet,', number of patients and numbers with severe pneumonitis treated to at least a given EUD'],'fontsize',16);
            xlabel('gEUD doses (Gy)'); ylabel('log_1_0a'); % set(gca,'fontsize',16);
        end
        
        function fAtlasFig_DVH(OCobj,dose_step,vol_step,fontsize)
            if isempty(OCobj.mAtlasTotal_DVH)
                disp('Empty member "mAtlasTotal_DVH", cannot display its figure.'); return;
            end
            dosestep = dose_step;
            volstep = vol_step;
            
%             OCgrp = [OCobj.mGrp];
%             nfx_flg = ones(length(OCgrp),1);
%             
%             if nfx>0 % nfx=-1 -> all fx
%             
%                 grp_flgs = [OCgrp.mFxNum];
%                 nfx_flg = grp_flgs==nfx;
%                 
%                 if sum(nfx_flg)==0
%                     disp(['Invalid nfx parameter']);
%                     return;
%                 end
%             end
%             
            doses=OCobj.mBinsDose; x=mod(doses,dosestep)==0;
            %vols=100.*(OCobj.mBinsVol./max(OCobj.mBinsVol)); y=mod(vols,volstep)==0;
            vols=[0:5:100]; y=mod(vols,volstep)==0;
            [xx,yy]=ndgrid(doses(x),vols(y)); xx=num2cell(xx); yy=num2cell(yy);
            strComp=(OCobj.mAtlasComp_DVH(x,y)); strTotal=(OCobj.mAtlasTotal_DVH(x,y));
            strAtlas=arrayfun(@(a,b) strcat('$$\frac{',num2str(a),'}{',num2str(b),'}$$'),strComp,strTotal,'UniformOutpu',false);
            cellfun(@(a,b,c) text(a,b,c,'fontsize',fontsize,'HorizontalAlignment','center','Interpreter','latex'),xx,yy,strAtlas);
            set(gca,'XLim',[xx{1,1}-(dosestep/2),xx{end,1}+(dosestep/2)]); 
            set(gca,'YLim',[yy{1,1}-(volstep),yy{1,end}+(volstep)]);
            set(gca,'XTick',[0:dosestep:xx{end,1}]);
          
            pos=get(gcf,'Position'); set(gcf,'Position',[pos(1:2), 900/3*4, 900]);
            set(gca,'xminortick','on','yminortick','on');
            set(gca,'box','on');
%             title([OCobj.xlsSheet,', number of patients and numbers with severe pneumonitis treated to at least a given EUD'],'fontsize',16);
            xlabel('BED Dose [Gy_3]','FontSize',24); 
            ylabel('Volume [%]','FontSize',24); % set(gca,'fontsize',16);
            set(gca,'FontSize',16);
        end
        
        function fAtlasRotatedFig_EUD(OCobj,fonts,ticks)
            if ~exist('fonts','var')
                fonts = 10;
            end
            if ~exist('ticks','var')
                ticks = 10;
            end
            
            if isempty(OCobj.mAtlasTotal)
                disp('Empty member "mAtlasTotal", cannot display its figure.'); return;
            end
            % check columns with informaiton
            % complication part
            f = diff(OCobj.mAtlasComp);
            for k = 1:size(OCobj.mAtlasComp,1) % upper side
                if any(f(k,:))
                    colcompu = k;
                    break;
                end
            end
            for k = size(OCobj.mAtlasComp,1)-1 : -1 : 1 % lower side
                if any(f(k,:))
                    colcompl = k+1;
                    break;
                end
            end
            % total part
            f = diff(OCobj.mAtlasTotal);
            for k = 1:size(OCobj.mAtlasTotal,1) % upper side
                if any(f(k,:))
                    coltotalu = k;
                    break;
                end
            end
            for k = size(OCobj.mAtlasTotal,1)-1 : -1 : 1 % right side
                if any(f(k,:))
                    coltotall = k+1;
                    break;
                end
            end
            % combination of the comp and total
            colu = min( colcompu, coltotalu );
            coll = max( colcompl, coltotall );
            % prepare to write the table in a figure
%             doses = zeros(size(OCobj.mBinsDose));
            doses=OCobj.mBinsDose(colu:coll);
            [xx,yy]=ndgrid(doses,1:length(OCobj.mLymanN));
            xx=num2cell(xx); 
            
            
            yy=num2cell(yy);
            
            
            strComp=OCobj.mAtlasComp(colu:coll,:); strTotal=OCobj.mAtlasTotal(colu:coll,:);
            strAtlas=arrayfun(@(a,b) strcat(num2str(a),'/',num2str(b)),strComp,strTotal,'UniformOutpu',false);
%             figure(1);
            clf reset;
            cellfun(@(a,b,c) text(b,a,c,'fontsize',fonts),xx,yy,strAtlas);
            set(gca,'YLim',[xx{1,1}-1,xx{end,1}+2]);
            set(gca,'XLim',[1-1,length(OCobj.mLymanN)+1]);
            set(gca,'XTick',1:2:length(OCobj.mLymanN)); set(gca,'XTickLabel',OCobj.mLymanN(1:2:end));
            set(gca,'fontsize',ticks);
            set(gca,'xminortick','on','yminortick','on');
            set(gca,'Position',[0.04,0.04,0.955,0.955]);
            set(gca,'box','on');
%             pos=get(gcf,'Position'); set(gcf,'Position',[pos(1:2), 900/3*4, 900]);
%             title([OCobj.xlsSheet,', number of patients and numbers with severe pneumonitis treated to at least a given EUD'],'fontsize',16);
            xlabel('gEUD doses (Gy)'); ylabel('log_1_0a'); % set(gca,'fontsize',16);
        end
        function fAtlasCompactFig_EUD(OCobj,fonts,ticks,stepsz)
            if ~exist('fonts','var')
                fonts = 10;
            end
            if ~exist('ticks','var')
                ticks = 10;
            end
            
            if isempty(OCobj.mAtlasTotal)
                disp('Empty member "mAtlasTotal", cannot display its figure.'); return;
            end
            
            % check duplicated columns for each n to determine the talbe size
            dupcol = logical(diff(OCobj.mAtlasComp));
            dupcol = dupcol | logical(diff(OCobj.mAtlasTotal));
            shiftpos = zeros(2,length(OCobj.mLymanN));
            for n = 1:length(OCobj.mLymanN)
                f = find(dupcol(:,n));
                shiftpos(1,n) = f(1);
                shiftpos(2,n) = f(end);
            end
            f = diff(shiftpos); fy = max(f)+1;

            % prepare data
            strComp = zeros(fy+1,length(OCobj.mLymanN));
            strTotal = zeros(fy+1,length(OCobj.mLymanN));
            y = size(OCobj.mAtlasTotal,1);
            for n = 1:length(OCobj.mLymanN)
                x = min(shiftpos(1,n)+fy, y);
                strComp(1:1+x-shiftpos(1,n),n) = OCobj.mAtlasComp(shiftpos(1,n):x,n);
                strTotal(1:1+x-shiftpos(1,n),n) = OCobj.mAtlasTotal(shiftpos(1,n):x,n);
            end
            strAtlas = arrayfun(@(a,b) strcat(num2str(a),'/',num2str(b)),strComp,strTotal,'UniformOutput',false);
%             strshift = num2cell(shiftpos(1,:))';
            strshift = arrayfun(@(a) num2str(OCobj.mBinsDose(a)),shiftpos(1,:),'UniformOutput',false);

            % plot the table
            clf reset;
            step_flg = [0:fy];
            step_flg = mod(step_flg,stepsz)==0;
            [x,y] = ndgrid(0:stepsz:fy,1:length(OCobj.mLymanN));
            xx = num2cell(x); 
            yy = num2cell(y);
            
            cellfun(@(a,b,c) text(b,a,c,'fontsize',fonts),xx,yy,strAtlas(step_flg,:));
            set(gca,'YLim',[xx{1,1}-1,xx{end,1}+2]);
            set(gca,'XLim',[1-1,length(OCobj.mLymanN)+1]);
            set(gca,'XTick',1:1:length(OCobj.mLymanN)); set(gca,'XTickLabel',OCobj.mLymanN(end:-1:1));
            set(gca,'fontsize',ticks);

%             axis manual;
            cellfun(@(b,x) text(b,-4,x,'fontsize',fonts),yy(1,:),strshift);
            set(gca,'Position',[0.04,0.1,0.956,0.89]);
            set(gca,'xminortick','off','yminortick','on');
            set(gca,'box','on');
            xlabel('log_1_0a'); ylabel('gEUD doses (Gy)'); % set(gca,'fontsize',16);
        end        
        function [fig mn_llhd mx_llhd best_loga]= fLogisticRegressionLikelyhoodExactFig_a_EUD(OCobj,loga,strMarker,lw)
            if ~exist('strMarker','var')
                strMarker = 'b';
            elseif ~ischar(strMarker)
                warning('classOutcomeAnalysis:EndPointObservedFig_EUD:strMarker','strMarker is not a string, omit');
                strMarker = 'b';
            end
            if ~exist('lw','var')
                lw = 1;
            elseif ~isnumeric(lw)
                warning('classOutcomeAnalysis:EndPointObservedFig_EUD:linewidth','lw is not a number, omit');
                lw = 1;
            end
            st = [OCobj.mLogisticRegressionMat];
            dpf = [st.dev]; % deviations
            st =[st.stats];
            df = [st.dfe]; % degree of freedom
            dpf = dpf./df; % deviations per degree of freedom
            [~,loc] = min(dpf); % the minimum deviation corresponds to the maximum likelihood
            
            if ~exist('loga','var') || ischar(loga)
                loga = -OCobj.mLymanN(loc);
            end
            disp(['the log10(a) in coefficient searching is: ',num2str(loga)]);
            best_loga = loga;
            [~,loc] = min(abs(OCobj.mLymanN+loga));
            %disp('the corresponding coefficients, sd, and 95% CI are:');
            %disp(num2str([st(loc).beta, st(loc).se, st(loc).beta-1.96*st(loc).se, st(loc).beta+1.96*st(loc).se]));

            loglikelyhood = -0.5*dpf;
            [mx,loc] = max(loglikelyhood); % the maximum loglikelyhood
            mx_llhd=mx;
            disp([]);
            disp([]);
            disp(['Max LLHD: ',num2str(mx_llhd,4),' at log10a: ',num2str(-OCobj.mLymanN(loc))]);
            disp([]);
            disp([]);
            [mn_llhd,~] = min(loglikelyhood);
            loglikelyhood68 = repmat(mx-0.5* 1 /df(loc),size(OCobj.mLymanN));
            loglikelyhood95 = repmat(mx-0.5* (1.96*2) /df(loc),size(OCobj.mLymanN));
            hold on;
            fig=plot(-OCobj.mLymanN, loglikelyhood,strMarker,'LineWidth',lw);
            plot(-OCobj.mLymanN, loglikelyhood68,strcat(strMarker(1),'--'),'LineWidth',1);
            %plot(-OCobj.mLymanN, loglikelyhood95,strcat(strMarker(1),'-'),'LineWidth',1);
            if strcmpi(strMarker(2),'s')
                plot(-OCobj.mLymanN(loc), loglikelyhood(loc),strcat(strMarker(1),'o'),'LineWidth',lw*2);
            else
                plot(-OCobj.mLymanN(loc), loglikelyhood(loc),strcat(strMarker(1),'s-'),'LineWidth',lw*2);
            end
            hold off;
            set(gca,'xminortick','on','yminortick','on');
            set(gca,'box','on'); 
            set(gca,'fontsize',16);
            xlabel('log_1_0(a)','fontsize',20);
            ylabel('Log-likelihood / df','fontsize',20);
            
%             % disp the result from atlas
%             if ~isempty(OCobj.mLogisticRegressionMatBin)
%                 st = [OCobj.mLogisticRegressionMatBin];
%                 disp(['the coefficients for the atlas are: ',num2str(st(loc).stats.beta')]);
%             end
        end
        
        function [fig mx_llhd min_aic aic_loga]= fLogisticRegressionAicFig_a_EUD(OCobj,loga,n_param,strMarker,lw)
            if ~exist('strMarker','var')
                strMarker = 'b';
            elseif ~ischar(strMarker)
                warning('classOutcomeAnalysis:EndPointObservedFig_EUD:strMarker','strMarker is not a string, omit');
                strMarker = 'b';
            end
            if ~exist('lw','var')
                lw = 1;
            elseif ~isnumeric(lw)
                warning('classOutcomeAnalysis:EndPointObservedFig_EUD:linewidth','lw is not a number, omit');
                lw = 1;
            end
            st = [OCobj.mLogisticRegressionMat];

            dpf = [st.dev]; % deviations
            st =[st.stats];
            
            df = [st.dfe]; % degree of freedom
            %dpf = dpf./df; % deviations per degree of freedom
            [~,loc] = min(dpf); % the minimum deviation corresponds to the maximum likelihood
            disp(['best log10(a) of Logistic Regression of exact gEUD is: ',num2str(-OCobj.mLymanN(loc))]);
            if ~exist('loga','var') || ischar(loga)
                loga = -OCobj.mLymanN(loc);
            end
            disp(['the log10(a) in coefficient searching is: ',num2str(loga)]);
            [~,loc] = min(abs(OCobj.mLymanN+loga));
            disp('the corresponding coefficients, sd, and 95% CI are:');
            disp(num2str([st(loc).beta, st(loc).se, st(loc).beta-1.96*st(loc).se, st(loc).beta+1.96*st(loc).se]));

            loglikelyhood = -0.5*dpf;
            mx_llhd = max(loglikelyhood./df);
            
            aic = -2*loglikelyhood + 2*(n_param);
            [mn,loc] = min(aic); % the minimum AIC
            mx_aic_llhd = loglikelyhood(loc);
            
            
            min_aic=mn;
            aic_loga = -OCobj.mLymanN(loc);
            %[mn_llhd,~] = min(loglikelyhood);
            
            loglikelyhood68 = repmat(mx_aic_llhd-0.5,size(OCobj.mLymanN));
            aic68 = -2*(loglikelyhood68) + 2*(n_param);
            
            %loglikelyhood95 = repmat(mx-0.5* (1.96*2) /df(loc),size(OCobj.mLymanN));
            hold on;
            fig=plot(-OCobj.mLymanN, aic,strMarker,'LineWidth',lw);
            plot(-OCobj.mLymanN, aic68,strcat(strMarker(1),'--'),'LineWidth',1);
            %plot(-OCobj.mLymanN, loglikelyhood95,strcat(strMarker(1),'-'),'LineWidth',1);
            if strcmpi(strMarker(2),'s')
                plot(-OCobj.mLymanN(loc), aic(loc),strcat(strMarker(1),'o'),'LineWidth',lw*2);
            else
                plot(-OCobj.mLymanN(loc), aic(loc),strcat(strMarker(1),'s-'),'LineWidth',lw*2);
            end
            set(gca,'xminortick','on','yminortick','on');
            set(gca,'box','on'); 
            set(gca,'fontsize',16);
            xlabel('log_1_0(a)','fontsize',20);
            ylabel('AIC','fontsize',20);
            
%             % disp the result from atlas
%             if ~isempty(OCobj.mLogisticRegressionMatBin)
%                 st = [OCobj.mLogisticRegressionMatBin];
%                 disp(['the coefficients for the atlas are: ',num2str(st(loc).stats.beta')]);
%             end
        end
        
    
          function fig= fLogisticRegressionTvalueExactFig_a_EUD(OCobj,strMarker,lw)
            if ~exist('strMarker','var')
                strMarker = 'b';
            elseif ~ischar(strMarker)
                warning('classOutcomeAnalysis:EndPointObservedFig_EUD:strMarker','strMarker is not a string, omit');
                strMarker = 'b';
            end
            if ~exist('lw','var')
                lw = 1;
            elseif ~isnumeric(lw)
                warning('classOutcomeAnalysis:EndPointObservedFig_EUD:linewidth','lw is not a number, omit');
                lw = 1;
            end
            st = [OCobj.mLogisticRegressionMat];
            st =[st.stats];
            tstats = [st.t];
            tstats = tstats(2,:); % the p-value corresponding to gEUD
            fig=semilogy(-OCobj.mLymanN, tstats,strMarker,'LineWidth',lw);
            hold on;grid on;
%            semilogy(-OCobj.mLymanN, repmat(0.05,size(OCobj.mLymanN)), 'r--','LineWidth',2);
            semilogy(-OCobj.mLymanN, repmat(1.96,size(OCobj.mLymanN),1), 'r--','LineWidth',2);
            semilogy(-OCobj.mLymanN, repmat(-1.96,size(OCobj.mLymanN),1), 'r--','LineWidth',2);
            hold off;
            set(gca,'xminortick','on','yminortick','on');
            set(gca,'box','on');
            set(gca,'fontsize',16);
            xlabel('log_1_0(a)','fontsize',20);
            ylabel('t-statistic','fontsize',20);
            %fig = gca;
            
          end
        
        
          
        
        function [fig min_p]= fLogisticRegressionPvalueExactFig_a_EUD(OCobj,strMarker,lw)
            if ~exist('strMarker','var')
                strMarker = 'b';
            elseif ~ischar(strMarker)
                warning('classOutcomeAnalysis:EndPointObservedFig_EUD:strMarker','strMarker is not a string, omit');
                strMarker = 'b';
            end
            if ~exist('lw','var')
                lw = 1;
            elseif ~isnumeric(lw)
                warning('classOutcomeAnalysis:EndPointObservedFig_EUD:linewidth','lw is not a number, omit');
                lw = 1;
            end
            st = [OCobj.mLogisticRegressionMat];
            st =[st.stats];
            pvalue = [st.p];
            pvalue = pvalue(2,:); % the p-value corresponding to gEUD
            [min_p,loc] = min(pvalue); % the location of mininum p-value
            fig=semilogy(-OCobj.mLymanN, pvalue,strMarker,'LineWidth',lw);
            disp(['Minimum p-value of ',num2str(min_p),' for log10a = ',num2str(-OCobj.mLymanN(loc))]);
            hold on;grid on;
            semilogy(-OCobj.mLymanN(loc), pvalue(loc),strMarker,'LineWidth',lw+2);
            %semilogy(-OCobj.mLymanN, repmat(0.05,size(OCobj.mLymanN)), strcat(strMarker(1),'--'),'LineWidth',1);
            semilogy(-OCobj.mLymanN, repmat(0.05,size(OCobj.mLymanN)), 'r--','LineWidth',2);
            hold off;
            set(gca,'xminortick','on','yminortick','on');
            set(gca,'box','on');
            set(gca,'fontsize',16);
            xlabel('log_1_0(a)','fontsize',20);
            ylabel('p-value','fontsize',20);
            %fig = gca;
            
        end
        
           function [fig min_p t_stats]= fLogisticRegressionPvalueExactFig_MM_a_EUD(OCobj,strMarker,lw)
            if ~exist('strMarker','var')
                strMarker = 'b';
            elseif ~ischar(strMarker)
                warning('classOutcomeAnalysis:EndPointObservedFig_EUD:strMarker','strMarker is not a string, omit');
                strMarker = 'b';
            end
            if ~exist('lw','var')
                lw = 1;
            elseif ~isnumeric(lw)
                warning('classOutcomeAnalysis:EndPointObservedFig_EUD:linewidth','lw is not a number, omit');
                lw = 1;
            end
            log_reg_data = OCobj.mLogisticRegressionGridMixtureModel;
            %log_reg_data = OCobj.mLogisticRegressionGrid;
            
            %% grid of b0,b1,n values, 
            % find best model
            llhds = log_reg_data.loglikelihood;
            b1s = log_reg_data.b1;
            b0s = log_reg_data.b0;
            
            pvalue = inf(length(OCobj.mLymanN),1);
            t_stats = inf(length(OCobj.mLymanN),1);
            % for each n
            
            for i=1:length(OCobj.mLymanN)
                 cur_llhds = llhds(:,:,i);
%                 
                 [max_llhd,loc] = max(cur_llhds(:)); 
                 [b0_ind,b1_ind] = ind2sub(size(cur_llhds),loc); 
% 
                %low95 = max_llhd -0.5*(chi2inv(0.95,3));
                low95 = max_llhd -0.5*(chi2inv(0.95,1));
      

                ctr_b0_b1 = contourc(b0s,b1s,cur_llhds,[low95 low95]);
                b1_95cl = ctr_b0_b1(2,2:end);%first entry is ?
                ci_hi = max(b1_95cl);
                ci_low = min(b1_95cl);
                
                
                b1_se = abs(max(b1_95cl)-min(b1_95cl))/1.96;
                b1_tstat = b1s(b1_ind)/b1_se;
                t_stats(i)=b1_tstat;
                b1_pval = 2 * normcdf(-abs(b1_tstat));  
                %b1_pval = 2 * tcdf(-abs(b1_tstat),122);  

                pvalue(i)=b1_pval;

        

                
%                 %fix b0
%                 b1_llhds = cur_llhds(b0_ind,:); 
%                 % should be same as max_llhd    
%                 b1_mx = b1_llhds(b1_ind); 
%                 
%                 %
%                 low95 = b1_mx -0.5*(1.96*2);
%                 
%                 % get CI 
%                 b1_95cl = b1s(b1_llhds>=low95);
%                 % get standard error
%                 b1_se = abs(max(b1_95cl)-min(b1_95cl))/(2*1.96);
%             
%                 b1_tstat = b1s(b1_ind)/b1_se;
%                 b1_pval = 2 * normcdf(-abs(b1_tstat));
%             
%                 pvalue(i)=b1_pval;
            end
%              [max_llhd,loc] = max(llhds(:));
%             [b0_ind,b1_ind,n_ind] = ind2sub(size(llhds),loc);
%             % want to find 95% CI for b1 to get se
%             b1_llhds = llhds(b0_ind,:,n_ind);
%             
%             b1_mx = b1_llhds(b1_ind); % should be same as max_llhd
%             
%             low95 = b1_mx -0.5*(chi2inv(0.95,1));
%             b1s = log_reg_data.b1;
%             
%             % if contains 0, not significant
%             b1_95cl = b1s(b1_llhds>=low95)
%             b1_se = abs(max(b1_95cl)-min(b1_95cl))/(2*1.96);
%             
%             b1_tstat = b1s(b1_ind)/b1_se;
%             b1_pval = 2 * normcdf(-abs(b1_tstat));
            
            %pvalue = log_reg_data.pvals;
            [min_p,loc] = min(pvalue); % the location of mininum p-value
            fig=semilogy(-OCobj.mLymanN, pvalue,strMarker,'LineWidth',lw);hold on;
            disp(['Minimum p-value of ',num2str(min_p),' for log10a = ',num2str(-OCobj.mLymanN(loc))]);
            hold on;grid on;
            semilogy(-OCobj.mLymanN(loc), pvalue(loc),strMarker,'LineWidth',lw+2);
            %semilogy(-OCobj.mLymanN, repmat(0.05,size(OCobj.mLymanN)), strcat(strMarker(1),'--'),'LineWidth',1);
            semilogy(-OCobj.mLymanN, repmat(0.05,size(OCobj.mLymanN)), 'r--','LineWidth',2);
            hold off;
            set(gca,'xminortick','on','yminortick','on');
            set(gca,'box','on');
            set(gca,'fontsize',16);
            xlabel('log_1_0(a)','fontsize',20);
            ylabel('p-value','fontsize',20);
            %fig = gca;
            
        end
       
        
        function fig=fLKBPvalueFig_a_EUD(OCobj,strMarker,lw)
            if ~exist('strMarker','var')
                strMarker = 'b';
            elseif ~ischar(strMarker)
                warning('classOutcomeAnalysis:EndPointObservedFig_EUD:strMarker','strMarker is not a string, omit');
                strMarker = 'b';
            end
            if ~exist('lw','var')
                lw = 1;
            elseif ~isnumeric(lw)
                warning('classOutcomeAnalysis:EndPointObservedFig_EUD:linewidth','lw is not a number, omit');
                lw = 1;
            end
            %st = [OCobj.mLogisticRegressionMat];
            %st =[st.stats];
            %pvalue = [st.p];
            %pvalue = pvalue(2,:); % the p-value corresponding to gEUD
            
                        
            pvalue = [OCobj.mLymanGrid.pvals];
            [~,loc] = min(pvalue); % the location of mininum p-value
            
            disp(['$$$$$$']);
            disp(['Min p-value: ',num2str(pvalue(loc)),...
                ' at log10(a) = ',num2str(-OCobj.mLymanN(loc))]);
            disp(['Mean dose p-value: ',num2str(pvalue(11))]);
            % to find loga10(a) for p-value <0.05
            % interpl1(pvalue, -OCobj.mLymanN,0.05)
            %  interp1(pvalue(9:15)*1e8,-OCobj.mLymanN(9:15),1)
            fig=semilogy(-OCobj.mLymanN, pvalue,strMarker,'LineWidth',lw);
            hold on;
            %semilogy(-OCobj.mLymanN(loc), pvalue(loc),strMarker,'LineWidth',lw+2);
            semilogy(-OCobj.mLymanN, repmat(0.05,size(OCobj.mLymanN)), 'r--','LineWidth',1);
            set(gca,'xminortick','on','yminortick','on');
            set(gca,'box','on');
            set(gca,'fontsize',14);
            xlabel('log_1_0(a)','fontsize',20);
            ylabel('p-value','fontsize',20);
        end 
        function fLKBPvalueBinFig_a_EUD(OCobj,strMarker,lw)
            if ~exist('strMarker','var')
                strMarker = 'b';
            elseif ~ischar(strMarker)
                warning('classOutcomeAnalysis:EndPointObservedFig_EUD:strMarker','strMarker is not a string, omit');
                strMarker = 'b';
            end
            if ~exist('lw','var')
                lw = 1;
            elseif ~isnumeric(lw)
                warning('classOutcomeAnalysis:EndPointObservedFig_EUD:linewidth','lw is not a number, omit');
                lw = 1;
            end
            %st = [OCobj.mLogisticRegressionMat];
            %st =[st.stats];
            %pvalue = [st.p];
            %pvalue = pvalue(2,:); % the p-value corresponding to gEUD
            
                        
            pvalue = [OCobj.mLymanGridBin.pvals];
            [~,loc] = min(pvalue); % the location of mininum p-value
            semilogy(-OCobj.mLymanN, pvalue,strMarker,'LineWidth',lw);
            hold on;
            semilogy(-OCobj.mLymanN(loc), pvalue(loc),strMarker,'LineWidth',lw+2);
            semilogy(-OCobj.mLymanN, repmat(0.05,size(OCobj.mLymanN)), strcat(strMarker(1),'--'),'LineWidth',1);
            set(gca,'xminortick','on','yminortick','on');
            set(gca,'box','on');
            xlabel('log_1_0a'); ylabel('p-value');
        end 
        
        function [loga,pval] = fLogisticRegressionRespondingCurveExactFig_a_EUD(OCobj,loga,strMarker,lw)
            % parameters
            if ~exist('strMarker','var')
                strMarker = 'b';
            elseif ~ischar(strMarker)
                warning('classOutcomeAnalysis:EndPointObservedFig_EUD:strMarker','strMarker is not a string, omit');
                strMarker = 'b';
            end
            if ~exist('lw','var')
                lw = 1;
            elseif ~isnumeric(lw)
                warning('classOutcomeAnalysis:EndPointObservedFig_EUD:linewidth','lw is not a number, omit');
                lw = 1;
            end
            st = [OCobj.mLogisticRegressionMat];
            dpf = [st.dev]; % deviations
            [~,loc] = min(dpf); % the minimum deviation corresponds to the maximum likelihood
            disp(['best log10(a) of Logistic Regression of exact gEUD is: ',num2str(-OCobj.mLymanN(loc))]);
            if ~exist('loga','var') || ischar(loga)
                loga = -OCobj.mLymanN(loc);
            end
            [~,loc] = min(abs(OCobj.mLymanN+loga)); % the n whose corresponding responding function will be ploted
            disp(['the log10(a) in responding curve is: ',num2str(loga)]);

            % responding curve
            st = OCobj.mLogisticRegressionMat(loc); % the fitting result of that n
            euds = [OCobj.mGrp.mEUD]; euds = euds(loc,:); % the gEUDs of that n
            %here
            
            pvalue = st.stats.p;
            pval = pvalue(2); % the p-value corresponding to gEUD
            
            %tmp
            doses = (0:max(euds))'; % doses (gEUD) whose RP probability will be computed
            %doses = (0:90)';
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
            
                
            disp(['***']);
            disp(['gEUD(log10a = ',num2str(loga),...
                ') < ',num2str(interp1(rpb,doses,.2)), 10,...
                ' for RP < 20%']);
            disp(['***']);
            
        end
        function fLogisticRegressionGoodnessOfFitFig(OCobj,loga)
            % compute the p-value curves
            nmn = 5; nmx = 30;
            nGrp = [nmn:nmx, OCobj.mNumInGrp]';
            nGrp = [nmn:nmx]';
            pValExact = zeros(length(nGrp),3);
            pValBin = zeros(length(nGrp),3);
            for n = 1:length(nGrp)
                % assign number of groups
                OCobj.mLogisticRegressionHosmerLemeshow.n = nGrp(n);
                OCobj.mLogisticRegressionHosmerLemeshowBin.n = nGrp(n);
                OCobj.mLogisticRegressionGTest.n = nGrp(n);
                OCobj.mLogisticRegressionGTestBin.n = nGrp(n);
                OCobj.mLogisticRegressionPearson.n = nGrp(n);
                OCobj.mLogisticRegressionPearsonBin.n = nGrp(n);

                % compute the goodness of fit
                OCobj = OCobj.LogitHosmerLemeshowTestAnalysisExact_EUD(loga);
                OCobj = OCobj.LogitHosmerLemeshowTestAnalysisBin_EUD(loga);
                OCobj = OCobj.LogitGTestAnalysisExact_EUD(loga);
                OCobj = OCobj.LogitGTestAnalysisBin_EUD(loga);
                OCobj = OCobj.LogitPearsonTestAnalysisExact_EUD(loga);
                OCobj = OCobj.LogitPearsonTestAnalysisBin_EUD(loga);

                % save
                pValExact(n,1) = OCobj.mLogisticRegressionHosmerLemeshow.p_value;
                pValExact(n,2) = OCobj.mLogisticRegressionGTest.p_value;
                pValExact(n,3) = OCobj.mLogisticRegressionPearson.p_value;
                pValBin(n,1) = OCobj.mLogisticRegressionHosmerLemeshowBin.p_value;
                pValBin(n,2) = OCobj.mLogisticRegressionGTestBin.p_value;
                pValBin(n,3) = OCobj.mLogisticRegressionPearsonBin.p_value;
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

        function [loga,pval,cur_plot] = fLogisticRegressionGridRespondingCurveExactFig_a_EUD(OCobj,loga,strMarker,lw,is_mm)
            % parameters
            if ~exist('strMarker','var')
                strMarker = 'b';
            elseif ~ischar(strMarker)
                warning('classOutcomeAnalysis:EndPointObservedFig_EUD:strMarker','strMarker is not a string, omit');
                strMarker = 'b';
            end
            if ~exist('lw','var')
                lw = 1;
            elseif ~isnumeric(lw)
                warning('classOutcomeAnalysis:EndPointObservedFig_EUD:linewidth','lw is not a number, omit');
                lw = 1;
            end
            
             if ~exist('is_mm','var')
                is_mm=false;
            elseif isequal(is_mm,'MM') || isequal(is_mm,'mm')
                is_mm=true;
            else
               warning('classOutcomeAnalysis:fLogisticRegressionGridExactFig_a_loglielihood_EUD:mixture model','is_mm is not reconized, omit');
               is_mm=false;
            end

            % find the best "a" 
            
            if is_mm,
                log_reg_data = OCobj.mLogisticRegressionGridMixtureModel;
            else
                log_reg_data = OCobj.mLogisticRegressionGrid;
            end
            
            [mx,loc] = max(log_reg_data.loglikelihood(:));
            [~,~,loc] = ind2sub(size(log_reg_data.loglikelihood),loc);
            disp(['the maximum log likelihood and its "log10(a)" in Logit of exact gEUD are: ',num2str([mx, -OCobj.mLymanN(loc)])]);
            if ~exist('loga','var') || ischar(loga)
                loga = -OCobj.mLymanN(loc);
            end
            disp(['the plotted response curve is at log10(a) = ',num2str(loga)]);

            % coefficients for the Lyman model
            [~,n] = min(abs(OCobj.mLymanN+loga)); % the n whose corresponding responding function will be ploted
            ll = log_reg_data.loglikelihood(:,:,n); % log likelihood of log10(a) = loga
            mx = max(ll(:));
            [xx,yy] = find(ll == mx); % the coefficients
            b0 = log_reg_data.b0(xx);
            b1 = log_reg_data.b1(yy);
            disp(['the b0 and b1 for the response function are: ',num2str([b0,b1])]);
            pval = log_reg_data.pvals(n);
            
            % curves
            euds = [OCobj.mGrp.mEUD]; euds = euds(n,:); % the gEUDs of that n
            doses = (0:max(euds))'; % doses (gEUD) whose RP probability will be computed
             
                  
            rpb = exp(b0+b1*doses);
            rpb = rpb./(1+rpb); % logistic probability
            % 95% CI
            low95 = mx-0.5* (1.96*2); % CI level
            cxy = contourc(log_reg_data.b0,log_reg_data.b1,ll',[low95,low95]); % coordinates of CI contours at CI level
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
            cur_plot=plot(doses,rpb,strMarker,'LineWidth',lw*2); % responding function
            plot(doses,CImx,strcat(strMarker,'--'),'LineWidth',lw); % responding function
            plot(doses,CImn,strcat(strMarker,'--'),'LineWidth',lw); % responding function
       
            
        end
        function fig=fSpearmanCorrelationCoefficient_a_EUD(OCobj,strMarker,lw)
            % parameters
            if ~exist('strMarker','var')
                strMarker = 'b-';
            elseif ~ischar(strMarker)
                warning('classOutcomeAnalysis:fSpearmanCorrelationCoefficient_a_EUD:strMarker','strMarker is not a string, omit');
                strMarker = 'b-';
            end
            if ~exist('lw','var')
                lw = 1;
            elseif ~isnumeric(lw)
                warning('classOutcomeAnalysis:fSpearmanCorrelationCoefficient_a_EUD:linewidth','lw is not a number, omit');
                lw = 1;
            end
            
            comps = [OCobj.mGrp.mFlgCensor]';
            comps = ~comps;
            euds = [OCobj.mGrp.mEUD]';
            rho = corr(euds,comps,'type','Spearman');
            %rs = sign(rho).*(rho.^2);
            rs = rho.^2;
            %rs = rho;
                        
            
            [mx,loc] = max(rs);
          %  [~,~,loc] = ind2sub(size(rs),loc);
            disp(['the maximum R^2 and its "log10(a)" : ',num2str([mx, -OCobj.mLymanN(loc)])]);
            
              % Calculation of confidence levels for R
            % See Altman Gardner 1988
            [mx_rho,~] = max(rho);
            z = 0.5*log((1.+mx_rho)/(1-mx_rho));
            z68 = z - sqrt(1/(OCobj.mNumInGrp-3));
            z95 = z - (1.96/sqrt(OCobj.mNumInGrp-3));
            
            rho68 = (exp(2*z68)-1)/(exp(2*z68)+1);
            rs68 = repmat(rho68^2,size(rs));
            rho95 = (exp(2*z95)-1)/(exp(2*z95)+1);
            rs95 = repmat(rho95^2,size(rs));
            
            hold on;
            fig=plot(-OCobj.mLymanN,rs,strMarker,'LineWidth',lw);
            plot(-OCobj.mLymanN, rs68,strcat(strMarker(1),'--'),'LineWidth',1);
            plot(-OCobj.mLymanN, rs95,strcat(strMarker(1),'-'),'LineWidth',1);
            plot(-OCobj.mLymanN(loc),mx,strMarker,'LineWidth',lw+2);
            set(gca,'xminortick','on','yminortick','on');
            set(gca,'box','on');
            set(gca,'fontsize',14);
            xlabel('log_1_0(a)','fontsize',18);
            ylabel('Spearman R^2','fontsize',18);
            
            %ylabel('Spearman R');
            
        end
       function fig=fWilcoxonRankSum_a_EUD(OCobj,strMarker,lw)
            % parameters
            if ~exist('strMarker','var')
                strMarker = 'b-';
            elseif ~ischar(strMarker)
                warning('classOutcomeAnalysis:fSpearmanCorrelationCoefficient_a_EUD:strMarker','strMarker is not a string, omit');
                strMarker = 'b-';
            end
            if ~exist('lw','var')
                lw = 1;
            elseif ~isnumeric(lw)
                warning('classOutcomeAnalysis:fSpearmanCorrelationCoefficient_a_EUD:linewidth','lw is not a number, omit');
                lw = 1;
            end
            
            comps = [OCobj.mGrp.mFlgCensor]';
            comps = ~comps;
            euds = [OCobj.mGrp.mEUD]';
            comp_euds = euds(comps,:)';
            cens_euds = euds(~comps,:)';
            pvals = ones(length(OCobj.mLymanN),1);
            for i=1:length(OCobj.mLymanN);
                pvals(i) = ranksum(cens_euds(i,:),comp_euds(i,:));                        
            end
            
            [mn,loc] = min(pvals);
          %  [~,~,loc] = ind2sub(size(rs),loc);
            disp(['the minimum RankSum pval and its "log10(a)" : ',num2str([mn, -OCobj.mLymanN(loc)])]);
            

            fig=semilogy(-OCobj.mLymanN,pvals,strMarker,'LineWidth',lw);hold on;
            semilogy(-OCobj.mLymanN(loc),mn,strMarker,'LineWidth',lw+2);
            semilogy(-OCobj.mLymanN, repmat(0.05,size(OCobj.mLymanN)), strcat(strMarker(1),'--'),'LineWidth',1);
            
            set(gca,'xminortick','on','yminortick','on');
            set(gca,'box','on');
            set(gca,'fontsize',14);
            xlabel('log_1_0(a)','fontsize',18);
            ylabel('Wilcoxon Rank-Sum','fontsize',18);
            
            %ylabel('Spearman R');
            
        end 
        
        function fig=fLogisticRegressionGridExactFig_a_loglikelhood_EUD(OCobj,strMarker,lw,is_mm)
            % parameters
            if ~exist('strMarker','var')
                strMarker = 'b-';
            elseif ~ischar(strMarker)
                warning('classOutcomeAnalysis:EndPointObservedFig_EUD:strMarker','strMarker is not a string, omit');
                strMarker = 'b-';
            end
            if ~exist('lw','var')
                lw = 1;
            elseif ~isnumeric(lw)
                warning('classOutcomeAnalysis:EndPointObservedFig_EUD:linewidth','lw is not a number, omit');
                lw = 1;
            end
            if ~exist('is_mm','var')
                is_mm=false;
            elseif isequal(is_mm,'MM') || isequal(is_mm,'mm')
                is_mm=true;
            else
               warning('classOutcomeAnalysis:fLogisticRegressionGridExactFig_a_loglielihood_EUD:mixture model','is_mm is not reconized, omit');
               is_mm=false;
            end

            % find the best "a" 
            if is_mm,
                [mx,loc] = max(OCobj.mLogisticRegressionGridMixtureModel.loglikelihood(:));
                [~,~,loc] = ind2sub(size(OCobj.mLogisticRegressionGridMixtureModel.loglikelihood),loc);
                llmx = max(max(OCobj.mLogisticRegressionGridMixtureModel.loglikelihood));
            else
                [mx,loc] = max(OCobj.mLogisticRegressionGrid.loglikelihood(:));
                [~,~,loc] = ind2sub(size(OCobj.mLogisticRegressionGrid.loglikelihood),loc);
                llmx = max(max(OCobj.mLogisticRegressionGrid.loglikelihood));
            end
                        disp(['the maximum log likelihood and its "log10(a)" in Logistic Regression of exact gEUD are: ',num2str([mx, -OCobj.mLymanN(loc)])]);
            %             llmx = zeros(size(OCobj.mLymanN)); %log likelihood's max for each log10(a)
%             for kk = 1:length(llmx)
%                 ll = OCobj.mLogisticRegressionGrid.loglikelihood(:,:,kk);
%                 llmx(kk) = max(ll(:));
%             end
            llmx = llmx(:) / (OCobj.mNumInGrp-2);
            mx = mx / (OCobj.mNumInGrp-2);
            loglikelyhood68 = repmat(mx-0.5* 1 /(OCobj.mNumInGrp-2),size(llmx));
            loglikelyhood95 = repmat(mx-0.5* (1.96*2) /(OCobj.mNumInGrp-2),size(llmx));
            hold on;
            fig=plot(-OCobj.mLymanN,llmx,strMarker,'LineWidth',lw);
            plot(-OCobj.mLymanN, loglikelyhood68,strcat(strMarker(1),'--'),'LineWidth',1);
            plot(-OCobj.mLymanN, loglikelyhood95,strcat(strMarker(1),'-'),'LineWidth',1);
            if length(strMarker)==1
                strMarker = [strMarker,'*'];
            elseif strMarker(2) == '-' || strMarker(2) == '.'
                strMarker(2) = '*';
            end
            plot(-OCobj.mLymanN(loc),mx,strMarker,'LineWidth',lw+2);
       
        end
        function fLogisticRegressionGridExactFig_b0_loglikelihood_EUD(OCobj,strMarker,lw)
            % parameters
            if ~exist('strMarker','var')
                strMarker = 'b';
            elseif ~ischar(strMarker)
                warning('classOutcomeAnalysis:EndPointObservedFig_EUD:strMarker','strMarker is not a string, omit');
                strMarker = 'b';
            end
            if ~exist('lw','var')
                lw = 1;
            elseif ~isnumeric(lw)
                warning('classOutcomeAnalysis:EndPointObservedFig_EUD:linewidth','lw is not a number, omit');
                lw = 1;
            end

            % find the best "b0"
            [mx,loc] = max(OCobj.mLogisticRegressionGrid.loglikelihood(:));
            [loc,~,~] = ind2sub(size(OCobj.mLogisticRegressionGrid.loglikelihood),loc);
            disp(['the maximum log likelihood and its b0 in Logistic Regression of exact gEUD are: ',num2str([mx, OCobj.mLogisticRegressionGrid.b0(loc)])]);

            % plot the maximum log likelihood curve w.r.t. b0
            llmx = zeros(size(OCobj.mLogisticRegressionGrid.b0)); %log likelihood's max for each log10(a)
            for kk = 1:length(llmx)
                ll = OCobj.mLogisticRegressionGrid.loglikelihood(kk,:,:);
                llmx(kk) = max(ll(:));
            end
            llmx = llmx / (OCobj.mNumInGrp-2);
            mx = mx / (OCobj.mNumInGrp-2);
            loglikelyhood68 = repmat(mx-0.5* 1 /(OCobj.mNumInGrp-2),size(llmx));
            loglikelyhood95 = repmat(mx-0.5* (1.96*2) /(OCobj.mNumInGrp-2),size(llmx));
            hold on;
            plot(OCobj.mLogisticRegressionGrid.b0,llmx,strMarker,'LineWidth',lw);
            plot(OCobj.mLogisticRegressionGrid.b0, loglikelyhood68,strcat(strMarker(1),'--'),'LineWidth',1);
            plot(OCobj.mLogisticRegressionGrid.b0, loglikelyhood95,strcat(strMarker(1),'-'),'LineWidth',1);
            if length(strMarker)==1
                strMarker = [strMarker,'*'];
            elseif strMarker(2) == '-' || strMarker(2) == '.'
                strMarker(2) = '*';
            end
            plot(OCobj.mLogisticRegressionGrid.b0(loc),mx,strMarker,'LineWidth',lw+2);
            hold off;
            set(gca,'xminortick','on','yminortick','on');
            set(gca,'box','on');
            xlabel('b0'); ylabel('log likelihood per degree of freedom');
        end
        function fLogisticRegressionGridExactFig_b1_loglikelihood_EUD(OCobj,strMarker,lw)
            % parameters
            if ~exist('strMarker','var')
                strMarker = 'b';
            elseif ~ischar(strMarker)
                warning('classOutcomeAnalysis:EndPointObservedFig_EUD:strMarker','strMarker is not a string, omit');
                strMarker = 'b';
            end
            if ~exist('lw','var')
                lw = 1;
            elseif ~isnumeric(lw)
                warning('classOutcomeAnalysis:EndPointObservedFig_EUD:linewidth','lw is not a number, omit');
                lw = 1;
            end

            % find the best "TD50"
            [mx,loc] = max(OCobj.mLogisticRegressionGrid.loglikelihood(:));
            [~,loc,~] = ind2sub(size(OCobj.mLogisticRegressionGrid.loglikelihood),loc);
            disp(['the maximum log likelihood and its "b1" in Logistic Regression of exact gEUD are: ',num2str([mx, OCobj.mLogisticRegressionGrid.b1(loc)])]);

            % plot the maximum log likelihood curve w.r.t. b1
            llmx = zeros(size(OCobj.mLogisticRegressionGrid.b1)); %log likelihood's max for each log10(a)
            for kk = 1:length(llmx)
                ll = OCobj.mLogisticRegressionGrid.loglikelihood(:,kk,:);
                llmx(kk) = max(ll(:));
            end
            llmx = llmx / (OCobj.mNumInGrp-2);
            mx = mx / (OCobj.mNumInGrp-2);
            loglikelyhood68 = repmat(mx-0.5* 1 /(OCobj.mNumInGrp-2),size(llmx));
            loglikelyhood95 = repmat(mx-0.5* (1.96*2) /(OCobj.mNumInGrp-2),size(llmx));
            hold on;
            plot(OCobj.mLogisticRegressionGrid.b1,llmx,strMarker,'LineWidth',lw);
            plot(OCobj.mLogisticRegressionGrid.b1, loglikelyhood68,strcat(strMarker(1),'--'),'LineWidth',1);
            plot(OCobj.mLogisticRegressionGrid.b1, loglikelyhood95,strcat(strMarker(1),'-'),'LineWidth',1);
            if length(strMarker)==1
                strMarker = [strMarker,'*'];
            elseif strMarker(2) == '-' || strMarker(2) == '.'
                strMarker(2) = '*';
            end
            plot(OCobj.mLogisticRegressionGrid.b1(loc),mx,strMarker,'LineWidth',lw+2);
            hold off;
            set(gca,'xminortick','on','yminortick','on');
            set(gca,'box','on');
            xlabel('b1'); ylabel('log likelihood per degree of freedom');
        end

        function fLogisticRegressionGridExactFig_b0_b1_EUD(OCobj,loga,zdepth)
            if ~exist('zdepth','var')
                zdepth = 10;
            end
            
            % find the best "a"
            [mx,loc] = max(OCobj.mLogisticRegressionGrid.loglikelihood(:));
            [~,~,loc] = ind2sub(size(OCobj.mLogisticRegressionGrid.loglikelihood),loc);
            disp(['the "log10(a)" for the maximum likelihood in Logistic Regression of exact gEUD is: ',num2str(-OCobj.mLymanN(loc))]);
            if ~exist('loga','var') || ischar(loga)
                loga = -OCobj.mLymanN(loc);
            end
            [~,loc] = min(abs(OCobj.mLymanN+loga)); % the n whose corresponding responding function will be ploted
            disp(['the "log10(a)" for the map of likelihood in Logistic Regression is: ',num2str(-OCobj.mLymanN(loc))]);

            
            % b0 and b1 of the best "log10(a)"
            ll = OCobj.mLogisticRegressionGrid.loglikelihood(:,:,loc);
            mx = max(ll(:));
            [dd,mm] = find(ll == mx,1);
            disp(['b0 & b1 are: ',num2str([OCobj.mLogisticRegressionGrid.b0(dd),OCobj.mLogisticRegressionGrid.b1(mm)])]);
            
            % map of the grid of the best "a"
            if 0 % show the mesh surface
                ll( ll < mx-zdepth ) = -inf;
                surf(OCobj.mLogisticRegressionGrid.b1,OCobj.mLogisticRegressionGrid.b0,ll,'EdgeColor','none');
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
                [c,h] = contour(OCobj.mLogisticRegressionGrid.b0,OCobj.mLogisticRegressionGrid.b1,ll',[low99,low95,low68]);
                hold on; plot(OCobj.mLogisticRegressionGrid.b0(dd),OCobj.mLogisticRegressionGrid.b1(mm),'k*'); hold off;
%                 text_handle = clabel(c,h,'LabelSpacing',100000);%'String',{'low 68%';'low 95%';'low 99.7%'},
%                 set(text_handle,'BackgroundColor',[1 1 .6], 'Edgecolor',[.7 .7 .7], 'String',{'low 68%';'low 95%';'low 99.7%'});
                set(gca,'xminortick','on','yminortick','on');
                set(gca,'box','on');
                xlabel('b0'); ylabel('b1 (Gy^-^1)');
            end

%             % the coresponding parameters computed from atlas
%             if ~isempty(OCobj.mLogisticRegressionGridBin)
%                 ll = OCobj.mLogisticRegressionGridBin.loglikelihood(:,:,loc); % log likelihood of log10(a) = loga
%                 mx = max(ll(:)); % the best point for this log10(a)
%                 [dd,mm] = find(ll == mx); % the coefficients
%                 disp(['b0 & b1 for the same log10(a) in the atlas are: ', num2str([OCobj.mLogisticRegressionGridBin.b0(dd),OCobj.mLogisticRegressionGridBin.b1(mm)])]);
%             end
        end
        function fLogisticRegressionGridExactFig_b0_a_EUD(OCobj,zdepth)
            if ~exist('zdepth','var')
                zdepth = 10;
            end
            
            % find the best "m"
            [mx,loc] = max(OCobj.mLogisticRegressionGrid.loglikelihood(:));
            [~,loc,~] = ind2sub(size(OCobj.mLogisticRegressionGrid.loglikelihood),loc);
            disp(['the "b1" for the maximum likelihood in Logistic Regression of exact gEUD is: ',num2str(OCobj.mLogisticRegressionGrid.b1(loc))]);
            
            % TD50 and log10(a) corresponding to the best m
            ll = squeeze(OCobj.mLogisticRegressionGrid.loglikelihood(:,loc,:));
            [dd,aa] = find(ll == mx,1);
            disp(['b0 & log10(a) are: ',num2str([OCobj.mLogisticRegressionGrid.b0(dd),-OCobj.mLymanN(aa)])]);
            
            % map of the grid of the best "m"
            if 0 % show the mesh surface
                ll( ll < mx-zdepth ) = -inf;
                surf(-OCobj.mLymanN,OCobj.mLogisticRegressionGrid.b0,ll,'EdgeColor','none');
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
                contour(-OCobj.mLymanN,OCobj.mLogisticRegressionGrid.b0,ll,[low99,low95,low68]);
                hold on; plot(-OCobj.mLymanN(aa),OCobj.mLogisticRegressionGrid.b0(dd),'k*'); hold off;
                set(gca,'xminortick','on','yminortick','on');
                set(gca,'box','on');
                xlabel('log_1_0(a)'); ylabel('b0');
            end
        end
        function fLogisticRegressionGridExactFig_b1_a_EUD(OCobj,zdepth)
            if ~exist('zdepth','var')
                zdepth = 10;
            end
            
            % find the best "m"
            [mx,loc] = max(OCobj.mLogisticRegressionGrid.loglikelihood(:));
            [loc,~,~] = ind2sub(size(OCobj.mLogisticRegressionGrid.loglikelihood),loc);
            disp(['the b0 for the maximum likelihood in Logistic Regression of exact gEUD is: ',num2str(OCobj.mLogisticRegressionGrid.b0(loc))]);
            
            % TD50 and log10(a) corresponding to the best m
            ll = squeeze(OCobj.mLogisticRegressionGrid.loglikelihood(loc,:,:));
            [mm,aa] = find(ll == mx,1);
            disp(['b1 & log10(a) are: ',num2str([OCobj.mLogisticRegressionGrid.b1(mm),-OCobj.mLymanN(aa)])]);
            
            % map of the grid of the best "m"
            if 0 % show the mesh surface
                ll( ll < mx-zdepth ) = -inf;
                surf(-OCobj.mLymanN,OCobj.mLogisticRegressionGrid.b1,ll,'EdgeColor','none');
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
                contour(-OCobj.mLymanN,OCobj.mLogisticRegressionGrid.b1,ll,[low99,low95,low68]);
                hold on; plot(-OCobj.mLymanN(aa),OCobj.mLogisticRegressionGrid.b1(mm),'k*'); hold off;
                set(gca,'xminortick','on','yminortick','on');
                set(gca,'box','on');
                xlabel('log_1_0(a)'); ylabel('b1 (Gy^-^1)');
            end
        end

        function fLymanGridExactFig_TD50_m_EUD(OCobj,loga,zdepth)
%             if ~exist('strMarker','var')
%                 strMarker = 'b';
%             elseif ~ischar(strMarker)
%                 warning('classOutcomeAnalysis:EndPointObservedFig_EUD:strMarker','strMarker is not a string, omit');
%                 strMarker = 'b';
%             end
%             if ~exist('lw','var')
%                 lw = 1;
%             elseif ~isnumeric(lw)
%                 warning('classOutcomeAnalysis:EndPointObservedFig_EUD:linewidth','lw is not a number, omit');
%                 lw = 1;
%             end

            if ~exist('zdepth','var')
                zdepth = 10;
            end
            
            % find the best "a"
            [mx,loc] = max(OCobj.mLymanGrid.loglikelihood(:));
            [~,~,loc] = ind2sub(size(OCobj.mLymanGrid.loglikelihood),loc);
%             st = [OCobj.mLogisticRegressionMat];
%             dpf = [st.dev]; % deviations
%             st =[st.stats];
%             df = [st.dfe]; % degree of freedom
%             dpf = dpf./df; % deviations per degree of freedom
%             loglikelyhood = -0.5*dpf;
%             [~,loc] = max(loglikelyhood); % the maximum loglikelyhood
            disp(['the "log10(a)" for the maximum likelihood in Lyman model of exact gEUD is: ',num2str(-OCobj.mLymanN(loc))]);
            if ~exist('loga','var') || ischar(loga)
                loga = -OCobj.mLymanN(loc);
            end
            [~,loc] = min(abs(OCobj.mLymanN+loga)); % the n whose corresponding responding function will be ploted
            disp(['the "log10(a)" for the map of likelihood in Lyman model is: ',num2str(-OCobj.mLymanN(loc))]);

            
            % TD50 and gamma of the best "log10(a)"
%             lymangrid = OCobj.mLymanGrid;
            ll = OCobj.mLymanGrid.loglikelihood(:,:,loc);
            mx = max(ll(:));
            [dd,mm] = find(ll == mx,1);
            disp(['TD50 & m are: ',num2str([OCobj.mLymanGrid.TD50(dd),OCobj.mLymanGrid.m(mm)])]);
            
            % map of the grid of the best "a"
            if 0 % show the mesh surface
                ll( ll < mx-zdepth ) = -inf;
                surf(OCobj.mLymanGrid.m,OCobj.mLymanGrid.TD50,ll,'EdgeColor','none');
                ylabel('TD50 (Gy)'); xlabel('m'); zlabel('Log likelihood');
            else % show the contours
                % compute the 68% and 95% CIs
             
                %% Fan's implementation (wrong?)
%                 low68 = mx-0.5* (1);
%                 low95 = mx-0.5* (1.96^2);
%                 low99 = mx-0.5* (3^2);

                %% Wilks
                %L >= L(max) - chi2inv(alpha,n(parameters)/2;
                low68 = mx -0.5*(chi2inv(0.68,2));
                low95 = mx -0.5*(chi2inv(0.95,2));
                low99 = mx -0.5*(chi2inv(0.99,2));

                % display the contours
                [c,h] = contour(OCobj.mLymanGrid.TD50,OCobj.mLymanGrid.m,ll',[low99,low95,low68],'LineWidth',2);
                hold on; plot(OCobj.mLymanGrid.TD50(dd),OCobj.mLymanGrid.m(mm),'k+','LineWidth',2,'MarkerSize',12); hold off;
%                 text_handle = clabel(c,h,'LabelSpacing',100000);%'String',{'low 68%';'low 95%';'low 99.7%'},
%                 set(text_handle,'BackgroundColor',[1 1 .6], 'Edgecolor',[.7 .7 .7], 'String',{'low 68%';'low 95%';'low 99.7%'});
                set(gca,'xminortick','on','yminortick','on');
                set(gca,'box','on');
                set(gca,'FontSize',16)
                xlabel('TD_{50} (Gy)','FontSize',22); ylabel('m','FontSize',22);
            end

%             % the coresponding parameters computed from atlas
%             if ~isempty(OCobj.mLymanGridBin)
%                 ll = OCobj.mLymanGridBin.loglikelihood(:,:,loc); % log likelihood of log10(a) = loga
%                 mx = max(ll(:)); % the best point for this log10(a)
%                 [dd,mm] = find(ll == mx); % the coefficients
%                 disp(['TD50 & m for the same log10(a) in the atlas are: ', num2str([OCobj.mLymanGridBin.TD50(dd),OCobj.mLymanGridBin.m(mm)])]);
%             end
        end
        
        
        function fig=fLymanGridExactFig_IsoSurface(OCobj,vw_arg)
            
            
            if nargin < 2
                cur_view=3;
            else
                cur_view=zeros(3,1)';
                cur_view(vw_arg)=1;
            end
            
            [mx,loc] = max(OCobj.mLymanGrid.loglikelihood(:));
            %[~,~,loc] = ind2sub(size(OCobj.mLymanGrid.loglikelihood),loc);
            %low95 = mx-0.5* (1.96^2);
            low95 = mx-0.5*chi2inv(0.95,3);
            
            %tmp
            ll=OCobj.mLymanGrid.loglikelihood;
            [x,y,z] =meshgrid(1:length(OCobj.mLymanGrid.m),...
                            1:length(OCobj.mLymanGrid.TD50),...
                            1:length(OCobj.mLymanN));
       
            p=patch(isosurface(x,y,z,ll,low95));
            
            %p=patch(isosurface(flipud(OCobj.mLymanGrid.m),OCobj.mLymanGrid
            %.TD50,-OCobj.mLymanN,ll,low95));
            
            %p=patch(isosurface(OCobj.mLymanGrid.m,OCobj.mLymanGrid.TD50,-O
            %Cobj.mLymanN,ll,low95));
            
            %p=patch(isosurface(OCobj.mLymanGrid.loglikelihood,low95));
            isonormals(OCobj.mLymanGrid.loglikelihood,p);
            set(p,'FaceColor','red','EdgeColor','none');
            alpha(p,0.3);
            view(3);
            camlight;
            lighting gouraud;
            %isosurface(OCobj.mLymanGrid.loglikelihood,low95);
             grid on;
             hold on; 
             [~,dd_ind]= min(abs(OCobj.mLymanGrid.TD50-19.2));
            [~,mm_ind] = min(abs(OCobj.mLymanGrid.m-0.38));
            [~,log10a_ind] = min(abs(-OCobj.mLymanN-(-0.16)));
            plot3(mm_ind,dd_ind,log10a_ind,'marker','o',...
                'markerfacecolor','k','color','k',...
                'markersize',10);
 %             plot3(0.38,19.2,-0.16,'marker','o',...
%                'markerfacecolor','k','color','k',...
 %                'markersize',10);
          
            % RE meta analysis results 95% CI
            % a: 0.817 (0.20-1.43)
            % log10(a) = -0.0862  (-0.699 - 0.155)
            % td50: 21.88 (15.92-27.84)
            % m: 0.287 (0.131-0.443)
            
            %% for ellispoid radii
            %% use width from 95% CI pointing toward pooled result
            
             [~,re_dd_ind]= min(abs(OCobj.mLymanGrid.TD50-21.88));
             [~,re_dd_rad_ind]= min(abs(OCobj.mLymanGrid.TD50-15.92));
            
             [~,re_mm_ind] = min(abs(OCobj.mLymanGrid.m-0.287));
             [~,re_mm_rad_ind]= min(abs(OCobj.mLymanGrid.m-0.443));
            
             [~,re_log10a_ind] = min(abs(-OCobj.mLymanN-(-0.0862)));
              [~,re_log10a_rad_ind]= min(abs(-OCobj.mLymanN-(-0.699)));
            plot3(re_mm_ind,re_dd_ind,re_log10a_ind,'marker','o','markerfacecolor','k','color','k','markersize',10);
            %plot3(0.287,21.88,-0.0862,'marker','o','markerfacecolor','k','color','k','markersize',10);
            
            [elx,ely,elz]=ellipsoid(re_mm_ind,re_dd_ind,re_log10a_ind,...
                abs(re_mm_ind-re_mm_rad_ind),...
                abs(re_dd_ind-re_dd_rad_ind),...
                abs(re_log10a_ind-re_log10a_rad_ind),100);
% [elx,ely,elz]=ellipsoid(0.287,21.88,-0.0862,...
%                 abs(0.287-0.443),...
%                 abs(21.88-15.92),...
%                 abs(-0.0862-(-0.699)),100);
            el_surf=surf(elx,ely,elz);
            el_p = patch(surf2patch(el_surf));
            delete(el_surf);
            set(el_p,'FaceColor','green','EdgeColor','none');
            alpha(el_p,0.3);
            %shading faceted;
             camlight;
            lighting gouraud;
            
            min_z = min(get(gca,'ZTick'));
 
            view(cur_view);  
           
            %set(gca,'FontSize',16);
            if isequal(cur_view,3),
                set(gca,'XTickLabel',num2str(OCobj.mLymanGrid.m(get(gca,'XTick')))); 
                xlabel('m (X)','FontSize',20);
                set(gca,'YTickLabel',num2str(OCobj.mLymanGrid.TD50(get(gca,'YTick'))));
                ylabel('TD_{50} [Gy] (Y)','FontSize',20);
                 zlim([min_z size(z,3)]);
                set(gca,'ZTickLabel',num2str(-OCobj.mLymanN(get(gca,'ZTick'))));
                zlabel('log_{10}(a) (Z)','FontSize',20);
                set(gca,'Zdir','reverse');
            elseif isequal(cur_view,[1 0 0])
                set(gca,'YTickLabel',num2str(OCobj.mLymanGrid.TD50(get(gca,'YTick'))));
                ylabel('TD_{50} [Gy] (Y)','FontSize',20);
                 zlim([min_z size(z,3)]);
                set(gca,'ZTickLabel',num2str(-OCobj.mLymanN(get(gca,'ZTick'))));
                zlabel('log_{10}(a) (Z)','FontSize',20);
                set(gca,'Zdir','reverse');
                
            elseif isequal(cur_view,[0 1 0])
                set(gca,'XTickLabel',num2str(OCobj.mLymanGrid.m(get(gca,'XTick')))); 
                xlabel('m (X)','FontSize',20);
                set(gca,'Xdir','reverse');
                 zlim([min_z size(z,3)]);
                set(gca,'ZTickLabel',num2str(-OCobj.mLymanN(get(gca,'ZTick'))));
                zlabel('log_{10}(a) (Z)','FontSize',20);
                set(gca,'Zdir','reverse');
           elseif isequal(cur_view,[0 0 1])
                set(gca,'XTickLabel',num2str(OCobj.mLymanGrid.m(get(gca,'XTick')))); 
                xlabel('m (X)','FontSize',20);
                set(gca,'YTickLabel',num2str(OCobj.mLymanGrid.TD50(get(gca,'YTick'))));
                ylabel('TD_{50} [Gy] (Y)','FontSize',20);
              
            end
               
               
            fig=gcf;
        end
        
        function fLymanGridExactFig_Profile_TD50_a_EUD(OCobj)

            llhds= OCobj.mLymanGrid.loglikelihood;
            %loop over m
            m_size = size(llhds,2);
            td50s = inf(m_size,1);
            as = inf(m_size,1);
            for i=1:m_size
               td50_a_llhds = squeeze(llhds(:,i,:));
               
               %find best llhd
               [~,loc] = max(td50_a_llhds(:));
               [td50_ind,a_ind] = ind2sub(size(td50_a_llhds),loc);
               td50s(i) = OCobj.mLymanGrid.TD50(td50_ind);
               as(i) = -OCobj.mLymanN(a_ind);
               
            end
            
            [ax,h1,h2]=plotyy(OCobj.mLymanGrid.m,as,OCobj.mLymanGrid.m,td50s);
            set(h1,'LineWidth',2);
            set(h2,'LineWidth',2);
            set(ax,'FontSize',16);
            axes(ax(1));ylabel('log_{10}(a)','FontSize',22);
            set(ax(1),'ytick',-1:0.2:1)
            axes(ax(2));ylabel('TD_{50}','FontSize',22);
            xlabel('m','FontSize',22);
            set(ax,'xlim',[0 1.2]);
         end
         
         function fLymanGridExactFig_Profile_td50_m_EUD(OCobj)

            llhds= OCobj.mLymanGrid.loglikelihood;
            %loop over m
            a_size = size(llhds,3);
            td50s = inf(a_size,1);
            ms = inf(a_size,1);
            for i=1:a_size
               td50_m_llhds = squeeze(llhds(:,:,i));
               
               %find best llhd
               [~,loc] = max(td50_m_llhds(:));
               [td50_ind,m_ind] = ind2sub(size(td50_m_llhds),loc);
               td50s(i) = OCobj.mLymanGrid.TD50(td50_ind);
               ms(i) = OCobj.mLymanGrid.m(m_ind);
               
            end
            
            [ax,h1,h2]=plotyy(-OCobj.mLymanN,td50s,-OCobj.mLymanN,ms);
            set(h1,'LineWidth',2);
            set(h2,'LineWidth',2);
            set(ax,'FontSize',16);
            axes(ax(1));ylabel('TD_{50}','FontSize',22);
            axes(ax(2));ylabel('m','FontSize',22);
            xlabel('log_{10}(a)','FontSize',22);

         end
         function fLymanGridExactFig_Profile_m_a_EUD(OCobj)
            
            llhds= OCobj.mLymanGrid.loglikelihood;
            %loop over m
            td50_size = size(llhds,1);
            ms = inf(td50_size,1);
            as = inf(td50_size,1);
          
            for i=1:td50_size
               m_a_llhds = squeeze(llhds(i,:,:));
               
               %find best llhd
               [~,loc] = max(m_a_llhds(:));
               [m_ind,a_ind] = ind2sub(size(m_a_llhds),loc);
               ms(i) = OCobj.mLymanGrid.m(m_ind);
               as(i) = -OCobj.mLymanN(a_ind);
               
            end
            
            [ax,h1,h2]=plotyy(OCobj.mLymanGrid.TD50,as,OCobj.mLymanGrid.TD50,ms);
            set(h1,'LineWidth',2);
            set(h2,'LineWidth',2);
            set(ax,'FontSize',16);
            axes(ax(1));ylabel('log_{10}(a)','FontSize',22);
            set(ax(1),'ytick',-1:0.2:1)
            axes(ax(2));ylabel('m','FontSize',22);
            xlabel('TD_{50}','FontSize',22);
            %set(ax,'xlim',[0 1.2]);
           end
         
        function fLymanGridExactFig_TD50_a_EUD(OCobj,zdepth)
            if ~exist('zdepth','var')
                zdepth = 10;
            end
            
            % find the best "m"
            [mx,loc] = max(OCobj.mLymanGrid.loglikelihood(:));
            [~,loc,~] = ind2sub(size(OCobj.mLymanGrid.loglikelihood),loc);
            disp(['the "m" for the maximum likelihood in Lyman model of exact gEUD is: ',num2str(OCobj.mLymanGrid.m(loc))]);
            
            % TD50 and log10(a) corresponding to the best m
            ll = squeeze(OCobj.mLymanGrid.loglikelihood(:,loc,:));
            [dd,aa] = find(ll == mx,1);
            disp(['TD50 & a are: ',num2str([OCobj.mLymanGrid.TD50(dd),-OCobj.mLymanN(aa)])]);
            
            % map of the grid of the best "m"
            if 0 % show the mesh surface
                ll( ll < mx-zdepth ) = -inf;
                surf(-OCobj.mLymanN,OCobj.mLymanGrid.TD50,ll,'EdgeColor','none');
                ylabel('TD50 (Gy)'); xlabel('log_1_0(a)'); zlabel('Log likelihood');
            else % show the contours
                % compute the 68% and 95% CIs
%                low68 = mx-0.5* (1);
%                low95 = mx-0.5* (1.96^2);
%                low99 = mx-0.5* (3^2);
                       %% Wilks
                %L >= L(max) - chi2inv(alpha,n(parameters)/2;
                low68 = mx -0.5*(chi2inv(0.68,2));
                low95 = mx -0.5*(chi2inv(0.95,2));
                low99 = mx -0.5*(chi2inv(0.99,2));

                % display the contours
                [C,~]=contour(-OCobj.mLymanN,OCobj.mLymanGrid.TD50,ll,[low99,low95,low68],'LineWidth',2);
                hold on; plot(-OCobj.mLymanN(aa),OCobj.mLymanGrid.TD50(dd),'k+','LineWidth',2,'MarkerSize',12); hold off;
                set(gca,'xminortick','on','yminortick','on');
                set(gca,'box','on');
                set(gca,'FontSize',16);
                xlabel('log_1_0(a)','FontSize',22); ylabel('TD50 (Gy)','FontSize',22);
            end
        end
        function fLymanGridExactFig_m_a_EUD(OCobj,zdepth)
            if ~exist('zdepth','var')
                zdepth = 10;
            end
            
            % find the best "m"
            [mx,loc] = max(OCobj.mLymanGrid.loglikelihood(:));
            [loc,~,~] = ind2sub(size(OCobj.mLymanGrid.loglikelihood),loc);
            disp(['the "TD50" for the maximum likelihood in Lyman model of exact gEUD is: ',num2str(OCobj.mLymanGrid.TD50(loc))]);
            
            % TD50 and log10(a) corresponding to the best m
            ll = squeeze(OCobj.mLymanGrid.loglikelihood(loc,:,:));
            [mm,aa] = find(ll == mx,1);
            disp(['m & a are: ',num2str([OCobj.mLymanGrid.m(mm),-OCobj.mLymanN(aa)])]);
            
            % map of the grid of the best "m"
            if 0 % show the mesh surface
                ll( ll < mx-zdepth ) = -inf;
                surf(-OCobj.mLymanN,OCobj.mLymanGrid.m,ll,'EdgeColor','none');
                ylabel('m'); xlabel('log_1_0(a)'); zlabel('Log likelihood');
            else % show the contours
                % compute the 68% and 95% CIs
%                 low68 = mx-0.5* (1);
%                 low95 = mx-0.5* (1.96^2);
%                 low99 = mx-0.5* (3^2);

          %% Wilks
                %L >= L(max) - chi2inv(alpha,n(parameters)/2;
                low68 = mx -0.5*(chi2inv(0.68,2));
                low95 = mx -0.5*(chi2inv(0.95,2));
                low99 = mx -0.5*(chi2inv(0.99,2));
%                 low68 = mx-0.5* (2.30);
%                 low95 = mx-0.5* (6.17);
%                 low99 = mx-0.5* (9.21);
                % display the contours
                contour(-OCobj.mLymanN,OCobj.mLymanGrid.m,ll,[low99,low95,low68],'LineWidth',2);
                hold on; plot(-OCobj.mLymanN(aa),OCobj.mLymanGrid.m(mm),'k+','LineWidth',2,'MarkerSize',12); hold off;
                set(gca,'xminortick','on','yminortick','on');
                set(gca,'box','on');
                set(gca,'FontSize',16)
                xlabel('log_1_0(a)','FontSize',22); ylabel('m','FontSize',22);
            end
        end
        
        function fig=fLymanGridExactFig_a_loglikelihood(OCobj,strMarker,lw)
            % parameters
            if ~exist('strMarker','var')
                strMarker = 'b-';
            elseif ~ischar(strMarker)
                warning('classOutcomeAnalysis:EndPointObservedFig_EUD:strMarker','strMarker is not a string, omit');
                strMarker = 'b-';
            end
            if ~exist('lw','var')
                lw = 1;
            elseif ~isnumeric(lw)
                warning('classOutcomeAnalysis:EndPointObservedFig_EUD:linewidth','lw is not a number, omit');
                lw = 1;
            end

            % find the best "a"
            [mx,loc] = max(OCobj.mLymanGrid.loglikelihood(:));
            [~,~,loc] = ind2sub(size(OCobj.mLymanGrid.loglikelihood),loc);
            disp(['Max llhd and log10(a): ',num2str([mx, -OCobj.mLymanN(loc)])]);
                        
            % plot the maximum log likelihood curve w.r.t. log10(a)
            llmx = zeros(size(OCobj.mLymanN)); %log likelihood's max for each log10(a)
            for kk = 1:length(llmx)
                ll = OCobj.mLymanGrid.loglikelihood(:,:,kk);
                llmx(kk) = max(ll(:));
            end
            llmx = llmx / (OCobj.mNumInGrp-2);
            mx = mx / (OCobj.mNumInGrp-2);
            
            lowlog68 = mx-0.5*1/(OCobj.mNumInGrp-2);
            lowlog95 = mx-0.5*(1.96*2)/(OCobj.mNumInGrp-2);
             
            [~, CI68loga] = ConfidenceInterval(OCobj.mLymanN,flipdim(llmx,1), lowlog68);
            [~, CI95loga] = ConfidenceInterval(OCobj.mLymanN,flipdim(llmx,1), lowlog95);
            
% 
%             disp(['68% log(a) CIs: ',...
%                 num2str(-OCobj.mLymanN(loc)),...
%                 ' [',...
%                 num2str(CI68loga(1)),...
%                 ', ',...
%                 num2str(CI68loga(2)),...
%                 ']']);
%             disp(['95% log(a) CIs: ',...
%                 num2str(-OCobj.mLymanN(loc)),...
%                 ' [',...
%                 num2str(CI95loga(1)),...
%                 ', ',...
%                 num2str(CI95loga(2)),...
%                 ']']);
            
         disp(['68% a CIs: ',...
                num2str(10.^(-OCobj.mLymanN(loc)),3),...
                ' [',...
                num2str(10.^(CI68loga(1)),3),...
                ', ',...
                num2str(10.^(CI68loga(2)),3),...
                ']']);
            disp(['95% a CIs: ',...
                num2str(10.^(-OCobj.mLymanN(loc)),3),...
                ' [',...
                num2str(10.^(CI95loga(1)),3),...
                ', ',...
                num2str(10.^(CI95loga(2)),3),...
                ']']);
            
            loglikelyhood68 = repmat(lowlog68,size(llmx));
            loglikelyhood95 = repmat(lowlog95,size(llmx));
            hold on;
            fig=plot(-OCobj.mLymanN,llmx,strMarker,'LineWidth',lw);
            plot(-OCobj.mLymanN, loglikelyhood68,strcat(strMarker(1),'--'),'LineWidth',1);
            plot(-OCobj.mLymanN, loglikelyhood95,strcat(strMarker(1),'-'),'LineWidth',1);
            if length(strMarker)==1
                strMarker = [strMarker,'*'];
            elseif strMarker(2) == '-' || strMarker(2) == '.'
                strMarker(2) = '*';
            end
            plot(-OCobj.mLymanN(loc),mx,strMarker,'LineWidth',lw+2,'MarkerSize',10);
            ylim([-0.43 -0.33]);
            set(gca,'xminortick','on','yminortick','on');
            set(gca,'box','on');
            %set(gca,'fontsize',14);
            xlabel('log_1_0(a)');
            ylabel('log likelihood per degree of freedom');
        end
        function fig=fLymanGridExactFig_TD50_loglikelihood(OCobj,strMarker,lw)
            % parameters
            if ~exist('strMarker','var')
                strMarker = 'b';
            elseif ~ischar(strMarker)
                warning('classOutcomeAnalysis:EndPointObservedFig_EUD:strMarker','strMarker is not a string, omit');
                strMarker = 'b';
            end
            if ~exist('lw','var')
                lw = 1;
            elseif ~isnumeric(lw)
                warning('classOutcomeAnalysis:EndPointObservedFig_EUD:linewidth','lw is not a number, omit');
                lw = 1;
            end

            % find the best "TD50"
            [mx,loc] = max(OCobj.mLymanGrid.loglikelihood(:));
            [loc,~,~] = ind2sub(size(OCobj.mLymanGrid.loglikelihood),loc);
            %disp(['the maximum log likelihood and its corresponding "TD50" in Lyman model of exact gEUD are: ',num2str([mx, OCobj.mLymanGrid.TD50(loc)])]);
            disp(['Max LLHD and TD50: ',num2str([mx, OCobj.mLymanGrid.TD50(loc)])]);
            % plot the maximum log likelihood curve w.r.t. TD50
            llmx = zeros(size(OCobj.mLymanGrid.TD50)); %log likelihood's max for each log10(a)
            for kk = 1:length(llmx)
                ll = OCobj.mLymanGrid.loglikelihood(kk,:,:);
                llmx(kk) = max(ll(:));
            end
            llmx = llmx / (OCobj.mNumInGrp-2);
            mx = mx / (OCobj.mNumInGrp-2);

            lowlog68 = mx-0.5*1/(OCobj.mNumInGrp-2);
            lowlog95 = mx-0.5*(1.96*2)/(OCobj.mNumInGrp-2);
            
            [~, CI68TD50] = ConfidenceInterval(OCobj.mLymanGrid.TD50,llmx, lowlog68);
            [~, CI95TD50] = ConfidenceInterval(OCobj.mLymanGrid.TD50,llmx, lowlog95);
            
            disp(['68% TD50 CIs: ',...
                num2str(OCobj.mLymanGrid.TD50(loc)),...
                ' [',...
                num2str(CI68TD50(1)),...
                ', ',...
                num2str(CI68TD50(2)),...
                ']']);
            disp(['95% TD50 CIs: ',...
                num2str(OCobj.mLymanGrid.TD50(loc)),...
                ' [',...
                num2str(CI95TD50(1)),...
                ', ',...
                num2str(CI95TD50(2)),...
                ']']);
            
            
            loglikelyhood68 = repmat(lowlog68,size(llmx));
            loglikelyhood95 = repmat(lowlog95,size(llmx));
            hold on;
            fig=plot(OCobj.mLymanGrid.TD50,llmx,strMarker,'LineWidth',lw);
            plot(OCobj.mLymanGrid.TD50, loglikelyhood68,strcat(strMarker(1),'--'),'LineWidth',1);
            plot(OCobj.mLymanGrid.TD50, loglikelyhood95,strcat(strMarker(1),'-'),'LineWidth',1);
            if length(strMarker)==1
                strMarker = [strMarker,'*'];
            elseif strMarker(2) == '-' || strMarker(2) == '.'
                strMarker(2) = '*';
            end
            plot(OCobj.mLymanGrid.TD50(loc),mx,strMarker,'LineWidth',lw+2,'MarkerSize',10);
            set(gca,'xminortick','on','yminortick','on');
            set(gca,'box','on');
            xlabel('TD50'); ylabel('log likelihood per degree of freedom');
        end
        function fig=fLymanGridExactFig_m_loglikelihood(OCobj,strMarker,lw)
            % parameters
            if ~exist('strMarker','var')
                strMarker = 'b';
            elseif ~ischar(strMarker)
                warning('classOutcomeAnalysis:EndPointObservedFig_EUD:strMarker','strMarker is not a string, omit');
                strMarker = 'b';
            end
            if ~exist('lw','var')
                lw = 1;
            elseif ~isnumeric(lw)
                warning('classOutcomeAnalysis:EndPointObservedFig_EUD:linewidth','lw is not a number, omit');
                lw = 1;
            end

            % find the best "TD50"
            [mx,loc] = max(OCobj.mLymanGrid.loglikelihood(:));
            [~,loc,~] = ind2sub(size(OCobj.mLymanGrid.loglikelihood),loc);
            disp(['the maximum log likelihood and its corresponding "m" in Lyman model of exact gEUD are: ',num2str([mx, OCobj.mLymanGrid.m(loc)])]);

            % plot the maximum log likelihood curve w.r.t. m
            llmx = zeros(size(OCobj.mLymanGrid.m)); %log likelihood's max for each log10(a)
            for kk = 1:length(llmx)
                ll = OCobj.mLymanGrid.loglikelihood(:,kk,:);
                llmx(kk) = max(ll(:));
            end
            llmx = llmx / (OCobj.mNumInGrp-2);
            mx = mx / (OCobj.mNumInGrp-2);
            
               
            lowlog68 = mx-0.5*1/(OCobj.mNumInGrp-2);
            lowlog95 = mx-0.5*(1.96*2)/(OCobj.mNumInGrp-2);
            
            [~, CI68M] = ConfidenceInterval(OCobj.mLymanGrid.m,llmx, lowlog68);
            [~, CI95M] = ConfidenceInterval(OCobj.mLymanGrid.m,llmx, lowlog95);
            
            disp(['68% m CIs: ',...
                num2str(OCobj.mLymanGrid.m(loc)),...
                ' [',...
                num2str(CI68M(1)),...
                ', ',...
                num2str(CI68M(2)),...
                ']']);
            disp(['95% m CIs: ',...
                num2str(OCobj.mLymanGrid.m(loc)),...
                ' [',...
                num2str(CI95M(1)),...
                ', ',...
                num2str(CI95M(2)),...
                ']']);
            
            
            loglikelyhood68 = repmat(lowlog68,size(llmx));
            loglikelyhood95 = repmat(lowlog95,size(llmx));
                        
           
            
            hold on;
            fig=plot(OCobj.mLymanGrid.m,llmx,strMarker,'LineWidth',lw);
            plot(OCobj.mLymanGrid.m, loglikelyhood68,strcat(strMarker(1),'--'),'LineWidth',1);
            plot(OCobj.mLymanGrid.m, loglikelyhood95,strcat(strMarker(1),'-'),'LineWidth',1);
            if length(strMarker)==1
                strMarker = [strMarker,'*'];
            elseif strMarker(2) == '-' || strMarker(2) == '.'
                strMarker(2) = '*';
            end
            plot(OCobj.mLymanGrid.m(loc),mx,strMarker,'LineWidth',lw+2,'MarkerSize',10);
            set(gca,'xminortick','on','yminortick','on');
            set(gca,'box','on');
            xlabel('m'); ylabel('log likelihood per degree of freedom');
        end
        
        function fig=fLymanGridExactFig_TD50_loglikelihoodAtLoga_EUD(OCobj,loga,strMarker,lw)
            % parameters
            if ~exist('strMarker','var')
                strMarker = 'b';
            elseif ~ischar(strMarker)
                warning('classOutcomeAnalysis:EndPointObservedFig_EUD:strMarker','strMarker is not a string, omit');
                strMarker = 'b';
            end
            if ~exist('lw','var')
                lw = 1;
            elseif ~isnumeric(lw)
                warning('classOutcomeAnalysis:EndPointObservedFig_EUD:linewidth','lw is not a number, omit');
                lw = 1;
            end

            % find the best "a"
            [~,loc] = max(OCobj.mLymanGrid.loglikelihood(:));
            [~,~,loc] = ind2sub(size(OCobj.mLymanGrid.loglikelihood),loc);
            %disp(['the "log10(a)" for the maximum likelihood in Lyman model of exact gEUD is: ',num2str(-OCobj.mLymanN(loc))]);
            if ~exist('loga','var') || ischar(loga)
                loga = -OCobj.mLymanN(loc);
            end
            [~,loc] = min(abs(OCobj.mLymanN+loga)); % the n whose corresponding responding function will be ploted
            disp(['log10(a): ',num2str(-OCobj.mLymanN(loc))]);

            % plot the maximum log likelihood curve w.r.t. TD50 under a spicific loga
            llmx = max(OCobj.mLymanGrid.loglikelihood(:,:,loc),[],2);
            llmx = llmx / (OCobj.mNumInGrp-2);
            [mx,loc] = max(llmx);
            
            disp(['Best TD50: ',num2str(OCobj.mLymanGrid.TD50(loc))]);
           
            lowlog68 = mx-0.5*1/(OCobj.mNumInGrp-2);
            lowlog95 = mx-0.5*(1.96*2)/(OCobj.mNumInGrp-2);
                     
            [~, CI68M] = ConfidenceInterval(OCobj.mLymanGrid.TD50,llmx, lowlog68);
            [~, CI95M] = ConfidenceInterval(OCobj.mLymanGrid.TD50,llmx, lowlog95);
            
            disp(['68% TD50 CIs: ',...
                num2str(OCobj.mLymanGrid.TD50(loc)),...
                ' [',...
                num2str(CI68M(1)),...
                ', ',...
                num2str(CI68M(2)),...
                ']']);
            disp(['95% TD50 CIs: ',...
                num2str(OCobj.mLymanGrid.TD50(loc)),...
                ' [',...
                num2str(CI95M(1)),...
                ', ',...
                num2str(CI95M(2)),...
                ']']);
            
            
            loglikelyhood68 = repmat(lowlog68,size(llmx));
            loglikelyhood95 = repmat(lowlog95,size(llmx));
            
            hold on;
            fig=plot(OCobj.mLymanGrid.TD50,llmx,strMarker,'LineWidth',lw);
            plot(OCobj.mLymanGrid.TD50, loglikelyhood68,strcat(strMarker(1),'--'),'LineWidth',1);
            plot(OCobj.mLymanGrid.TD50, loglikelyhood95,strcat(strMarker(1),'-'),'LineWidth',1);
            if length(strMarker)==1
                strMarker = [strMarker,'*'];
            elseif strMarker(2) == '-' || strMarker(2) == '.'
                strMarker(2) = '*';
            end
            plot(OCobj.mLymanGrid.TD50(loc),mx,strMarker,'LineWidth',lw+2);
            set(gca,'xminortick','on','yminortick','on');
            set(gca,'box','on');
            xlabel('TD50'); ylabel('log likelihood per degree of freedom');
        end
        function fig=fLymanGridExactFig_m_loglikelihoodAtLoga_EUD(OCobj,loga,strMarker,lw)
            % parameters
            if ~exist('strMarker','var')
                strMarker = 'b';
            elseif ~ischar(strMarker)
                warning('classOutcomeAnalysis:EndPointObservedFig_EUD:strMarker','strMarker is not a string, omit');
                strMarker = 'b';
            end
            if ~exist('lw','var')
                lw = 1;
            elseif ~isnumeric(lw)
                warning('classOutcomeAnalysis:EndPointObservedFig_EUD:linewidth','lw is not a number, omit');
                lw = 1;
            end

            % find the best "a"
            [~,loc] = max(OCobj.mLymanGrid.loglikelihood(:));
            [~,~,loc] = ind2sub(size(OCobj.mLymanGrid.loglikelihood),loc);
            disp(['the "log10(a)" for the maximum likelihood in Lyman model of exact gEUD is: ',num2str(-OCobj.mLymanN(loc))]);
            if ~exist('loga','var') || ischar(loga)
                loga = -OCobj.mLymanN(loc);
            end
            [~,loc] = min(abs(OCobj.mLymanN+loga)); % the n whose corresponding responding function will be ploted
            disp(['the "log10(a)" for the coefficients is: ',num2str(-OCobj.mLymanN(loc))]);

            % plot the maximum log likelihood curve w.r.t. m under a spicific loga
            llmx = max(OCobj.mLymanGrid.loglikelihood(:,:,loc),[],1);
            llmx = llmx / (OCobj.mNumInGrp-2);
            [mx,loc] = max(llmx);
            
             lowlog68 = mx-0.5*1/(OCobj.mNumInGrp-2);
            lowlog95 = mx-0.5*(1.96*2)/(OCobj.mNumInGrp-2);
            
            [~, CI68M] = ConfidenceInterval(OCobj.mLymanGrid.m,llmx, lowlog68);
            [~, CI95M] = ConfidenceInterval(OCobj.mLymanGrid.m,llmx, lowlog95);
            
            disp(['68% m CIs: ',...
                num2str(OCobj.mLymanGrid.m(loc)),...
                ' [',...
                num2str(CI68M(1)),...
                ', ',...
                num2str(CI68M(2)),...
                ']']);
            disp(['95% m CIs: ',...
                num2str(OCobj.mLymanGrid.m(loc)),...
                ' [',...
                num2str(CI95M(1)),...
                ', ',...
                num2str(CI95M(2)),...
                ']']);
            
            
            
            
            disp(['the best m is: ',num2str(OCobj.mLymanGrid.m(loc))]);
            loglikelyhood68 = repmat(mx-0.5* 1 /(OCobj.mNumInGrp-2),size(llmx));
            loglikelyhood95 = repmat(mx-0.5* (1.96*2) /(OCobj.mNumInGrp-2),size(llmx));
            hold on;
            fig=plot(OCobj.mLymanGrid.m,llmx,strMarker,'LineWidth',lw);
            plot(OCobj.mLymanGrid.m, loglikelyhood68,strcat(strMarker(1),'--'),'LineWidth',1);
            plot(OCobj.mLymanGrid.m, loglikelyhood95,strcat(strMarker(1),'-'),'LineWidth',1);
            if length(strMarker)==1
                strMarker = [strMarker,'*'];
            elseif strMarker(2) == '-' || strMarker(2) == '.'
                strMarker(2) = '*';
            end
            plot(OCobj.mLymanGrid.m(loc),mx,strMarker,'LineWidth',lw+2);
            set(gca,'xminortick','on','yminortick','on');
            set(gca,'box','on');
            xlabel('m'); ylabel('log likelihood per degree of freedom');
        end

        function fLymanGridGoodnessOfFitFig(OCobj,loga)
            % compute the p-value curves
            nmn = 5; nmx = 30;
            nGrp = [nmn:nmx, OCobj.mNumInGrp]';
            nGrp = [nmn:nmx]';
            pValExact = zeros(length(nGrp),3);
            pValBin = zeros(length(nGrp),3);
            for n = 1:length(nGrp)
                % assign number of groups
                OCobj.mLymanHosmerLemeshow.n = nGrp(n);
                OCobj.mLymanHosmerLemeshowBin.n = nGrp(n);
                OCobj.mLymanGTest.n = nGrp(n);
                OCobj.mLymanGTestBin.n = nGrp(n);
                OCobj.mLymanPearson.n = nGrp(n);
                OCobj.mLymanPearsonBin.n = nGrp(n);

                % compute the goodness of fit
                OCobj = OCobj.LymanHosmerLemeshowTestAnalysisExact_EUD(loga);
                OCobj = OCobj.LymanHosmerLemeshowTestAnalysisBin_EUD(loga);
                OCobj = OCobj.LymanGTestAnalysisExact_EUD(loga);
                OCobj = OCobj.LymanGTestAnalysisBin_EUD(loga);
                OCobj = OCobj.LymanPearsonTestAnalysisExact_EUD(loga);
                OCobj = OCobj.LymanPearsonTestAnalysisBin_EUD(loga);

                % save
                pValExact(n,1) = OCobj.mLymanHosmerLemeshow.p_value;
                pValExact(n,2) = OCobj.mLymanGTest.p_value;
                pValExact(n,3) = OCobj.mLymanPearson.p_value;
                pValBin(n,1) = OCobj.mLymanHosmerLemeshowBin.p_value;
                pValBin(n,2) = OCobj.mLymanGTestBin.p_value;
                pValBin(n,3) = OCobj.mLymanPearsonBin.p_value;
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

        function loga = fLymanGridResponseExactFig_a_EUD(OCobj,loga,strMarker,lw)
            % parameters
            if ~exist('strMarker','var')
                strMarker = 'b';
            elseif ~ischar(strMarker)
                warning('classOutcomeAnalysis:EndPointObservedFig_EUD:strMarker','strMarker is not a string, omit');
                strMarker = 'b';
            end
            if ~exist('lw','var')
                lw = 1;
            elseif ~isnumeric(lw)
                warning('classOutcomeAnalysis:EndPointObservedFig_EUD:linewidth','lw is not a number, omit');
                lw = 1;
            end
            [mx,loc] = max(OCobj.mLymanGrid.loglikelihood(:));
            [~,~,loc] = ind2sub(size(OCobj.mLymanGrid.loglikelihood),loc);
            disp(['the maximum log likelihood and its "log10(a)" in Lyman model of exact gEUD are: ',num2str([mx, -OCobj.mLymanN(loc)])]);
            if ~exist('loga','var') || ischar(loga)
                loga = -OCobj.mLymanN(loc);
            end
            disp(['the plotted response curve is at log10(a) = ',num2str(loga)]);

            % coefficients for the Lyman model
            [~,n] = min(abs(OCobj.mLymanN+loga)); % the n whose corresponding responding function will be ploted
            ll = OCobj.mLymanGrid.loglikelihood(:,:,n); % log likelihood of log10(a) = loga
            mx = max(ll(:));
            disp(['LLHD(log(a) = ',num2str(loga), ') = ',num2str(mx)]);
               
            [xx,yy] = find(ll == mx); % the coefficients
            TD50 = OCobj.mLymanGrid.TD50(xx);
            m = OCobj.mLymanGrid.m(yy);
            disp(['the TD50 and m for the response function are: ',num2str([TD50,m])]);

            % curves
            euds = [OCobj.mGrp.mEUD]; euds = euds(n,:); % the gEUDs of that n
            doses = (0:max(euds))'; % doses (gEUD) whose RP probability will be computed
            rpb = normcdf((doses-TD50)/(m*TD50),0,1); % Lyman probability
            % 95% CI
            low95 = mx-0.5* (1.96*2); % CI level
            cxy = contourc(OCobj.mLymanGrid.TD50,OCobj.mLymanGrid.m,ll',[low95,low95]); % coordinates of CI contours at CI level
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
            disp(['gEUD(a = ',num2str(10^(loga)),...
                '< ',num2str(interp1(rpb,doses,.2)), 10,...
                ' for RP < 20%']);
                
            plot(doses,rpb,strMarker,'LineWidth',lw*2); % responding function
            plot(doses,CImx,strMarker,'LineWidth',lw); % responding function
            plot(doses,CImn,strMarker,'LineWidth',lw); % responding function
            set(gca,'xminortick','on','yminortick','on');
            set(gca,'box','on');
            xlabel('gEUD'); ylabel('RP probability');
        end
        
        
        function loga=fLymanResponseFig_m_EUD(OCobj,m)
            % parameters
            
            % get best fit TD50 and a for each m
              llhds= OCobj.mLymanGrid.loglikelihood;
            %loop over m
            %m_size = size(llhds,2);
            
                        
                % find index for current m
            [~,m_ind] = min(abs(OCobj.mLymanGrid.m - m));
               
                %find best llhd, given m 
            td50_a_llhds = squeeze(llhds(:,m_ind,:));
            [mx,loc] = max(td50_a_llhds(:));
            %mx = max(ll(:));% used to find CIs
            
            [td50_ind,a_ind] = ind2sub(size(td50_a_llhds),loc);
               
             TD50 = OCobj.mLymanGrid.TD50(td50_ind);
             loga = -OCobj.mLymanN(a_ind);
             str = ['m = ',num2str(m),10,...
                 'TD50 = ',num2str(TD50),10,...
                 'log_{10}(a) = ',num2str(loga)];
             
            % coefficients for the Lyman model
             [~,n] = min(abs(OCobj.mLymanN+loga)); % the n whose corresponding responding function will be ploted
            % curves
            euds = [OCobj.mGrp.mEUD]; euds = euds(n,:); % the gEUDs of that n

            doses = (0:max(euds))'; % doses (gEUD) whose RP probability will be computed
            rpb = normcdf((doses-TD50)/(m*TD50),0,1); % Lyman probability

            %             % 95% CI
%             low68 = mx-0.5* (1);
%             low95 = mx-0.5*(1.96*2); % CI level
%             cxy = contourc(OCobj.mLymanGrid.TD50,-OCobj.mLymanN,td50_a_llhds',[low68,low68]); % coordinates of CI contours at CI level
%             CImx = zeros(length(doses),1); % upper end of CI for response at each gEUD
%             CImn = zeros(length(doses),1); % lower end of CI for response at each gEUD
%             TD50 = cxy(1,2:end);
%             m = cxy(2,2:end); 
%             for d = 1:length(doses)
%                 CI = normcdf((doses(d)-TD50)./(m.*TD50),0,1); % probability at each iso point
%                 CImx(d) = max(CI);
%                 CImn(d) = min(CI);
          %  end
            % plot
            
            hold on;
                
            plot(doses,rpb,'-','LineWidth',2); % responding function
%             plot(doses,CImx,'-.','LineWidth',1); % responding function
%             plot(doses,CImn,'-.','LineWidth',1); % responding function
            
            text(0.15,0.85,str,'FontSize',20,'Units','normalized');
            set(gca,'xminortick','on','yminortick','on');
            set(gca,'box','on');
            xlabel('gEUD'); ylabel('RP probability');
        end
        function [medianeud,betainv84,betainv16] = fComplicationObservedFig_EUD(OCobj,loga,numIntervals,strMarker,lw)
            % parameters
            if ~exist('strMarker','var')
                strMarker = 'b';
            elseif ~ischar(strMarker)
                warning('classOutcomeAnalysis:EndPointObservedFig_EUD:strMarker','strMarker is not a string, omit');
                strMarker = 'b';
            end
            if ~exist('lw','var')
                lw = 1;
            elseif ~isnumeric(lw)
                warning('classOutcomeAnalysis:EndPointObservedFig_EUD:linewidth','lw is not a number, omit');
                lw = 1;
            end
            if ~exist('loga','var')
                st = [OCobj.mLogisticRegressionMatBin];
                dpf = [st.dev]; % deviations
                [~,loc] = min(dpf); % the minimum deviation corresponds to the maximum likelihood
                loga = -OCobj.mLymanN(loc);
            end
            if ~exist('numIntervals','var')
                numIntervals = 4;
            end
            
            % parameters
            if ~exist('loga','var') || ischar(loga)
                st = [OCobj.mLogisticRegressionMat];
                dpf = [st.dev]; % deviations
                [~,loc] = min(dpf); % the minimum deviation corresponds to the maximum likelihood
                loga = -OCobj.mLymanN(loc);
            end
            [~,n] = min(abs(OCobj.mLymanN + loga)); % the n whose gEUDs will be ploted
            euds = [OCobj.mGrp.mEUD]; euds = euds(n,:); % the gEUDs of that n
            if ~exist('numIntervals','var')
                numIntervals = 4;
            end

            % grouping patients
            flg=[OCobj.mGrp.mFlgCensor]; % censor flags of patients
            [medianeud,~,~,binlow,binhigh,numcomp,numtotal,betainv84,betainv16] = EventObserved(flg,euds,numIntervals);
            prob = numcomp./numtotal;
            % plot
            errorbar(medianeud,prob,max(0,prob-betainv16),max(0,betainv84-prob),strMarker,'LineWidth',lw);
            errorbar_x(medianeud,prob,(medianeud-binlow),(binhigh-medianeud),strMarker);
            
            %tmp
            ylim([0 0.7]);
            xlim([0 max(binhigh)*1.1]);
            
            %xlims = xlim;
            %xmax = xlims(2);
            %xlim([0 xmax]);
            %ylim([0 0.7]);
            
            set(gca,'xminortick','on','yminortick','on');
            set(gca,'box','on');
            set(gca,'FontSize',16);
            xlabel('gEUD [Gy]','FontSize',22);
            ylabel('Complication rate observed','FontSize',22);
        end
        
         function fProbabilityFig_DVH(OCobj,cur_title)
            % map of the probability of have complication rate of at least value (e.g. 20%)
            if isempty(OCobj.mBetaCumulativeMat)
                disp('Empty member (BetaCumulativeMat), can not display its figure.'); return;
            end

            % prepare
            cm = colormap(jet(300)); cm=cm(1:256,:); %cm(end,:) = 0.5;
            imgmsk = OCobj.mAtlasTotal_DVH > 0; imgmsk=imgmsk';
            eudmx = OCobj.mBinsDose(end); eudmn = OCobj.mBinsDose(1);

            % image data 
            img=1-OCobj.mBetaCumulativeMat; img=img';
            mx = ceil(max(img(imgmsk))*256);
            img(~imgmsk) = NaN;
            colormap(cm(1:mx,:));
            contourf(img); axis xy;
            h_cb=colorbar;
            %set(h_cb,'ylim',[0 1]);
            %caxis([0 1]);

            %             set(gca,'Position',[0.05,0.05,0.8,0.9]);
            
            set(gca,'YTick',[1:2:21]); set(gca,'YTickLabel',0:10:100);
            %xlim = get(gca,'XLim'); % x limit
            %ylim = get(gca,'YLim'); % x limit
            %xticklabel = 0:5:eudmx; xtick = diff(xlim)/(eudmx-eudmn)*xticklabel+xlim(1);
            %set(gca,'XTick',xtick);
            %set(gca,'XTickLabel',xticklabel);
             
            
            set(gca,'xminortick','on','yminortick','on');
            set(gca,'box','on');
            set(gca,'fontsize',18);
            xlabel('BED Dose [Gy_3]','fontsize',24);
            ylabel('Volume [%]','fontsize',24);
            %if isempty(cur_title),
             %             title('Probability that RPS rate \geq 20%','FontSize',18);              
            %else
            %title([cur_title,10,'Probability that RPS rate \geq 20%'],'FontSize',18);              
            %end           
        end
        
        
        
        function fProbabilityFig_EUD(OCobj,cur_title)
            % map of the probability of have complication rate of at least value (e.g. 20%)
            if isempty(OCobj.mBetaCumulativeMat)
                disp('Empty member (BetaCumulativeMat), can not display its figure.'); return;
            end

            % prepare
            cm = colormap(jet(300)); cm=cm(1:256,:); %cm(end,:) = 0.5;
            imgmsk = OCobj.mAtlasTotal > 0; imgmsk=imgmsk';
            eudmx = OCobj.mBinsDose(end); eudmn = OCobj.mBinsDose(1);

            % image data 
            img=1-OCobj.mBetaCumulativeMat; img=img';
            mx = ceil(max(img(imgmsk))*256);
            img(~imgmsk) = NaN;
            colormap(cm(1:mx,:));
            contourf(img); axis xy;
            colorbar;
       
%             set(gca,'Position',[0.05,0.05,0.8,0.9]);
            
            set(gca,'YTick',1:length(OCobj.mLymanN)); set(gca,'YTickLabel',-OCobj.mLymanN);
            xlim = get(gca,'XLim'); % x limit
            %ylim = get(gca,'YLim'); % x limit
            xticklabel = 0:5:eudmx; xtick = diff(xlim)/(eudmx-eudmn)*xticklabel+xlim(1);
            set(gca,'XTick',xtick);
            set(gca,'XTickLabel',xticklabel);
             
            
            set(gca,'xminortick','on','yminortick','on');
            set(gca,'box','on');
            set(gca,'fontsize',14);
            ylabel('log_1_0(a)','fontsize',18);
            xlabel('gEUD [Gy]','fontsize',18);
            if isempty(cur_title),
                          title('Probability that RPS rate \geq 20%','FontSize',18);              
            else
            title([cur_title,10,'Probability that RPS rate \geq 20%'],'FontSize',18);              
            end            
%         xlim=get(gca,'XLim');
%         xtick=linspace(xlim(1),xlim(2),xticknum);
%         set(gca,'XTick',xtick);
%         set(gca,'XTickLabel',doseticks);
%         xlabel('EUD dose (Gy)');
% 
%         set(gca,'YTick',1:length(CGobjs(k).mLymanN)); set(gca,'YTickLabel',CGobjs(k).mLymanN);

%             img=1-OCobj.BetaCumulativeMat'; img(1,end)=1; img(2,end)=0;
%             contourf(img); %contourf(flipud(rot90(OCobj.BetaCumulativeMat,1)));
%             xsteps=0:OCobj.GyStep:max(OCobj.mEUD(:)); xtickvec=10:10:xsteps(end); xtickstep=(size(OCobj.BetaCumulativeMat,1)/xsteps(end)*10);
%             set(gca,'XTick',xtickstep:xtickstep:size(OCobj.BetaCumulativeMat,1)); set(gca,'XTickLabel',xtickvec);
%             ytickvec=OCobj.log10n; ytickstep=size(OCobj.BetaCumulativeMat,2)/length(ytickvec);
%             set(gca,'YTick',ytickstep:ytickstep:size(OCobj.BetaCumulativeMat,2)); set(gca,'YTickLabel',ytickvec);
%             %             pmin=min(OCobj.BetaCumulativeMat(:)); pmax=max(OCobj.BetaCumulativeMat(:)); %bstep=(pmax-pmin)/10*256;
%             colorbar; %colorbar('YTickLabel',round([pmin:(pmax-pmin)/9:pmax]*10)/10);
%             title([OCobj.xlsSheet,', the probability that observed complications arise from true rate > 20%']);
        end
        function fLow68pctConfidenceFig_EUD(OCobj)
            % map of lower 68% confidence limit on complication probability
            if isempty(OCobj.mBetaCumulativeMat)
                disp('Empty member (BetaCumulativeMat), can not display its figure.'); return;
            end

            % prepare
            cm = colormap(jet(300)); cm=cm(1:256,:); %cm(end,:) = 0.5;
            imgmsk = OCobj.mAtlasTotal > 0; imgmsk=imgmsk';
            eudmx = OCobj.mBinsDose(end); eudmn = OCobj.mBinsDose(1);

            img=OCobj.mBetaInverseMat';
            mx = ceil(max(img(imgmsk))*256);
            img(~imgmsk) = NaN;
            colormap(cm(1:mx,:));
            colormap(cm);
            contourf(img); axis xy;
            colorbar;

            caxis([0 0.5]);        
%             set(gca,'Position',[0.05,0.05,0.8,0.9]);
           
            set(gca,'YTick',1:length(OCobj.mLymanN)); set(gca,'YTickLabel',-OCobj.mLymanN);
            xlim = get(gca,'XLim'); % x limit
            xticklabel = 0:5:eudmx; xtick = diff(xlim)/(eudmx-eudmn)*xticklabel+xlim(1);
            set(gca,'XTick',xtick);
            set(gca,'XTickLabel',xticklabel);
            
            set(gca,'xminortick','on','yminortick','on');
            set(gca,'box','on');
            set(gca,'fontsize',14);
            ylabel('log_1_0(a)','fontsize',18);
            xlabel('gEUD [Gy]','fontsize',18);
            title('Low 68% CL on RP rate','FontSize',18);            
            
%             if isempty(OCobj.BetaInverseMat)
%                 disp('Empty member (BetaInverseMat), can not display its figure.'); return;
%             end
%             img=OCobj.BetaInverseMat'; img(1,end)=1; img(2,end)=0;
%             contourf(img); %contourf(flipud(rot90(OCobj.BetaInverseMat,1)));
%             xsteps=0:OCobj.GyStep:max(OCobj.mEUD(:)); xtickvec=10:10:xsteps(end); xtickstep=(size(OCobj.BetaInverseMat,1)/xsteps(end)*10);
%             set(gca,'XTick',xtickstep:xtickstep:size(OCobj.BetaInverseMat,1)); set(gca,'XTickLabel',xtickvec);
%             ytickvec=OCobj.log10n; ytickstep=size(OCobj.BetaInverseMat,2)/length(ytickvec);
%             set(gca,'YTick',ytickstep:ytickstep:size(OCobj.BetaInverseMat,2)); set(gca,'YTickLabel',ytickvec);
%             colorbar;
%             title([OCobj.xlsSheet,', lower 68% confidence limit on complication probability']);
        end
    function fLow68pctConfidenceFig_DVH(OCobj)
            % map of lower 68% confidence limit on complication probability
            if isempty(OCobj.mBetaCumulativeMat)
                disp('Empty member (BetaCumulativeMat), can not display its figure.'); return;
            end

            % prepare
            cm = colormap(jet(300)); cm=cm(1:256,:); %cm(end,:) = 0.5;
            imgmsk = OCobj.mAtlasTotal_DVH > 0; imgmsk=imgmsk';
            eudmx = OCobj.mBinsDose(end); eudmn = OCobj.mBinsDose(1);

            img=OCobj.mBetaInverseMat';
            mx = ceil(max(img(imgmsk))*256);
            img(~imgmsk) = NaN;
            colormap(cm(1:mx,:));
            colormap(cm);
            contourf(img); axis xy;
            colorbar;

            %caxis([0 0.5]);        
%             set(gca,'Position',[0.05,0.05,0.8,0.9]);
            set(gca,'YTick',[1:2:21]); set(gca,'YTickLabel',0:10:100);
            
%             set(gca,'YTick',1:length(OCobj.mLymanN)); set(gca,'YTickLabel',-OCobj.mLymanN);
%             xlim = get(gca,'XLim'); % x limit
%             xticklabel = 0:5:eudmx; xtick = diff(xlim)/(eudmx-eudmn)*xticklabel+xlim(1);
%             set(gca,'XTick',xtick);
%             set(gca,'XTickLabel',xticklabel);
%             
            set(gca,'xminortick','on','yminortick','on');
            set(gca,'box','on');
            set(gca,'fontsize',18);
            xlabel('BED Dose [Gy_3]','fontsize',24);
            ylabel('Volume [%]','fontsize',24);
            %title('Low 68% CL on RP rate','FontSize',18);            
            

        end
    end
    
    methods(Static)
        function OCobj = fCombineAtlas_EUD(CGobj1,CGobj2)
            % integrity check
            % mLymanN
            if any(abs(CGobj1.mLymanN-CGobj2.mLymanN)>1e-6)
                error('The log10(n) is not equal in the two CGobjs so they can not be combined');
            end
            % dose bins in atlas computation
            if max(CGobj1.mBinsDose) >= max(CGobj2.mBinsDose)
                dosebins1 = CGobj1.mBinsDose;
                dosebins2 = CGobj2.mBinsDose;
            else
                dosebins1 = CGobj2.mBinsDose;
                dosebins2 = CGobj1.mBinsDose;
            end
            d = abs( dosebins1(1:size(dosebins2,1)) - dosebins2 );
            if any(d>1e-6)
                error('The dose bins in atlas computation were not consistant so they are not combined');
            end
            
            % combination
            OCobj = CGobj1;
            OCobj = OCobj.AddPatient(CGobj2.mGrp);
            
            OCobj.mBinsDose = dosebins1;
            size1 = size(CGobj1.mAtlasTotal); size2 = size(CGobj2.mAtlasTotal);
            OCobj.mAtlasTotal = zeros(max(size1,size2));
            OCobj.mAtlasTotal(1:size1(1),1:size1(2)) = CGobj1.mAtlasTotal;
            OCobj.mAtlasTotal(1:size2(1),1:size2(2)) = OCobj.mAtlasTotal(1:size2(1),1:size2(2))+CGobj2.mAtlasTotal;
            OCobj.mAtlasComp = zeros(max(size1,size2));
            OCobj.mAtlasComp(1:size1(1),1:size1(2)) = CGobj1.mAtlasComp;
            OCobj.mAtlasComp(1:size2(1),1:size2(2)) = OCobj.mAtlasComp(1:size2(1),1:size2(2))+CGobj2.mAtlasComp;
            
%             % compute probabilities
%             OCobj=OCobj.BetaCumulativeProbability_EUD();
%             OCobj=OCobj.BetaInverseProbability_EUD();
%             OCobj=OCobj.LogitAnalysis_EUD();
        end
        function [Q,Qp, I2] = fLogisticRegressionDataConsistency_EUD(OCobjs) % check data set consistency for combination
            % Q - the Q-value of heterogeneuity between parameters of data
            % sets, b0, b1,..., a
            % Qp - the corresponding p-value
            % I2 - I^2 the inconsistency between parameters of data sets,
            % b0, b1,..., a
            
            % prepare
            num = length(OCobjs(:));
            Q = zeros(num,1); Qp = zeros(num,1); I2 = zeros(num,1);
            b0 = zeros(num,1); b1 = zeros(num,1); a = zeros(num,1);
            b0sd = zeros(num,1); b1sd = zeros(num,1); asd = zeros(num,1);
            % best parameters and standard deviations
            for k = 1:num
                % best parameters
                mx = max(OCobjs(k).mLogisticRegressionGrid.loglikelihood(:));
%                 mx = mx/(OCobjs(k).mNumInGrp-2);
%                 [x,y,z] = ind2sub(size(OCobjs(k).mLogisticRegressionGrid.loglikelihood),loc); % locations of the parameters of the best fit
%                 b0(k) = OCobjs(k).mLogisticRegressionGridBetaRange{1}(x);
%                 b1(k) = OCobjs(k).mLogisticRegressionGridBetaRange{2}(y);
%                 a(k) = -OCobjs(k).mLymanN(z);
                % extract the 68% confidence intervals
                likely1 = max(max(OCobjs(k).mLogisticRegressionGrid.loglikelihood,[],3),[],2);
%                 likely1 = likely1/(OCobjs(k).mNumInGrp-2);
                likely2 = max(max(OCobjs(k).mLogisticRegressionGrid.loglikelihood,[],3),[],1);
%                 likely2 = likely2/(OCobjs(k).mNumInGrp-2);
                likely3 = max(max(OCobjs(k).mLogisticRegressionGrid.loglikelihood,[],2),[],1);
%                 likely3 = likely3/(OCobjs(k).mNumInGrp-2);
%                 CI68 = mx - 0.5*1/(OCobjs(k).mNumInGrp-2); % 68% confidence interval formula, *1 because it is 1 sigma away from the peak
                CI68 = mx - 0.5*1; % 68% confidence interval formula, *1 because 68% is 1 sigma away from the peak
                % 68% point of b0
                [b0(k), CI] = ConfidenceInterval(OCobjs(k).mLogisticRegressionGrid.b0,likely1, CI68);
                b0sd(k) = min(CI(2)-b0(k), b0(k)-CI(1));
                % 68% point of b1
                [b1(k), CI] = ConfidenceInterval(OCobjs(k).mLogisticRegressionGrid.b1,likely2, CI68);
                b1sd(k) = min(CI(2)-b1(k), b1(k)-CI(1));
                % 68% point of a
                [a(k), CI] = ConfidenceInterval(OCobjs(k).mLymanN,likely3, CI68);
                asd(k) = min(CI(2)-a(k), a(k)-CI(1));
            end
            % the heterogenuneity between data sets
            a = -a; % transfer n to a
            [Q(1),Qp(1), I2(1)] = DataSetInconsistency(b0, b0sd); % b0
            [Q(2),Qp(2), I2(2)] = DataSetInconsistency(b1, b1sd); % b1
            [Q(3),Qp(3), I2(3)] = DataSetInconsistency(a, asd); % a
        end
        function [Q,Qp, I2,I2up,I2down] = fLKBDataConsistency_EUD(OCobjs) % check data set consistency for combination
            % Q - the Q-value of heterogeneuity between parameters of data
            % sets,TD50,m, a
            % Qp - the corresponding p-value
            % I2 - I^2 the inconsistency between parameters of data sets,
            % b0, b1,..., a
            
            % prepare
            num = length(OCobjs(:));
            Q = zeros(3,1); Qp = zeros(3,1); 
            I2 = zeros(3,1); I2up = zeros(3,1); I2down = zeros(3,1); 
            b0 = zeros(num,1); b1 = zeros(num,1); a = zeros(num,1);
            b0sd = zeros(num,1); b1sd = zeros(num,1); asd = zeros(num,1);
            % best parameters and standard deviations
            % 68% CI found for each parameter, SD taken from those.
            % TODO, use CIs to get I2 range?
            for k = 1:num
                % best parameters
                mx = max(OCobjs(k).mLymanGrid.loglikelihood(:));
                
                %OCobj.mLymanGrid.TD50
                
                % HERE
                likely1 = max(max(OCobjs(k).mLymanGrid.loglikelihood,[],3),[],2);
                likely2 = max(max(OCobjs(k).mLymanGrid.loglikelihood,[],3),[],1);
                likely3 = max(max(OCobjs(k).mLymanGrid.loglikelihood,[],2),[],1);

                %CI68 = mx - 0.5*1; % 68% confidence interval formula, *1
                %because 68% is 1 sigma away from the peak
                %tmp
                CI68 = mx - 0.5*1; % 68% confidence interval formula, *1 because 68% is 1 sigma away from the peak
                % 68% point of b0
                [b0(k), CI] = ConfidenceInterval(OCobjs(k).mLymanGrid.TD50,likely1, CI68);
                b0sd(k) = min(CI(2)-b0(k), b0(k)-CI(1));
                % 68% point of b1
                [b1(k), CI] = ConfidenceInterval(OCobjs(k).mLymanGrid.m,likely2, CI68);
                b1sd(k) = min(CI(2)-b1(k), b1(k)-CI(1));
                % 68% point of a
                [a(k), CI] = ConfidenceInterval(OCobjs(k).mLymanN,likely3, CI68);
                asd(k) = min(CI(2)-a(k), a(k)-CI(1));
            end
            % the heterogenuneity between data sets
            a = -a; % transfer n to a
            [Q(1),Qp(1), I2(1),I2up(1),I2down(1)] = DataSetInconsistency(b0, b0sd); % b0
            [Q(2),Qp(2), I2(2),I2up(2),I2down(2)] = DataSetInconsistency(b1, b1sd); % b1
            [Q(3),Qp(3), I2(3),I2up(3),I2down(3)] = DataSetInconsistency(a, asd); % a
        end
    end
    
    methods % unused
        
        function OCobj = fLogisticRegressionHosmerLemeshowTestExact_EUD(OCobj,loga)
            % loga determination
            st = [OCobj.mLogisticRegressionMat];
            dpf = [st.dev]; % deviations
            [~,loc] = min(dpf); % the minimum deviation corresponds to the maximum likelihood
            disp(['The best loga of Logistic Regression is: ',num2str(-OCobj.mLymanN(loc))]);
            if ~exist('loga','var') || ischar(loga)
                loga = -OCobj.mLymanN(loc);
            end
            [~,loc] = min(abs(OCobj.mLymanN - (-loga))); % find the loga whose goodness of fit is computed
            disp(['The loga for the goodness of fit computation is: ',num2str(-OCobj.mLymanN(loc))]);

            % group number
            if isequal(OCobj.mLogisticRegressionHosmerLemeshow.n, 0)
                if OCobj.mNumInGrp<100
                    OCobj.mLogisticRegressionHosmerLemeshow.n = 5;
                else
                    OCobj.mLogisticRegressionHosmerLemeshow.n = 10;
                end
            end
%             OCobj.mLogisticRegressionHosmerLemeshow.n = OCobj.mNumInGrp;

            % prepare goodness of fit computation
            flg=[OCobj.mGrp.mFlgCensor]; % censor flags of patients
            euds = [OCobj.mGrp.mEUD]; euds = euds(loc,:); % the gEUDs of that n
            st = OCobj.mLogisticRegressionMat(loc); % the fitting result of that n
            disp(['The b0 and b1 are: ',num2str(st.b')]);
            % goodness of fit
            [medianVal,~,~,~,~,numComp,numTotal,BetaInv84,BetaInv16] = EventObserved(flg,euds,OCobj.mLogisticRegressionHosmerLemeshow.n); % observations
            [rpb] = glmval(st.b, medianVal,'logit',st.stats); % predictions
            numE= rpb.*numTotal; % expectations
            OCobj.mLogisticRegressionHosmerLemeshow.Chi2 = sum( (numComp-numE).^2 ./ (numTotal.*rpb.*(1-rpb)) );
            OCobj.mLogisticRegressionHosmerLemeshow.df = OCobj.mLogisticRegressionHosmerLemeshow.n - (length(st.b)+1); % (length(st.b)+1) because the loga is not counted in st.b
            OCobj.mLogisticRegressionHosmerLemeshow.p_value = 1 - chi2cdf( OCobj.mLogisticRegressionHosmerLemeshow.Chi2, OCobj.mLogisticRegressionHosmerLemeshow.df );
            disp(OCobj.mLogisticRegressionHosmerLemeshow);
        end
        function OCobj = fLogisticRegressionHosmerLemeshowTestBin_EUD(OCobj,loga)
            % loga determination
            st = [OCobj.mLogisticRegressionMatBin];
            dpf = [st.dev]; % deviations
            [~,loc] = min(dpf); % the minimum deviation corresponds to the maximum likelihood
            disp(['The best loga of Logistic Regression is: ',num2str(-OCobj.mLymanN(loc))]);
            if ~exist('loga','var') || ischar(loga)
                loga = -OCobj.mLymanN(loc);
            end
            [~,loc] = min(abs(OCobj.mLymanN - (-loga))); % find the loga whose goodness of fit is computed
            disp(['The loga for the goodness of fit computation is: ',num2str(-OCobj.mLymanN(loc))]);

            % group number
            if isequal(OCobj.mLogisticRegressionHosmerLemeshowBin.n, 0)
                if OCobj.mNumInGrp<100
                    OCobj.mLogisticRegressionHosmerLemeshowBin.n = 5;
                else
                    OCobj.mLogisticRegressionHosmerLemeshowBin.n = 10;
                end
            end
%             OCobj.mLogisticRegressionHosmerLemeshowBin.n = OCobj.mNumInGrp;

            % prepare goodness of fit computation
            flg=[OCobj.mGrp.mFlgCensor]; % censor flags of patients
            euds = [OCobj.mGrp.mEUD]; euds = euds(loc,:); % the gEUDs of that n
            st = OCobj.mLogisticRegressionMatBin(loc); % the fitting result of that n
            disp(['The b0 and b1 are: ',num2str(st.b')]);
            % goodness of fit
            [medianVal,~,~,~,~,numComp,numTotal,BetaInv84,BetaInv16] = EventObserved(flg,euds,OCobj.mLogisticRegressionHosmerLemeshowBin.n); % observations
            [rpb] = glmval(st.b, medianVal,'logit',st.stats); % predictions
            numE= rpb.*numTotal; % expectations
            OCobj.mLogisticRegressionHosmerLemeshowBin.Chi2 = sum( (numComp-numE).^2 ./ (numTotal.*rpb.*(1-rpb)) );
            OCobj.mLogisticRegressionHosmerLemeshowBin.df = OCobj.mLogisticRegressionHosmerLemeshowBin.n - (length(st.b)+1); % (length(st.b)+1) because the loga is not counted in st.b
            OCobj.mLogisticRegressionHosmerLemeshowBin.p_value = 1 - chi2cdf( OCobj.mLogisticRegressionHosmerLemeshowBin.Chi2, OCobj.mLogisticRegressionHosmerLemeshowBin.df );
            disp(OCobj.mLogisticRegressionHosmerLemeshowBin);
        end
        function OCobj = fLogisticRegressionGTestExact_EUD(OCobj,loga)
            % loga determination
            st = [OCobj.mLogisticRegressionMat];
            dpf = [st.dev]; % deviations
            [~,loc] = min(dpf); % the minimum deviation corresponds to the maximum likelihood
            disp(['The best loga of Logistic Regression is: ',num2str(-OCobj.mLymanN(loc))]);
            if ~exist('loga','var') || ischar(loga)
                loga = -OCobj.mLymanN(loc);
            end
            [~,loc] = min(abs(OCobj.mLymanN - (-loga))); % find the loga whose goodness of fit is computed
            disp(['The loga for the goodness of fit computation is: ',num2str(-OCobj.mLymanN(loc))]);

            % group number
            if isequal(OCobj.mLogisticRegressionGTest.n, 0)
                if OCobj.mNumInGrp<100
                    OCobj.mLogisticRegressionGTest.n = 5;
                else
                    OCobj.mLogisticRegressionGTest.n = 10;
                end
            end
%             OCobj.mLogisticRegressionGTest.n = OCobj.mNumInGrp;

            % prepare goodness of fit computation
            flg=[OCobj.mGrp.mFlgCensor]; % censor flags of patients
            euds = [OCobj.mGrp.mEUD]; euds = euds(loc,:); % the gEUDs of that n
            st = OCobj.mLogisticRegressionMat(loc); % the fitting result of that n
            disp(['The b0 and b1 are: ',num2str(st.b')]);
            % goodness of fit
            [medianVal,~,~,~,~,numComp,numTotal,BetaInv84,BetaInv16] = EventObserved(flg,euds,OCobj.mLogisticRegressionGTest.n); % observations
            [rpb] = glmval(st.b, medianVal,'logit',st.stats); % predictions
            numE= rpb.*numTotal; % expectations
            f = find(numComp);
            OCobj.mLogisticRegressionGTest.Chi2 = 2 * sum( numComp(f) .* log(numComp(f)./numE(f)) );
            OCobj.mLogisticRegressionGTest.df = OCobj.mLogisticRegressionGTest.n - (length(st.b)+1); % (length(st.b)+1) because the loga is not counted in st.b
            OCobj.mLogisticRegressionGTest.p_value = 1 - chi2cdf( OCobj.mLogisticRegressionGTest.Chi2, OCobj.mLogisticRegressionGTest.df );
            disp(OCobj.mLogisticRegressionGTest);
        end
        function OCobj = fLogisticRegressionGTestBin_EUD(OCobj,loga)
            % loga determination
            st = [OCobj.mLogisticRegressionMatBin];
            dpf = [st.dev]; % deviations
            [~,loc] = min(dpf); % the minimum deviation corresponds to the maximum likelihood
            disp(['The best loga of Logistic Regression is: ',num2str(-OCobj.mLymanN(loc))]);
            if ~exist('loga','var') || ischar(loga)
                loga = -OCobj.mLymanN(loc);
            end
            [~,loc] = min(abs(OCobj.mLymanN - (-loga))); % find the loga whose goodness of fit is computed
            disp(['The loga for the goodness of fit computation is: ',num2str(-OCobj.mLymanN(loc))]);

            % group number
            if isequal(OCobj.mLogisticRegressionGTestBin.n, 0)
                if OCobj.mNumInGrp<100
                    OCobj.mLogisticRegressionGTestBin.n = 5;
                else
                    OCobj.mLogisticRegressionGTestBin.n = 10;
                end
            end
%             OCobj.mLogisticRegressionGTestBin.n = OCobj.mNumInGrp;

            % prepare goodness of fit computation
            flg=[OCobj.mGrp.mFlgCensor]; % censor flags of patients
            euds = [OCobj.mGrp.mEUD]; euds = euds(loc,:); % the gEUDs of that n
            st = OCobj.mLogisticRegressionMatBin(loc); % the fitting result of that n
            disp(['The b0 and b1 are: ',num2str(st.b')]);
            % goodness of fit
            [medianVal,~,~,~,~,numComp,numTotal,BetaInv84,BetaInv16] = EventObserved(flg,euds,OCobj.mLogisticRegressionGTestBin.n); % observations
            [rpb] = glmval(st.b, medianVal,'logit',st.stats); % predictions
            numE= rpb.*numTotal; % expectations
            f = find(numComp);
            OCobj.mLogisticRegressionGTestBin.Chi2 = 2 * sum( numComp(f) .* log(numComp(f)./numE(f)) );
            OCobj.mLogisticRegressionGTestBin.df = OCobj.mLogisticRegressionGTestBin.n - (length(st.b)+1); % (length(st.b)+1) because the loga is not counted in st.b
            OCobj.mLogisticRegressionGTestBin.p_value = 1 - chi2cdf( OCobj.mLogisticRegressionGTestBin.Chi2, OCobj.mLogisticRegressionGTestBin.df );
            disp(OCobj.mLogisticRegressionGTestBin);
        end
        function OCobj = fLogisticRegressionPearsonTestExact_EUD(OCobj,loga)
            % loga determination
            st = [OCobj.mLogisticRegressionMat];
            dpf = [st.dev]; % deviations
            [~,loc] = min(dpf); % the minimum deviation corresponds to the maximum likelihood
            disp(['The best loga of Logistic Regression is: ',num2str(-OCobj.mLymanN(loc))]);
            if ~exist('loga','var') || ischar(loga)
                loga = -OCobj.mLymanN(loc);
            end
            [~,loc] = min(abs(OCobj.mLymanN - (-loga))); % find the loga whose goodness of fit is computed
            disp(['The loga for the goodness of fit computation is: ',num2str(-OCobj.mLymanN(loc))]);

            % group number
            if isequal(OCobj.mLogisticRegressionPearson.n, 0)
                if OCobj.mNumInGrp<100
                    OCobj.mLogisticRegressionPearson.n = 5;
                else
                    OCobj.mLogisticRegressionPearson.n = 10;
                end
            end

            % prepare goodness of fit computation
            flg=[OCobj.mGrp.mFlgCensor]; % censor flags of patients
            euds = [OCobj.mGrp.mEUD]; euds = euds(loc,:); % the gEUDs of that n
            st = OCobj.mLogisticRegressionMat(loc); % the fitting result of that n
            disp(['The b0 and b1 are: ',num2str(st.b')]);
            % goodness of fit
            [medianVal,~,~,~,~,numComp,numTotal,BetaInv84,BetaInv16] = EventObserved(flg,euds,OCobj.mLogisticRegressionPearson.n); % observations
            [rpb] = glmval(st.b, medianVal,'logit',st.stats); % predictions
            numE= rpb.*numTotal; % expectations
            OCobj.mLogisticRegressionPearson.Chi2 = sum( (numComp-numE).^2 ./ numE );
            OCobj.mLogisticRegressionPearson.df = OCobj.mLogisticRegressionPearson.n - (length(st.b)+1); % (length(st.b)+1) because the loga is not counted in st.b
            OCobj.mLogisticRegressionPearson.p_value = 1 - chi2cdf( OCobj.mLogisticRegressionPearson.Chi2, OCobj.mLogisticRegressionPearson.df );
            disp(OCobj.mLogisticRegressionPearson);
        end
        function OCobj = fLogisticRegressionPearsonTestBin_EUD(OCobj,loga)
            % loga determination
            st = [OCobj.mLogisticRegressionMatBin];
            dpf = [st.dev]; % deviations
            [~,loc] = min(dpf); % the minimum deviation corresponds to the maximum likelihood
            disp(['The best loga of Logistic Regression is: ',num2str(-OCobj.mLymanN(loc))]);
            if ~exist('loga','var') || ischar(loga)
                loga = -OCobj.mLymanN(loc);
            end
            [~,loc] = min(abs(OCobj.mLymanN - (-loga))); % find the loga whose goodness of fit is computed
            disp(['The loga for the goodness of fit computation is: ',num2str(-OCobj.mLymanN(loc))]);

            % group number
            if isequal(OCobj.mLogisticRegressionPearsonBin.n, 0)
                if OCobj.mNumInGrp<100
                    OCobj.mLogisticRegressionPearsonBin.n = 5;
                else
                    OCobj.mLogisticRegressionPearsonBin.n = 10;
                end
            end

            % prepare goodness of fit computation
            flg=[OCobj.mGrp.mFlgCensor]; % censor flags of patients
            euds = [OCobj.mGrp.mEUD]; euds = euds(loc,:); % the gEUDs of that n
            st = OCobj.mLogisticRegressionMatBin(loc); % the fitting result of that n
            disp(['The b0 and b1 are: ',num2str(st.b')]);
            % goodness of fit
            [medianVal,~,~,~,~,numComp,numTotal,BetaInv84,BetaInv16] = EventObserved(flg,euds,OCobj.mLogisticRegressionPearsonBin.n); % observations
            [rpb] = glmval(st.b, medianVal,'logit',st.stats); % predictions
            numE= rpb.*numTotal; % expectations
            OCobj.mLogisticRegressionPearsonBin.Chi2 = sum( (numComp-numE).^2 ./ numE );
            OCobj.mLogisticRegressionPearsonBin.df = OCobj.mLogisticRegressionPearsonBin.n - (length(st.b)+1); % (length(st.b)+1) because the loga is not counted in st.b
            OCobj.mLogisticRegressionPearsonBin.p_value = 1 - chi2cdf( OCobj.mLogisticRegressionPearsonBin.Chi2, OCobj.mLogisticRegressionPearsonBin.df );
            disp(OCobj.mLogisticRegressionPearsonBin);
        end

        function fLogisticRegressionLikelyhoodBinFig_a_EUD(OCobj,loga,strMarker,lw)
            if ~exist('strMarker','var')
                strMarker = 'b';
            elseif ~ischar(strMarker)
                warning('classOutcomeAnalysis:EndPointObservedFig_EUD:strMarker','strMarker is not a string, omit');
                strMarker = 'b';
            end
            if ~exist('lw','var')
                lw = 1;
            elseif ~isnumeric(lw)
                warning('classOutcomeAnalysis:EndPointObservedFig_EUD:linewidth','lw is not a number, omit');
                lw = 1;
            end
            st = [OCobj.mLogisticRegressionMatBin];
            dpf = [st.dev]; % deviations
            st =[st.stats];
            df = [st.dfe]; % degree of freedom
            dpf = dpf./df; % deviations per degree of freedom
            [~,loc] = min(dpf); % the minimum deviation corresponds to the maximum likelihood
            disp(['best log10(a) of Logistic Regression of gEUD Atlas is: ',num2str(-OCobj.mLymanN(loc))]);
            if ~exist('loga','var') || ischar(loga)
                loga = -OCobj.mLymanN(loc);
            end
            disp(['the log10(a) in coefficient searching is: ',num2str(loga)]);
            [~,loc] = min(abs(OCobj.mLymanN+loga));
            disp('the corresponding coefficients, sd, and 95% CI are:');
            disp(num2str([st(loc).beta, st(loc).se, st(loc).beta-1.96*st(loc).se, st(loc).beta+1.96*st(loc).se]));

            loglikelyhood = -0.5*dpf;
            [mx,loc] = max(loglikelyhood); % the maximum loglikelyhood
            loglikelyhood68 = repmat(mx-0.5* 1/df(loc),size(OCobj.mLymanN));
            loglikelyhood95 = repmat(mx-0.5* (1.96*2) /df(loc),size(OCobj.mLymanN));
            hold on;
            plot(-OCobj.mLymanN, loglikelyhood,strMarker,'LineWidth',lw);
%             plot(-OCobj.mLymanN, loglikelyhood68,strcat(strMarker(1),'--'),'LineWidth',1);
%             plot(-OCobj.mLymanN, loglikelyhood95,strcat(strMarker(1),'-'),'LineWidth',1);
%             plot(-OCobj.mLymanN(loc), loglikelyhood(loc),strMarker,'LineWidth',lw+2);
            set(gca,'xminortick','on','yminortick','on');
            set(gca,'box','on');
            xlabel('log_1_0a'); ylabel('loglikely per degree of freedom');
        end
        function fLogisticRegressionPvalueBinFig_a_EUD(OCobj,strMarker,lw)
            if ~exist('strMarker','var')
                strMarker = 'b';
            elseif ~ischar(strMarker)
                warning('classOutcomeAnalysis:EndPointObservedFig_EUD:strMarker','strMarker is not a string, omit');
                strMarker = 'b';
            end
            if ~exist('lw','var')
                lw = 1;
            elseif ~isnumeric(lw)
                warning('classOutcomeAnalysis:EndPointObservedFig_EUD:linewidth','lw is not a number, omit');
                lw = 1;
            end
            st = [OCobj.mLogisticRegressionMatBin];
            st =[st.stats];
            pvalue = [st.p];
            pvalue = pvalue(2,:); % the p-value corresponding to gEUD
            semilogy(-OCobj.mLymanN, pvalue,strMarker,'LineWidth',lw);
            hold on;
            semilogy(-OCobj.mLymanN, repmat(0.05,size(OCobj.mLymanN)), strcat(strMarker(1),'--'),'LineWidth',1);
            set(gca,'xminortick','on','yminortick','on');
            set(gca,'box','on');
            xlabel('log_1_0a'); ylabel('p-value');
        end
        function loga = fLogisticRegressionRespondingCurveBinFig_a_EUD(OCobj,loga,strMarker,lw)
            % parameters
            if ~exist('strMarker','var')
                strMarker = 'b';
            elseif ~ischar(strMarker)
                warning('classOutcomeAnalysis:EndPointObservedFig_EUD:strMarker','strMarker is not a string, omit');
                strMarker = 'b';
            end
            if ~exist('lw','var')
                lw = 1;
            elseif ~isnumeric(lw)
                warning('classOutcomeAnalysis:EndPointObservedFig_EUD:linewidth','lw is not a number, omit');
                lw = 1;
            end
            st = [OCobj.mLogisticRegressionMatBin];
            dpf = [st.dev]; % deviations
            [~,loc] = min(dpf); % the minimum deviation corresponds to the maximum likelihood
            disp(['best log10(a) of Logistic Regression of gEUD Atlas is: ',num2str(-OCobj.mLymanN(loc))]);
            if ~exist('loga','var') || ischar(loga)
                loga = -OCobj.mLymanN(loc);
            end
            disp(['the log10(a) in responding curve is: ',num2str(loga)]);

            % responding curve
            [~,loc] = min(abs(OCobj.mLymanN+loga)); % the n whose corresponding responding function will be ploted
            st = OCobj.mLogisticRegressionMatBin(loc); % the fitting result of that n
            euds = [OCobj.mGrp.mEUD]; euds = euds(loc,:); % the gEUDs of that n
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
        function loga = fLogisticRegressionGridRespondingCurveBinFig_a_EUD(OCobj,loga,strMarker,lw)
            % parameters
            if ~exist('strMarker','var')
                strMarker = 'b';
            elseif ~ischar(strMarker)
                warning('classOutcomeAnalysis:EndPointObservedFig_EUD:strMarker','strMarker is not a string, omit');
                strMarker = 'b';
            end
            if ~exist('lw','var')
                lw = 1;
            elseif ~isnumeric(lw)
                warning('classOutcomeAnalysis:EndPointObservedFig_EUD:linewidth','lw is not a number, omit');
                lw = 1;
            end
            [mx,loc] = max(OCobj.mLogisticRegressionGridBin.loglikelihood(:));
            [~,~,loc] = ind2sub(size(OCobj.mLogisticRegressionGridBin.loglikelihood),loc);
            disp(['the maximum log likelihood and its "log10(a)" in Logit of Bin gEUD are: ',num2str([mx, -OCobj.mLymanN(loc)])]);
            if ~exist('loga','var') || ischar(loga)
                loga = -OCobj.mLymanN(loc);
            end
            disp(['the plotted response curve is at log10(a) = ',num2str(loga)]);

            % coefficients for the Lyman model
            [~,n] = min(abs(OCobj.mLymanN+loga)); % the n whose corresponding responding function will be ploted
            ll = OCobj.mLogisticRegressionGridBin.loglikelihood(:,:,n); % log likelihood of log10(a) = loga
            mx = max(ll(:));
            [xx,yy] = find(ll == mx); % the coefficients
            b0 = OCobj.mLogisticRegressionGridBin.b0(xx);
            b1 = OCobj.mLogisticRegressionGridBin.b1(yy);
            disp(['the b0 and b1 for the response function are: ',num2str([b0,b1])]);

            % curves
            euds = [OCobj.mGrp.mEUD]; euds = euds(n,:); % the gEUDs of that n
            doses = (0:max(euds))'; % doses (gEUD) whose RP probability will be computed
            rpb = exp(b0+b1*doses);
            rpb = rpb./(1+rpb); % logistic probability
            % 95% CI
            low95 = mx- 0.5* (1.96*2); % CI level
            cxy = contourc(OCobj.mLogisticRegressionGridBin.b0,OCobj.mLogisticRegressionGridBin.b1,ll',[low95,low95]); % coordinates of CI contours at CI level
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
        function fLogisticRegressionGridBinFig_a_loglikelihood_EUD(OCobj,strMarker,lw)
            % parameters
            if ~exist('strMarker','var')
                strMarker = 'b';
            elseif ~ischar(strMarker)
                warning('classOutcomeAnalysis:EndPointObservedFig_EUD:strMarker','strMarker is not a string, omit');
                strMarker = 'b';
            end
            if ~exist('lw','var')
                lw = 1;
            elseif ~isnumeric(lw)
                warning('classOutcomeAnalysis:EndPointObservedFig_EUD:linewidth','lw is not a number, omit');
                lw = 1;
            end

            % find the best "a"
            [mx,loc] = max(OCobj.mLogisticRegressionGridBin.loglikelihood(:));
            [~,~,loc] = ind2sub(size(OCobj.mLogisticRegressionGridBin.loglikelihood),loc);
            disp(['the maximum log likelihood and its "log10(a)" in Logistic Regression of gEUD Atlas are: ',num2str([mx, -OCobj.mLymanN(loc)])]);

            % plot the maximum log likelihood curve w.r.t. log10(a)
            llmx = zeros(size(OCobj.mLymanN)); %log likelihood's max for each log10(a)
            for kk = 1:length(llmx)
                ll = OCobj.mLogisticRegressionGridBin.loglikelihood(:,:,kk);
                llmx(kk) = max(ll(:));
            end
            llmx = llmx / (OCobj.mNumInGrp-2);
            mx = mx / (OCobj.mNumInGrp-2);
            loglikelyhood68 = repmat(mx-0.5* 1 /(OCobj.mNumInGrp-2),size(llmx));
            loglikelyhood95 = repmat(mx-0.5* (1.96*2) /(OCobj.mNumInGrp-2),size(llmx));
            hold on;
            plot(-OCobj.mLymanN,llmx,strMarker,'LineWidth',lw);
            plot(-OCobj.mLymanN(loc),mx,strMarker,'LineWidth',lw+2);
%             plot(-OCobj.mLymanN, loglikelyhood68,strcat(strMarker(1),'--'),'LineWidth',1);
%             plot(-OCobj.mLymanN, loglikelyhood95,strcat(strMarker(1),'-'),'LineWidth',1);
            set(gca,'xminortick','on','yminortick','on');
            set(gca,'box','on');
            xlabel('log_1_0(a)'); ylabel('log likelihood per degree of freedom');
        end
        function fLogisticRegressionGridBinFig_b0_loglikelihood_EUD(OCobj,strMarker,lw)
            % parameters
            if ~exist('strMarker','var')
                strMarker = 'b';
            elseif ~ischar(strMarker)
                warning('classOutcomeAnalysis:EndPointObservedFig_EUD:strMarker','strMarker is not a string, omit');
                strMarker = 'b';
            end
            if ~exist('lw','var')
                lw = 1;
            elseif ~isnumeric(lw)
                warning('classOutcomeAnalysis:EndPointObservedFig_EUD:linewidth','lw is not a number, omit');
                lw = 1;
            end

            % find the best "b0"
            [mx,loc] = max(OCobj.mLogisticRegressionGridBin.loglikelihood(:));
            [loc,~,~] = ind2sub(size(OCobj.mLogisticRegressionGridBin.loglikelihood),loc);
            disp(['the maximum log likelihood and its b0 in Logistic Regression of gEUD Atlas are: ',num2str([mx, OCobj.mLogisticRegressionGridBin.b0(loc)])]);

            % plot the maximum log likelihood curve w.r.t. b0
            llmx = zeros(size(OCobj.mLogisticRegressionGridBin.b0)); %log likelihood's max for each log10(a)
            for kk = 1:length(llmx)
                ll = OCobj.mLogisticRegressionGridBin.loglikelihood(kk,:,:);
                llmx(kk) = max(ll(:));
            end
            llmx = llmx / (OCobj.mNumInGrp-2);
            mx = mx / (OCobj.mNumInGrp-2);
            loglikelyhood68 = repmat(mx-0.5* 1 /(OCobj.mNumInGrp-2),size(llmx));
            loglikelyhood95 = repmat(mx-0.5* (1.96*2) /(OCobj.mNumInGrp-2),size(llmx));
            hold on;
            plot(OCobj.mLogisticRegressionGridBin.b0,llmx,strMarker,'LineWidth',lw);
%             plot(OCobj.mLogisticRegressionGridBin.b0, loglikelyhood68,strcat(strMarker(1),'--'),'LineWidth',1);
%             plot(OCobj.mLogisticRegressionGridBin.b0, loglikelyhood95,strcat(strMarker(1),'-'),'LineWidth',1);
            if length(strMarker)==1
                strMarker = [strMarker,'*'];
            elseif strMarker(2) == '-' || strMarker(2) == '.'
                strMarker(2) = '*';
            end
            plot(OCobj.mLogisticRegressionGridBin.b0(loc),mx,strMarker,'LineWidth',lw+2);
            hold off;
            set(gca,'xminortick','on','yminortick','on');
            set(gca,'box','on');
            xlabel('b0'); ylabel('log likelihood per degree of freedom');
        end
        function fLogisticRegressionGridBinFig_b1_loglikelihood_EUD(OCobj,strMarker,lw)
            % parameters
            if ~exist('strMarker','var')
                strMarker = 'b';
            elseif ~ischar(strMarker)
                warning('classOutcomeAnalysis:EndPointObservedFig_EUD:strMarker','strMarker is not a string, omit');
                strMarker = 'b';
            end
            if ~exist('lw','var')
                lw = 1;
            elseif ~isnumeric(lw)
                warning('classOutcomeAnalysis:EndPointObservedFig_EUD:linewidth','lw is not a number, omit');
                lw = 1;
            end

            % find the best "b1"
            [mx,loc] = max(OCobj.mLogisticRegressionGridBin.loglikelihood(:));
            [~,loc,~] = ind2sub(size(OCobj.mLogisticRegressionGridBin.loglikelihood),loc);
            disp(['the maximum log likelihood and its "b1" in Logistic Regression of gEUD Atlas are: ',num2str([mx, OCobj.mLogisticRegressionGridBin.b1(loc)])]);

            % plot the maximum log likelihood curve w.r.t. b1
            llmx = zeros(size(OCobj.mLogisticRegressionGridBin.b1)); %log likelihood's max for each log10(a)
            for kk = 1:length(llmx)
                ll = OCobj.mLogisticRegressionGridBin.loglikelihood(:,kk,:);
                llmx(kk) = max(ll(:));
            end
            llmx = llmx / (OCobj.mNumInGrp-2);
            mx = mx / (OCobj.mNumInGrp-2);
            loglikelyhood68 = repmat(mx-0.5* 1 /(OCobj.mNumInGrp-2),size(llmx));
            loglikelyhood95 = repmat(mx-0.5* (1.96*2) /(OCobj.mNumInGrp-2),size(llmx));
            hold on;
            plot(OCobj.mLogisticRegressionGridBin.b1,llmx,strMarker,'LineWidth',lw);
%             plot(OCobj.mLogisticRegressionGridBin.b1, loglikelyhood68,strcat(strMarker(1),'--'),'LineWidth',1);
%             plot(OCobj.mLogisticRegressionGridBin.b1, loglikelyhood95,strcat(strMarker(1),'-'),'LineWidth',1);
            if length(strMarker)==1
                strMarker = [strMarker,'*'];
            elseif strMarker(2) == '-' || strMarker(2) == '.'
                strMarker(2) = '*';
            end
            plot(OCobj.mLogisticRegressionGridBin.b1(loc),mx,strMarker,'LineWidth',lw+2);
            hold off;
            set(gca,'xminortick','on','yminortick','on');
            set(gca,'box','on');
            xlabel('b1'); ylabel('log likelihood per degree of freedom');
        end
        function fLogisticRegressionGridBinFig_b0_b1_EUD(OCobj,loga,zdepth)
            if ~exist('zdepth','var')
                zdepth = 10;
            end
            
            % find the best "a"
            [mx,loc] = max(OCobj.mLogisticRegressionGridBin.loglikelihood(:));
            [~,~,loc] = ind2sub(size(OCobj.mLogisticRegressionGridBin.loglikelihood),loc);
            disp(['the "log10(a)" for the maximum likelihood in Logistic Regression of gEUD Atlas is: ',num2str(-OCobj.mLymanN(loc))]);
            if ~exist('loga','var') || ischar(loga)
                loga = -OCobj.mLymanN(loc);
            end
            [~,loc] = min(abs(OCobj.mLymanN+loga)); % the n whose corresponding responding function will be ploted
            disp(['the "log10(a)" for the map of likelihood in Logistic Regression is: ',num2str(-OCobj.mLymanN(loc))]);

            
            % b0 and b1 of the best "log10(a)"
            ll = OCobj.mLogisticRegressionGridBin.loglikelihood(:,:,loc);
            mx = max(ll(:));
            [dd,mm] = find(ll == mx,1);
            disp(['b0 & b1 are: ',num2str([OCobj.mLogisticRegressionGridBin.b0(dd),OCobj.mLogisticRegressionGridBin.b1(mm)])]);
            
            % map of the grid of the best "a"
            if 0 % show the mesh surface
                ll( ll < mx-zdepth ) = -inf;
                surf(OCobj.mLogisticRegressionGridBin.b1,OCobj.mLogisticRegressionGridBin.b0,ll,'EdgeColor','none');
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
                [c,h] = contour(OCobj.mLogisticRegressionGridBin.b0,OCobj.mLogisticRegressionGridBin.b1,ll',[low99,low95,low68]);
                hold on; plot(OCobj.mLogisticRegressionGridBin.b0(dd),OCobj.mLogisticRegressionGridBin.b1(mm),'k*'); hold off;
%                 text_handle = clabel(c,h,'LabelSpacing',100000);%'String',{'low 68%';'low 95%';'low 99.7%'},
%                 set(text_handle,'BackgroundColor',[1 1 .6], 'Edgecolor',[.7 .7 .7], 'String',{'low 68%';'low 95%';'low 99.7%'});
                set(gca,'xminortick','on','yminortick','on');
                set(gca,'box','on');
                xlabel('b0'); ylabel('b1 (Gy^-^1)');
            end
        end
        function fLogisticRegressionGridBinFig_b0_a_EUD(OCobj,zdepth)
            if ~exist('zdepth','var')
                zdepth = 10;
            end
            
            % find the best "m"
            [mx,loc] = max(OCobj.mLogisticRegressionGridBin.loglikelihood(:));
            [~,loc,~] = ind2sub(size(OCobj.mLogisticRegressionGridBin.loglikelihood),loc);
            disp(['the "b1" for the maximum likelihood in Logistic Regression of gEUD Atlas is: ',num2str(OCobj.mLogisticRegressionGridBin.b1(loc))]);
            
%             % reconstruct the log likelihood matrix using image reconstruction
%             ll = -inf(size(OCobj.mLogisticRegressionGridBin.loglikelihood));
%             [mx,loc1] = max(OCobj.mLogisticRegressionGridBin.loglikelihood(:));
%             ll(loc1) = mx;
%             OCobj.mLogisticRegressionGridBin.loglikelihood = imreconstruct(ll,OCobj.mLogisticRegressionGridBin.loglikelihood,26);
            
            % TD50 and log10(a) corresponding to the best m
            ll = squeeze(OCobj.mLogisticRegressionGridBin.loglikelihood(:,loc,:));
            [dd,aa] = find(ll == mx,1);
            disp(['b0 & log10(a) are: ',num2str([OCobj.mLogisticRegressionGridBin.b0(dd),-OCobj.mLymanN(aa)])]);
            
            % map of the grid of the best "m"
            if 0 % show the mesh surface
                ll( ll < mx-zdepth ) = -inf;
                surf(-OCobj.mLymanN,OCobj.mLogisticRegressionGridBin.b0,ll,'EdgeColor','none');
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
                contour(-OCobj.mLymanN,OCobj.mLogisticRegressionGridBin.b0,ll,[low99,low95,low68]);
                hold on; plot(-OCobj.mLymanN(aa),OCobj.mLogisticRegressionGridBin.b0(dd),'k*'); hold off;
                set(gca,'xminortick','on','yminortick','on');
                set(gca,'box','on');
                xlabel('log_1_0(a)'); ylabel('b0');
            end
        end
        function fLogisticRegressionGridBinFig_b1_a_EUD(OCobj,zdepth)
            if ~exist('zdepth','var')
                zdepth = 10;
            end
            
            % find the best "m"
            [mx,loc] = max(OCobj.mLogisticRegressionGridBin.loglikelihood(:));
            [loc,~,~] = ind2sub(size(OCobj.mLogisticRegressionGridBin.loglikelihood),loc);
            disp(['the b0 for the maximum likelihood in Logistic Regression of gEUD Atlas is: ',num2str(OCobj.mLogisticRegressionGridBin.b0(loc))]);
            
            % TD50 and log10(a) corresponding to the best m
            ll = squeeze(OCobj.mLogisticRegressionGridBin.loglikelihood(loc,:,:));
            [mm,aa] = find(ll == mx,1);
            disp(['b1 & log10(a) are: ',num2str([OCobj.mLogisticRegressionGridBin.b1(mm),-OCobj.mLymanN(aa)])]);
            
            % map of the grid of the best "m"
            if 0 % show the mesh surface
                ll( ll < mx-zdepth ) = -inf;
                surf(-OCobj.mLymanN,OCobj.mLogisticRegressionGridBin.b1,ll,'EdgeColor','none');
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
                contour(-OCobj.mLymanN,OCobj.mLogisticRegressionGridBin.b1,ll,[low99,low95,low68]);
                hold on; plot(-OCobj.mLymanN(aa),OCobj.mLogisticRegressionGridBin.b1(mm),'k*'); hold off;
                set(gca,'xminortick','on','yminortick','on');
                set(gca,'box','on');
                xlabel('log_1_0(a)'); ylabel('b1 (Gy^-^1)');
            end
        end

        function OCobj = fLymanAnalysisFittingExact_EUD(OCobj)
            % preparation
            fttp = fittype('normcdf((eud-TD50)/(m*TD50),0,1)', 'independent','eud', 'coefficients',{'TD50','m'});
            fttp = fittype('normcdf(eud,TD50,m)', 'independent','eud', 'coefficients',{'TD50','m'});
            ftopt = fitoptions('method','LinearLeastSquares','algorithm','Levenberg-Marquardt','Display','on');
            euds = [OCobj.mGrp.mEUD]';
            flg = double([OCobj.mGrp.mFlgCensor]');

%             warning('off','MATLAB:singularMatrix');
            % for each mLymanN, compute the parameters of Lyman model, TD50 and m
            kk = 1; % the first fit
            [dose,indx] = sort(euds(:,kk),'ascend');
            comp = ~flg(indx); f = find(comp); % rearrange the complication info so it corresponds to the dose
            comp = cumsum(comp)/OCobj.mNumInGrp; % cumulative complication for dose
            % found the best grid as start point
                llh = OCobj.mLymanGrid.loglikelihood(:,:,kk); %log likeli-hood from grid computation
                mn = max(llh(:)); [dd,mm]=find(llh==mn);
            [fittedobj,goodness,output] = fit(dose,comp,fttp,fitoptions(ftopt,'StartPoint',[OCobj.mLymanGrid.TD50(dd(1)),OCobj.mLymanGrid.m(mm(1))]+10));
            [fittedobj,goodness,output] = fit(euds(:,:,kk),flg,fttp,fitoptions(ftopt,'StartPoint',[OCobj.mLymanGrid.TD50(dd(1)),OCobj.mLymanGrid.m(mm(1))]+10));
            [fittedobj,goodness,output] = fit(euds(:,:,kk),flg,fttp,'StartPoint',[OCobj.mLymanGrid.TD50(dd(1)),OCobj.mLymanGrid.m(mm(1))]+10);
            fitobj = cell(size(OCobj.mLymanN));
            fitobj{kk} = fittedobj;
            goodness = repmat(goodness,size(OCobj.mLymanN));
            output = repmat(output,size(OCobj.mLymanN));
            for kk = 2:size(OCobj.mLymanN,1)
                llh = OCobj.mLymanGrid.loglikelihood(:,:,kk); %log likeli-hood from grid computation
                mn = max(llh(:)); [dd,mm]=find(llh==mn);

                [fittedobj,goodness(kk),output(kk)] = fit(euds(:,kk),flg,fttp,fitoptions(ftopt,'StartPoint',[OCobj.mLymanGrid.TD50(dd(1)),OCobj.mLymanGrid.m(mm(1))]));
                fittedobj
                fitobj{kk} = fittedobj;
            end
%             warning('on','MATLAB:singularMatrix');
            OCobj.mLyman = struct('FitObj',fitobj,'Goodness',goodness,'Output',output);
        end
        function OCobj = fLymanHosmerLemeshowTestAnalysisExact_EUD(OCobj,loga)
            % loga determination
            [~,loc] = max(OCobj.mLymanGrid.loglikelihood(:));
            [~,~,loc] = ind2sub(size(OCobj.mLymanGrid.loglikelihood),loc);
            disp(['the "log10(a)" for the maximum likelihood in Lyman model is: ',num2str(-OCobj.mLymanN(loc))]);
            if ~exist('loga','var') || ischar(loga)
                loga = -OCobj.mLymanN(loc);
            end
            [~,loc] = min(abs(OCobj.mLymanN+loga)); 
            disp(['The loga for the goodness of fit computation is: ',num2str(-OCobj.mLymanN(loc))]);
            ll = OCobj.mLymanGrid.loglikelihood(:,:,loc);
            mx = max(ll(:));
            [dd,mm] = find(ll == mx,1);
            TD50 = OCobj.mLymanGrid.TD50(dd);
            m = OCobj.mLymanGrid.m(mm);
            disp(['TD50 & m are: ',num2str([TD50,m])]);

            % group number
            if isequal(OCobj.mLymanHosmerLemeshow.n, 0)
                if OCobj.mNumInGrp<100
                    OCobj.mLymanHosmerLemeshow.n = 5;
                else
                    OCobj.mLymanHosmerLemeshow.n = 10;
                end
            end
%             OCobj.mLymanHosmerLemeshow.n = OCobj.mNumInGrp;

            % prepare goodness of fit computation
            flg=[OCobj.mGrp.mFlgCensor]; % censor flags of patients
            euds = [OCobj.mGrp.mEUD]; euds = euds(loc,:); % the gEUDs of that n
            % goodness of fit
            [medianVal,~,~,~,~,numComp,numTotal,~,~] = EventObserved(flg,euds,OCobj.mLymanHosmerLemeshow.n); % observations
            rpb = normcdf((medianVal-TD50)/(m*TD50),0,1); % Lyman probability
            numE= rpb.*numTotal; % expectations
            OCobj.mLymanHosmerLemeshow.Chi2 = sum( (numComp-numE).^2 ./ (numTotal.*rpb.*(1-rpb)) );
            OCobj.mLymanHosmerLemeshow.df = OCobj.mLymanHosmerLemeshow.n - 2 - 1;
            OCobj.mLymanHosmerLemeshow.p_value = 1 - chi2cdf( OCobj.mLymanHosmerLemeshow.Chi2, OCobj.mLymanHosmerLemeshow.df );
            disp(OCobj.mLymanHosmerLemeshow);
        end
        function OCobj = fLymanHosmerLemeshowTestAnalysisBin_EUD(OCobj,loga)
            % loga determination
            [~,loc] = max(OCobj.mLymanGridBin.loglikelihood(:));
            [~,~,loc] = ind2sub(size(OCobj.mLymanGridBin.loglikelihood),loc);
            disp(['the "log10(a)" for the maximum likelihood in Lyman model is: ',num2str(-OCobj.mLymanN(loc))]);
            if ~exist('loga','var') || ischar(loga)
                loga = -OCobj.mLymanN(loc);
            end
            [~,loc] = min(abs(OCobj.mLymanN+loga)); % the n whose corresponding responding function will be ploted
            disp(['The loga for the goodness of fit computation is: ',num2str(-OCobj.mLymanN(loc))]);
            ll = OCobj.mLymanGridBin.loglikelihood(:,:,loc);
            mx = max(ll(:));
            [dd,mm] = find(ll == mx,1);
            TD50 = OCobj.mLymanGridBin.TD50(dd);
            m = OCobj.mLymanGridBin.m(mm);
            disp(['TD50 & m are: ',num2str([TD50,m])]);

            % group number
            if isequal(OCobj.mLymanHosmerLemeshowBin.n, 0)
                if OCobj.mNumInGrp<100
                    OCobj.mLymanHosmerLemeshowBin.n = 5;
                else
                    OCobj.mLymanHosmerLemeshowBin.n = 10;
                end
            end
%             OCobj.mLymanHosmerLemeshowBin.n = OCobj.mNumInGrp;

            % prepare goodness of fit computation
            flg=[OCobj.mGrp.mFlgCensor]; % censor flags of patients
            euds = [OCobj.mGrp.mEUD]; euds = euds(loc,:); % the gEUDs of that n
            % goodness of fit
            [medianVal,~,~,~,~,numComp,numTotal,BetaInv84,BetaInv16] = EventObserved(flg,euds,OCobj.mLymanHosmerLemeshowBin.n); % observations
            rpb = normcdf((medianVal-TD50)/(m*TD50),0,1); % Lyman probability
            numE= rpb.*numTotal; % expectations
            OCobj.mLymanHosmerLemeshowBin.Chi2 = sum( (numComp-numE).^2 ./ (numTotal.*rpb.*(1-rpb)) );
            OCobj.mLymanHosmerLemeshowBin.df = OCobj.mLymanHosmerLemeshowBin.n - 2 - 1;
            OCobj.mLymanHosmerLemeshowBin.p_value = 1 - chi2cdf( OCobj.mLymanHosmerLemeshowBin.Chi2, OCobj.mLymanHosmerLemeshowBin.df );
            disp(OCobj.mLymanHosmerLemeshowBin);
        end
        function OCobj = fLymanGTestAnalysisExact_EUD(OCobj,loga)
            % loga determination
            [~,loc] = max(OCobj.mLymanGrid.loglikelihood(:));
            [~,~,loc] = ind2sub(size(OCobj.mLymanGrid.loglikelihood),loc);
            disp(['the "log10(a)" for the maximum likelihood in Lyman model is: ',num2str(-OCobj.mLymanN(loc))]);
            if ~exist('loga','var') || ischar(loga)
                loga = -OCobj.mLymanN(loc);
            end
            [~,loc] = min(abs(OCobj.mLymanN+loga)); % the n whose corresponding responding function will be ploted
            disp(['The loga for the goodness of fit computation is: ',num2str(-OCobj.mLymanN(loc))]);
            ll = OCobj.mLymanGrid.loglikelihood(:,:,loc);
            mx = max(ll(:));
            [dd,mm] = find(ll == mx,1);
            TD50 = OCobj.mLymanGrid.TD50(dd);
            m = OCobj.mLymanGrid.m(mm);
            disp(['TD50 & m are: ',num2str([TD50,m])]);

            % group number
            if isequal(OCobj.mLymanGTest.n, 0)
                if OCobj.mNumInGrp<100
                    OCobj.mLymanGTest.n = 5;
                else
                    OCobj.mLymanGTest.n = 10;
                end
            end
%             OCobj.mLymanGTest.n = OCobj.mNumInGrp;

            % prepare goodness of fit computation
            flg=[OCobj.mGrp.mFlgCensor]; % censor flags of patients
            euds = [OCobj.mGrp.mEUD]; euds = euds(loc,:); % the gEUDs of that n
            % goodness of fit
            [medianVal,~,~,~,~,numComp,numTotal,BetaInv84,BetaInv16] = EventObserved(flg,euds,OCobj.mLymanGTest.n); % observations
            rpb = normcdf((medianVal-TD50)/(m*TD50),0,1); % Lyman probability
            numE= rpb.*numTotal; % expectations
            f = find(numComp);
            OCobj.mLymanGTest.Chi2 = 2 * sum( numComp(f) .* log(numComp(f)./numE(f)) );
            OCobj.mLymanGTest.df = OCobj.mLymanGTest.n - 2 - 1;
            OCobj.mLymanGTest.p_value = 1 - chi2cdf( OCobj.mLymanGTest.Chi2, OCobj.mLymanGTest.df );
            disp(OCobj.mLymanGTest);
        end
        function OCobj = fLymanGTestAnalysisBin_EUD(OCobj,loga)
            % loga determination
            [~,loc] = max(OCobj.mLymanGridBin.loglikelihood(:));
            [~,~,loc] = ind2sub(size(OCobj.mLymanGridBin.loglikelihood),loc);
            disp(['the "log10(a)" for the maximum likelihood in Lyman model is: ',num2str(-OCobj.mLymanN(loc))]);
            if ~exist('loga','var') || ischar(loga)
                loga = -OCobj.mLymanN(loc);
            end
            [~,loc] = min(abs(OCobj.mLymanN+loga)); % the n whose corresponding responding function will be ploted
            disp(['The loga for the goodness of fit computation is: ',num2str(-OCobj.mLymanN(loc))]);
            ll = OCobj.mLymanGridBin.loglikelihood(:,:,loc);
            mx = max(ll(:));
            [dd,mm] = find(ll == mx,1);
            TD50 = OCobj.mLymanGridBin.TD50(dd);
            m = OCobj.mLymanGridBin.m(mm);
            disp(['TD50 & m are: ',num2str([TD50,m])]);

            % group number
            if isequal(OCobj.mLymanGTestBin.n, 0)
                if OCobj.mNumInGrp<100
                    OCobj.mLymanGTestBin.n = 5;
                else
                    OCobj.mLymanGTestBin.n = 10;
                end
            end
%             OCobj.mLymanGTestBin.n = OCobj.mNumInGrp;

            % prepare goodness of fit computation
            flg=[OCobj.mGrp.mFlgCensor]; % censor flags of patients
            euds = [OCobj.mGrp.mEUD]; euds = euds(loc,:); % the gEUDs of that n
            % goodness of fit
            [medianVal,~,~,~,~,numComp,numTotal,BetaInv84,BetaInv16] = EventObserved(flg,euds,OCobj.mLymanGTestBin.n); % observations
            rpb = normcdf((medianVal-TD50)/(m*TD50),0,1); % Lyman probability
            numE= rpb.*numTotal; % expectations
            f = find(numComp);
            OCobj.mLymanGTestBin.Chi2 = 2 * sum( numComp(f) .* log(numComp(f)./numE(f)) );
            OCobj.mLymanGTestBin.df = OCobj.mLymanGTestBin.n - 2 - 1;
            OCobj.mLymanGTestBin.p_value = 1 - chi2cdf( OCobj.mLymanGTestBin.Chi2, OCobj.mLymanGTestBin.df );
            disp(OCobj.mLymanGTestBin);
        end
        function OCobj = fLymanPearsonTestAnalysisExact_EUD(OCobj,loga)
            % loga determination
            [~,loc] = max(OCobj.mLymanGrid.loglikelihood(:));
            [~,~,loc] = ind2sub(size(OCobj.mLymanGrid.loglikelihood),loc);
            disp(['the "log10(a)" for the maximum likelihood in Lyman model is: ',num2str(-OCobj.mLymanN(loc))]);
            if ~exist('loga','var') || ischar(loga)
                loga = -OCobj.mLymanN(loc);
            end
            [~,loc] = min(abs(OCobj.mLymanN+loga)); % the n whose corresponding responding function will be ploted
            disp(['The loga for the goodness of fit computation is: ',num2str(-OCobj.mLymanN(loc))]);
            ll = OCobj.mLymanGrid.loglikelihood(:,:,loc);
            mx = max(ll(:));
            [dd,mm] = find(ll == mx,1);
            TD50 = OCobj.mLymanGrid.TD50(dd);
            m = OCobj.mLymanGrid.m(mm);
            disp(['TD50 & m are: ',num2str([TD50,m])]);

            % group number
            if isequal(OCobj.mLymanPearson.n, 0)
                if OCobj.mNumInGrp<100
                    OCobj.mLymanPearson.n = 5;
                else
                    OCobj.mLymanPearson.n = 10;
                end
            end

            % prepare goodness of fit computation
            flg=[OCobj.mGrp.mFlgCensor]; % censor flags of patients
            euds = [OCobj.mGrp.mEUD]; euds = euds(loc,:); % the gEUDs of that n
            % goodness of fit
            [medianVal,~,~,~,~,numComp,numTotal,BetaInv84,BetaInv16] = EventObserved(flg,euds,OCobj.mLymanPearson.n); % observations
            rpb = normcdf((medianVal-TD50)/(m*TD50),0,1); % Lyman probability
            numE= rpb.*numTotal; % expectations
            OCobj.mLymanPearson.Chi2 = sum( (numComp-numE).^2 ./ numE );
            OCobj.mLymanPearson.df = OCobj.mLymanPearson.n - 2 - 1;
            OCobj.mLymanPearson.p_value = 1 - chi2cdf( OCobj.mLymanPearson.Chi2, OCobj.mLymanPearson.df );
            disp(OCobj.mLymanPearson);
        end
        function OCobj = fLymanPearsonTestAnalysisBin_EUD(OCobj,loga)
            % loga determination
            [~,loc] = max(OCobj.mLymanGridBin.loglikelihood(:));
            [~,~,loc] = ind2sub(size(OCobj.mLymanGridBin.loglikelihood),loc);
            disp(['the "log10(a)" for the maximum likelihood in Lyman model is: ',num2str(-OCobj.mLymanN(loc))]);
            if ~exist('loga','var') || ischar(loga)
                loga = -OCobj.mLymanN(loc);
            end
            [~,loc] = min(abs(OCobj.mLymanN+loga)); % the n whose corresponding responding function will be ploted
            disp(['The loga for the goodness of fit computation is: ',num2str(-OCobj.mLymanN(loc))]);
            ll = OCobj.mLymanGridBin.loglikelihood(:,:,loc);
            mx = max(ll(:));
            [dd,mm] = find(ll == mx,1);
            TD50 = OCobj.mLymanGridBin.TD50(dd);
            m = OCobj.mLymanGridBin.m(mm);
            disp(['TD50 & m are: ',num2str([TD50,m])]);

            % group number
            if isequal(OCobj.mLymanPearsonBin.n, 0)
                if OCobj.mNumInGrp<100
                    OCobj.mLymanPearsonBin.n = 5;
                else
                    OCobj.mLymanPearsonBin.n = 10;
                end
            end

            % prepare goodness of fit computation
            flg=[OCobj.mGrp.mFlgCensor]; % censor flags of patients
            euds = [OCobj.mGrp.mEUD]; euds = euds(loc,:); % the gEUDs of that n
            % goodness of fit
            [medianVal,~,~,~,~,numComp,numTotal,BetaInv84,BetaInv16] = EventObserved(flg,euds,OCobj.mLymanPearsonBin.n); % observations
            rpb = normcdf((medianVal-TD50)/(m*TD50),0,1); % Lyman probability
            numE= rpb.*numTotal; % expectations
            OCobj.mLymanPearsonBin.Chi2 = sum( (numComp-numE).^2 ./ numE );
            OCobj.mLymanPearsonBin.df = OCobj.mLymanPearsonBin.n - 2 - 1;
            OCobj.mLymanPearsonBin.p_value = 1 - chi2cdf( OCobj.mLymanPearsonBin.Chi2, OCobj.mLymanPearsonBin.df );
            disp(OCobj.mLymanPearsonBin);
        end
        function fLymanGridBinFig_TD50_m_EUD(OCobj,loga,zdepth)
            if ~exist('zdepth','var')
                zdepth = 10;
            end
            
            % find the best "a"
            [mx,loc] = max(OCobj.mLymanGridBin.loglikelihood(:));
            [~,~,loc] = ind2sub(size(OCobj.mLymanGridBin.loglikelihood),loc);
            disp(['the "log10(a)" for the maximum likelihood in Lyman model of gEUD Atlas is: ',num2str(-OCobj.mLymanN(loc))]);
            if ~exist('loga','var') || ischar(loga)
                loga = -OCobj.mLymanN(loc);
            end
            [~,loc] = min(abs(OCobj.mLymanN+loga)); % the n whose corresponding responding function will be ploted
            disp(['the "log10(a)" for the map of likelihood in Lyman model is: ',num2str(-OCobj.mLymanN(loc))]);
            
            % TD50 and gamma of the best "log10(a)"
            ll = OCobj.mLymanGridBin.loglikelihood(:,:,loc);
            mx = max(ll(:));
            [dd,mm] = find(ll == mx,1);
            disp(['TD50 & m are: ',num2str([OCobj.mLymanGridBin.TD50(dd),OCobj.mLymanGridBin.m(mm)])]);
            
            % map of the grid of the best "a"
            if 0 % show the mesh surface
                ll( ll < mx-zdepth ) = -inf;
                surf(OCobj.mLymanGridBin.m,OCobj.mLymanGridBin.TD50,ll,'EdgeColor','none');
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
                contour(OCobj.mLymanGridBin.TD50,OCobj.mLymanGridBin.m,ll',[low99,low95,low68]);
                hold on; plot(OCobj.mLymanGridBin.TD50(dd),OCobj.mLymanGridBin.m(mm),'k*'); hold off;
                set(gca,'xminortick','on','yminortick','on');
                set(gca,'box','on');
                xlabel('TD50 (Gy)'); ylabel('m');
            end
        end
        function fLymanGridBinFig_TD50_a_EUD(OCobj,zdepth)
            if ~exist('zdepth','var')
                zdepth = 10;
            end
            
            % find the best "m"
            [mx,loc] = max(OCobj.mLymanGridBin.loglikelihood(:));
            [~,loc,~] = ind2sub(size(OCobj.mLymanGridBin.loglikelihood),loc);
            disp(['the "m" for the maximum likelihood in Lyman model of gEUD Atlas is: ',num2str(OCobj.mLymanGridBin.m(loc))]);
            
            % TD50 and log10(a) corresponding to the best m
            ll = squeeze(OCobj.mLymanGridBin.loglikelihood(:,loc,:));
            [dd,aa] = find(ll == mx,1);
            disp(['TD50 & a are: ',num2str([OCobj.mLymanGridBin.TD50(dd),-OCobj.mLymanN(aa)])]);
            
            % map of the grid of the best "m"
            if 0 % show the mesh surface
                ll( ll < mx-zdepth ) = -inf;
                surf(-OCobj.mLymanN,OCobj.mLymanGridBin.TD50,ll,'EdgeColor','none');
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
                contour(-OCobj.mLymanN,OCobj.mLymanGridBin.TD50,ll,[low99,low95,low68]);
                hold on; plot(-OCobj.mLymanN(aa),OCobj.mLymanGridBin.TD50(dd),'k*'); hold off;
                set(gca,'xminortick','on','yminortick','on');
                set(gca,'box','on');
                xlabel('log_1_0(a)'); ylabel('TD50 (Gy)');
            end
        end
        function fLymanGridBinFig_m_a_EUD(OCobj,zdepth)
            if ~exist('zdepth','var')
                zdepth = 10;
            end
            
            % find the best "m"
            [mx,loc] = max(OCobj.mLymanGridBin.loglikelihood(:));
            [loc,~,~] = ind2sub(size(OCobj.mLymanGridBin.loglikelihood),loc);
            disp(['the "TD50" for the maximum likelihood in Lyman model of gEUD Atlas is: ',num2str(OCobj.mLymanGridBin.TD50(loc))]);
            
            % TD50 and log10(a) corresponding to the best m
            ll = squeeze(OCobj.mLymanGridBin.loglikelihood(loc,:,:));
            [mm,aa] = find(ll == mx,1);
            disp(['m & a are: ',num2str([OCobj.mLymanGridBin.m(mm),-OCobj.mLymanN(aa)])]);
            
            % map of the grid of the best "m"
            if 0 % show the mesh surface
                ll( ll < mx-zdepth ) = -inf;
                surf(-OCobj.mLymanN,OCobj.mLymanGridBin.m,ll,'EdgeColor','none');
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
                contour(-OCobj.mLymanN,OCobj.mLymanGridBin.m,ll,[low99,low95,low68]);
                hold on; plot(-OCobj.mLymanN(aa),OCobj.mLymanGridBin.m(mm),'k*'); hold off;
                set(gca,'xminortick','on','yminortick','on');
                set(gca,'box','on');
                xlabel('log_1_0(a)'); ylabel('m');
            end
        end
        function fLymanGridBinFig_a_loglikelihood(OCobj,strMarker,lw)
            % parameters
            if ~exist('strMarker','var')
                strMarker = 'b';
            elseif ~ischar(strMarker)
                warning('classOutcomeAnalysis:EndPointObservedFig_EUD:strMarker','strMarker is not a string, omit');
                strMarker = 'b';
            end
            if ~exist('lw','var')
                lw = 1;
            elseif ~isnumeric(lw)
                warning('classOutcomeAnalysis:EndPointObservedFig_EUD:linewidth','lw is not a number, omit');
                lw = 1;
            end

            % find the best "a"
            [mx,loc] = max(OCobj.mLymanGridBin.loglikelihood(:));
            [~,~,loc] = ind2sub(size(OCobj.mLymanGridBin.loglikelihood),loc);
            disp(['the maximum log likelihood and its corresponding "log10(a)" in Lyman model of gEUD Atlas are: ',num2str([mx, -OCobj.mLymanN(loc)])]);

            % plot the maximum log likelihood curve w.r.t. log10(a)
            llmx = zeros(size(OCobj.mLymanN)); %log likelihood's max for each log10(a)
            for kk = 1:length(llmx)
                ll = OCobj.mLymanGridBin.loglikelihood(:,:,kk);
                llmx(kk) = max(ll(:));
            end
            llmx = llmx / (OCobj.mNumInGrp-2);
            mx = mx / (OCobj.mNumInGrp-2);
          
            
            lowlog68 = mx-0.5*1/(OCobj.mNumInGrp-2);
            lowlog95 = mx-0.5*(1.96*2)/(OCobj.mNumInGrp-2);
            
            [~, CI68loga] = ConfidenceInterval(OCobj.mLymanN,flipdim(llmx,1), lowlog68);
            [~, CI95loga] = ConfidenceInterval(OCobj.mLymanN,flipdim(llmx,1), lowlog95);
            

            disp(['68% log(a) CIs: ',...
                num2str(-OCobj.mLymanN(loc)),...
                ' [',...
                num2str(CI68loga(1)),...
                ', ',...
                num2str(CI68loga(2)),...
                ']']);
            disp(['95% log(a) CIs: ',...
                num2str(-OCobj.mLymanN(loc)),...
                ' [',...
                num2str(CI95loga(1)),...
                ', ',...
                num2str(CI95loga(2)),...
                ']']);
            
            %loglikelyhood68 = repmat(lowlog68),size(llmx));
            %loglikelyhood95 = repmat(lowlog95),size(llmx));
            
            hold on;
            plot(-OCobj.mLymanN,llmx,strMarker,'LineWidth',lw);
            plot(-OCobj.mLymanN(loc),mx,strMarker,'LineWidth',lw+2);
%             plot(-OCobj.mLymanN, loglikelyhood68,strcat(strMarker(1),'--'),'LineWidth',1);
%             plot(-OCobj.mLymanN, loglikelyhood95,strcat(strMarker(1),'-'),'LineWidth',1);
            set(gca,'xminortick','on','yminortick','on');
            set(gca,'box','on');
            xlabel('log_1_0(a)'); ylabel('log likelihood per degree of freedom');
        end
        function fLymanGridBinFig_TD50_loglikelihood(OCobj,strMarker,lw)
            % parameters
            if ~exist('strMarker','var')
                strMarker = 'b';
            elseif ~ischar(strMarker)
                warning('classOutcomeAnalysis:EndPointObservedFig_EUD:strMarker','strMarker is not a string, omit');
                strMarker = 'b';
            end
            if ~exist('lw','var')
                lw = 1;
            elseif ~isnumeric(lw)
                warning('classOutcomeAnalysis:EndPointObservedFig_EUD:linewidth','lw is not a number, omit');
                lw = 1;
            end

            % find the best "TD50"
            [mx,loc] = max(OCobj.mLymanGridBin.loglikelihood(:));
            [loc,~,~] = ind2sub(size(OCobj.mLymanGridBin.loglikelihood),loc);
            disp(['the maximum log likelihood and its corresponding "TD50" in Lyman model of gEUD Atlas are: ',num2str([mx, OCobj.mLymanGridBin.TD50(loc)])]);

            % plot the maximum log likelihood curve w.r.t. TD50
            llmx = zeros(size(OCobj.mLymanGridBin.TD50)); %log likelihood's max for each log10(a)
            for kk = 1:length(llmx)
                ll = OCobj.mLymanGridBin.loglikelihood(kk,:,:);
                llmx(kk) = max(ll(:));
            end
            llmx = llmx / (OCobj.mNumInGrp-2);
            mx = mx / (OCobj.mNumInGrp-2);
             
            lowlog68 = mx-0.5*1/(OCobj.mNumInGrp-2);
            lowlog95 = mx-0.5*(1.96*2)/(OCobj.mNumInGrp-2);
            
            [~, CI68TD50] = ConfidenceInterval(OCobj.mLymanGrid.TD50,llmx, lowlog68);
            [~, CI95TD50] = ConfidenceInterval(OCobj.mLymanGrid.TD50,llmx, lowlog95);
            
            disp(['68% TD50 CIs: ',...
                num2str(OCobj.mLymanGrid.TD50(loc)),...
                ' [',...
                num2str(CI68TD50(1)),...
                ', ',...
                num2str(CI68TD50(2)),...
                ']']);
            disp(['95% TD50 CIs: ',...
                num2str(OCobj.mLymanGrid.TD50(loc)),...
                ' [',...
                num2str(CI95TD50(1)),...
                ', ',...
                num2str(CI95TD50(2)),...
                ']']);
            
            
            loglikelyhood68 = repmat(lowlog68,size(llmx));
            loglikelyhood95 = repmat(lowlog95,size(llmx));
            hold on;
            plot(OCobj.mLymanGridBin.TD50,llmx,strMarker,'LineWidth',lw);
            plot(OCobj.mLymanGridBin.TD50(loc),mx,strMarker,'LineWidth',lw+2);
%             plot(OCobj.mLymanGrid.TD50, loglikelyhood68,strcat(strMarker(1),'--'),'LineWidth',1);
%             plot(OCobj.mLymanGrid.TD50, loglikelyhood95,strcat(strMarker(1),'-'),'LineWidth',1);
            set(gca,'xminortick','on','yminortick','on');
            set(gca,'box','on');
            xlabel('TD50'); ylabel('log likelihood per degree of freedom');
        end
        function fLymanGridBinFig_m_loglikelihood(OCobj,strMarker,lw)
            % parameters
            if ~exist('strMarker','var')
                strMarker = 'b';
            elseif ~ischar(strMarker)
                warning('classOutcomeAnalysis:EndPointObservedFig_EUD:strMarker','strMarker is not a string, omit');
                strMarker = 'b';
            end
            if ~exist('lw','var')
                lw = 1;
            elseif ~isnumeric(lw)
                warning('classOutcomeAnalysis:EndPointObservedFig_EUD:linewidth','lw is not a number, omit');
                lw = 1;
            end

            % find the best "m"
            [mx,loc] = max(OCobj.mLymanGridBin.loglikelihood(:));
            [~,loc,~] = ind2sub(size(OCobj.mLymanGridBin.loglikelihood),loc);
            disp(['the maximum log likelihood and its corresponding "m" in Lyman model of gEUD Atlas are: ',num2str([mx, OCobj.mLymanGridBin.m(loc)])]);

            % plot the maximum log likelihood curve w.r.t. TD50
            llmx = zeros(size(OCobj.mLymanGridBin.m)); %log likelihood's max for each log10(a)
            for kk = 1:length(llmx)
                ll = OCobj.mLymanGridBin.loglikelihood(:,kk,:);
                llmx(kk) = max(ll(:));
            end
            llmx = llmx / (OCobj.mNumInGrp-2);
            mx = mx / (OCobj.mNumInGrp-2);

              
            lowlog68 = mx-0.5*1/(OCobj.mNumInGrp-2);
            lowlog95 = mx-0.5*(1.96*2)/(OCobj.mNumInGrp-2);
            
            [~, CI68M] = ConfidenceInterval(OCobj.mLymanGrid.m,llmx, lowlog68);
            [~, CI95M] = ConfidenceInterval(OCobj.mLymanGrid.m,llmx, lowlog95);
            
            disp(['68% m CIs: ',...
                num2str(OCobj.mLymanGrid.m(loc)),...
                ' [',...
                num2str(CI68M(1)),...
                ', ',...
                num2str(CI68M(2)),...
                ']']);
            disp(['95% m CIs: ',...
                num2str(OCobj.mLymanGrid.m(loc)),...
                ' [',...
                num2str(CI95M(1)),...
                ', ',...
                num2str(CI95M(2)),...
                ']']);
            
            
            loglikelyhood68 = repmat(lowlog68,size(llmx));
            loglikelyhood95 = repmat(lowlog95,size(llmx));
            
            
            hold on;
            plot(OCobj.mLymanGridBin.m,llmx,strMarker,'LineWidth',lw);
            plot(OCobj.mLymanGridBin.m(loc),mx,strMarker,'LineWidth',lw+2);
%             plot(OCobj.mLymanGrid.m, loglikelyhood68,strcat(strMarker(1),'--'),'LineWidth',1);
%             plot(OCobj.mLymanGrid.m, loglikelyhood95,strcat(strMarker(1),'-'),'LineWidth',1);
            set(gca,'xminortick','on','yminortick','on');
            set(gca,'box','on');
            xlabel('m'); ylabel('log likelihood per degree of freedom');
        end
        function loga = fLymanGridResponseBinFig_a_EUD(OCobj,loga,strMarker,lw)
            % parameters
            if ~exist('strMarker','var')
                strMarker = 'b';
            elseif ~ischar(strMarker)
                warning('classOutcomeAnalysis:EndPointObservedFig_EUD:strMarker','strMarker is not a string, omit');
                strMarker = 'b';
            end
            if ~exist('lw','var')
                lw = 1;
            elseif ~isnumeric(lw)
                warning('classOutcomeAnalysis:EndPointObservedFig_EUD:linewidth','lw is not a number, omit');
                lw = 1;
            end
            [mx,loc] = max(OCobj.mLymanGridBin.loglikelihood(:));
            [~,~,loc] = ind2sub(size(OCobj.mLymanGridBin.loglikelihood),loc);
            disp(['the maximum log likelihood and its "log10(a)" in Lyman model of Bin gEUD are: ',num2str([mx, -OCobj.mLymanN(loc)])]);
            if ~exist('loga','var') || ischar(loga)
                loga = -OCobj.mLymanN(loc);
            end
            disp(['the plotted response curve is at log10(a) = ',num2str(loga)]);

            % coefficients for the Lyman model
            [~,n] = min(abs(OCobj.mLymanN+loga)); % the n whose corresponding responding function will be ploted
            ll = OCobj.mLymanGridBin.loglikelihood(:,:,n); % log likelihood of log10(a) = loga
            mx = max(ll(:));
            [xx,yy] = find(ll == mx); % the coefficients
            TD50 = OCobj.mLymanGridBin.TD50(xx);
            m = OCobj.mLymanGridBin.m(yy);
            disp(['the TD50 and m for the response function are: ',num2str([TD50,m])]);

            % curves
            euds = [OCobj.mGrp.mEUD]; euds = euds(n,:); % the gEUDs of that n
            doses = (0:max(euds))'; % doses (gEUD) whose RP probability will be computed
            rpb = normcdf((doses-TD50)/(m*TD50),0,1); % Lyman probability
            % 95% CI
            low95 = mx-0.5* (1.96*2); % CI level
            cxy = contourc(OCobj.mLymanGridBin.TD50,OCobj.mLymanGridBin.m,ll',[low95,low95]); % coordinates of CI contours at CI level
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

        function fWriteXls_EUD(OCobj) % write results from EUD into a spread sheet
            warning('off','MATLAB:xlswrite:AddSheet');
            if OCobj.flgColOutput
                % EUD
                xlswrite(OCobj.xlsFile_output,OCobj.log10n',strcat(OCobj.xlsSheet,'_EUD'),'B1');
                if ~isempty(OCobj.PatientCode)
                    xlswrite(OCobj.xlsFile_output,OCobj.PatientCode(OCobj.PatientRows),strcat(OCobj.xlsSheet,'_EUD'),'A2');
                end
                if ~isempty(OCobj.mEUD)
                    xlswrite(OCobj.xlsFile_output,OCobj.mEUD(OCobj.PatientRows,:),strcat(OCobj.xlsSheet,'_EUD'),'B2');
                end
                
                % total patients at log10n and dose
                xlswrite(OCobj.xlsFile_output,OCobj.log10n', strcat(OCobj.xlsSheet,'_Total'),'B1');
                if ~isempty(OCobj.mAtlasTotal);
%                     doses = (0:size(OCobj.mAtlasTotal,1)-1)*OCobj.GyStep;
                    xlswrite(OCobj.xlsFile_output,OCobj.AtlasBins,strcat(OCobj.xlsSheet,'_Total'),'A2');
                    xlswrite(OCobj.xlsFile_output,OCobj.mAtlasTotal,strcat(OCobj.xlsSheet,'_Total'),'B2');
                end
                
                % complication patients at log10n and dose
                xlswrite(OCobj.xlsFile_output,OCobj.log10n', strcat(OCobj.xlsSheet,'_Comp'),'B1');
                if ~isempty(OCobj.mAtlasComp)
%                     doses = (0:size(OCobj.mAtlasComp,1)-1)*OCobj.GyStep;
                    xlswrite(OCobj.xlsFile_output,OCobj.AtlasBins,strcat(OCobj.xlsSheet,'_Comp'),'A2');
                    xlswrite(OCobj.xlsFile_output,OCobj.mAtlasComp,strcat(OCobj.xlsSheet,'_Comp'),'B2');
                end
                
                % Beta probability
                if ~isempty(OCobj.BetaCumulativeMat)
%                     doses = (0:size(OCobj.mAtlasComp,1)-1)*OCobj.GyStep;
                    for k=1:length(OCobj.mBetaCumulativeThreshold)
                        xlswrite(OCobj.xlsFile_output, OCobj.log10n', strcat(OCobj.xlsSheet,'_prob_',num2str(OCobj.mBetaCumulativeThreshold(k))),'B1');
                        xlswrite(OCobj.xlsFile_output,OCobj.AtlasBins,strcat(OCobj.xlsSheet,'_prob_',num2str(OCobj.mBetaCumulativeThreshold(k))),'A2');
                        xlswrite(OCobj.xlsFile_output,OCobj.BetaCumulativeMat(:,:,k),strcat(OCobj.xlsSheet,'_prob_',num2str(OCobj.mBetaCumulativeThreshold(k))),'B2');
                    end
                end
                if ~isempty(OCobj.BetaInverseMat)
%                     doses = (0:size(OCobj.mAtlasComp,1)-1)*OCobj.GyStep;
                    for k=1:length(OCobj.mBetaInverseThreshold)
                        xlswrite(OCobj.xlsFile_output,OCobj.log10n', strcat(OCobj.xlsSheet,'_Low_',num2str(OCobj.mBetaInverseThreshold(k))),'B1');
                        xlswrite(OCobj.xlsFile_output,OCobj.AtlasBins,strcat(OCobj.xlsSheet,'_Low_',num2str(OCobj.mBetaInverseThreshold(k))),'A2');
                        xlswrite(OCobj.xlsFile_output,OCobj.BetaCumulativeMat(:,:,k),strcat(OCobj.xlsSheet,'_Low_',num2str(OCobj.mBetaInverseThreshold(k))),'B2');
                    end
                end
            else
                % EUD
                xlswrite(OCobj.xlsFile_output, OCobj.log10n, strcat(OCobj.xlsSheet,'_EUD'),'A2');
                if ~isempty(OCobj.PatientCode)
                    xlswrite(OCobj.xlsFile_output,OCobj.PatientCode(OCobj.PatientRows)',strcat(OCobj.xlsSheet,'_EUD'),'B1');
                end
                if ~isempty(OCobj.mEUD)
                    xlswrite(OCobj.xlsFile_output,OCobj.mEUD(OCobj.PatientRows,:)',strcat(OCobj.xlsSheet,'_EUD'),'B2');
                end
                
                % total patients at log10n and dose
                xlswrite(OCobj.xlsFile_output, OCobj.log10n, strcat(OCobj.xlsSheet,'_Total'),'A2');
                if ~isempty(OCobj.mAtlasTotal);
%                     doses = (0:size(OCobj.mAtlasTotal,1)-1)*OCobj.GyStep;
                    xlswrite(OCobj.xlsFile_output,OCobj.AtlasBins',strcat(OCobj.xlsSheet,'_Total'),'B1');
                    xlswrite(OCobj.xlsFile_output,OCobj.mAtlasTotal',strcat(OCobj.xlsSheet,'_Total'),'B2');
                end
                
                % complication patients at log10n and dose
                xlswrite(OCobj.xlsFile_output, OCobj.log10n, strcat(OCobj.xlsSheet,'_Comp'),'A2');
                if ~isempty(OCobj.mAtlasComp)
%                     doses = (0:size(OCobj.mAtlasComp,1)-1)*OCobj.GyStep;
                    xlswrite(OCobj.xlsFile_output,OCobj.AtlasBins',strcat(OCobj.xlsSheet,'_Comp'),'B1');
                    xlswrite(OCobj.xlsFile_output,OCobj.mAtlasComp',strcat(OCobj.xlsSheet,'_Comp'),'B2');
                end
                
                % Beta probability
                if ~isempty(OCobj.BetaCumulativeMat)
%                     doses = (0:size(OCobj.mAtlasComp,1)-1)*OCobj.GyStep;
                    for k=1:length(OCobj.mBetaCumulativeThreshold)
                        xlswrite(OCobj.xlsFile_output, OCobj.log10n, strcat(OCobj.xlsSheet,'_prob_',num2str(OCobj.mBetaCumulativeThreshold(k))),'A2');
                        xlswrite(OCobj.xlsFile_output,OCobj.AtlasBins',strcat(OCobj.xlsSheet,'_prob_',num2str(OCobj.mBetaCumulativeThreshold(k))),'B1');
                        xlswrite(OCobj.xlsFile_output,OCobj.BetaCumulativeMat(:,:,k)',strcat(OCobj.xlsSheet,'_prob_',num2str(OCobj.mBetaCumulativeThreshold(k))),'B2');
                    end
                end
                if ~isempty(OCobj.BetaInverseMat)
%                     doses = (0:size(OCobj.mAtlasComp,1)-1)*OCobj.GyStep;
                    for k=1:length(OCobj.mBetaInverseThreshold)
                        xlswrite(OCobj.xlsFile_output, OCobj.log10n, strcat(OCobj.xlsSheet,'_Low_',num2str(OCobj.mBetaInverseThreshold(k))),'A2');
                        xlswrite(OCobj.xlsFile_output,OCobj.AtlasBins',strcat(OCobj.xlsSheet,'_Low_',num2str(OCobj.mBetaInverseThreshold(k))),'B1');
                        xlswrite(OCobj.xlsFile_output,OCobj.BetaInverseMat(:,:,k)',strcat(OCobj.xlsSheet,'_Low_',num2str(OCobj.mBetaInverseThreshold(k))),'B2');
                    end
                end
            end
            warning('on','MATLAB:xlswrite:AddSheet');
        end
    end
end