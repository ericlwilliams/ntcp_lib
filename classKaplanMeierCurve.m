classdef classKaplanMeierCurve
    properties
        mSurvivalTime = cell(0); % original survival data
        mFlgCensor =cell(0); % censor information of the survival data
        
        mSurvivalTimeSorted
        mSurvivalCurve
        mSurvivalCurveVariance
        mCensorStatistics
        mSurvivalTimeTable
        mCumulativeHazard
        
        mSurvivalTimeCombined
        mSurvivalCurveCombined
        mSurvivalCurveVarianceCombined
        mCensorStatisticsCombined
        mSurvivalTimeTableCombined
        mCumulativeHazardCombined
        
        mFlgCompareWithVariations
        mChi2
        mpValue
        mCurveArea
        mHR
    end

    methods
        function SAobj = classKaplanMeierCurve(mSurvivalTime,mFlgCensor,mFlgCompareWithVariations)
            if exist('mSurvivalTime','var')
                SAobj.mSurvivalTime=mSurvivalTime(:);
                if exist('mFlgCensor','var')
                    SAobj.mFlgCensor=mFlgCensor(:);
                else
                    SAobj.mFlgCensor=cell(size(SAobj.mSurvivalTime));
                    for k=1:length(SAobj.mSurvivalTime)
                        SAobj.mFlgCensor{k}=false(size(SAobj.mSurvivalTime{k}));
                    end
                end
            end
            if exist('mFlgCompareWithVariations','var')
                SAobj.mFlgCompareWithVariations=mFlgCompareWithVariations;
            else
                SAobj.mFlgCompareWithVariations=true;
            end
        end
        function SAobj = set.mSurvivalTime(SAobj,survivaltime)
            if ~isa(survivaltime,'cell')
                error('The mSurvivalTime in classKaplanMeierCurve should be a cell array');
            end
            f = cellfun('isempty', survivaltime);
            if any(f)
                warning(['The mSurvivalTime is assigned a cell array with empty cell(s): ', num2str((find(f(:)))')]);
            end
%             if all(cellfun('isempty',survivaltime))
%                 error('The mSurvivalTime in classKaplanMeierCurve is assigned by empty cell(s)');
%             end
            SAobj.mSurvivalTime = survivaltime(:);
        end
        function SAobj = set.mFlgCensor(SAobj,flgcensor)
            if ~isa(flgcensor,'cell')
                error('The mSurvivalTime in classKaplanMeierCurve should be a cell array');
            end
            f = cellfun('isempty', flgcensor);
            if any(f)
                warning(['The mFlgCensor is assigned a cell array with empty cell(s): ', num2str((find(f(:)))')]);
            end
%             if all(cellfun('isempty',flgcensor))
%                 error('The mFlgCensor in classKaplanMeierCurve is assigned by empty cell(s)');
%             end
            SAobj.mFlgCensor = flgcensor(:);
        end
        
        function SAobj=fCalculateSurvivalCurve(SAobj)
            SAobj.mSurvivalTime=SAobj.mSurvivalTime(:);
            SAobj.mFlgCensor=SAobj.mFlgCensor(:);
            num=size(SAobj.mSurvivalTime,1);
            SAobj.mCumulativeHazard=cell(num,1);
            SAobj.mCensorStatistics=cell(num,1);
            SAobj.mSurvivalCurveVariance=cell(num,1);
            for k=1:num
                [SAobj.mSurvivalTimeSorted{k,1}, SAobj.mSurvivalCurve{k,1}, SAobj.mCensorStatistics{k,1}, ...
                    SAobj.mSurvivalTimeTable{k,1}, SAobj.mCumulativeHazard{k,1}, SAobj.mSurvivalCurveVariance{k,1}]...
                    =KaplanMeierSurvivalCurve(SAobj.mSurvivalTime{k},SAobj.mFlgCensor{k});
            end
        end
        
        function SAobj=fCombineSurvivalTime(SAobj)
            SAobj.mSurvivalTime=SAobj.mSurvivalTime(:);
            SAobj.mFlgCensor=SAobj.mFlgCensor(:);

            num=size(SAobj.mSurvivalTime,1);
            SAobj.mSurvivalTimeCombined=[]; FlagCensorCombined=logical([]);
            for k=1:num
                SAobj.mSurvivalTimeCombined=[SAobj.mSurvivalTimeCombined; SAobj.mSurvivalTime{k}];
                FlagCensorCombined=[FlagCensorCombined;SAobj.mFlgCensor{k}];
            end
            if isempty(SAobj.mSurvivalTimeCombined)
                disp('no survival data to compute');
                return;
            end
            [SAobj.mSurvivalTimeCombined, SAobj.mSurvivalCurveCombined, SAobj.mCensorStatisticsCombined, ...
                SAobj.mSurvivalTimeTableCombined, SAobj.mCumulativeHazardCombined, SAobj.mSurvivalCurveVarianceCombined]...
                =KaplanMeierSurvivalCurve(SAobj.mSurvivalTimeCombined, FlagCensorCombined);
        end
        
        function SAobj=fCompareSurvivalByLogrank(SAobj)
            if SAobj.mFlgCompareWithVariations && size(SAobj.mSurvivalTime,1)==2
                SAobj=fLogrankTestByExpectationAndVariation(SAobj);
            else
                SAobj=fLogrankTestByExpectation(SAobj);
            end
        end

        function SAobj=fLogrankTestByExpectation(SAobj)
            % % calculate the last time point for comparison
            %     mSurvivalTimeTable=SAobj.SurvivalTimeTible;
            %     mSurvivalTimeTableCombined=SAobj.mSurvivalTimeTableCombined;
            %     f=cellfun(@(x) x(end,1), mSurvivalTimeTable); f=min(f); % the shorted survival curve's last time point
            
            numgroups=size(SAobj.mSurvivalTimeTable,1); % number of groups for comparison
            numevent=size(SAobj.mSurvivalTimeTableCombined,1); % number of total events
            
            % computer es and os for chi2 computation
            tbn=zeros(numevent+1,numgroups); % table for ni, the (surplus) one row is for computation convenience
            tbe=zeros(numevent+1,numgroups); % table for events
            
            %     tbn(1,:)=cellfun(@(x) x(1,3),SAobj.mSurvivalTimeTable); % fill ni's of first row with patient at risk in each group
            %     tbn(2:end,end)=mSurvivalTimeTable{end}(:,1); % fill the last column with total ni in each row
            
            for k=1:numgroups % fill the row group by group
                [flg,loc]=ismember(SAobj.mSurvivalTimeTableCombined(:,1),SAobj.mSurvivalTimeTable{k}(:,1));
                for n=numevent:-1:1 % fill the ni and event(i) row by row (timepoint by timepoint)
                    if flg(n) % it occurs in the k-th group
                        tbn(n,k)=SAobj.mSurvivalTimeTable{k}(loc(n),3); % ni, (n+1) reflects the first surplus row in the table
                        tbe(n,k)=SAobj.mSurvivalTimeTable{k}(loc(n),2); % event number
                    else % no events in the k-th group, copy the patient at risk from the next row
                        tbn(n,k)=tbn(n+1,k); % no event on the timepoint, record the ni
                    end
                end
            end
            tbn(end,:)=[]; tbe(end,:)=[]; % the surplus last row is not used any more
            if any(sum(tbn,2)-SAobj.mSurvivalTimeTableCombined(:,3)) || any(sum(tbe,2)-SAobj.mSurvivalTimeTableCombined(:,2)) % check if the table is correctly filled
                error('table is not correct');
            end
            
            % search the proper stop time point for chi2 computation, which is the row where no more event happen after, or one group survivors run out, whichever comes first
            n=numevent;
            while length(find(tbn(n,:)))<numgroups && n>0 % for matlab, it is faster if the time point when one group runs out is the stop time point
                n=n-1;
            end
            
            % compute the chi2 and p-value
            dtbe=repmat(sum(tbe,2),[1,numgroups]); % denominator table of events
            tbn=dtbe.*tbn./repmat(SAobj.mSurvivalTimeTableCombined(:,3),[1,numgroups]);
            tbn=tbn(1:n,:); tbe=tbe(1:n,:); % only the first 1:n rows are applied in computation
            tbn=sum(tbn); % the expected events in each group
            tbe=sum(tbe); % the events in each group
            if abs(sum(tbn)-sum(tbe))>1e-6
                error('table is not correct');
            end
            SAobj.mChi2 = sum (((tbe-tbn).^2)./tbn);
            SAobj.mpValue = 1 - cdf('chi2',SAobj.mChi2,numgroups-1);
            
            % compute the survival order of groups (which group's survival curve is above which)
            mCurveArea=zeros(numgroups,1); endtime=SAobj.mSurvivalTimeTableCombined(n,1);
            for k=1:numgroups
                f=find(endtime>SAobj.mSurvivalTimeSorted{k}); % f(end) is the last time point for group k
                mCurveArea(k)=sum(SAobj.mSurvivalCurve{k}(1:f(end)).*diff([SAobj.mSurvivalTimeSorted{k}(1:f(end));endtime]));
            end
            SAobj.mCurveArea=mCurveArea;
        end
        
        function SAobj=fLogrankTestByExpectationAndVariation(SAobj)
            
            numgroups=size(SAobj.mSurvivalTimeTable,1); % number of groups for comparison
            numevent=size(SAobj.mSurvivalTimeTableCombined,1); % number of total events
            
            % computer es and os for chi2 computation
            tbn=zeros(numevent+1,numgroups); % table for ni, the (surplus) one row is for computation convenience
            tbe=zeros(numevent+1,numgroups); % table for events
            
            for k=1:numgroups % fill the row group by group
                [flg,loc]=ismember(SAobj.mSurvivalTimeTableCombined(:,1),SAobj.mSurvivalTimeTable{k}(:,1));
                for n=numevent:-1:1 % fill the ni and event(i) row by row (timepoint by timepoint)
                    if flg(n) % it occurs in the k-th group
                        tbn(n,k)=SAobj.mSurvivalTimeTable{k}(loc(n),3); % ni, (n+1) reflects the first surplus row in the table
                        tbe(n,k)=SAobj.mSurvivalTimeTable{k}(loc(n),2); % event number
                    else % no events in the k-th group, copy the patient at risk from the next row
                        tbn(n,k)=tbn(n+1,k); % no event on the timepoint, record the ni
                    end
                end
            end
            tbn(end,:)=[]; tbe(end,:)=[]; % the surplus last row is not used any more
            if any(sum(tbn,2)-SAobj.mSurvivalTimeTableCombined(:,3)) || any(sum(tbe,2)-SAobj.mSurvivalTimeTableCombined(:,2)) % check if the table is correctly filled
                error('table is not correct');
            end
            
            % search the proper stop time point for chi2 computation, which is the row where no more event happen after, or one group survivors run out, whichever comes first
            neffectiverows=numevent;
            while length(find(tbn(neffectiverows,:)))<numgroups && neffectiverows>0 % for matlab, it is faster if the time point when one group runs out is the stop time point
                neffectiverows=neffectiverows-1;
            end
            
            % compute the variations and covariations
            sttc=SAobj.mSurvivalTimeTableCombined(1:neffectiverows,:); % survival time table combined, for coding convenience
            tbn=tbn(1:neffectiverows,:); tbe=tbe(1:neffectiverows,:); % keep effective rows only for chi2 computation
            if numgroups==2
                covmat=0;
                for k=1:neffectiverows
                    covmat=covmat+ sttc(k,2)*tbn(k,1)/sttc(k,3)*(1-tbn(k,1)/sttc(k,3))*(sttc(k,3)-sttc(k,2))/(sttc(k,3)-1);
                end
            else
                covmat=zeros(numgroups,numgroups);
                for k=1:neffectiverows % compute through each time point
                    for m=1:numgroups % compute cov by searching each group vs each other group
                        covmn = (- sttc(k,2) * tbn(k,m) / sttc(k,3)^2 * (sttc(k,3)-sttc(k,2))/(sttc(k,3)-1)) * tbn(k,:); % cov of group m vs all other groups
                        covmn(m) = tbn(k,m)/sttc(k,3) * (1-tbn(k,m)/sttc(k,3)) * sttc(k,2) * (sttc(k,3)-sttc(k,2))/(sttc(k,3)-1); % cov of group m vs. group m is its variation, compute it specifically
                        covmat(m,:)=covmat(m,:)+covmn;
                    end
                end
            end
            
            % compute the chi2 and p-value
            dtbe=repmat(sum(tbe,2),[1,numgroups]); % duplicate of total events at each time point for each group
            tbexp=dtbe.*tbn./repmat(sttc(:,3),[1,numgroups]); % expected events at each time point for each group
            tbexp=sum(tbexp); % the total expected events in each group
            tbevents=sum(tbe); % the total actual events in each group
            if abs(sum(tbexp)-sum(tbevents))>1e-6
                error('table is not correct');
            end
            
            SAobj.mChi2=(abs(tbevents(1)-tbexp(1)))^2/sum(sttc(:,2).*tbn(:,1).*tbn(:,2)./sttc(:,3).^2);
            SAobj.mpValue=1-cdf('chi2',SAobj.mChi2,1);
            
            %     SAobj.mChi2=(abs(tbevents(1)-tbexp(1))-0.5)/sqrt(covmat(1,1));
            %     SAobj.mpValue=1-cdf('norm',SAobj.mChi2,0,1);
            
            %     SAobj.mChi2 = (tbevents-tbexp) * (covmat \ (tbevents-tbexp)');
            %     SAobj.mpValue = 1 - cdf('chi2',SAobj.mChi2,numgroups-1);
            
            % compute the survival order of groups (which group's survival curve is above which)
            mCurveArea=zeros(numgroups,1); endtime=SAobj.mSurvivalTimeTableCombined(neffectiverows,1);
            for k=1:numgroups
                f=find(endtime>SAobj.mSurvivalTimeSorted{k}); % f(end) is the last time point for group k
                mCurveArea(k)=sum(SAobj.mSurvivalCurve{k}(1:f(end)).*diff([SAobj.mSurvivalTimeSorted{k}(1:f(end));endtime]));
            end
            SAobj.mCurveArea=mCurveArea;
        end
    end
end


function [mSurvivalTimeSorted,mSurvivalCurve,mCensorStatistics,mSurvivalTimeTable,mCumulativeHazard,CurveVar]=KaplanMeierSurvivalCurve(mSurvivalTime,mFlgCensor)
% input: mSurvivalTime -- survival nx1 vector of survival time

% prepare
    if ~exist('mFlgCensor','var') % default is no censor
        mFlgCensor=false(length(mSurvivalTime));
    end
    mCensorStatistics=[];
    
% sort data by survival time, and generate histogram table
    [mSurvivalTimeTable,f]=sort(mSurvivalTime); mFlgCensor=mFlgCensor(f); % mSurvivalTimeTable - table of survive
    if any(mFlgCensor)
        CensorTimeTable=tabulate(mSurvivalTimeTable(mFlgCensor)); 
        CensorTimeTable(CensorTimeTable(:,2)==0,:)=[]; % CensorTimeTable - table of censors (the censor time)
   
        %CompTimeTable=tabulate(mSurvivalTimeTable(~mFlgCensor)); 
        %CompTimeTable(CompTimeTable(:,2)==0,:)=[]; % CensorTimeTable - table of censors (the censor time) 
    end
    mSurvivalTimeTable=tabulate(mSurvivalTimeTable); 
    mSurvivalTimeTable(mSurvivalTimeTable(:,2)==0,:)=[]; % survival table. for interger survival times not observed time added by tabulate should be removed
    
% computer the number at risk and the number of events, note to remove censored data from the number of events
    f=flipud(mSurvivalTimeTable(:,2)); 
    f=cumsum(f); 
    mSurvivalTimeTable(:,3)=flipud(f); % column 3 is the number at risk
    
    
% computer the number at risk and the number of events, note to remove censored data from the number of events
    f=flipud(mSurvivalTimeTable(:,2)); f=cumsum(f); mSurvivalTimeTable(:,3)=flipud(f); % column 3 is the number at risk
    if any(mFlgCensor)
        [f,g]=ismember(mSurvivalTimeTable(:,1),CensorTimeTable(:,1)); % search censored events
        g(g==0)=[]; % remove non-matching events
        mCensorStatistics=[find(f)+1,CensorTimeTable(g,2)]; % find(f)+1 - add the location one to offset the shift of survival curve
        mSurvivalTimeTable(f,2)=mSurvivalTimeTable(f,2)-CensorTimeTable(g,2); % remove cencored number from event number
    end

% compute the survival curve
    mSurvivalTimeSorted=[0;mSurvivalTimeTable(:,1)]; % the survival time points. 0 indicate the beginning of experiment
    mSurvivalCurve=[1; cumprod( 1 - mSurvivalTimeTable(:,2)./mSurvivalTimeTable(:,3) )];
    mFlgCensor=[false;mFlgCensor];
    mCumulativeHazard=[0;cumsum(mSurvivalTimeTable(:,2)./mSurvivalTimeTable(:,3))];
    % variation of KM curve
    CurveVar = mSurvivalTimeTable(:,2) ./ (mSurvivalTimeTable(:,3).*(mSurvivalTimeTable(:,3)-mSurvivalTimeTable(:,2)));
    CurveVar = cumsum(CurveVar);
    CurveVar = (mSurvivalCurve.^2) .* [0;CurveVar];
    CurveVar = sqrt(CurveVar);
end

