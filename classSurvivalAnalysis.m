classdef classSurvivalAnalysis
    properties
        SurvivalTime = cell(0); % original survival data
        flgCensor =cell(0); % censor information of the survival data
        
        SurvivalTimeSorted
        SurvivalCurve
        SurvivalCurveVariance
        CensorStatistics
        SurvivalTimeTable
        CumulativeHazard
        
        SurvivalTimeCombined
        SurvivalCurveCombined
        SurvivalCurveVarianceCombined
        CensorStatisticsCombined
        SurvivalTimeTableCombined
        CumulativeHazardCombined
        
        FlagCompareWithVariations
        Chi2
        pValue
        CurveArea
    end

    methods
        function SAobj = classSurvivalAnalysis(SurvivalTime,flgCensor,FlagCompareWithVariations)
            if exist('SurvivalTime','var')
                SAobj.SurvivalTime=SurvivalTime(:);
                if exist('flgCensor','var')
                    SAobj.flgCensor=flgCensor(:);
                else
                    SAobj.flgCensor=cell(size(SAobj.SurvivalTime));
                    for k=1:length(SAobj.SurvivalTime)
                        SAobj.flgCensor{k}=false(size(SAobj.SurvivalTime{k}));
                    end
                end
            end
            if exist('FlagCompareWithVariations','var')
                SAobj.FlagCompareWithVariations=FlagCompareWithVariations;
            else
                SAobj.FlagCompareWithVariations=true;
            end
        end
        function SAobj = set.SurvivalTime(SAobj,survivaltime)
            if ~isa(survivaltime,'cell')
                error('The SurvivalTime in classSurvivalAnalysis should be a cell array');
            end
%             if all(cellfun('isempty',survivaltime))
%                 error('The SurvivalTime in classSurvivalAnalysis is assigned by empty cell(s)');
%             end
            SAobj.SurvivalTime = survivaltime(:);
        end
        function SAobj = set.flgCensor(SAobj,flgcensor)
            if ~isa(flgcensor,'cell')
                error('The SurvivalTime in classSurvivalAnalysis should be a cell array');
            end
%             if all(cellfun('isempty',flgcensor))
%                 error('The flgCensor in classSurvivalAnalysis is assigned by empty cell(s)');
%             end
            SAobj.flgCensor = flgcensor(:);
        end
        
        function SAobj=CalculateSurvivalCurve(SAobj)
            SAobj.SurvivalTime=SAobj.SurvivalTime(:);
            SAobj.flgCensor=SAobj.flgCensor(:);
            num=size(SAobj.SurvivalTime,1);
            SAobj.CumulativeHazard=cell(num,1);
            SAobj.CensorStatistics=cell(num,1);
            SAobj.SurvivalCurveVariance=cell(num,1);
            for k=1:num
                [SAobj.SurvivalTimeSorted{k,1}, SAobj.SurvivalCurve{k,1}, SAobj.CensorStatistics{k,1}, ...
                    SAobj.SurvivalTimeTable{k,1}, SAobj.CumulativeHazard{k,1}, SAobj.SurvivalCurveVariance{k,1}]...
                    =KaplanMeierSurvivalCurve(SAobj.SurvivalTime{k},SAobj.flgCensor{k});
            end
        end
        
        function SAobj=CombineSurvivalTime(SAobj)
            SAobj.SurvivalTime=SAobj.SurvivalTime(:);
            SAobj.flgCensor=SAobj.flgCensor(:);

            num=size(SAobj.SurvivalTime,1);
            SAobj.SurvivalTimeCombined=[]; FlagCensorCombined=logical([]);
            for k=1:num
                SAobj.SurvivalTimeCombined=[SAobj.SurvivalTimeCombined; SAobj.SurvivalTime{k}];
                FlagCensorCombined=[FlagCensorCombined;SAobj.flgCensor{k}];
            end
            if isempty(SAobj.SurvivalTimeCombined)
                disp('no survival data to compute');
                return;
            end
            [SAobj.SurvivalTimeCombined, SAobj.SurvivalCurveCombined, SAobj.CensorStatisticsCombined, ...
                SAobj.SurvivalTimeTableCombined, SAobj.CumulativeHazardCombined, SAobj.SurvivalCurveVarianceCombined]...
                =KaplanMeierSurvivalCurve(SAobj.SurvivalTimeCombined, FlagCensorCombined);
        end
        
        function SAobj=CompareSurvivalByLogrank(SAobj)
            if SAobj.FlagCompareWithVariations && size(SAobj.SurvivalTime,1)==2
                SAobj=LogrankTestByExpectationAndVariation(SAobj);
            else
                SAobj=LogrankTestByExpectation(SAobj);
            end
        end

        function SAobj=LogrankTestByExpectation(SAobj)
            % % calculate the last time point for comparison
            %     SurvivalTimeTable=SAobj.SurvivalTimeTible;
            %     SurvivalTimeTableCombined=SAobj.SurvivalTimeTableCombined;
            %     f=cellfun(@(x) x(end,1), SurvivalTimeTable); f=min(f); % the shorted survival curve's last time point
            
            numgroups=size(SAobj.SurvivalTimeTable,1); % number of groups for comparison
            numevent=size(SAobj.SurvivalTimeTableCombined,1); % number of total events
            
            % computer es and os for chi2 computation
            tbn=zeros(numevent+1,numgroups); % table for ni, the (surplus) one row is for computation convenience
            tbe=zeros(numevent+1,numgroups); % table for events
            
            %     tbn(1,:)=cellfun(@(x) x(1,3),SAobj.SurvivalTimeTable); % fill ni's of first row with patient at risk in each group
            %     tbn(2:end,end)=SurvivalTimeTable{end}(:,1); % fill the last column with total ni in each row
            
            for k=1:numgroups % fill the row group by group
                [flg,loc]=ismember(SAobj.SurvivalTimeTableCombined(:,1),SAobj.SurvivalTimeTable{k}(:,1));
                for n=numevent:-1:1 % fill the ni and event(i) row by row (timepoint by timepoint)
                    if flg(n) % it occurs in the k-th group
                        tbn(n,k)=SAobj.SurvivalTimeTable{k}(loc(n),3); % ni, (n+1) reflects the first surplus row in the table
                        tbe(n,k)=SAobj.SurvivalTimeTable{k}(loc(n),2); % event number
                    else % no events in the k-th group, copy the patient at risk from the next row
                        tbn(n,k)=tbn(n+1,k); % no event on the timepoint, record the ni
                    end
                end
            end
            tbn(end,:)=[]; tbe(end,:)=[]; % the surplus last row is not used any more
            if any(sum(tbn,2)-SAobj.SurvivalTimeTableCombined(:,3)) || any(sum(tbe,2)-SAobj.SurvivalTimeTableCombined(:,2)) % check if the table is correctly filled
                error('table is not correct');
            end
            
            % search the proper stop time point for chi2 computation, which is the row where no more event happen after, or one group survivors run out, whichever comes first
            n=numevent;
            while length(find(tbn(n,:)))<numgroups && n>0 % for matlab, it is faster if the time point when one group runs out is the stop time point
                n=n-1;
            end
            
            % compute the chi2 and p-value
            dtbe=repmat(sum(tbe,2),[1,numgroups]); % denominator table of events
            tbn=dtbe.*tbn./repmat(SAobj.SurvivalTimeTableCombined(:,3),[1,numgroups]);
            tbn=tbn(1:n,:); tbe=tbe(1:n,:); % only the first 1:n rows are applied in computation
            tbn=sum(tbn); % the expected events in each group
            tbe=sum(tbe); % the events in each group
            if abs(sum(tbn)-sum(tbe))>1e-6
                error('table is not correct');
            end
            SAobj.Chi2 = sum (((tbe-tbn).^2)./tbn);
            SAobj.pValue = 1 - cdf('chi2',SAobj.Chi2,numgroups-1);
            
            % compute the survival order of groups (which group's survival curve is above which)
            CurveArea=zeros(numgroups,1); endtime=SAobj.SurvivalTimeTableCombined(n,1);
            for k=1:numgroups
                f=find(endtime>SAobj.SurvivalTimeSorted{k}); % f(end) is the last time point for group k
                CurveArea(k)=sum(SAobj.SurvivalCurve{k}(1:f(end)).*diff([SAobj.SurvivalTimeSorted{k}(1:f(end));endtime]));
            end
            SAobj.CurveArea=CurveArea;
        end
        
        function SAobj=LogrankTestByExpectationAndVariation(SAobj)
            
            numgroups=size(SAobj.SurvivalTimeTable,1); % number of groups for comparison
            numevent=size(SAobj.SurvivalTimeTableCombined,1); % number of total events
            
            % computer es and os for chi2 computation
            tbn=zeros(numevent+1,numgroups); % table for ni, the (surplus) one row is for computation convenience
            tbe=zeros(numevent+1,numgroups); % table for events
            
            for k=1:numgroups % fill the row group by group
                [flg,loc]=ismember(SAobj.SurvivalTimeTableCombined(:,1),SAobj.SurvivalTimeTable{k}(:,1));
                for n=numevent:-1:1 % fill the ni and event(i) row by row (timepoint by timepoint)
                    if flg(n) % it occurs in the k-th group
                        tbn(n,k)=SAobj.SurvivalTimeTable{k}(loc(n),3); % ni, (n+1) reflects the first surplus row in the table
                        tbe(n,k)=SAobj.SurvivalTimeTable{k}(loc(n),2); % event number
                    else % no events in the k-th group, copy the patient at risk from the next row
                        tbn(n,k)=tbn(n+1,k); % no event on the timepoint, record the ni
                    end
                end
            end
            tbn(end,:)=[]; tbe(end,:)=[]; % the surplus last row is not used any more
            if any(sum(tbn,2)-SAobj.SurvivalTimeTableCombined(:,3)) || any(sum(tbe,2)-SAobj.SurvivalTimeTableCombined(:,2)) % check if the table is correctly filled
                error('table is not correct');
            end
            
            % search the proper stop time point for chi2 computation, which is the row where no more event happen after, or one group survivors run out, whichever comes first
            neffectiverows=numevent;
            while length(find(tbn(neffectiverows,:)))<numgroups && neffectiverows>0 % for matlab, it is faster if the time point when one group runs out is the stop time point
                neffectiverows=neffectiverows-1;
            end
            
            % compute the variations and covariations
            sttc=SAobj.SurvivalTimeTableCombined(1:neffectiverows,:); % survival time table combined, for coding convenience
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
            
            SAobj.Chi2=(abs(tbevents(1)-tbexp(1)))^2/sum(sttc(:,2).*tbn(:,1).*tbn(:,2)./sttc(:,3).^2);
            SAobj.pValue=1-cdf('chi2',SAobj.Chi2,1);
            
            %     SAobj.Chi2=(abs(tbevents(1)-tbexp(1))-0.5)/sqrt(covmat(1,1));
            %     SAobj.pValue=1-cdf('norm',SAobj.Chi2,0,1);
            
            %     SAobj.Chi2 = (tbevents-tbexp) * (covmat \ (tbevents-tbexp)');
            %     SAobj.pValue = 1 - cdf('chi2',SAobj.Chi2,numgroups-1);
            
            % compute the survival order of groups (which group's survival curve is above which)
            CurveArea=zeros(numgroups,1); endtime=SAobj.SurvivalTimeTableCombined(neffectiverows,1);
            for k=1:numgroups
                f=find(endtime>SAobj.SurvivalTimeSorted{k}); % f(end) is the last time point for group k
                CurveArea(k)=sum(SAobj.SurvivalCurve{k}(1:f(end)).*diff([SAobj.SurvivalTimeSorted{k}(1:f(end));endtime]));
            end
            SAobj.CurveArea=CurveArea;
        end
    end
end


function [SurvivalTimeSorted,SurvivalCurve,CensorStatistics,SurvivalTimeTable,CumulativeHazard,CurveVar]=KaplanMeierSurvivalCurve(SurvivalTime,flgCensor)
% input: SurvivalTime -- survival nx1 vector of survival time

% prepare
    if ~exist('flgCensor','var') % default is no censor
        flgCensor=false(length(SurvivalTime));
    end
    CensorStatistics=[];
    
% sort data by survival time, and generate histogram table
    [SurvivalTimeTable,f]=sort(SurvivalTime); flgCensor=flgCensor(f); % SurvivalTimeTable - table of survive
    if any(flgCensor)
        CensorTimeTable=tabulate(SurvivalTimeTable(flgCensor)); CensorTimeTable(CensorTimeTable(:,2)==0,:)=[]; % CensorTimeTable - table of censors (the censor time)
    end
    SurvivalTimeTable=tabulate(SurvivalTimeTable); SurvivalTimeTable(SurvivalTimeTable(:,2)==0,:)=[]; % survival table. for interger survival times not observed time added by tabulate should be removed
    
% computer the number at risk and the number of events, note to remove censored data from the number of events
    f=flipud(SurvivalTimeTable(:,2)); f=cumsum(f); SurvivalTimeTable(:,3)=flipud(f); % column 3 is the number at risk
    if any(flgCensor)
        [f,g]=ismember(SurvivalTimeTable(:,1),CensorTimeTable(:,1)); % search censored events
        g(g==0)=[]; % remove non-matching events
        CensorStatistics=[find(f)+1,CensorTimeTable(g,2)]; % find(f)+1 - add the location one to offset the shift of survival curve
        SurvivalTimeTable(f,2)=SurvivalTimeTable(f,2)-CensorTimeTable(g,2); % remove cencored number from event number
    end

% compute the survival curve
    SurvivalTimeSorted=[0;SurvivalTimeTable(:,1)]; % the survival time points. 0 indicate the beginning of experiment
    SurvivalCurve=[1; cumprod( 1 - SurvivalTimeTable(:,2)./SurvivalTimeTable(:,3) )];
    flgCensor=[false;flgCensor];
    CumulativeHazard=[0;cumsum(SurvivalTimeTable(:,2)./SurvivalTimeTable(:,3))];
    % variation of KM curve
    CurveVar = SurvivalTimeTable(:,2) ./ (SurvivalTimeTable(:,3).*(SurvivalTimeTable(:,3)-SurvivalTimeTable(:,2)));
    CurveVar = cumsum(CurveVar);
    CurveVar = (SurvivalCurve.^2) .* [0;CurveVar];
    CurveVar = sqrt(CurveVar);
end

