classdef classComplicationInfoFromXls
    properties
        xlsRaw % raw data of xls sheet
        
        PatientCol % patient code column name
        PatientID
        PatientRows % flags of which row has patient code
        
        FxNumCol % name of the fraction column
        FxNum; % number of fractions
        FxNumRows % 1 - has fraction number, 0 - no fraction number
        
        BaselineCol % column name of the start point of computation of complicaiton
        BaselineDate % the date of the start point (in number)
        BaselineRows % the rows with baseline date
        LastFollowupCol % the column of last follow ups
        LastFollowupDate % time from start point to last follow ups
        LastFollowupRows % rows with last followup date

        ComplicationCol % patient complication code column
        ComplicationGrade % patient complication grade
        ComplicationThreshold % threshold of severe complications
        
        ComplicationDateCol % Grade x occurred date
        ComplicationDate % times before patient complication occured from start point, usually in months
        
        ComplicationRows % flags of complication grade
        CensorFlag % 1 - censored
    end

    methods
        function PTobj=classComplicationInfoFromXls()
        end
        function PTobj=PatientCodeInXlsraw(PTobj) % extract the column according to column name
            [colcontent,colflg]=ColumnInXlsraw(PTobj.xlsRaw,PTobj.PatientCol); % data in the column
            PTobj.PatientID=colcontent;
            f=cellfun(@(x) any(isnan(x)), colcontent); % nan cells shall be excluded
            colflg(f)=false; PTobj.PatientRows=colflg; % flags for the data in the column
            % check integrity
            f=unique(colcontent(colflg));
            if length(find(colflg))~=length(f)
                error('IDs in the spread sheet are not unique');
            end
        end
        function PTobj=FractionNumbersInXlsraw(PTobj)
            [colcontent,colflg]=ColumnInXlsraw(PTobj.xlsRaw,PTobj.FxNumCol); % data in the column
            f=cellfun(@(x) any(ischar(x))||any(isnan(x)), colcontent); colflg(f)=false; % nan cells shall be excluded
            PTobj.FxNumRows=colflg; % flags for the data in the column
            PTobj.FxNum=zeros(size(colcontent)); PTobj.FxNum(colflg)=cell2mat(colcontent(colflg));
            if any(PTobj.FxNum(colflg)==0)
                disp('There are fraction numbers equal to zero');
            end
        end
        function PTobj=ComplicationInXlsraw(PTobj)
            [colcontent,colflg]=ColumnInXlsraw(PTobj.xlsRaw,PTobj.ComplicationCol); % data in the column
            f=cellfun(@(x) any(ischar(x))||any(isnan(x)), colcontent); colflg(f)=false; % nan and character cells shall be excluded
            PTobj.ComplicationRows=colflg; % flags for the data in the column
            PTobj.ComplicationGrade=zeros(size(colcontent)); PTobj.ComplicationGrade(colflg)=cell2mat(colcontent(colflg));
            PTobj.CensorFlag=true(size(colcontent)); PTobj.CensorFlag(colflg)=PTobj.ComplicationGrade(colflg)<PTobj.ComplicationThreshold;
        end
        function PTobj=BaselineInXlsraw(PTobj)
            % content of the column
            [colcontent,colflg]=ColumnInXlsraw(PTobj.xlsRaw,PTobj.BaselineCol); % data in the column
            f=cellfun(@(x) any(isnan(x)), colcontent); colflg(f)=false; % nan and character cells shall be excluded
            PTobj.BaselineRows=colflg; % flags for the data in the column
            PTobj.BaselineDate=-inf(size(colcontent)); PTobj.BaselineDate(colflg)=cellfun(@(x) datenum(x,'mm/dd/yyyy'),colcontent(colflg));
            if any(isinf(PTobj.BaselineDate(colflg)))
                disp('There are missing baseline dates');
            end
        end
        function PTobj=ComplicationDaysInXlsraw(PTobj)
            % complication flag
            [colcontent,colflg]=ColumnInXlsraw(PTobj.xlsRaw,PTobj.ComplicationDateCol); % data in the column
            f=cellfun(@(x) any(isnan(x)), colcontent); colflg(f)=false; % nan and character cells shall be excluded
            PTobj.ComplicationRows=colflg; % flags for the data in the column
            PTobj.ComplicationDate=inf(size(colcontent)); PTobj.ComplicationDate(colflg)=cellfun(@(x) datenum(x,'mm/dd/yyyy'),colcontent(colflg));
            if any(isinf(PTobj.ComplicationDate(colflg)))
                disp('There are missing complication dates');
            end
%             PTobj.ComplicationDays(colflg)=PTobj.ComplicationDays(colflg)-PTobj.BaselineDate(colflg);
        end
        function PTobj=LastFollowupDaysInXlsraw(PTobj)
            [colcontent,colflg]=ColumnInXlsraw(PTobj.xlsRaw,PTobj.LastFollowupCol); % data in the column
            f=cellfun(@(x) any(isnan(x)), colcontent); colflg(f)=false; % nan and character cells shall be excluded
            PTobj.LastFollowupRows=colflg;
            PTobj.LastFollowupDate=inf(size(colcontent)); PTobj.LastFollowupDate(colflg)=cellfun(@(x) datenum(x,'mm/dd/yyyy'),colcontent(colflg));
            if any(isinf(PTobj.LastFollowupDate(colflg)))
                disp('There are missing last followup dates');
            end
%             PTobj.LastFollowupDays(colflg)=PTobj.LastFollowupDays(colflg)-PTobj.BaselineDate(colflg);
        end
        function PTobj=SelectPatientsInXlsraw(PTobj,idx)
            PTobj.xlsRaw=PTobj.xlsRaw(idx,:); % raw data of xls sheet
            PTobj.PatientID=PTobj.PatientID(idx,:);
            PTobj.FxNum=PTobj.FxNum(idx,:); % number of fractions
            try
                PTobj.BaselineDate=PTobj.BaselineDate(idx,:); % the date of the start point (in number)
            catch
            end
            PTobj.LastFollowupDate=PTobj.LastFollowupDate(idx,:); % time from start point to last follow ups
            try
                PTobj.ComplicationGrade=PTobj.ComplicationGrade(idx,:); % patient complication grade
            catch
            end
            PTobj.ComplicationDate=PTobj.ComplicationDate(idx,:); % times before patient complication occured from start point, usually in months
            PTobj.ComplicationRows=PTobj.ComplicationRows(idx,:); % complication cases
            PTobj.CensorFlag=PTobj.CensorFlag(idx,:); % 1 - censored
        end
    end
end


function [colcontent,colflg]=ColumnInXlsraw(xlsraw,colname)
% extract the column starting from the row under column name
    f=cellfun(@(x) strcmpi(x,colname), xlsraw);
    [m,n]=find(f); % location of the column
    if length(m)>1
        warning('There are more than one cells containing the same column name, pick the first one.');
    elseif isempty(m)
        error('no column name found, can not continue.');
    end
    colcontent = xlsraw(:,n); colflg=true(size(xlsraw,1),1); colflg(1:m(1))=false;
end