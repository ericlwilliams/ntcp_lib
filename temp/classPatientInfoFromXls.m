classdef classPatientInfoFromXls
    properties
        xlsRaw % raw data of xls sheet
        
        PatientColName % patient code column name
        PatientID
        PatientRows % flags of which row has patient code

        ComplicationStartColName % column name of the start point of computation of complicaiton
        ComplicationStartDate % the date of the start point (in number)
        
        ComplicationColName % patient complication code column
        ComplicationGrade % patient complication grade
        ComplicationThreshold % threshold of severe complications
        ComplicationDateColName % Grade x occurred date
        ComplicationDays % times before patient complication occured from start point, usually in months
        ComplicationRows % flags of complication grade
        ComplicationFlag % 1 - the complication is severe (grade >= ComplicationThreshold)
        
        FxNumColName % name of the fraction column
        FxNum; % number of fractions
        FxNumRows % 1 - has fraction number, 0 - no fraction number
    end

    methods
        function PTobj=PatientCodeInXlsraw(PTobj) % extract the column according to column name
            [colcontent,colflg]=ColumnInXlsraw(PTobj.xlsRaw,PTobj.PatientColName); % data in the column
            PTobj.PatientID=colcontent;
            f=cellfun(@(x) any(isnan(x)), colcontent); % nan cells shall be excluded
            colflg(f)=false; PTobj.PatientRows=colflg; % flags for the data in the column
        end
        function PTobj=FractionNumbers(PTobj)
            [colcontent,colflg]=ColumnInXlsraw(PTobj.xlsRaw,PTobj.FxNumColName); % data in the column
            PTobj.FxNum=colcontent;
            f=cellfun(@(x) any(isnan(x)), colcontent); % nan cells shall be excluded
            colflg(f)=false; PTobj.FxNumRows=colflg; % flags for the data in the column
        end
        function PTobj=ComplicationInXlsraw(PTobj)
            [colcontent,colflg]=ColumnInXlsraw(PTobj.xlsRaw,PTobj.ComplicationColName); % data in the column
            f=cellfun(@(x) any(ischar(x))||any(isnan(x)), colcontent); colflg(f)=false; % nan and character cells shall be excluded
            colflg=colflg&PTobj.PatientRows; % exclude non-patient data
            PTobj.ComplicationRows=colflg; % flags for the data in the column
            PTobj.ComplicationGrade=zeros(size(colcontent)); PTobj.ComplicationGrade(colflg)=cell2mat(colcontent(colflg));
            PTobj.ComplicationFlag=false(size(colcontent)); PTobj.ComplicationFlag(colflg)=PTobj.ComplicationGrade(colflg)>=PTobj.ComplicationThreshold;
        end
        function PTobj=ComplicationStartPoint(PTobj)
            % content of the column
            [colcontent,colflg]=ColumnInXlsraw(PTobj.xlsRaw,PTobj.ComplicationStartColName); % data in the column
            f=cellfun(@(x) any(ischar(x))||any(isnan(x)), colcontent); colflg(f)=false; % nan and character cells shall be excluded
            colflg=colflg&PTobj.PatientRows; % exclude non-patient data
            PTobj.ComplicationRows=colflg; % flags for the data in the column
            PTobj.ComplicationStartDate=zeros(size(colcontent)); PTobj.ComplicationStartDate(colflg)=cellfun(@(x) datenum(x,'mm/dd/yyyy'),colcontent(colflg));
        end
        function PTobj=ComplicationDaysInXlsraw(PTobj)
            % complication flag
            [colcontent,colflg]=ColumnInXlsraw(PTobj.xlsRaw,PTobj.ComplicationDateColName); % data in the column
            f=cellfun(@(x) any(ischar(x))||any(isnan(x)), colcontent); colflg(f)=false; % nan and character cells shall be excluded
            colflg=colflg&PTobj.PatientRows; % exclude non-patient data
            PTobj.ComplicationRows=colflg; % flags for the data in the column
            PTobj.ComplicationDays=zeros(size(colcontent)); PTobj.ComplicationDays(colflg)=cellfun(@(x) datenum(x,'mm/dd/yyyy'),colcontent(colflg));
            PTobj.ComplicationDays(colflg)=PTobj.ComplicationDays(colflg)-PTobj.ComplicationStartDate(colflg);
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