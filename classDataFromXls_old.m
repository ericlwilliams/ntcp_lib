classdef classDataFromXls_old
    properties
        xlsRaw % raw data of xls sheet
        
        ColName % column name
        ColData % data in the column
        flgDataRows % flags of which row has the data
        
        RowName % row name
        RowData % data in the row
        flgDataCols % flags for which column has the data
    end

    methods
        function Xobj = classDataFromXls()
        end
        function Xobj = ExtractColData(Xobj) % extract the data in the column according to column name
            % find the location of the column
            f = cellfun(@(x) strcmpi(x,Xobj.ColName), Xobj.xlsRaw);
            [m,n] = find(f); % location of the column
            if length(m)>1
                warning(['There are ',num2str(length(m)),' cells containing the same column name "',Xobj.ColName,'", pick the first one.']);
            elseif isempty(m)
                error(['Column name "',Xobj.ColName,'" not found, can not continue.']);
            end
            % extract the column data
            Xobj.ColData = Xobj.xlsRaw(:,n(1));
            Xobj.flgDataRows = true(size(Xobj.xlsRaw,1),1);
            Xobj.flgDataRows(1:m(1)) = false;
            % exculde nan and empty cells
            f = cellfun(@(x) any(isnan(x)), Xobj.ColData);
            f = f | cellfun('isempty', Xobj.ColData);
            Xobj.flgDataRows(f) = false;
        end
        function Xobj = ExtractRowData(Xobj)
            % find the location of the row
            f = cellfun(@(x) strcmpi(x,Xobj.RowName), Xobj.xlsRaw);
            [m,n] = find(f); % location of the column
            if length(m)>1
                warning(['There are ',num2str(length(m)),' cells containing the same row name "',Xobj.RowName,'", pick the first one.']);
            elseif isempty(m)
                error(['Row name "',Xobj.RowName,'" not found, can not continue.']);
            end
            % extract the row data
            Xobj.RowData = Xobj.xlsRaw(m(1),:)';
            Xobj.flgDataCols = true(size(Xobj.xlsRaw,2),1);
            Xobj.flgDataCols(1:n(1),:) = false;
            % exculde nan cells
            f = cellfun(@(x) any(isnan(x)), Xobj.RowData);
            Xobj.flgDataCols(f) = false;
        end
    end
end