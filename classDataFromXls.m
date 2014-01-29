classdef classDataFromXls
    properties
        mXlsRaw % raw data of xls sheet
        
        mLabel % column name
        mData % data in the column/row
        mFlg % flags of which column/row has the data
    end

    methods
        function Xobj = classDataFromXls()
        end
        function Xobj = fExtractColData(Xobj) % extract the data in the column according to column name
            % find the location of the column
            f = cellfun(@(x) strcmpi(x,Xobj.mLabel), Xobj.mXlsRaw);
            [m,n] = find(f); % location of the column
            if length(m)>1
                for k = 2:length(m)
                    if ~isequal(Xobj.mXlsRaw(:,n(1)),Xobj.mXlsRaw(:,n(k)))
                        f = cellfun(@(x,y) isequal(x,y), Xobj.mXlsRaw(:,n(1)), Xobj.mXlsRaw(:,n(k)));
                        warning(['There are ',num2str(length(m)),' cells containing the same column name "',Xobj.mLabel,'", but the data are not consistant, pick the first column.']);
                        disp(['Inconsistence compared with the first column: column ', num2str(n(k))]);
                        disp(['Inconsistence comparted with the first column: row ', num2str(find(~f'))]);
                    end
                end
            elseif isempty(m)
                error(['Column name "',Xobj.mLabel,'" not found, can not continue.']);
            end
            % extract the column data
            Xobj.mData = Xobj.mXlsRaw(:,n(1));
            Xobj.mFlg = true(size(Xobj.mXlsRaw,1),1);
            Xobj.mFlg(1:m(1)) = false; % mark cells preceeding data to be removed
            % exculde nan and empty cells
            f = cellfun(@(x) any(isnan(x)), Xobj.mData);
            f = f | cellfun('isempty', Xobj.mData);
            Xobj.mFlg(f) = false;
        end
        function Xobj = fExtractRowData(Xobj)
            % find the location of the row
            f = cellfun(@(x) strcmpi(x,Xobj.mLabel), Xobj.mXlsRaw);
            [m,n] = find(f); % location of the column
            if length(m)>1
                for k = 2:length(m)
                    if ~isequal(Xobj.mXlsRaw(m(1),:), Xobj.mXlsRaw(m(k),:))
                        f = cellfun(@(x,y) isequal(x,y), Xobj.mXlsRaw(m(1),:), Xobj.mXlsRaw(m(k),:));
                        warning(['There are ',num2str(length(m)),' cells containing the same row name "',Xobj.mLabel,'", but the data are not consistant, pick the first row.']);
                        disp(['Inconsistence compared with the first row: row ', num2str(n(k))]);
                        disp(['Inconsistence comparted with the first row: column ', num2str(find(~f))]);
                    end
                end
            elseif isempty(m)
                error(['Row name "',Xobj.mLabel,'" not found, can not continue.']);
            end
            % extract the row data
            Xobj.mData = Xobj.mXlsRaw(m(1),:)';
            Xobj.mFlg = true(size(Xobj.mXlsRaw,2),1);
            Xobj.mFlg(1:n(1),:) = false;
            % exculde nan cells
            f = cellfun(@(x) any(isnan(x)), Xobj.mData);
            Xobj.mFlg(f) = false;
        end
    end
end