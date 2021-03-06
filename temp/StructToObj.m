function obj=StructToObj(strct)
    if ~isfield(strct,'ClassName'); % not an object, return with an empty variable
        obj=[];
        return;
    end
    
    strct1 = strct(:); strctnum=size(strct1,1);
    props = fieldnames(strct1(1)); % properties of an object
    f = cellfun(@(x) strcmp(x,'ClassName'), props); props(f)=[]; % remove the field "ClassName" which is not a property of an object
    
    % first object
        k=1;
        obj = eval(strct1(k).ClassName); % set up an object of the class
        objprops = properties(obj); % properties of the object
        f = ismember(props,objprops); props = props(f); % remove properties that the object does not have
        propsnum=numel(props);
        flgobj = false(propsnum,1); % flag for the object members
        for p = 1:propsnum
            if isstruct(strct1(k).(props{p}))
                o = StructToObj(strct1(k).(props{p}));
                if isempty(o) % not an object, just a structure
                    obj.(props{p})=strct1(k).(props{p});
                else
                    obj.(props{p})=o;
                    flgobj(p)=true;
                end
            else
                obj.(props{p}) = strct1(k).(props{p});
            end
        end
        
    % rest objects
        obj = repmat(obj,size(strct));
        % objects
        f = find(flgobj); fl=length(f);
        for k=2:strctnum
            for p=1:fl
                o = StructToObj(strct1(k).(props{f(p)}));
                obj(k).(props{f(p)}) = o;
            end
        end
        % non object properties
        f = find(~flgobj); fl=length(f);
        for k = 2:strctnum
            for p = 1:fl
                obj(k).(props{f(p)}) = strct1(k).(props{f(p)});
            end
        end
end