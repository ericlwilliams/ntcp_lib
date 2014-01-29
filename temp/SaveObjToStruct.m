function strct=ObjToStruct(obj)
    obj1=obj(:); objnum=size(obj1,1);
    props = properties(obj(1));
    propsnum=numel(props);
    flgobj=false(propsnum,1); % record the object members
    % first object
        k=1;
        for p = 1:propsnum
            if isobject(obj1(k).(props{p}))
                s = SaveObjToStruct(obj1(k).(props{p}));
                strct.(props{p}) = s;
                flgobj(p)=true;
            else
                strct.(props{p}) = obj1(k).(props{p});
            end
        end
    % rest objects
        strct = repmat(strct,[size(obj)]);
        % object properties
        f = find(flgobj); fl = length(f);
        for k = 2:objnum
            for p = 1:fl
                s = SaveObjToStruct(obj1(k).(props{f(p)}));
                strct(k).(props{f(p)}) = s;
            end
        end
        % non object properties
        f = find(~flgobj); fl = length(f);
        for k = 2:objnum
            for p = 1:fl
                strct(k).(props{f(p)}) = obj1(k).(props{f(p)});
            end
        end
end