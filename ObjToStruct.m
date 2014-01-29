function strct=ObjToStruct(obj)
    obj1=obj(:); objnum=size(obj1,1);
    props = properties(obj); propsnum=numel(props);
    
    flgobj=false(propsnum,1); % record the object members
    
    strct.ClassName = class(obj); % keep the class name for restoration
    % check if the object is an empty object
    if objnum == 0 % empty obj, construct an empty object
        obj1 = eval(strct.ClassName); % set up an object of the class
    end
    
    % first object
        k=1;
        for p = 1:propsnum
            if isobject(obj1(k).(props{p}))
                s = ObjToStruct(obj1(k).(props{p}));
                strct.(props{p}) = s;
                flgobj(p)=true;
            else
                strct.(props{p}) = obj1(k).(props{p});
            end
        end
    % rest objects
        strct = repmat(strct,size(obj));
        % object properties
        f = find(flgobj); fl = length(f);
        for k = 2:objnum
            for p = 1:fl
                s = ObjToStruct(obj1(k).(props{f(p)}));
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