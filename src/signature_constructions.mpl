bound_info := proc(x, bound, eps)
$ifdef LOG_TIME
    INIT_START_LOG_TIME("bound_info",0)
$endif
local i1, i2, j1, j2;
    # This is a bounded bounduality
    if(nops(bound) = 2) then
        i1 := simplify(op(bound[1])[1]);
        i2 := simplify(op(bound[1])[2]);
        j1 := simplify(op(bound[2])[1]);
        j2 := simplify(op(bound[2])[2]);
        if evalb(i1 = x) then
            if evalb(j1 = x) then
$ifdef LOG_TIME
                END_LOG_TIME("bound_info",0)
$endif
                return [min(i2, j2)+eps, max(i2, j2)-eps];
            else
$ifdef LOG_TIME
                END_LOG_TIME("bound_info",0)
$endif
                return [min(i2, j1)+eps, max(i2, j1)-eps];
            end if;
        else
            if evalb(j1 = x) then
$ifdef LOG_TIME
                END_LOG_TIME("bound_info",0)
$endif
                return [min(i1, j2)+eps, max(i1, j2)-eps];
            else
$ifdef LOG_TIME
                END_LOG_TIME("bound_info",0)
$endif
                return [min(i1, j1)+eps, max(i1, j1)-eps];
            end if;
        end if;
        # This is an equality or unbounded bounduality
    else
        i1 := simplify(op(bound[1])[1]);
        j1 := simplify(op(bound[1])[2]);
        if type(bound[1], `=`) then
            if evalb(i1 = x) and evalb(j1 = x) then
$ifdef LOG_TIME
                END_LOG_TIME("bound_info",0)
$endif
                return [-infinity, infinity];
            end if;
            if evalb(i1 = x) then
$ifdef LOG_TIME
                END_LOG_TIME("bound_info",0)
$endif
                return [j1, j1];
            else
$ifdef LOG_TIME
                END_LOG_TIME("bound_info",0)
$endif
                return [i1, i1];
            end if;
        end if;
        if type(bound[1], `<=`) then
            if evalb(i1 = x) then
$ifdef LOG_TIME
                END_LOG_TIME("bound_info",0)
$endif
                return [-infinity, j1-eps];
            else
$ifdef LOG_TIME
                END_LOG_TIME("bound_info",0)
$endif
                return [i1+eps, infinity];
            end if;
        end if;
    end if;
end proc;
