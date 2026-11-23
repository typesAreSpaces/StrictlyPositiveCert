# Shared cooperative deadline for the public certificate pipeline.
local univariate_search_remaining := proc(deadline, stage)
local remaining := deadline-time[real]();
    if remaining<=0 then
        error "certificate search exhausted max_search_time during %1",stage;
    end if;
    return remaining;
end proc;

# A direction is sum(weights[i]*gens[i]); weights are fixed rational squares.
# A nonpositive f cannot be repaired where P>=0, including zeros of P.
local univariate_scale_regions := proc(f, P, cells, x)
local chart := univariate_sign_chart(P,x), nonnegative_cells;
    nonnegative_cells := univariate_intersect_cells(cells,univariate_select_cells(chart,"nonnegative"));
    if not univariate_positive_on_cells(f,nonnegative_cells,x) then return FAIL; end if;
    return [univariate_intersect_cells(cells,univariate_select_cells(chart,"positive")),
            univariate_intersect_cells(cells,univariate_select_cells(chart,"negative"))];
end proc;

local univariate_scale_status := proc(f, P, lambda, cells, regions, x)
local h := expand(f-lambda*P), big, small;
    if univariate_positive_on_cells(h,cells,x) then return "ok"; end if;
    big := not univariate_positive_on_cells(h,regions[1],x);
    small := not univariate_positive_on_cells(h,regions[2],x);
    if big and small then return "infeasible";
    elif big then return "big";
    elif small then return "small";
    else return "infeasible";
    end if;
end proc;

# Return [positive rational scale or FAIL, number of trials used].
local univariate_scale_direction := proc(f, P0, cells, x, max_trials, max_updates, deadline)
local P, regions, lambda, lo := 0, hi := 0, status, trial;
    if max_trials=0 then return [FAIL,0]; end if;
    P := timelimit(univariate_search_remaining(deadline,"scalar direction"),expand(P0));
    regions := timelimit(univariate_search_remaining(deadline,"scalar sign regions"),
                         univariate_scale_regions(f,P,cells,x));
    if regions=FAIL then return [FAIL,0]; end if;
    if P=0 then return [1,0]; end if;
    # Scale the initial trial by coefficient size; tiny constant directions
    # should not need dozens of doublings merely to reach a useful scale.
    lambda := timelimit(univariate_search_remaining(deadline,"scalar normalization"),
                        max(1,op(map(abs,[coeffs(expand(f),x)])))/max(op(map(abs,[coeffs(P,x)]))));
    for trial from 1 to min(max_trials,max_updates+1) do
        status := timelimit(univariate_search_remaining(deadline,"scalar trial"),
                            univariate_scale_status(f,P,lambda,cells,regions,x));
        if status="ok" then return [lambda,trial]; end if;
        if status="infeasible" then return [FAIL,trial]; end if;
        if status="small" then lo:=lambda; else hi:=lambda; end if;
        if lo=0 then lambda:=hi/2;
        elif hi=0 then lambda:=2*lo;
        else lambda:=(lo+hi)/2;
        end if;
    end do;
    return [FAIL,min(max_trials,max_updates+1)];
end proc;

# Returns SOS multipliers with a strictly positive residual on {B>=0},
# or FAIL. The residual need not be globally nonnegative at this stage.
local univariate_reduction := proc(f, gens, x, B,
                                  {max_trials::nonnegint := 64,
                                   max_updates::nonnegint := 32,
                                   deadline := time[real]()+300})
local cells, i, q, r, weights, result, trials := 0, seen := {}, square;
local try_direction;
    cells := timelimit(univariate_search_remaining(deadline,"reduction sign cells"),
                       univariate_select_cells(univariate_sign_chart(B,x),"nonnegative"));
    if timelimit(univariate_search_remaining(deadline,"reduction positivity"),
                 univariate_positive_on_cells(f,cells,x)) then
        return map(g->0,gens);
    end if;

    for i from 1 to nops(gens) do
        if gens[i]=0 then next; end if;
        q := timelimit(univariate_search_remaining(deadline,"polynomial division"),quo(f,gens[i],x));
        if not timelimit(univariate_search_remaining(deadline,"quotient nonnegativity"),
                         univariate_nonnegative_on_R(q,x)) then next; end if;
        r := timelimit(univariate_search_remaining(deadline,"polynomial remainder"),rem(f,gens[i],x));
        if timelimit(univariate_search_remaining(deadline,"remainder positivity"),
                     univariate_positive_on_cells(r,cells,x)) then
            weights := map(g->0,gens); weights[i]:=q;
            DEBUG(__FILE__, __LINE__, lprint(">> Univariate division reduction",i));
            return weights;
        end if;
    end do;

    try_direction := proc(w)
    local j, direction, scaled;
        if trials>=max_trials then return FAIL; end if;
        direction := timelimit(univariate_search_remaining(deadline,"direction construction"),
                               expand(add(w[j]*gens[j],j=1..nops(gens))));
        if member(direction,seen) then return FAIL; end if;
        seen := seen union {direction};
        scaled := univariate_scale_direction(f,direction,cells,x,max_trials-trials,max_updates,deadline);
        trials := trials+scaled[2];
        if scaled[1]=FAIL then return FAIL; end if;
        return map(weight->scaled[1]*weight,w);
    end proc;

    # Constants first, then low-degree polynomial squares. Joint directions
    # allow cancellation between generators; individual ones allow sparsity.
    for square in [1,x^2,(x-1)^2,(x+1)^2] do
        if trials>=max_trials then return FAIL; end if;
        result := try_direction(map(g->square,gens));
        if result<>FAIL then return result; end if;
        for i from 1 to nops(gens) do
            if trials>=max_trials then return FAIL; end if;
            weights:=map(g->0,gens); weights[i]:=square;
            result:=try_direction(weights);
            if result<>FAIL then return result; end if;
        end do;
    end do;
    return FAIL;
end proc;

local restore_generator_multipliers := proc(H, indices, original_count)
local restored, i;
    restored := [H[1],seq(0,i=1..original_count)];
    for i from 1 to nops(indices) do restored[indices[i]+1]:=H[i+1]; end do;
    return restored;
end proc;
