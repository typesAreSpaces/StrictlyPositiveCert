# Exact sign cells: [left, right, include_left, include_right].
# Keep algebraic endpoints exact; isolating-interval midpoints are only samples
# between roots, never substitutes for the roots themselves.
local univariate_sign := proc(value)
local v := evala(value);
    if is(v > 0) = true then return 1; end if;
    if is(v < 0) = true then return -1; end if;
    if is(v = 0) = true then return 0; end if;
    error "could not determine an exact algebraic sign: %1", value;
end proc;

local univariate_compare := proc(a, b)
    if evalb(a = b) then return 0; end if;
    if a = -infinity or b = infinity then return -1; end if;
    if a = infinity or b = -infinity then return 1; end if;
    return univariate_sign(a-b);
end proc;

local univariate_root_data := proc(p, x)
local q, boxes, roots, samples, i, width := 1, separated;
    if degree(p,x) <= 0 then return [[],[0]]; end if;
    q := primpart(quo(p,gcd(p,diff(p,x)),x),x);
    # Separate even neighboring boxes sharing a rational root as an endpoint.
    while true do
        boxes := realroot(q,width);
        separated := true;
        for i from 1 to nops(boxes)-1 do
            if boxes[i][2] >= boxes[i+1][1] then
                separated := false; break;
            end if;
        end do;
        if separated then break; end if;
        width := width/2;
    end do;
    if boxes = [] then return [[],[0]]; end if;
    roots := [seq(`if`(boxes[i][1]=boxes[i][2], boxes[i][1],
                       RootOf(q,x,boxes[i][1]..boxes[i][2])),i=1..nops(boxes))];
    samples := [boxes[1][1]-1,
                seq((boxes[i][2]+boxes[i+1][1])/2,i=1..nops(boxes)-1),
                boxes[-1][2]+1];
    return [roots,samples];
end proc;

local univariate_sign_chart := proc(p0, x)
local p := expand(p0), roots, samples, chart := [], i, a, b;
    if p = 0 then return [[[-infinity,infinity,false,false],0]]; end if;
    roots,samples := op(univariate_root_data(p,x));
    for i from 1 to nops(samples) do
        a := `if`(i=1,-infinity,roots[i-1]);
        b := `if`(i>nops(roots),infinity,roots[i]);
        chart := [op(chart),[[a,b,false,false],univariate_sign(subs(x=samples[i],p))]];
        if i <= nops(roots) then
            chart := [op(chart),[[b,b,true,true],0]];
        end if;
    end do;
    return chart;
end proc;

local univariate_matches := proc(s, relation)
    if relation = "positive" then return evalb(s>0);
    elif relation = "nonnegative" then return evalb(s>=0);
    elif relation = "negative" then return evalb(s<0);
    elif relation = "nonpositive" then return evalb(s<=0);
    else error "unknown sign relation: %1",relation;
    end if;
end proc;

local univariate_select_cells := proc(chart, relation)
local result := [], entry, cell, last;
    for entry in chart do
        if not univariate_matches(entry[2],relation) then next; end if;
        cell := entry[1];
        if result <> [] then
            last := result[-1];
            if evalb(last[2]=cell[1]) and (last[4] or cell[3]) then
                result[-1] := [last[1],cell[2],last[3],cell[4]];
                next;
            end if;
        end if;
        result := [op(result),cell];
    end do;
    return result;
end proc;

local univariate_intersect_cells := proc(A, B)
local result := [], i := 1, j := 1, a, b, left, right, lc, rc, c, width_sign;
    while i<=nops(A) and j<=nops(B) do
        a,b := A[i],B[j];
        c := univariate_compare(a[1],b[1]);
        if c>0 then left,lc := a[1],a[3];
        elif c<0 then left,lc := b[1],b[3];
        else left,lc := a[1],evalb(a[3] and b[3]); end if;
        c := univariate_compare(a[2],b[2]);
        if c<0 then right,rc := a[2],a[4];
        elif c>0 then right,rc := b[2],b[4];
        else right,rc := a[2],evalb(a[4] and b[4]); end if;
        width_sign := univariate_compare(left,right);
        if width_sign<0 then
            result := [op(result),[left,right,lc,rc]];
        elif width_sign=0 and lc and rc then
            # Equal algebraic numbers can have different defining polynomials.
            result := [op(result),[left,left,true,true]];
        end if;
        if c<=0 then i:=i+1; end if;
        if c>=0 then j:=j+1; end if;
    end do;
    return result;
end proc;

local univariate_cells_to_sas := proc(cells, x)
local result := [], cell, constraints;
    for cell in cells do
        if evalb(cell[1]=cell[2]) then
            result := [op(result),[x=cell[1]]]; next;
        end if;
        constraints := [];
        if cell[1] <> -infinity then
            constraints := [`if`(cell[3],cell[1]<=x,cell[1]<x)];
        end if;
        if cell[2] <> infinity then
            constraints := [op(constraints),`if`(cell[4],x<=cell[2],x<cell[2])];
        end if;
        if constraints=[] then constraints := [x=x]; end if;
        result := [op(result),constraints];
    end do;
    return result;
end proc;

# Conjunction of univariate polynomial sign conditions; no general CAD call.
local univariate_sign_set := proc(polys, relations, x)
local cells := [[-infinity,infinity,false,false]], i;
    for i from 1 to nops(polys) do
        cells := univariate_intersect_cells(cells,
                                            univariate_select_cells(univariate_sign_chart(polys[i],x),relations[i]));
        if cells=[] then break; end if;
    end do;
    return univariate_cells_to_sas(cells,x);
end proc;

local univariate_odd_part := proc(p, x)
local decomposition := sqrfree(p,x), item, result;
    result := decomposition[1];
    for item in decomposition[2] do
        if type(item[2],odd) then result := result*item[1]; end if;
    end do;
    return expand(result);
end proc;

# Primitive squarefree part with rational coefficients; same real roots as p.
local univariate_squarefree := proc(p, x)
    if degree(p,x) <= 0 then return p; end if;
    return primpart(quo(p,gcd(p,diff(p,x)),x),x);
end proc;

# realroot isolates by Descartes' rule with exact rationals. On high-degree
# residuals with large coefficients it is orders of magnitude cheaper than
# building a Sturm sequence, which suffers from coefficient growth.
local univariate_has_real_root := proc(p, x)
    return evalb(degree(p,x) > 0 and realroot(univariate_squarefree(p,x)) <> []);
end proc;

local univariate_nonnegative_on_R := proc(p0, x)
local p := expand(p0);
    if p=0 then return true; end if;
    if lcoeff(p,x)<0 then return false; end if;
    return not univariate_has_real_root(univariate_odd_part(p,x),x);
end proc;

# Decide whether some root of the squarefree q lies in an open cell interior.
# A nondegenerate isolating interval holds one root in its open interior;
# when it straddles a cell endpoint, either that endpoint is the root (it is
# then outside the interior) or a finer isolation separates them. Return FAIL
# if max_refinements finer isolations still leave a straddling interval.
local univariate_roots_in_interiors := proc(q, cells, x, max_refinements::nonnegint := 40)
local intervals, interval, cell, l, u, a, b, undecided, k, refinement, inside;
    intervals := realroot(q);
    for refinement from 0 to max_refinements do
        undecided := false;
        for interval in intervals do
            l,u := op(interval);
            for cell in cells do
                a,b := cell[1],cell[2];
                if evalb(a=b) then next; end if;
                if l = u then
                    inside := evalb(univariate_compare(a,l)<0 and univariate_compare(l,b)<0);
                elif univariate_compare(u,a) <= 0 or univariate_compare(b,l) <= 0 then
                    inside := false;
                elif univariate_compare(a,l) <= 0 and univariate_compare(u,b) <= 0 then
                    inside := true;
                else
                    # The interval straddles a or b. An endpoint inside (l,u)
                    # that is a root of q is the isolated root itself.
                    inside := FAIL;
                    for k in [a,b] do
                        if k <> -infinity and k <> infinity and
                        univariate_compare(l,k) < 0 and univariate_compare(k,u) < 0 and
                        univariate_sign(subs(x=k,q)) = 0 then
                            inside := false;
                        end if;
                    end do;
                end if;
                if inside = true then return true; end if;
                if inside = FAIL then undecided := true; end if;
            end do;
        end do;
        if not undecided then return false; end if;
        intervals := realroot(q,1/2^refinement);
    end do;
    return FAIL;
end proc;

# Sturm sequence of a squarefree p. Remainders are divided by positive
# rationals only, which preserves every sign variation.
local univariate_sturm_sequence := proc(p, x)
local sequence := [p, diff(p,x)], r;
    while degree(sequence[-1],x) > 0 do
        r := -rem(sequence[-2],sequence[-1],x);
        if r = 0 then break; end if;
        sequence := [op(sequence), r/abs(content(r,x))];
    end do;
    return sequence;
end proc;

# Number of distinct roots in (a,b] of the squarefree polynomial that
# generated the sequence. Signs are exact, also at algebraic endpoints;
# Maple's sturm drops a root at a RootOf endpoint b, counting only (a,b).
local univariate_sturm_count := proc(sequence, x, a, b)
local variations := proc(t)
    local signs, count := 0, i;
        if t = infinity then
            signs := map(s -> sign(lcoeff(s,x)), sequence);
        elif t = -infinity then
            signs := map(s -> sign(lcoeff(s,x))*(-1)^degree(s,x), sequence);
        else
            signs := map(s -> univariate_sign(subs(x=t,s)), sequence);
        end if;
        signs := remove(s -> s = 0, signs);
        for i from 1 to nops(signs)-1 do
            if signs[i] <> signs[i+1] then count := count+1; end if;
        end do;
        return count;
    end proc;
    return variations(a)-variations(b);
end proc;

# Non-strict positivity permits even-multiplicity roots, hence locates only
# roots of the odd part in interiors. Root isolation uses realroot; Sturm
# sequences remain only as a fallback for unresolved endpoint proximity.
local univariate_positive_on_cells := proc(p0, cells, x, strict::truefalse := true)
local p := expand(p0), q, qs, cell, a, b, sample, count, s, located, sequence;
    if cells=[] then return true; end if;
    if p=0 then return not strict; end if;
    q := `if`(strict,p,univariate_odd_part(p,x));
    for cell in cells do
        if cell[3] then
            s := univariate_sign(subs(x=cell[1],p));
            if s<0 or (strict and s=0) then return false; end if;
        end if;
        if cell[4] then
            s := univariate_sign(subs(x=cell[2],p));
            if s<0 or (strict and s=0) then return false; end if;
        end if;
    end do;
    if degree(q,x)>0 then
        qs := univariate_squarefree(q,x);
        located := univariate_roots_in_interiors(qs,cells,x);
        if located = true then return false; end if;
        if located = FAIL then
            sequence := univariate_sturm_sequence(qs,x);
            for cell in cells do
                a,b := cell[1],cell[2];
                if evalb(a=b) then next; end if;
                # The count covers (a,b]; discount a root at b, an endpoint.
                count := univariate_sturm_count(sequence,x,a,b);
                if b<>infinity and univariate_sign(subs(x=b,qs))=0 then count:=count-1; end if;
                if count>0 then return false; end if;
            end do;
        end if;
    end if;
    # q has no root in any interior, so one sample fixes its sign there.
    for cell in cells do
        a,b := cell[1],cell[2];
        if evalb(a=b) then next; end if;
        sample := samplePoint(a,b);
        if univariate_sign(subs(x=sample,q))<=0 then return false; end if;
    end do;
    return true;
end proc;
