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
local p := expand(p0), data, roots, samples, chart := [], i, a, b;
    if p = 0 then return [[[-infinity,infinity,false,false],0]]; end if;
    data := univariate_root_data(p,x); roots,samples := op(data);
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

local univariate_nonnegative_on_R := proc(p0, x)
local p := expand(p0), oddpart;
    if p=0 then return true; end if;
    if lcoeff(p,x)<0 then return false; end if;
    oddpart := univariate_odd_part(p,x);
    return evalb(degree(oddpart,x)<=0 or sturm(oddpart,x,-infinity,infinity)=0);
end proc;

# A single sequence is reused across every component. Non-strict positivity
# permits even-multiplicity roots, hence counts only the odd part in interiors.
local univariate_positive_on_cells := proc(p0, cells, x, strict::truefalse := true)
local p := expand(p0), q, sequence := FAIL, cell, a, b, sample, count, s;
    if cells=[] then return true; end if;
    if p=0 then return not strict; end if;
    q := `if`(strict,p,univariate_odd_part(p,x));
    for cell in cells do
        a,b := cell[1],cell[2];
        if cell[3] then
            s := univariate_sign(subs(x=a,p));
            if s<0 or (strict and s=0) then return false; end if;
        end if;
        if cell[4] then
            s := univariate_sign(subs(x=b,p));
            if s<0 or (strict and s=0) then return false; end if;
        end if;
        if evalb(a=b) then next; end if;
        if degree(q,x)>0 then
            if sequence=FAIL then sequence := sturmseq(q,x); end if;
            count := sturm(sequence,x,a,b);
            # sturm counts (a,b]; excluded endpoints and non-strict tests
            # permit a root at b. An included strict endpoint was checked above.
            if b<>infinity and (not strict or not cell[4]) and
            univariate_sign(subs(x=b,q))=0 then count:=count-1; end if;
            if count>0 then return false; end if;
        end if;
        sample := samplePoint(a,b);
        if univariate_sign(subs(x=sample,q))<=0 then return false; end if;
    end do;
    return true;
end proc;
