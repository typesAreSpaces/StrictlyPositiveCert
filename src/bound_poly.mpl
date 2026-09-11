# Returns a list of the form [g, h1, h2, index1, index2]
# where
# - g is the bound poly
# - h1, h2 are the sums of squares which make the linear combination for g
# - index1, index2 are the indexes in basis which make the linear combination for g
local bound_poly := proc(basis, x)
$ifdef LOG_TIME
    INIT_START_LOG_TIME("bound_poly",0)
$endif
local degrees, fst_coeffs, snd_coeffs, h1, h2;
local i, j;
    degrees := map(poly -> degree(poly, x), basis);
    DEBUG(__FILE__, __LINE__, ENABLE_DEBUGGING, lprint(">> degrees", degrees));

    fst_coeffs := map[indices](
        i -> coeff(basis[i], x^degrees[i]), basis);
    DEBUG(__FILE__, __LINE__, ENABLE_DEBUGGING, lprint(">> fst_coeffs", fst_coeffs));

    snd_coeffs := map[indices](
        i ->
        if degrees[i] = 1 then
            subs(x = 0, basis[i])/abs(fst_coeffs[i])
        else
            coeff(basis[i]/abs(fst_coeffs[i]), x^(degrees[i] - 1))
        end if, basis);
    DEBUG(__FILE__, __LINE__, ENABLE_DEBUGGING, lprint(">> snd_coeffs", snd_coeffs));

    for i from 1 to nops(basis) do
        if type(degrees[i], even) and fst_coeffs[i] < 0 then
$ifdef LOG_TIME
            END_LOG_TIME("bound_poly",0)
$endif
            return [basis[i], 1, 0, i, i];
        end if;
    end do;

# At this point, every element in the basis has odd degree
    for i from 1 to nops(basis) - 1 do
        for j from i + 1 to nops(basis) do
            if fst_coeffs[i]*fst_coeffs[j] > 0 then
                next;
            end if;
            if degrees[i] = degrees[j] then
                h1 := x^2/abs(fst_coeffs[i]);
                h2 := (x + sign(fst_coeffs[i])*(1/2*snd_coeffs[i] + 1/2*snd_coeffs[j] + 1))^2
                /abs(fst_coeffs[j]);
            elif degrees[j] < degrees[i] then
                h1 := 1/abs(fst_coeffs[i]);
                h2 := x^(degrees[i] - degrees[j] - 2)
                *(x + sign(fst_coeffs[i])*(1/2*snd_coeffs[i] + 1/2*snd_coeffs[j] + 1))^2
                /abs(fst_coeffs[j]);
            else
                h1 := x^(degrees[j] - degrees[i] - 2)
                *(x - sign(fst_coeffs[i])*(1/2*snd_coeffs[i] + 1/2*snd_coeffs[j] + 1))^2
                /abs(fst_coeffs[i]);
                h2 := 1/abs(fst_coeffs[j]);
            end if;
$ifdef LOG_TIME
            END_LOG_TIME("bound_poly",0)
$endif
DEBUG(__FILE__, __LINE__, ENABLE_DEBUGGING, lprint(">> h1", h1));
DEBUG(__FILE__, __LINE__, ENABLE_DEBUGGING, lprint(">> h2", h2));
            return [expand(h1*basis[i] + h2*basis[j]), h1, h2, i, j];
        end do;
    end do;
end proc;
