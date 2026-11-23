# Returns a list of the form [g, h1, h2, index1, index2]
# where
# - g is the bound poly
# - h1, h2 are the sums of squares which make the linear combination for g
# - index1, index2 are the indexes in gens which make the linear combination for g
local bound_poly := proc(gens, x)
local degrees, fst_coeffs, snd_coeffs, h1, h2;
local i, j;
$ifdef LOG_TIME
    INIT_START_LOG_TIME("bound_poly")
    try
$endif
    degrees := map(poly -> degree(poly, x), gens);
    DEBUG(__FILE__, __LINE__, lprint(">> degrees", degrees));

    fst_coeffs := map[indices](
        i -> coeff(gens[i], x^degrees[i]), gens);
    DEBUG(__FILE__, __LINE__, lprint(">> fst_coeffs", fst_coeffs));

    snd_coeffs := map[indices](
        i ->
        if degrees[i] = 1 then
            subs(x = 0, gens[i])/abs(fst_coeffs[i])
        else
            coeff(gens[i]/abs(fst_coeffs[i]), x^(degrees[i] - 1))
        end if, gens);
    DEBUG(__FILE__, __LINE__, lprint(">> snd_coeffs", snd_coeffs));

    for i from 1 to nops(gens) do
        if type(degrees[i], even) and fst_coeffs[i] < 0 then
            return [gens[i], 1, 0, i, i];
        end if;
    end do;

# At this point, every element in the gens has odd degree
# or even degree with positive leading term
    for i from 1 to nops(gens) - 1 do
        if type(degrees[i], 'even') then
            next;
        end if;
        for j from i + 1 to nops(gens) do
            if type(degrees[j], 'even') then
                next;
            end if;
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
            DEBUG(__FILE__, __LINE__, lprint(">> h1", h1));
            DEBUG(__FILE__, __LINE__, lprint(">> h2", h2));
            return [expand(h1*gens[i] + h2*gens[j]), h1, h2, i, j];
        end do;
    end do;
    return NULL;
$ifdef LOG_TIME
    finally
        END_LOG_TIME("bound_poly");
    end try;
$endif
end proc;
