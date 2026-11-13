spCertificates := proc(f, basis, x)
$ifdef LOG_TIME
    INIT_START_LOG_TIME("spCertificates",0)
$endif
    if SemiAlgebraic([f < 0], [x]) = [] then
        return [f, op(map(0, basis))];
    end if;
local g, H2, f2, H3, f3, H4, certificates, eps_LS := -1;
local i;
local quickcheck, quickcertificates;
local fpos, fnonneg, bound_multiplier, remainder;
local divided := false;
# Cumulative SOS factors for the input, first remainder, and final remainder.
# Each stage's certificate must be multiplied by its corresponding entry.
local sos_factors := [1, 1, 1];
    certificates := map(gen -> 0, basis);

    for i from 1 to nops(basis) do
        if basis[i] = -1 then
            certificates[i] := ((f - 1)/2)^2;
            return [((f + 1)/2)^2, op(certificates)];
        end if;
    end do;

$ifdef ENABLE_QUICK_CERT
    quickcheck, quickcertificates := quickCertificates(f, basis, x);
    if quickcheck then
      DEBUG(__FILE__, __LINE__, lprint(">> Quick certificates found", op(quickcertificates)));
      return quickcertificates;
    end if;
$endif

$ifdef FACTOR_SOS_PRE
    fnonneg, fpos := SplitNonnegativeFactor(f, x);
$else
    fnonneg, fpos := 1, f;
$endif
    sos_factors[1] := fnonneg;

    g := bound_poly(basis, x);
    DEBUG(__FILE__, __LINE__, lprint(">> Done with bound_poly"));
    DEBUG(__FILE__, __LINE__, lprint(">> bound_poly g", g));

    DEBUG(__FILE__, __LINE__, lprint(">> Poly f for averkov_lemma_7", fpos));
    DEBUG(__FILE__, __LINE__, lprint(">> Poly g[1] for averkov_lemma_7", g[1]));
    H2 := averkov_lemma_7(x, fpos, basis, g[1], N_GUESS_AVKL);
    DEBUG(__FILE__, __LINE__, lprint(">> Done with averkov_lemma_7"));

    DEBUG(__FILE__, __LINE__, lprint(">> H2", H2));
    f2 := fpos - dot_product(basis, H2);
$ifdef ENABLE_AVERKOV_CHECK
    lprint(">> 2. Checking correctness of averkov_lemma_7",
           SemiAlgebraic([g[1] >= 0, f2 <= 0], [x]));
$endif
    DEBUG(__FILE__, __LINE__, lprint(">> f2", f2));

$ifdef FACTOR_SOS_PRE
    fnonneg, fpos := SplitNonnegativeFactor(f2, x);
$else
    fnonneg, fpos := 1, f2;
$endif
    sos_factors[2] := sos_factors[1]*fnonneg;
    f2 := fpos;

$ifdef WEIFENG_OPTIMIZATION
    f3 := f2;
$else
# f2 = H3*g[1] + f3 with SOS H3 and f3 skips lower_bound_poly, and
# the final stage then closes with H4 = 0.
$ifdef ENABLE_STAGE_DIVISION
    divided, H3, f3 := divisionCertificate(f2, g[1], x);
    DEBUG(__FILE__, __LINE__, lprint(">> f2 division by g[1]", divided));
$endif
    if not divided then
        H3, eps_LS := lower_bound_poly(x, f2, g[1]);
        DEBUG(__FILE__, __LINE__, lprint(">> Done with Lower_bound_poly"));
        DEBUG(__FILE__, __LINE__, lprint(">> H3", H3));
        DEBUG(__FILE__, __LINE__, lprint(">> eps_LS", eps_LS));
        f3 := f2 - g[1]*H3;
    end if;
$endif
    DEBUG(__FILE__, __LINE__, lprint(">> f3", f3));
    DEBUG(__FILE__, __LINE__, lprint(">> g[1]", g[1]));
$ifdef FACTOR_SOS_PRE
    fnonneg, fpos := SplitNonnegativeFactor(f3, x);
$else
    fnonneg, fpos := 1, f3;
$endif
    sos_factors[3] := sos_factors[2]*fnonneg;
    f3 := fpos;
# f3 = H4*g[1] + remainder avoids the high-exponent Averkov multiplier.
$ifdef ENABLE_STAGE_DIVISION
    divided, H4, remainder := divisionCertificate(f3, g[1], x);
    DEBUG(__FILE__, __LINE__, lprint(">> f3 division by g[1]", divided));
$endif
    if not divided then
        H4 := averkov_extended_lemma(x, f3, g[1], N_GUESS_LS, eps_LS);
    end if;
    DEBUG(__FILE__, __LINE__, lprint(">> Done with Last_step"));
    DEBUG(__FILE__, __LINE__, lprint(">> H4", H4));

    certificates := map(h -> sos_factors[1]*h, H2);
$ifdef WEIFENG_OPTIMIZATION
    bound_multiplier := sos_factors[3]*H4;
$else
    bound_multiplier := sos_factors[2]*H3 + sos_factors[3]*H4;
$endif
    certificates[g[4]] := certificates[g[4]] + bound_multiplier*g[2];
    if g[3] <> 0 then
        certificates[g[5]] := certificates[g[5]] + bound_multiplier*g[3];
    end if;

    # Lift the final SOS as well, retaining the extracted factors explicitly.
    certificates := [sos_factors[3]*(f3 - g[1]*H4), op(certificates)];
$ifdef ENABLE_VERIFICATION
    lprint(">> This should be zero", expand(f - dot_product([1, op(basis)], certificates)));
$endif
    DEBUG(__FILE__, __LINE__, lprint(">> Certificates found", op(certificates)));
$ifdef LOG_TIME
    END_LOG_TIME("spCertificates",0);
$endif
    return certificates;
end proc;
