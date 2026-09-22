spCertificates := proc(f, basis, x)
$ifdef LOG_TIME
    INIT_START_LOG_TIME("spCertificates",0)
$endif
    if SemiAlgebraic([f < 0], [x]) = [] then
        return [f, op(map(0, basis))];
    end if;
local g, H2, f2, H3, f3, H4, certificates, eps_LS;
local i;
local quickcheck, quickcertificates;
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

    g := bound_poly(basis, x);
    DEBUG(__FILE__, __LINE__, lprint(">> Done with bound_poly"));
    DEBUG(__FILE__, __LINE__, lprint(">> bound_poly g", g));

    DEBUG(__FILE__, __LINE__, lprint(">> Poly f for averkov_lemma_7", f));
    DEBUG(__FILE__, __LINE__, lprint(">> Poly g[1] for averkov_lemma_7", g[1]));
    H2 := averkov_lemma_7(x, f, basis, g[1], N_GUESS_AVKL);
    DEBUG(__FILE__, __LINE__, lprint(">> Done with averkov_lemma_7"));

    DEBUG(__FILE__, __LINE__, lprint(">> H2", H2));
    f2 := f - dot_product(basis, H2);
$ifdef ENABLE_AVERKOV_CHECK
    lprint(">> 2. Checking correctness of averkov_lemma_7",
           SemiAlgebraic([g[1] >= 0, f2 <= 0], [x]));
$endif
    DEBUG(__FILE__, __LINE__, lprint(">> f2", f2));

$ifndef WEIFENG_OPTIMIZATION
    H3, eps_LS := lower_bound_poly(x, f2, g[1]);
    DEBUG(__FILE__, __LINE__, lprint(">> Done with Lower_bound_poly"));
    DEBUG(__FILE__, __LINE__, lprint(">> H3", H3));
$endif

$ifdef WEIFENG_OPTIMIZATION
    f3 := f2;
$else
    f3 := f2 - g[1]*H3;
$endif
    DEBUG(__FILE__, __LINE__, lprint(">> f3", f3));
    DEBUG(__FILE__, __LINE__, lprint(">> g[1]", g[1]));
    H4 := averkov_extended_lemma(x, f3, g[1], N_GUESS_LS, eps_LS);
    DEBUG(__FILE__, __LINE__, lprint(">> Done with Last_step"));
    DEBUG(__FILE__, __LINE__, lprint(">> H4", H4));

    certificates := H2;
$ifdef WEIFENG_OPTIMIZATION
    certificates[g[4]] := certificates[g[4]] + H4*g[2];
$else
    certificates[g[4]] := certificates[g[4]] + (H3+H4)*g[2];
$endif
    if g[3] <> 0 then
$ifdef WEIFENG_OPTIMIZATION
        certificates[g[5]] := certificates[g[5]] + H4*g[3];
$else
        certificates[g[5]] := certificates[g[5]] + (H3+H4)*g[3];
$endif
    end if;

    certificates := [f - dot_product(basis, certificates), op(certificates)];
$ifdef ENABLE_VERIFICATION
    lprint(">> This should be zero", expand(f - dot_product([1, op(basis)], certificates)));
$endif
    DEBUG(__FILE__, __LINE__, lprint(">> Certificates found", op(certificates)));
$ifdef LOG_TIME
    END_LOG_TIME("spCertificates",0);
$endif
    return certificates;
end proc;
