spCertificates := proc(f, basis, x,
    {max_k::nonnegint := 100, max_bisection::nonnegint := 200})
$ifdef LOG_TIME
    INIT_START_LOG_TIME("spCertificates",0)
$endif
    if SemiAlgebraic([f < 0], [x]) = [] then
        return [f, op(map(0, basis))];
    end if;
local g, H2, f2, monogenic, certificates;
local i;
local quickcheck, quickcertificates;
    for i from 1 to nops(basis) do
        if type(basis[i], rational) and basis[i] < 0 then
            certificates := map(gen -> 0, basis);
            certificates[i] := ((f - 1)/2)^2/(-basis[i]);
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

    # Global keyword names avoid substitution of the caller's parameter values.
    monogenic := monogenicCertificates(f2, g[1], x,
        ':-max_k'=max_k, ':-max_bisection'=max_bisection);
    DEBUG(__FILE__, __LINE__, lprint(">> Done with monogenicCertificates"));

    certificates := H2;
    certificates[g[4]] := certificates[g[4]] + monogenic[2]*g[2];
    if g[3] <> 0 then
        certificates[g[5]] := certificates[g[5]] + monogenic[2]*g[3];
    end if;

    certificates := [monogenic[1], op(certificates)];
$ifdef ENABLE_VERIFICATION
    lprint(">> This should be zero", expand(f - dot_product([1, op(basis)], certificates)));
$endif
    DEBUG(__FILE__, __LINE__, lprint(">> Certificates found", op(certificates)));
$ifdef LOG_TIME
    END_LOG_TIME("spCertificates",0);
$endif
    return certificates;
end proc;
