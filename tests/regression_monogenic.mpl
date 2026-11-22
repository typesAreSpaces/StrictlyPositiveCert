# Run from the repository root with the Makefile's -D values.
read "src/StrictlyPositiveCert.mpl":

check_monogenic := proc()
local x, y, i, item, f, g, basis, cert, term, squares, sigma, reconstructed;
local check_polynomials, caught, samples;
    check_polynomials := proc(f, basis, cert, x)
    local term;
        if nops(cert) <> nops(basis)+1 or
        expand(f-StrictlyPositiveCert:-dot_product([1,op(basis)],cert)) <> 0 then
            error "certificate reconstruction failed for %1", [f,basis];
        end if;
        for term in cert do
            if not type(term,polynom(rational,x)) or
            SemiAlgebraic([term < 0],[x]) <> [] then
                error "invalid SOS multiplier for %1", [f,basis];
            end if;
        end do;
    end proc;
    samples := [
        [10^9*(1-x^2)+1,1-x^2],
        [x+3/2,1-x^2],
        [x^2-1/2,(1-x^2)*(x^2-4)],
        [x^3+1/10,x*(1-x)],
        [2-x^3,1-x^2],
        [1-x^3,-x^2],
        [x+1,-x^2],
        [4/3-x^2,(1-x^2)^3],
        [x,-x^2-1],
        [x,-1],
        [x^3-2*x+1,-2/3],
        [-x^2,-2],
        [-1,-1],
        [-1,-2],
        [1,1-x^2],
        [0,-1],
        [3/2+2^260*x^4-x^2/2,1-x^2],
        [1+(1-x^2)*(x/2)^20+x^24/2^50,1-x^2]
               ];
    for item in samples do
        f,g := op(item);
        cert := StrictlyPositiveCert:-monogenicCertificates(f,g,x);
        check_polynomials(f,[g],cert,x);
        # Public integration retains the polynomial-list certificate format.
        cert := StrictlyPositiveCert:-spCertificates(f,[g],x);
        check_polynomials(f,[g],cert,x);
    end do;
    for item in [[1/10^6,1/1000],[1/10^12,1/10^9],[1,1]] do
        f := item[1]*(x+2); g := item[2]*(1-x^2);
        cert := StrictlyPositiveCert:-spCertificates(f,[g],x);
        check_polynomials(f,[g],cert,x);
    end do;
    for item in [[x+2,[x+1,1-x]], [x+1,[4-x^2,x]]] do
        cert := StrictlyPositiveCert:-spCertificates(item[1],item[2],x);
        check_polynomials(item[1],item[2],cert,x);
    end do;
    cert := StrictlyPositiveCert:-monogenicCertificates(y+3/2,1-y^2,y,max_k=0);
    check_polynomials(y+3/2,[1-y^2],cert,y);

    # Redundant constants must not reach bound construction. Restore each
    # removed entry as zero, preserving the original generator indices.
    for item in [[x+3/2,[2,1-x^2,0,3/2]],
                 [x+3/2,[x+1,2,1-x,0]],
                 [x^3,[1,0,-2/3,4]],
                 [x^2+1,[0,1,3/2]]] do
        f,basis := op(item);
        cert := StrictlyPositiveCert:-spCertificates(f,basis,x);
        check_polynomials(f,basis,cert,x);
        for i from 1 to nops(basis) do
            if type(basis[i],rational) and basis[i]>=0 and cert[i+1]<>0 then
                error "redundant constant received a nonzero multiplier: %1",item;
            end if;
        end do;
    end do;
    caught := false;
    try StrictlyPositiveCert:-spCertificates(x,[0,1,2],x);
    catch "no nonredundant generators remain": caught := true;
    end try;
    if not caught then error "constant-only basis did not fail explicitly"; end if;
    caught := false;
    try StrictlyPositiveCert:-spCertificates(2-x^3,[2,1-x^2],x,max_k=0);
    catch: caught := true;
    end try;
    if not caught then error "constant filtering did not forward the search budget"; end if;

    # Explicit square decompositions agree with the returned polynomial identity.
    for item in [[x+3/2,1-x^2], [2-x^3,1-x^2],
                 [3/2+2^260*x^4-x^2/2,1-x^2],
                 [x,-1], [x^3-2*x+1,-2/3], [1,-1], [-1,-1], [0,-2]] do
        cert := StrictlyPositiveCert:-monogenicCertificates(op(item),x,
                                                            output_sos=true,max_sos_iterations=66);
        reconstructed := [];
        for squares in cert do
            sigma := 0;
            for term in squares do
                if not type(term[1],rational) or term[1]<0 or
                not type(term[2],polynom(rational,x)) then
                    error "invalid weighted square";
                end if;
                sigma := sigma+term[1]*term[2]^2;
            end do;
            reconstructed := [op(reconstructed),expand(sigma)];
        end do;
        check_polynomials(item[1],[item[2]],reconstructed,x);
    end do;

    # Negative constants use the square identity without a search or reduction.
    cert := StrictlyPositiveCert:-monogenicCertificates(x,-1,x,
                              max_k=0,max_bisection=0,max_query_retries=0);
    check_polynomials(x,[-1],cert,x);
    cert := StrictlyPositiveCert:-monogenicCertificates(x^3,-2/3,x,
                              output_sos=true,max_k=0,max_bisection=0,
                              max_query_retries=0,max_sos_iterations=0);
    reconstructed := map(squares -> expand(add(term[1]*term[2]^2,term in squares)),cert);
    check_polynomials(x^3,[-2/3],reconstructed,x);

    # Invalid inputs and exhausted budgets must fail explicitly.
    for item in [[-1,1-x^2], [1,1+x^2], [x+1.5,1-x^2]] do
        caught := false;
        try StrictlyPositiveCert:-monogenicCertificates(op(item),x);
        catch: caught := true;
        end try;
        if not caught then error "invalid input was accepted: %1",item; end if;
    end do;
    caught := false;
    try StrictlyPositiveCert:-monogenicCertificates(2-x^3,1-x^2,x,max_k=0);
    catch: caught := true;
    end try;
    if not caught then error "max_k budget was ignored"; end if;
    caught := false;
    try StrictlyPositiveCert:-spCertificates(2-x^3,[1-x^2],x,max_k=0);
    catch: caught := true;
    end try;
    if not caught then error "public max_k budget was not forwarded"; end if;
    caught := false;
    try StrictlyPositiveCert:-monogenicCertificates(10^9*(1-x^2)+1,1-x^2,x,
                                                    max_k=0,max_bisection=0);
    catch: caught := true;
    end try;
    if not caught then error "max_bisection budget was ignored"; end if;
    caught := false;
    try StrictlyPositiveCert:-spCertificates(10^9*(1-x^2)+1,[1-x^2],x,
                                             max_k=0,max_bisection=0);
    catch: caught := true;
    end try;
    if not caught then error "public max_bisection budget was not forwarded"; end if;
    caught := false;
    try StrictlyPositiveCert:-monogenicCertificates(3/2+2^260*x^4-x^2/2,1-x^2,x,
                                                    output_sos=true,max_sos_iterations=64);
    catch: caught := true;
    end try;
    if not caught then error "SOS refinement budget was ignored"; end if;
    printf("PASS: monogenic and public certificates, square decompositions, edge cases, and budgets\n");
end proc:
check_monogenic():
