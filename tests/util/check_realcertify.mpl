_pwd := currentdir();
currentdir(FileTools[ParentDirectory](__FILE__));

$include "bitsize.mpl";
$include "quickdegree.mpl";

currentdir(_pwd);

#$define SHOW_CHECK_REALCERTIFY

extractCertificates := proc(out)
local i, j;
local certificates, sos, rlist, oneWithBasis, s, idx, idxi;
## Verifying output
    if (not(out[1] = false)) then
        certificates := [];
        s := 0;
        idx := 0;
        sos := out[4];
        rlist := out[3];
        oneWithBasis := [1, op(basis)];

        for i from 1 to nops(basis)+1 do
            certificates := [op(certificates), 0];
            idxi:=idx+rlist[i];

            for j from idx+1 to idxi do
                s := s + oneWithBasis[i]*sos[2*j-1]*sos[2*j]^2;
                certificates[i] := certificates[i] + sos[2*j-1]*sos[2*j]^2;
            od;

            idx:=idxi;
        od;
        return true, certificates, s;
    fi;
    return false, [0], 0;
end proc;

checkRealCertify := proc(f, basis, test_name)
    printf("\n>> Test\n%s\n", test_name);
    printf(">> input polynomial\n%s\n", convert(f, string));
    printf(">> basis\n%s\n", convert(basis, string));
local st := time();
local out, isvalid, certificates, s, h;
    out := multivsos_internal(f, glist=basis, gmp=true, relaxorder=2);
    #out := multivsos_internal(f, glist=basis, algo=1, gmp=true);
    #out := multivsos_internal(f, glist=basis, algo=2, gmp=true);
    st := time() - st;
    isvalid, certificates, s := extractCertificates(out);

    if isvalid then
      printf(">> Time taken\n%f\n", st);
      lprint(">> Sums of squares multipliers", certificates);
      printf(">> Degree size\n%s\n", convert(foldl((_x, _y) -> max(_x, _y), 0, op(map(h -> quickdegree(h, x), certificates))), string));
      #printf(">> Bitsize complexity of certificates: %s\n", convert(add(bitsizeP(h), h in certificates), string));

$ifdef SHOW_CHECK_REALCERTIFY
lprint("Certificates produced", certificates);
lprint("Test exact certificate", expand (f - s));
lprint("Test certificates produced are sums of squares");
$endif

      return 0;
    else
      printf(">> RealCertify couldnt find certificates\n%f\n", st);
      return 1;
    end if; 
end proc;
