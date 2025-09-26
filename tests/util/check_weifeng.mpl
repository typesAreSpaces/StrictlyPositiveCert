_pwd := currentdir();
currentdir(FileTools[ParentDirectory](__FILE__));

$include "bitsize.mpl";
$include "quickdegree.mpl";

currentdir(_pwd);

#$define SHOW_CHECK_WEIFENG

checkWeifeng := proc(f, basis, x, test_name)
    printf("\n>> Test\n%s\n", test_name);
    printf(">> input polynomial\n%s\n", convert(f, string));
    printf(">> basis\n%s\n", convert(basis, string));
local st := time();
local H := spCertificates(f, basis, x);
    printf(">> Time taken\n%f\n", time() - st);
    lprint(">> Sums of squares multipliers", H);
    printf(">> Degree size\n%s\n", convert(foldl((_x, _y) -> max(_x, _y), 0, op(map(h -> quickdegree(h, x), H))), string));
    #printf(">> Bitsize complexity of certificates: %s\n", convert(add(bitsizeP(h, x), h in H), string));

$ifdef SHOW_CHECK_WEIFENG
local g_1 := basis[1];
local g_2 := basis[2];
local s := H[1];
lprint("Test exact certificate", f - (s + g_1*H[2] + g_2*H[3]));
lprint("Test certificates produced are sums of squares");
lprint(map(poly -> SemiAlgebraic([poly < 0], [x]), H));
$endif

    return;
end proc;
