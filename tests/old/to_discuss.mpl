read "../check_weifeng.mpl";
read "../check_realcertify.mpl";

getenv(OMP_NUM_THREADS);
kernelopts(numcpus);

#printlevel := 5;

# -----------------------------------------------------------------------------------------------------
# To discuss with Weifeng

f := (x-1129/50)*(x-2617/100);
basis := [(x-321/100)*(x-221/25)*(x-1041/50)*(x-2707/100), -(x-321/100)*(x-221/25)*(x-257/25)^2*(x-1041/50)*(x-2707/100)];
# -----------------------------------------------------------------------------------------------------

#checkRealCertify(f, basis, x, "Test 1");
checkWeifeng(f, basis, x, "Test 1");
