restart:
# Search the experimental archive before any installed copy of the package.
libname := "../StrictlyPositiveCert.mla", libname:
with(StrictlyPositiveCert);

interface(echo=0);
getenv(OMP_NUM_THREADS);
kernelopts(numcpus);

printlevel := -1;
#printlevel := 8;

#spCertificates(x+2,[x^2,-x^5,x+1],x);
spCertificates(x+3/2,[1, 1-x^2],x);
