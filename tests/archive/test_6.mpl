with(SolveTools, SemiAlgebraic);
with(StrictlyPositiveCert):

basis := [-(x+3)*(x-3)*(x+1)*(x+1/2)];
lprint(spCertificates(-3/140*x+29/280, basis, x));
lprint(spCertificates(-3/140*x+19/140, basis, x));

