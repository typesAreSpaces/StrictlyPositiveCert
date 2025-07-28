with(SolveTools, SemiAlgebraic);
with(StrictlyPositiveCert):

basis := [-(x+1)^2*(x-1)^2];
lprint(spCertificates((x+1/2)*(x-1/2), basis, x));
