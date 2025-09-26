with(SolveTools, SemiAlgebraic);
with(StrictlyPositiveCert):

Gfix := [x+3, (x+2)*(x-1), (x+1)*(x-2), -x+3, -(x+3)*(x-3)];
sigma := 1/100+35483/1125000*(x+2)*(x-2);
lprint(spCertificates(sigma, Gfix, x));

findEps := proc(basis, T, x)
Gfix := [x+3, (x+2)*(x-1), (x+1)*(x-2), -x+3];
lprint(findEps(Gfix, SemiAlgebraic([-sigma >= 0], [x]), x));

