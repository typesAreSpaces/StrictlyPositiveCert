with(SolveTools, SemiAlgebraic);
with(StrictlyPositiveCert):

g1:=-(x+2)*(x+1)^3*x^2*(x-1)*(x-2);
g1_1:=-(x+2)*(x+1)^3*(x-1)*(x-2);
g2:=(x+1)*(x+1/2);

basis := [g1];
lprint(spCertificates(-(x+3)*(x-3), basis, x));

basis := [g1, g2];
lprint(spCertificates((x+3/4)*(x+1/4), basis, x));
lprint(spCertificates((x+5/8)*(x-1/2), [-(x+2)*(x+1)^3*(x-1)*(x-2)], x));

lprint(spCertificates(-(x-3/4), [-(x+2)*(x+1)^3*x^2], x));

lprint(spCertificates((x-1/2)*(x-3/4), [-(x+2)*(x+1)^3*x^2*(x-1)*(x-2)], x));
lprint(spCertificates((x-1/8)*(x-3/4), [-(x+2)*(x+1)^3*x^2*(x-1)*(x-2)], x));
