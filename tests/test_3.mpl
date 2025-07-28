with(SolveTools, SemiAlgebraic);
with(StrictlyPositiveCert):

basis := [-(x+2)*(x+1)*x^2*(x-1)*(x-2) + 4*x^2*(x+1)^2];
sigma := -13/1536*x + 53/1536 - 901/248984*(-13/1536*x + 53/1536)*(x + 1)*x^3*(-x^2 + x + 8);
lprint(spCertificates(sigma, basis, x));

tau := -13/1536*x^3 + 9/512*x^2 - 1/64*x + 1/8 + 901/248984*(-13/1536*x + 53/1536)*(x + 1)^2*x^6;
lprint(spCertificates(tau, basis, x));

basis := [-(x+3)*(x+2)*(x+1)*(x-1)*(x-2)*(x-3)];
lprint(spCertificates(66 - 27*x + x^3, basis, x));
lprint(spCertificates((x+6)*(66 - 27*x + x^3), basis, x));
lprint(spCertificates(x+6, basis, x));
lprint(spCertificates((x-13/10)*(x-17/10), basis, x));
lprint(spCertificates((x+6)*(x-13/10)*(x-17/10), basis, x));
lprint(spCertificates((x+6)*(x-13/10)*(x-17/10)*(-(x-6)), basis, x));
lprint(spCertificates((x+6)*(x-13/10)*(x-17/10)*(-(x-6))*(x-7)*(x-8), basis, x));
