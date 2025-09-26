with(SolveTools, SemiAlgebraic);
with(StrictlyPositiveCert):

basis := [-(x^2+1)*(x+2)*(x+1)*(x-1)*(x-2)];
sigma:=-1/10*x^2+3/5+7369/67212*(-1/10*x^2+3/5)*(x^2+1)*(x+2)*(x+1)*(x-1)*(x-2);
tau:=-1/10+7369/67212*(-1/10*x^2+3/5)*(x^2+1)^2;
lprint(spCertificates(sigma, basis, x));
lprint(spCertificates(tau, basis, x));

