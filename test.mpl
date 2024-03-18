with(StrictlyPositiveCert):

nat := [x, x*(x-1), (x-1)*(x-2), -(x-2)];

f := -81/16*x^4 + 369/16*x^3 - 63/2*x^2 + 19/2*x + 5;

lprint(spCertificates(f, nat, x));
