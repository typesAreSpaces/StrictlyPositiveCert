with(SolveTools, SemiAlgebraic);
with(StrictlyPositiveCert):

lprint(spCertificates((x + 1)*(x - 1), [-(x + 3)*(x + 2)*(x - 2)*(x - 3)], x));
