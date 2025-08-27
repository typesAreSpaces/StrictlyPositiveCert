with(SolveTools, SemiAlgebraic);
with(StrictlyPositiveCert):

#eps := 1/2;

#for k in [13, 15, 17, 19, 21] do
  #basis := [(1-x^2)^k]:
  #st := time():
  #lprint(spCertificates(1 + eps + x, basis, x)):
#end do;

#eps := 1/2:
#lprint(">> eps = 1/2");

#for k in [13, 15, 17, 19, 21] do
  #basis := [(1-x^2)^k]:
  #st := time():
  #lprint(spCertificates(1 + eps + x, basis, x));
  #print("It took ", time() - st, " seconds");
#end do;

#eps := 1/3:
#lprint(">> eps = 1/3");

#for k in [13, 15, 17, 19, 21] do
  #basis := [(1-x^2)^k]:
  #st := time():
  #lprint(spCertificates(1 + eps + x, basis, x));
  #print("It took ", time() - st, " seconds");
#end do;

#spCertificates(-(x - 5), [(x + 1)*(x - 2)*(x - 4), -(x - 1)], x);

#sigma := 2509/179200*x^2-33297/313600-303/1254400*x^4;
#tau := 1/56-303/627200*(x+3)*(x-3);

#f := (x+3)*(x-3);
#g := -(1/2)*(-4 - 21*x^2 + x^4);

#lprint(">> Certificate of sigma");
#map(ok -> lprint(ok), spCertificates(sigma, [f*g], x));
#lprint(">> Certificate of tau");
#map(ok -> lprint(ok), spCertificates(tau, [f*g], x));

#f := (-3 + x)*(3 + x);
#g := -(1/216)*(-6 + x)*(6 + x)*(79 + 41*x^2);

#sigma := 41/12096*x^2-257/3024-66559/14631321600*(x-6)*(x+6)*(41*x^2+79)+763
#/165183300*(41/12096*x^2-257/3024-66559/14631321600*(x-6)*(x+6)*(41*x^2+79))*(x
#+3)*(x-3)*(x-6)*(x+6)*(41*x^2+79);
#tau := 1/56-66559/67737600*(x+3)*(x-3)+1526/1529475*(41/12096*x^2-257/3024-\
#66559/14631321600*(x-6)*(x+6)*(41*x^2+79))*(x+3)^2*(x-3)^2;

#lprint(">> Certificate of sigma");
#map(ok -> lprint(ok, ","), spCertificates(sigma, [f*g], x));
#lprint(">> Certificate of tau");
#map(ok -> lprint(ok, ","), spCertificates(tau, [f*g], x));

#g1 := x + 3;
#g2 := (x+2)*(x+1);
#g3 := (x-1)*(x-2);
#g4 := -(x - 3);

#sigma := -1/240*x + 1/40;
#tau := -1/240*x^3 + 11/40 + 9/80*x;

#lprint(">> Certificate of sigma");
#map(ok -> lprint(ok, ","), spCertificates(sigma, [g1*g2*g3*g4], x));
#lprint(">> Certificate of tau");
#map(ok -> lprint(ok, ","), spCertificates(tau, [g1*g2*g3*g4], x));

#sigma := -3/40*x+1/20+2657411/299935000*(x+2)*(x+1);
#tau := -3/40*x+11/40-2657411/299935000*(x+3)*(-x+3);

#lprint(">> Certificate of sigma");
#map(ok -> lprint(ok, ","), spCertificates(sigma, [g1*g2*g4], x));
#lprint(">> Certificate of tau");
#map(ok -> lprint(ok, ","), spCertificates(tau, [g1*g2*g4], x));

#lprint(spCertificates(2-x^2, [10-x^8], x));

lprint(spCertificates((x + 1)*(x-1), [-1], x));
