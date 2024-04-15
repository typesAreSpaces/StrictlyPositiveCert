with(SolveTools, SemiAlgebraic);
with(StrictlyPositiveCert):

basis := [-(x+3)*(x+2)*(x-2)*(x-3)*(x-13)*(x+13)];

# -----------------------------------------------------------------------------
sigma := 2967271/21632698560*x^3+4191089/4506812200*x^2-168521967/7210899520
*x+4978748893/54081746400-87/10242755*x^4-71/50123439*(2967271/21632698560*x^3+
4191089/4506812200*x^2-168521967/7210899520*x+4978748893/54081746400-87/
10242755*x^4)*(x+3)*(x+2)*(2-71/50123439*(x+3)*(x+2)*(-x^4+5*x^3+163*x^2-845*x+
1014)+(1-71/50123439*(x+3)*(x+2)*(-x^4+5*x^3+163*x^2-845*x+1014))^160+(1-71/
50123439*(x+3)*(x+2)*(-x^4+5*x^3+163*x^2-845*x+1014))^159+(1-71/50123439*(x+3)*
(x+2)*(-x^4+5*x^3+163*x^2-845*x+1014))^158+(1-71/50123439*(x+3)*(x+2)*(-x^4+5*x
^3+163*x^2-845*x+1014))^157+(1-71/50123439*(x+3)*(x+2)*(-x^4+5*x^3+163*x^2-845*
x+1014))^156+(1-71/50123439*(x+3)*(x+2)*(-x^4+5*x^3+163*x^2-845*x+1014))^155+(1
-71/50123439*(x+3)*(x+2)*(-x^4+5*x^3+163*x^2-845*x+1014))^154+(1-71/50123439*(x
+3)*(x+2)*(-x^4+5*x^3+163*x^2-845*x+1014))^153+(1-71/50123439*(x+3)*(x+2)*(-x^4
+5*x^3+163*x^2-845*x+1014))^152+(1-71/50123439*(x+3)*(x+2)*(-x^4+5*x^3+163*x^2-\
845*x+1014))^151+(1-71/50123439*(x+3)*(x+2)*(-x^4+5*x^3+163*x^2-845*x+1014))^2+
(1-71/50123439*(x+3)*(x+2)*(-x^4+5*x^3+163*x^2-845*x+1014))^3+(1-71/50123439*(x
+3)*(x+2)*(-x^4+5*x^3+163*x^2-845*x+1014))^4+(1-71/50123439*(x+3)*(x+2)*(-x^4+5
*x^3+163*x^2-845*x+1014))^5+(1-71/50123439*(x+3)*(x+2)*(-x^4+5*x^3+163*x^2-845*
x+1014))^6+(1-71/50123439*(x+3)*(x+2)*(-x^4+5*x^3+163*x^2-845*x+1014))^7+(1-71/
50123439*(x+3)*(x+2)*(-x^4+5*x^3+163*x^2-845*x+1014))^8+(1-71/50123439*(x+3)*(x
+2)*(-x^4+5*x^3+163*x^2-845*x+1014))^9+(1-71/50123439*(x+3)*(x+2)*(-x^4+5*x^3+
163*x^2-845*x+1014))^10+(1-71/50123439*(x+3)*(x+2)*(-x^4+5*x^3+163*x^2-845*x+
1014))^11+(1-71/50123439*(x+3)*(x+2)*(-x^4+5*x^3+163*x^2-845*x+1014))^12+(1-71/
50123439*(x+3)*(x+2)*(-x^4+5*x^3+163*x^2-845*x+1014))^13+(1-71/50123439*(x+3)*(
x+2)*(-x^4+5*x^3+163*x^2-845*x+1014))^14+(1-71/50123439*(x+3)*(x+2)*(-x^4+5*x^3
+163*x^2-845*x+1014))^15+(1-71/50123439*(x+3)*(x+2)*(-x^4+5*x^3+163*x^2-845*x+
1014))^16+(1-71/50123439*(x+3)*(x+2)*(-x^4+5*x^3+163*x^2-845*x+1014))^17+(1-71/
50123439*(x+3)*(x+2)*(-x^4+5*x^3+163*x^2-845*x+1014))^18+(1-71/50123439*(x+3)*(
x+2)*(-x^4+5*x^3+163*x^2-845*x+1014))^19+(1-71/50123439*(x+3)*(x+2)*(-x^4+5*x^3
+163*x^2-845*x+1014))^20+(1-71/50123439*(x+3)*(x+2)*(-x^4+5*x^3+163*x^2-845*x+
1014))^21+(1-71/50123439*(x+3)*(x+2)*(-x^4+5*x^3+163*x^2-845*x+1014))^22+(1-71/
50123439*(x+3)*(x+2)*(-x^4+5*x^3+163*x^2-845*x+1014))^23+(1-71/50123439*(x+3)*(
x+2)*(-x^4+5*x^3+163*x^2-845*x+1014))^24+(1-71/50123439*(x+3)*(x+2)*(-x^4+5*x^3
+163*x^2-845*x+1014))^25+(1-71/50123439*(x+3)*(x+2)*(-x^4+5*x^3+163*x^2-845*x+
1014))^26+(1-71/50123439*(x+3)*(x+2)*(-x^4+5*x^3+163*x^2-845*x+1014))^27+(1-71/
50123439*(x+3)*(x+2)*(-x^4+5*x^3+163*x^2-845*x+1014))^28+(1-71/50123439*(x+3)*(
x+2)*(-x^4+5*x^3+163*x^2-845*x+1014))^29+(1-71/50123439*(x+3)*(x+2)*(-x^4+5*x^3
+163*x^2-845*x+1014))^30+(1-71/50123439*(x+3)*(x+2)*(-x^4+5*x^3+163*x^2-845*x+
1014))^31+(1-71/50123439*(x+3)*(x+2)*(-x^4+5*x^3+163*x^2-845*x+1014))^32+(1-71/
50123439*(x+3)*(x+2)*(-x^4+5*x^3+163*x^2-845*x+1014))^33+(1-71/50123439*(x+3)*(
x+2)*(-x^4+5*x^3+163*x^2-845*x+1014))^34+(1-71/50123439*(x+3)*(x+2)*(-x^4+5*x^3
+163*x^2-845*x+1014))^35+(1-71/50123439*(x+3)*(x+2)*(-x^4+5*x^3+163*x^2-845*x+
1014))^36+(1-71/50123439*(x+3)*(x+2)*(-x^4+5*x^3+163*x^2-845*x+1014))^37+(1-71/
50123439*(x+3)*(x+2)*(-x^4+5*x^3+163*x^2-845*x+1014))^38+(1-71/50123439*(x+3)*(
x+2)*(-x^4+5*x^3+163*x^2-845*x+1014))^39+(1-71/50123439*(x+3)*(x+2)*(-x^4+5*x^3
+163*x^2-845*x+1014))^40+(1-71/50123439*(x+3)*(x+2)*(-x^4+5*x^3+163*x^2-845*x+
1014))^41+(1-71/50123439*(x+3)*(x+2)*(-x^4+5*x^3+163*x^2-845*x+1014))^42+(1-71/
50123439*(x+3)*(x+2)*(-x^4+5*x^3+163*x^2-845*x+1014))^43+(1-71/50123439*(x+3)*(
x+2)*(-x^4+5*x^3+163*x^2-845*x+1014))^44+(1-71/50123439*(x+3)*(x+2)*(-x^4+5*x^3
+163*x^2-845*x+1014))^45+(1-71/50123439*(x+3)*(x+2)*(-x^4+5*x^3+163*x^2-845*x+
1014))^46+(1-71/50123439*(x+3)*(x+2)*(-x^4+5*x^3+163*x^2-845*x+1014))^47+(1-71/
50123439*(x+3)*(x+2)*(-x^4+5*x^3+163*x^2-845*x+1014))^48+(1-71/50123439*(x+3)*(
x+2)*(-x^4+5*x^3+163*x^2-845*x+1014))^49+(1-71/50123439*(x+3)*(x+2)*(-x^4+5*x^3
+163*x^2-845*x+1014))^50+(1-71/50123439*(x+3)*(x+2)*(-x^4+5*x^3+163*x^2-845*x+
1014))^51+(1-71/50123439*(x+3)*(x+2)*(-x^4+5*x^3+163*x^2-845*x+1014))^52+(1-71/
50123439*(x+3)*(x+2)*(-x^4+5*x^3+163*x^2-845*x+1014))^53+(1-71/50123439*(x+3)*(
x+2)*(-x^4+5*x^3+163*x^2-845*x+1014))^54+(1-71/50123439*(x+3)*(x+2)*(-x^4+5*x^3
+163*x^2-845*x+1014))^55+(1-71/50123439*(x+3)*(x+2)*(-x^4+5*x^3+163*x^2-845*x+
1014))^56+(1-71/50123439*(x+3)*(x+2)*(-x^4+5*x^3+163*x^2-845*x+1014))^57+(1-71/
50123439*(x+3)*(x+2)*(-x^4+5*x^3+163*x^2-845*x+1014))^58+(1-71/50123439*(x+3)*(
x+2)*(-x^4+5*x^3+163*x^2-845*x+1014))^59+(1-71/50123439*(x+3)*(x+2)*(-x^4+5*x^3
+163*x^2-845*x+1014))^60+(1-71/50123439*(x+3)*(x+2)*(-x^4+5*x^3+163*x^2-845*x+
1014))^61+(1-71/50123439*(x+3)*(x+2)*(-x^4+5*x^3+163*x^2-845*x+1014))^62+(1-71/
50123439*(x+3)*(x+2)*(-x^4+5*x^3+163*x^2-845*x+1014))^63+(1-71/50123439*(x+3)*(
x+2)*(-x^4+5*x^3+163*x^2-845*x+1014))^64+(1-71/50123439*(x+3)*(x+2)*(-x^4+5*x^3
+163*x^2-845*x+1014))^65+(1-71/50123439*(x+3)*(x+2)*(-x^4+5*x^3+163*x^2-845*x+
1014))^66+(1-71/50123439*(x+3)*(x+2)*(-x^4+5*x^3+163*x^2-845*x+1014))^67+(1-71/
50123439*(x+3)*(x+2)*(-x^4+5*x^3+163*x^2-845*x+1014))^68+(1-71/50123439*(x+3)*(
x+2)*(-x^4+5*x^3+163*x^2-845*x+1014))^69+(1-71/50123439*(x+3)*(x+2)*(-x^4+5*x^3
+163*x^2-845*x+1014))^70+(1-71/50123439*(x+3)*(x+2)*(-x^4+5*x^3+163*x^2-845*x+
1014))^71+(1-71/50123439*(x+3)*(x+2)*(-x^4+5*x^3+163*x^2-845*x+1014))^72+(1-71/
50123439*(x+3)*(x+2)*(-x^4+5*x^3+163*x^2-845*x+1014))^73+(1-71/50123439*(x+3)*(
x+2)*(-x^4+5*x^3+163*x^2-845*x+1014))^74+(1-71/50123439*(x+3)*(x+2)*(-x^4+5*x^3
+163*x^2-845*x+1014))^75+(1-71/50123439*(x+3)*(x+2)*(-x^4+5*x^3+163*x^2-845*x+
1014))^76+(1-71/50123439*(x+3)*(x+2)*(-x^4+5*x^3+163*x^2-845*x+1014))^77+(1-71/
50123439*(x+3)*(x+2)*(-x^4+5*x^3+163*x^2-845*x+1014))^78+(1-71/50123439*(x+3)*(
x+2)*(-x^4+5*x^3+163*x^2-845*x+1014))^79+(1-71/50123439*(x+3)*(x+2)*(-x^4+5*x^3
+163*x^2-845*x+1014))^80+(1-71/50123439*(x+3)*(x+2)*(-x^4+5*x^3+163*x^2-845*x+
1014))^81+(1-71/50123439*(x+3)*(x+2)*(-x^4+5*x^3+163*x^2-845*x+1014))^82+(1-71/
50123439*(x+3)*(x+2)*(-x^4+5*x^3+163*x^2-845*x+1014))^83+(1-71/50123439*(x+3)*(
x+2)*(-x^4+5*x^3+163*x^2-845*x+1014))^84+(1-71/50123439*(x+3)*(x+2)*(-x^4+5*x^3
+163*x^2-845*x+1014))^85+(1-71/50123439*(x+3)*(x+2)*(-x^4+5*x^3+163*x^2-845*x+
1014))^86+(1-71/50123439*(x+3)*(x+2)*(-x^4+5*x^3+163*x^2-845*x+1014))^87+(1-71/
50123439*(x+3)*(x+2)*(-x^4+5*x^3+163*x^2-845*x+1014))^88+(1-71/50123439*(x+3)*(
x+2)*(-x^4+5*x^3+163*x^2-845*x+1014))^89+(1-71/50123439*(x+3)*(x+2)*(-x^4+5*x^3
+163*x^2-845*x+1014))^90+(1-71/50123439*(x+3)*(x+2)*(-x^4+5*x^3+163*x^2-845*x+
1014))^91+(1-71/50123439*(x+3)*(x+2)*(-x^4+5*x^3+163*x^2-845*x+1014))^92+(1-71/
50123439*(x+3)*(x+2)*(-x^4+5*x^3+163*x^2-845*x+1014))^93+(1-71/50123439*(x+3)*(
x+2)*(-x^4+5*x^3+163*x^2-845*x+1014))^94+(1-71/50123439*(x+3)*(x+2)*(-x^4+5*x^3
+163*x^2-845*x+1014))^95+(1-71/50123439*(x+3)*(x+2)*(-x^4+5*x^3+163*x^2-845*x+
1014))^96+(1-71/50123439*(x+3)*(x+2)*(-x^4+5*x^3+163*x^2-845*x+1014))^97+(1-71/
50123439*(x+3)*(x+2)*(-x^4+5*x^3+163*x^2-845*x+1014))^98+(1-71/50123439*(x+3)*(
x+2)*(-x^4+5*x^3+163*x^2-845*x+1014))^99+(1-71/50123439*(x+3)*(x+2)*(-x^4+5*x^3
+163*x^2-845*x+1014))^100+(1-71/50123439*(x+3)*(x+2)*(-x^4+5*x^3+163*x^2-845*x+
1014))^101+(1-71/50123439*(x+3)*(x+2)*(-x^4+5*x^3+163*x^2-845*x+1014))^102+(1-\
71/50123439*(x+3)*(x+2)*(-x^4+5*x^3+163*x^2-845*x+1014))^103+(1-71/50123439*(x+
3)*(x+2)*(-x^4+5*x^3+163*x^2-845*x+1014))^104+(1-71/50123439*(x+3)*(x+2)*(-x^4+
5*x^3+163*x^2-845*x+1014))^105+(1-71/50123439*(x+3)*(x+2)*(-x^4+5*x^3+163*x^2-\
845*x+1014))^106+(1-71/50123439*(x+3)*(x+2)*(-x^4+5*x^3+163*x^2-845*x+1014))^
107+(1-71/50123439*(x+3)*(x+2)*(-x^4+5*x^3+163*x^2-845*x+1014))^108+(1-71/
50123439*(x+3)*(x+2)*(-x^4+5*x^3+163*x^2-845*x+1014))^109+(1-71/50123439*(x+3)*
(x+2)*(-x^4+5*x^3+163*x^2-845*x+1014))^110+(1-71/50123439*(x+3)*(x+2)*(-x^4+5*x
^3+163*x^2-845*x+1014))^111+(1-71/50123439*(x+3)*(x+2)*(-x^4+5*x^3+163*x^2-845*
x+1014))^112+(1-71/50123439*(x+3)*(x+2)*(-x^4+5*x^3+163*x^2-845*x+1014))^113+(1
-71/50123439*(x+3)*(x+2)*(-x^4+5*x^3+163*x^2-845*x+1014))^114+(1-71/50123439*(x
+3)*(x+2)*(-x^4+5*x^3+163*x^2-845*x+1014))^115+(1-71/50123439*(x+3)*(x+2)*(-x^4
+5*x^3+163*x^2-845*x+1014))^116+(1-71/50123439*(x+3)*(x+2)*(-x^4+5*x^3+163*x^2-\
845*x+1014))^117+(1-71/50123439*(x+3)*(x+2)*(-x^4+5*x^3+163*x^2-845*x+1014))^
118+(1-71/50123439*(x+3)*(x+2)*(-x^4+5*x^3+163*x^2-845*x+1014))^119+(1-71/
50123439*(x+3)*(x+2)*(-x^4+5*x^3+163*x^2-845*x+1014))^120+(1-71/50123439*(x+3)*
(x+2)*(-x^4+5*x^3+163*x^2-845*x+1014))^121+(1-71/50123439*(x+3)*(x+2)*(-x^4+5*x
^3+163*x^2-845*x+1014))^122+(1-71/50123439*(x+3)*(x+2)*(-x^4+5*x^3+163*x^2-845*
x+1014))^123+(1-71/50123439*(x+3)*(x+2)*(-x^4+5*x^3+163*x^2-845*x+1014))^124+(1
-71/50123439*(x+3)*(x+2)*(-x^4+5*x^3+163*x^2-845*x+1014))^125+(1-71/50123439*(x
+3)*(x+2)*(-x^4+5*x^3+163*x^2-845*x+1014))^126+(1-71/50123439*(x+3)*(x+2)*(-x^4
+5*x^3+163*x^2-845*x+1014))^127+(1-71/50123439*(x+3)*(x+2)*(-x^4+5*x^3+163*x^2-\
845*x+1014))^128+(1-71/50123439*(x+3)*(x+2)*(-x^4+5*x^3+163*x^2-845*x+1014))^
129+(1-71/50123439*(x+3)*(x+2)*(-x^4+5*x^3+163*x^2-845*x+1014))^130+(1-71/
50123439*(x+3)*(x+2)*(-x^4+5*x^3+163*x^2-845*x+1014))^131+(1-71/50123439*(x+3)*
(x+2)*(-x^4+5*x^3+163*x^2-845*x+1014))^132+(1-71/50123439*(x+3)*(x+2)*(-x^4+5*x
^3+163*x^2-845*x+1014))^133+(1-71/50123439*(x+3)*(x+2)*(-x^4+5*x^3+163*x^2-845*
x+1014))^134+(1-71/50123439*(x+3)*(x+2)*(-x^4+5*x^3+163*x^2-845*x+1014))^135+(1
-71/50123439*(x+3)*(x+2)*(-x^4+5*x^3+163*x^2-845*x+1014))^136+(1-71/50123439*(x
+3)*(x+2)*(-x^4+5*x^3+163*x^2-845*x+1014))^137+(1-71/50123439*(x+3)*(x+2)*(-x^4
+5*x^3+163*x^2-845*x+1014))^138+(1-71/50123439*(x+3)*(x+2)*(-x^4+5*x^3+163*x^2-\
845*x+1014))^139+(1-71/50123439*(x+3)*(x+2)*(-x^4+5*x^3+163*x^2-845*x+1014))^
140+(1-71/50123439*(x+3)*(x+2)*(-x^4+5*x^3+163*x^2-845*x+1014))^141+(1-71/
50123439*(x+3)*(x+2)*(-x^4+5*x^3+163*x^2-845*x+1014))^142+(1-71/50123439*(x+3)*
(x+2)*(-x^4+5*x^3+163*x^2-845*x+1014))^143+(1-71/50123439*(x+3)*(x+2)*(-x^4+5*x
^3+163*x^2-845*x+1014))^144+(1-71/50123439*(x+3)*(x+2)*(-x^4+5*x^3+163*x^2-845*
x+1014))^145+(1-71/50123439*(x+3)*(x+2)*(-x^4+5*x^3+163*x^2-845*x+1014))^146+(1
-71/50123439*(x+3)*(x+2)*(-x^4+5*x^3+163*x^2-845*x+1014))^147+(1-71/50123439*(x
+3)*(x+2)*(-x^4+5*x^3+163*x^2-845*x+1014))^148+(1-71/50123439*(x+3)*(x+2)*(-x^4
+5*x^3+163*x^2-845*x+1014))^149+(1-71/50123439*(x+3)*(x+2)*(-x^4+5*x^3+163*x^2-\
845*x+1014))^150)*(-x^4+5*x^3+163*x^2-845*x+1014);

tau := 1/10560*x+13/26400-87/10242755*(x+3)*(x+2)+71/50123439*(2967271/
21632698560*x^3+4191089/4506812200*x^2-168521967/7210899520*x+4978748893/
54081746400-87/10242755*x^4)*(x+3)^2*(x+2)^2*(2-71/50123439*(x+3)*(x+2)*(-x^4+5
*x^3+163*x^2-845*x+1014)+(1-71/50123439*(x+3)*(x+2)*(-x^4+5*x^3+163*x^2-845*x+
1014))^160+(1-71/50123439*(x+3)*(x+2)*(-x^4+5*x^3+163*x^2-845*x+1014))^159+(1-\
71/50123439*(x+3)*(x+2)*(-x^4+5*x^3+163*x^2-845*x+1014))^158+(1-71/50123439*(x+
3)*(x+2)*(-x^4+5*x^3+163*x^2-845*x+1014))^157+(1-71/50123439*(x+3)*(x+2)*(-x^4+
5*x^3+163*x^2-845*x+1014))^156+(1-71/50123439*(x+3)*(x+2)*(-x^4+5*x^3+163*x^2-\
845*x+1014))^155+(1-71/50123439*(x+3)*(x+2)*(-x^4+5*x^3+163*x^2-845*x+1014))^
154+(1-71/50123439*(x+3)*(x+2)*(-x^4+5*x^3+163*x^2-845*x+1014))^153+(1-71/
50123439*(x+3)*(x+2)*(-x^4+5*x^3+163*x^2-845*x+1014))^152+(1-71/50123439*(x+3)*
(x+2)*(-x^4+5*x^3+163*x^2-845*x+1014))^151+(1-71/50123439*(x+3)*(x+2)*(-x^4+5*x
^3+163*x^2-845*x+1014))^2+(1-71/50123439*(x+3)*(x+2)*(-x^4+5*x^3+163*x^2-845*x+
1014))^3+(1-71/50123439*(x+3)*(x+2)*(-x^4+5*x^3+163*x^2-845*x+1014))^4+(1-71/
50123439*(x+3)*(x+2)*(-x^4+5*x^3+163*x^2-845*x+1014))^5+(1-71/50123439*(x+3)*(x
+2)*(-x^4+5*x^3+163*x^2-845*x+1014))^6+(1-71/50123439*(x+3)*(x+2)*(-x^4+5*x^3+
163*x^2-845*x+1014))^7+(1-71/50123439*(x+3)*(x+2)*(-x^4+5*x^3+163*x^2-845*x+
1014))^8+(1-71/50123439*(x+3)*(x+2)*(-x^4+5*x^3+163*x^2-845*x+1014))^9+(1-71/
50123439*(x+3)*(x+2)*(-x^4+5*x^3+163*x^2-845*x+1014))^10+(1-71/50123439*(x+3)*(
x+2)*(-x^4+5*x^3+163*x^2-845*x+1014))^11+(1-71/50123439*(x+3)*(x+2)*(-x^4+5*x^3
+163*x^2-845*x+1014))^12+(1-71/50123439*(x+3)*(x+2)*(-x^4+5*x^3+163*x^2-845*x+
1014))^13+(1-71/50123439*(x+3)*(x+2)*(-x^4+5*x^3+163*x^2-845*x+1014))^14+(1-71/
50123439*(x+3)*(x+2)*(-x^4+5*x^3+163*x^2-845*x+1014))^15+(1-71/50123439*(x+3)*(
x+2)*(-x^4+5*x^3+163*x^2-845*x+1014))^16+(1-71/50123439*(x+3)*(x+2)*(-x^4+5*x^3
+163*x^2-845*x+1014))^17+(1-71/50123439*(x+3)*(x+2)*(-x^4+5*x^3+163*x^2-845*x+
1014))^18+(1-71/50123439*(x+3)*(x+2)*(-x^4+5*x^3+163*x^2-845*x+1014))^19+(1-71/
50123439*(x+3)*(x+2)*(-x^4+5*x^3+163*x^2-845*x+1014))^20+(1-71/50123439*(x+3)*(
x+2)*(-x^4+5*x^3+163*x^2-845*x+1014))^21+(1-71/50123439*(x+3)*(x+2)*(-x^4+5*x^3
+163*x^2-845*x+1014))^22+(1-71/50123439*(x+3)*(x+2)*(-x^4+5*x^3+163*x^2-845*x+
1014))^23+(1-71/50123439*(x+3)*(x+2)*(-x^4+5*x^3+163*x^2-845*x+1014))^24+(1-71/
50123439*(x+3)*(x+2)*(-x^4+5*x^3+163*x^2-845*x+1014))^25+(1-71/50123439*(x+3)*(
x+2)*(-x^4+5*x^3+163*x^2-845*x+1014))^26+(1-71/50123439*(x+3)*(x+2)*(-x^4+5*x^3
+163*x^2-845*x+1014))^27+(1-71/50123439*(x+3)*(x+2)*(-x^4+5*x^3+163*x^2-845*x+
1014))^28+(1-71/50123439*(x+3)*(x+2)*(-x^4+5*x^3+163*x^2-845*x+1014))^29+(1-71/
50123439*(x+3)*(x+2)*(-x^4+5*x^3+163*x^2-845*x+1014))^30+(1-71/50123439*(x+3)*(
x+2)*(-x^4+5*x^3+163*x^2-845*x+1014))^31+(1-71/50123439*(x+3)*(x+2)*(-x^4+5*x^3
+163*x^2-845*x+1014))^32+(1-71/50123439*(x+3)*(x+2)*(-x^4+5*x^3+163*x^2-845*x+
1014))^33+(1-71/50123439*(x+3)*(x+2)*(-x^4+5*x^3+163*x^2-845*x+1014))^34+(1-71/
50123439*(x+3)*(x+2)*(-x^4+5*x^3+163*x^2-845*x+1014))^35+(1-71/50123439*(x+3)*(
x+2)*(-x^4+5*x^3+163*x^2-845*x+1014))^36+(1-71/50123439*(x+3)*(x+2)*(-x^4+5*x^3
+163*x^2-845*x+1014))^37+(1-71/50123439*(x+3)*(x+2)*(-x^4+5*x^3+163*x^2-845*x+
1014))^38+(1-71/50123439*(x+3)*(x+2)*(-x^4+5*x^3+163*x^2-845*x+1014))^39+(1-71/
50123439*(x+3)*(x+2)*(-x^4+5*x^3+163*x^2-845*x+1014))^40+(1-71/50123439*(x+3)*(
x+2)*(-x^4+5*x^3+163*x^2-845*x+1014))^41+(1-71/50123439*(x+3)*(x+2)*(-x^4+5*x^3
+163*x^2-845*x+1014))^42+(1-71/50123439*(x+3)*(x+2)*(-x^4+5*x^3+163*x^2-845*x+
1014))^43+(1-71/50123439*(x+3)*(x+2)*(-x^4+5*x^3+163*x^2-845*x+1014))^44+(1-71/
50123439*(x+3)*(x+2)*(-x^4+5*x^3+163*x^2-845*x+1014))^45+(1-71/50123439*(x+3)*(
x+2)*(-x^4+5*x^3+163*x^2-845*x+1014))^46+(1-71/50123439*(x+3)*(x+2)*(-x^4+5*x^3
+163*x^2-845*x+1014))^47+(1-71/50123439*(x+3)*(x+2)*(-x^4+5*x^3+163*x^2-845*x+
1014))^48+(1-71/50123439*(x+3)*(x+2)*(-x^4+5*x^3+163*x^2-845*x+1014))^49+(1-71/
50123439*(x+3)*(x+2)*(-x^4+5*x^3+163*x^2-845*x+1014))^50+(1-71/50123439*(x+3)*(
x+2)*(-x^4+5*x^3+163*x^2-845*x+1014))^51+(1-71/50123439*(x+3)*(x+2)*(-x^4+5*x^3
+163*x^2-845*x+1014))^52+(1-71/50123439*(x+3)*(x+2)*(-x^4+5*x^3+163*x^2-845*x+
1014))^53+(1-71/50123439*(x+3)*(x+2)*(-x^4+5*x^3+163*x^2-845*x+1014))^54+(1-71/
50123439*(x+3)*(x+2)*(-x^4+5*x^3+163*x^2-845*x+1014))^55+(1-71/50123439*(x+3)*(
x+2)*(-x^4+5*x^3+163*x^2-845*x+1014))^56+(1-71/50123439*(x+3)*(x+2)*(-x^4+5*x^3
+163*x^2-845*x+1014))^57+(1-71/50123439*(x+3)*(x+2)*(-x^4+5*x^3+163*x^2-845*x+
1014))^58+(1-71/50123439*(x+3)*(x+2)*(-x^4+5*x^3+163*x^2-845*x+1014))^59+(1-71/
50123439*(x+3)*(x+2)*(-x^4+5*x^3+163*x^2-845*x+1014))^60+(1-71/50123439*(x+3)*(
x+2)*(-x^4+5*x^3+163*x^2-845*x+1014))^61+(1-71/50123439*(x+3)*(x+2)*(-x^4+5*x^3
+163*x^2-845*x+1014))^62+(1-71/50123439*(x+3)*(x+2)*(-x^4+5*x^3+163*x^2-845*x+
1014))^63+(1-71/50123439*(x+3)*(x+2)*(-x^4+5*x^3+163*x^2-845*x+1014))^64+(1-71/
50123439*(x+3)*(x+2)*(-x^4+5*x^3+163*x^2-845*x+1014))^65+(1-71/50123439*(x+3)*(
x+2)*(-x^4+5*x^3+163*x^2-845*x+1014))^66+(1-71/50123439*(x+3)*(x+2)*(-x^4+5*x^3
+163*x^2-845*x+1014))^67+(1-71/50123439*(x+3)*(x+2)*(-x^4+5*x^3+163*x^2-845*x+
1014))^68+(1-71/50123439*(x+3)*(x+2)*(-x^4+5*x^3+163*x^2-845*x+1014))^69+(1-71/
50123439*(x+3)*(x+2)*(-x^4+5*x^3+163*x^2-845*x+1014))^70+(1-71/50123439*(x+3)*(
x+2)*(-x^4+5*x^3+163*x^2-845*x+1014))^71+(1-71/50123439*(x+3)*(x+2)*(-x^4+5*x^3
+163*x^2-845*x+1014))^72+(1-71/50123439*(x+3)*(x+2)*(-x^4+5*x^3+163*x^2-845*x+
1014))^73+(1-71/50123439*(x+3)*(x+2)*(-x^4+5*x^3+163*x^2-845*x+1014))^74+(1-71/
50123439*(x+3)*(x+2)*(-x^4+5*x^3+163*x^2-845*x+1014))^75+(1-71/50123439*(x+3)*(
x+2)*(-x^4+5*x^3+163*x^2-845*x+1014))^76+(1-71/50123439*(x+3)*(x+2)*(-x^4+5*x^3
+163*x^2-845*x+1014))^77+(1-71/50123439*(x+3)*(x+2)*(-x^4+5*x^3+163*x^2-845*x+
1014))^78+(1-71/50123439*(x+3)*(x+2)*(-x^4+5*x^3+163*x^2-845*x+1014))^79+(1-71/
50123439*(x+3)*(x+2)*(-x^4+5*x^3+163*x^2-845*x+1014))^80+(1-71/50123439*(x+3)*(
x+2)*(-x^4+5*x^3+163*x^2-845*x+1014))^81+(1-71/50123439*(x+3)*(x+2)*(-x^4+5*x^3
+163*x^2-845*x+1014))^82+(1-71/50123439*(x+3)*(x+2)*(-x^4+5*x^3+163*x^2-845*x+
1014))^83+(1-71/50123439*(x+3)*(x+2)*(-x^4+5*x^3+163*x^2-845*x+1014))^84+(1-71/
50123439*(x+3)*(x+2)*(-x^4+5*x^3+163*x^2-845*x+1014))^85+(1-71/50123439*(x+3)*(
x+2)*(-x^4+5*x^3+163*x^2-845*x+1014))^86+(1-71/50123439*(x+3)*(x+2)*(-x^4+5*x^3
+163*x^2-845*x+1014))^87+(1-71/50123439*(x+3)*(x+2)*(-x^4+5*x^3+163*x^2-845*x+
1014))^88+(1-71/50123439*(x+3)*(x+2)*(-x^4+5*x^3+163*x^2-845*x+1014))^89+(1-71/
50123439*(x+3)*(x+2)*(-x^4+5*x^3+163*x^2-845*x+1014))^90+(1-71/50123439*(x+3)*(
x+2)*(-x^4+5*x^3+163*x^2-845*x+1014))^91+(1-71/50123439*(x+3)*(x+2)*(-x^4+5*x^3
+163*x^2-845*x+1014))^92+(1-71/50123439*(x+3)*(x+2)*(-x^4+5*x^3+163*x^2-845*x+
1014))^93+(1-71/50123439*(x+3)*(x+2)*(-x^4+5*x^3+163*x^2-845*x+1014))^94+(1-71/
50123439*(x+3)*(x+2)*(-x^4+5*x^3+163*x^2-845*x+1014))^95+(1-71/50123439*(x+3)*(
x+2)*(-x^4+5*x^3+163*x^2-845*x+1014))^96+(1-71/50123439*(x+3)*(x+2)*(-x^4+5*x^3
+163*x^2-845*x+1014))^97+(1-71/50123439*(x+3)*(x+2)*(-x^4+5*x^3+163*x^2-845*x+
1014))^98+(1-71/50123439*(x+3)*(x+2)*(-x^4+5*x^3+163*x^2-845*x+1014))^99+(1-71/
50123439*(x+3)*(x+2)*(-x^4+5*x^3+163*x^2-845*x+1014))^100+(1-71/50123439*(x+3)*
(x+2)*(-x^4+5*x^3+163*x^2-845*x+1014))^101+(1-71/50123439*(x+3)*(x+2)*(-x^4+5*x
^3+163*x^2-845*x+1014))^102+(1-71/50123439*(x+3)*(x+2)*(-x^4+5*x^3+163*x^2-845*
x+1014))^103+(1-71/50123439*(x+3)*(x+2)*(-x^4+5*x^3+163*x^2-845*x+1014))^104+(1
-71/50123439*(x+3)*(x+2)*(-x^4+5*x^3+163*x^2-845*x+1014))^105+(1-71/50123439*(x
+3)*(x+2)*(-x^4+5*x^3+163*x^2-845*x+1014))^106+(1-71/50123439*(x+3)*(x+2)*(-x^4
+5*x^3+163*x^2-845*x+1014))^107+(1-71/50123439*(x+3)*(x+2)*(-x^4+5*x^3+163*x^2-\
845*x+1014))^108+(1-71/50123439*(x+3)*(x+2)*(-x^4+5*x^3+163*x^2-845*x+1014))^
109+(1-71/50123439*(x+3)*(x+2)*(-x^4+5*x^3+163*x^2-845*x+1014))^110+(1-71/
50123439*(x+3)*(x+2)*(-x^4+5*x^3+163*x^2-845*x+1014))^111+(1-71/50123439*(x+3)*
(x+2)*(-x^4+5*x^3+163*x^2-845*x+1014))^112+(1-71/50123439*(x+3)*(x+2)*(-x^4+5*x
^3+163*x^2-845*x+1014))^113+(1-71/50123439*(x+3)*(x+2)*(-x^4+5*x^3+163*x^2-845*
x+1014))^114+(1-71/50123439*(x+3)*(x+2)*(-x^4+5*x^3+163*x^2-845*x+1014))^115+(1
-71/50123439*(x+3)*(x+2)*(-x^4+5*x^3+163*x^2-845*x+1014))^116+(1-71/50123439*(x
+3)*(x+2)*(-x^4+5*x^3+163*x^2-845*x+1014))^117+(1-71/50123439*(x+3)*(x+2)*(-x^4
+5*x^3+163*x^2-845*x+1014))^118+(1-71/50123439*(x+3)*(x+2)*(-x^4+5*x^3+163*x^2-\
845*x+1014))^119+(1-71/50123439*(x+3)*(x+2)*(-x^4+5*x^3+163*x^2-845*x+1014))^
120+(1-71/50123439*(x+3)*(x+2)*(-x^4+5*x^3+163*x^2-845*x+1014))^121+(1-71/
50123439*(x+3)*(x+2)*(-x^4+5*x^3+163*x^2-845*x+1014))^122+(1-71/50123439*(x+3)*
(x+2)*(-x^4+5*x^3+163*x^2-845*x+1014))^123+(1-71/50123439*(x+3)*(x+2)*(-x^4+5*x
^3+163*x^2-845*x+1014))^124+(1-71/50123439*(x+3)*(x+2)*(-x^4+5*x^3+163*x^2-845*
x+1014))^125+(1-71/50123439*(x+3)*(x+2)*(-x^4+5*x^3+163*x^2-845*x+1014))^126+(1
-71/50123439*(x+3)*(x+2)*(-x^4+5*x^3+163*x^2-845*x+1014))^127+(1-71/50123439*(x
+3)*(x+2)*(-x^4+5*x^3+163*x^2-845*x+1014))^128+(1-71/50123439*(x+3)*(x+2)*(-x^4
+5*x^3+163*x^2-845*x+1014))^129+(1-71/50123439*(x+3)*(x+2)*(-x^4+5*x^3+163*x^2-\
845*x+1014))^130+(1-71/50123439*(x+3)*(x+2)*(-x^4+5*x^3+163*x^2-845*x+1014))^
131+(1-71/50123439*(x+3)*(x+2)*(-x^4+5*x^3+163*x^2-845*x+1014))^132+(1-71/
50123439*(x+3)*(x+2)*(-x^4+5*x^3+163*x^2-845*x+1014))^133+(1-71/50123439*(x+3)*
(x+2)*(-x^4+5*x^3+163*x^2-845*x+1014))^134+(1-71/50123439*(x+3)*(x+2)*(-x^4+5*x
^3+163*x^2-845*x+1014))^135+(1-71/50123439*(x+3)*(x+2)*(-x^4+5*x^3+163*x^2-845*
x+1014))^136+(1-71/50123439*(x+3)*(x+2)*(-x^4+5*x^3+163*x^2-845*x+1014))^137+(1
-71/50123439*(x+3)*(x+2)*(-x^4+5*x^3+163*x^2-845*x+1014))^138+(1-71/50123439*(x
+3)*(x+2)*(-x^4+5*x^3+163*x^2-845*x+1014))^139+(1-71/50123439*(x+3)*(x+2)*(-x^4
+5*x^3+163*x^2-845*x+1014))^140+(1-71/50123439*(x+3)*(x+2)*(-x^4+5*x^3+163*x^2-\
845*x+1014))^141+(1-71/50123439*(x+3)*(x+2)*(-x^4+5*x^3+163*x^2-845*x+1014))^
142+(1-71/50123439*(x+3)*(x+2)*(-x^4+5*x^3+163*x^2-845*x+1014))^143+(1-71/
50123439*(x+3)*(x+2)*(-x^4+5*x^3+163*x^2-845*x+1014))^144+(1-71/50123439*(x+3)*
(x+2)*(-x^4+5*x^3+163*x^2-845*x+1014))^145+(1-71/50123439*(x+3)*(x+2)*(-x^4+5*x
^3+163*x^2-845*x+1014))^146+(1-71/50123439*(x+3)*(x+2)*(-x^4+5*x^3+163*x^2-845*
x+1014))^147+(1-71/50123439*(x+3)*(x+2)*(-x^4+5*x^3+163*x^2-845*x+1014))^148+(1
-71/50123439*(x+3)*(x+2)*(-x^4+5*x^3+163*x^2-845*x+1014))^149+(1-71/50123439*(x
+3)*(x+2)*(-x^4+5*x^3+163*x^2-845*x+1014))^150);
# -----------------------------------------------------------------------------

decomposeSOSFactors := proc(f, x);
    local sqrfree_f , todo, i;
    local _poly, _multiplicity;

    local sos := 1, rest;

    sqrfree_f := sqrfree(f, x);
    todo := op(sqrfree_f)[2];
    rest := op(sqrfree_f)[1];

    for i from 1 to nops(todo) do
      _poly := todo[i][1];
      _multiplicity := todo[i][2];
      if type(_multiplicity, even) then
          sos := sos*(_poly)^_multiplicity;
      else
        if evalb(SemiAlgebraic([_poly < 0], [x]) = []) then
          sos := sos*(_poly)^_multiplicity;
        else
          rest := rest*(_poly)^_multiplicity;
        end if;
      end if;
    end do;

    return simplify(sos), simplify(rest);
end proc;

print(">> Start");

#sigma_sos, sigma_rest := decomposeSOSFactors(sigma, x);
#lprint(">> sigma_sos", sigma_sos);
#lprint(">> sigma_rest", sigma_rest);

#sigma_real := -(-9957497786 + 2527829505*x - 100586136*x^2 - 14836355*x^3 + 918720*x^4);
#lprint(spCertificates(sigma_real, basis, x));

#_tau := sqrfree(tau);
#lprint(_tau);
#lprint(spCertificates(tau, basis, x));


# TODO
Gfix := [x+3, (x+2)*(x-1), (x+1)*(x-2), -x+3, -(x+3)*(x-3)];
Gfix := [x+3, (x+2)*(x-1), (x+1)*(x-2), -x+3];
sigma := 1/100+35483/1125000*(x+2)*(x-2);
#lprint(spCertificates(sigma, Gfix, x));

#findEps := proc(basis, T, x)
lprint(findEps(Gfix, SemiAlgebraic([-sigma >= 0], [x]), x));
