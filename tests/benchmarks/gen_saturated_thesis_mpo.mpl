local n, k, _k, external_benchmarks ; 

external_benchmarks := []; 

_k := 10;

for k from 1 to _k do
  #external_benchmarks := 
  #[op(external_benchmarks), [x + k + 1, -mul(x^2-n^2, n=1..k)]]:
  external_benchmarks := 
  [op(external_benchmarks), [-x + k + 1, -mul(x^2-n^2, n=1..k)]]:
end do;

lprint(">> benchmark");
lprint(external_benchmarks);

