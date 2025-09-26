$define QUICKDEGREE

$ifdef QUICKDEGREE
_quickdegree := proc(f, x);
if type(f, `+`) then
  return max(map(_f -> _quickdegree(_f, x), convert(f, list, `+`)));
end if;

if type(f, `*`) then
  return foldl((x1, x2) -> x1 + x2, 
    0, 
    op(map(_f -> _quickdegree(_f, x), convert(f, list, `*`))));
end if;

if type(f, `^`) then
  local _f, _exp;
  _f, _exp := op(f);
  return _exp*_quickdegree(_f, x);
end if;

if f = x then
  return 1;
else
  return 0;
end if;

end proc;
$endif

quickdegree := proc(f, x)
local output;
$ifdef QUICKDEGREE
output := _quickdegree(f, x);
$elsedef
output := degree(expand(f));
$endif
#print(">> Degree", output);
return output;
end proc;
