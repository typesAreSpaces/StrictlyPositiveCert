dir := FileTools:-JoinPath([currentdir(), "StrictlyPositiveCert.mla"]);
if not FileTools:-Exists(dir) then
    march('create', dir);
end if;
$include "src/StrictlyPositiveCert.mpl";
savelib('StrictlyPositiveCert', dir);
