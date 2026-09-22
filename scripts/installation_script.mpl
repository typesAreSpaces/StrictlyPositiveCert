dir := FileTools:-JoinPath([currentdir(), "StrictlyPositiveCert.mla"]);
march('create', dir);
$include "src/StrictlyPositiveCert.mpl";
savelib('StrictlyPositiveCert', dir);
