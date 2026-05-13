dir := FileTools:-JoinPath([currentdir(), "/StrictlyPositiveCert.mla"]);
march('create', dir);
read "src/StrictlyPositiveCert.mpl";
savelib('StrictlyPositiveCert', dir);
