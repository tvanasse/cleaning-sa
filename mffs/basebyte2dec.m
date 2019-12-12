% converts a vector x in base256 to a decimal value
% BAR 11/12/14
function d = basebyte2dec(vector)

digits = length(vector);
bases  = 256.^(0:digits-1)';

d = sum(vector.*bases);

