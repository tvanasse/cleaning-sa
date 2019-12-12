% value = dec2basebyte(d, n) or dec2base256 converts the nonnegative integer d to base 256. 
% using n  produces a representation with at least n digits.
% BAR 11/12/14
function vector = dec2basebyte(d,n)

min_n = ceil(log(d)/log(256)); % covert number to log256

if nargin < 2
    n = min_n;  % default value produces minimum number of digits
end

if n < min_n
    error('The number of digits used to represent the number is not large enough')
end

vector = NaN(n,1);
for i = 1:n
    vector(i) = mod(floor(d/256^(i-1)),256);
end