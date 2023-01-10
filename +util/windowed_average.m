function avg = windowed_average(f, k, type)

if ( nargin < 2 )
    k = 2;
end

if ( nargin < 3 )
    type = 'arith';
end

if ( k == 0 || k == 1)
    avg = f;
    return
end

if ( strcmpi(type, 'geom') )
    geom = true;
    op = @times;
else
    geom = false;
    op = @plus;
end

kl = floor(k/2);
kr = floor(k/2)-mod(k+1,2);

n = numel(f);
g = f(:);
g = [g(n-kl+1:n); g(:); g(1:kr)];
avg = g(1:n);
for j = 1:kl+kr
    avg = op(avg, g((1:n)+j));
end

if ( geom )
    avg = avg.^(1/k);
else
    avg = avg / k;
end
avg = reshape(avg, size(f));

end
