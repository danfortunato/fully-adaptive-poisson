function [step, bump] = makeMask(n, width, type)

if ( nargin < 1 )
    n = 1000;
end

if ( nargin < 2 )
    width = 30;
end

y = util.slepian(n, width/n);
bump = chebfun(y, [0 1], 'equi');
step = cumsum(bump) / sum(bump);

end
