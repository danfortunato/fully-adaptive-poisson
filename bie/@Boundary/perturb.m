function [C, xyfun] = perturb(S, width, dom)
%PERTURB   Construct a Boundary by perturbing another.

if ( nargin == 2 )
    dom = [0 2*pi 0 1];
end

if ( ~(isa(width, 'function_handle') || isa(width, 'chebfun')) )
    width = @(t) width+0*t;
end

nx = normal( chebfun(S.f, [0 2*pi], 'trig') ); 
nx = nx(:,1)+nx(:,2)*1i;
nx = nx./abs(nx);
xyfun = chebfun2(@(t,r) S.f(t) + r.*width(t).*nx(t), dom);

newz = cell(S.np,1);
for k = 1:S.np
    newz{k} = xyfun(S.s{k}, 1);
    newz{k} = newz{k}(:);
end

if ( strcmpi(S.rule, 'ptr') )
    newz = newz{1};
end

C = Boundary(newz);
%xyfun = chebfun2(@(t,r) (1-r).*S.f(t) + r.*obj.f(t), dom);

end
