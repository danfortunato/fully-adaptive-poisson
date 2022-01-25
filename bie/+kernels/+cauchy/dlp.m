function u = dlp(source, target, density)
%DLP   Evaluate Cauchy's integral formula.

N = source.N;
M = size(target,2);

u = zeros(M,N);
for j = 1:N
    r = target - source.x(:,j);
    r = r(1,:)+r(2,:)*1i;
    u(:,j) = 1 ./ r / (2*pi*1i) * source.w(j);
end

if ( nargin == 3 )
    u = u * density;
end

end
