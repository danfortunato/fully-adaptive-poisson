function A = DLP_nystrom(source)
%DLP_NYSTROM   Build Nystrom matrix for double-layer potential.

% Unpack things
N     = source.N * source.np;
x     = cell2mat(source.x);
w     = cell2mat(source.w);
nx    = cell2mat(source.normal);
kappa = cell2mat(source.curvature);

% Fill in double-layer potential
A = zeros(N);
for j = 1:N
    r = x - x(j,:);
    A(:,j) = r./sum(r.^2,2)*nx(j,:).'/(2*pi)*w(j);
end
% Limiting value when source == target
A(1:N+1:end) = -kappa/(4*pi).*w;
%[rows, cols] = find(isnan(A));
%A(isnan(A)) = -kappa(cols)/(4*pi).*w(cols);

end
