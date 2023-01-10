function D = FMM_eval_nystrom(source)

n = source.N;
np = source.np;
target = cell2mat(source.x);

D = zeros(n*np);
density = zeros(n*np, 1);

for k = 1:n*np
    density(k) = 1;
    D(:,k) = FMM_eval(source, target, 'dipstr', density);
    density(k) = 0;
end

kappa = cell2mat(source.curvature);
w = cell2mat(source.w);
D(1:n*np+1:end) = -kappa/(4*pi).*w;

end
