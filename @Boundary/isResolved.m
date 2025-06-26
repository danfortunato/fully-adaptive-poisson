function resolved = isResolved(Gamma, f, tol)

n = Gamma.N;
xy = cat(3, Gamma.x{:});
xx = squeeze(xy(:,1,:));
yy = squeeze(xy(:,2,:));
ff = feval(f, xx, yy);
hx = max(xx) - min(xx);
hy = max(yy) - min(yy);
hh = min(hx, hy);

V2C = legvals2legcoeffs(eye(n));

eta = 0;
resolved = false(length(Gamma), 1);
for k = 1:length(Gamma)
    vals = ff(:,k);
    coeffs = V2C * vals;
    coeffs = coeffs(1:n);
    err_cfs = sum(abs(coeffs(end-1:end))) / 2;
    h = hh(k);
    vmax = max(abs(vals(:)));
    resolved(k) = ( err_cfs * h^eta < tol * max(vmax, 1) );
end

end