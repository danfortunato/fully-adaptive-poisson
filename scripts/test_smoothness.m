n = 16;
Gamma = Boundary.circle(n, 'quadrature', 'panel');
Gamma = refine(Gamma, 2);
analytical = 0;

if ( analytical )
    zz = cell2mat(Gamma.z);
    normals = cell2mat(Gamma.normal);
    cn = normals(:,1) + 1i*normals(:,2);
    zz1 = zz - 0.15*cn;
else
    [x, y] = smoothStrip(Gamma, n, 0.8);
    %xleg = chebvals2legvals(reshape(x, n, Gamma.np)); xleg = xleg(:);
    %yleg = chebvals2legvals(reshape(y, n, Gamma.np)); yleg = yleg(:);
    xleg = x;
    yleg = y;
    zz1 = xleg + 1i*yleg;
end
z1 = mat2cell(zz1, repmat(n, Gamma.np, 1), 1);
Gamma1 = Boundary(z1);


% Test if curve is globally smooth
m = n;
vals = zeros(m*Gamma1.np, 1);
for k = 1:Gamma1.np
    fun = Gamma1.f.funs{k};
    tt = trigpts(m, fun.domain);
    vals((k-1)*m+1:k*m) = feval(fun, tt);
end

gamfun = chebfun(vals, [0 2*pi], 'trig');
plotcoeffs(gamfun)
shg

%% Test if panels have locally resolved the curve
tol = 1e-8;
f = {Gamma1.f, @(t) abs(Gamma1.df(t)), Gamma1.dff};
%f = {Gamma1.f, @(t) abs(Gamma1.df(t))};
resolved = zeros(Gamma.np, 1);
bent     = zeros(Gamma.np, 1);
for k = 1:Gamma1.np
    
    a = Gamma1.breaks(k);
    b = Gamma1.breaks(k+1);

    %%%%%%%% resolved = isResolved(f, a, b, N, tol);
    m = 2*n;
    pts = util.gauss(n, a, b);
    xx = linspace(a, b, m).';
    R = util.interpmat(util.gauss(n), linspace(-1, 1, m));
    nf = numel(f);
    vals  = zeros(n, nf);
    tvals = zeros(m, nf);
    for l = 1:nf
        vals(:,l)  = f{l}(pts);
        tvals(:,l) = f{l}(xx);
    end
    % Interpolation error
    err = tvals - R * vals;
    relerr = err ./ sum(abs(tvals));
    [min(norm(relerr,inf), norm(err,inf)),  10*m*max(1e-14, tol)/(b-a)]
    resolved(k) = min(norm(relerr,inf), norm(err,inf)) < 10*m*max(1e-14, tol)/(b-a);

    %%%%%%%% bent = isBent(f, a, b, N, tol);
    [pts, w, D] = util.gauss(n, a, b);
    zp  = D   * f{1}(pts);
    zpp = D^2 * f{1}(pts);
    % Bending energy (normalized by 2*pi)
    bending = sum(w .* imag(zpp./zp).^2 ./ abs(zp)) * (b-a); % / (2*pi);
    bent(k) = ( bending > 5/log10(1/tol) );

end
