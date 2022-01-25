n = 16;
Gamma = Boundary.circle(n, 'quadrature', 'panel');
%Gamma = Boundary.star(n, 'quadrature', 'panel');
%Gamma = refine(Gamma);
%Gamma = Boundary.multiscale_circle(n, 'quadrature', 'panel');
%Gamma = Boundary.gobbler(n, 'delta', 1e-1, 'align', -1, 'quadrature', 'panel');

syms x y
sol = @(x,y) -cos(5*x).*sin(5*y);
f = matlabFunction(simplify(laplacian(sol(x,y))));
bc = @(x,y) sol(x,y);

S = AdaptivePoissonSolver(Gamma, f);
u = S.solve(bc);

%% Plot
nbulk = 400;
[uu, xx, yy, inGamma] = S.sampleUniform(u, nbulk);

uu_sol = nan(nbulk);
uu_sol(inGamma) = sol(xx(inGamma),yy(inGamma));

max(max(abs(uu - uu_sol)))

figure(1)
%surf(xx, yy, ug - uu_sol1)
%surf(xx, yy, uu_sol1)
%surf(xx, yy, uu)
%surf(xx, yy, uu-uu_sol)
surf(xx, yy, log10(abs(uu-uu_sol)))
shading interp
shg

%% Error in strip
uu = cell(length(S.strip_dom), 1);
err = zeros(length(S.strip_dom), 1);
samp = linspace(-1, 1, S.n_sem+2).';
samp([1 end]) = [];
[tt, rr] = meshgrid(samp);
for k = 1:length(S.strip_dom)
    xfun = chebfun2(S.strip_dom(k).x);
    yfun = chebfun2(S.strip_dom(k).y);
    xx = xfun(tt,rr);
    yy = yfun(tt,rr);
    stripvals = util.bary2d(u.strip{k}, tt, rr);
    u2 = reshape(u.glue_e(xx(:), yy(:)), S.n_sem, S.n_sem);
    u3 = reshape(u.bc(xx(:), yy(:)), S.n_sem, S.n_sem);
    uu = stripvals + u2 + u3;
    err(k) = max(max(abs(uu - (sol(xx,yy)))));
    figure(1)
    view(3)
    hold on
    surf(xx, yy, log10(abs(uu - (sol(xx,yy))))), hold on
    shading interp
    drawnow
end
fprintf('   Error in strip = %g\n', max(err));