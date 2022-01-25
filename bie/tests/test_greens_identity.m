function pass = test_greens_identity()

tol = 1e-11;
pass = [];

% Let u solve homog PDE (Lap eq) in Omega.
% We test interior GRF formula: u = S u_n^- - D u^-   in Omega, and surface lim

fholom = @(z) z.^4 + exp(z);  fpholom = @(z) 4*z.^3 + exp(z);  % holomorphic
u = @(z) real(fholom(z));                                 % harmonic test func
ux = @(z) real(fpholom(z)); uy = @(z) -imag(fpholom(z));  % partials, sign!

%% First let's test using the periodic trapezoid rule
N = 200; % Number of quadrature points
S = Boundary.star(N);

% get bdry data u^-, u_n^- for test solution ...
sx  = cell2mat(S.x);      sx  = sx(:,1)  + sx(:,2)*1i;
snx = cell2mat(S.normal); snx = snx(:,1) + snx(:,2)*1i;
ub = u(sx);
unb = real(snx).*ux(sx) + imag(snx).*uy(sx);    % normal dot grad u

% ========= test far interior point via plain native quadrature eval...
tx = 0.5+0.2i; % targ: far inside pt (so plain quadr good)
Sunb = kernels.laplace.slp(S, 'target', tx, 'density', unb, 'side', 'i');
Dub  = kernels.laplace.dlp(S, 'target', tx, 'density', ub,  'side', 'i');
vt = Sunb - Dub;
pass(1) = abs(vt - u(tx)) < tol;

% ========= test close interior point
s0 = 0.3;
dist = 1e-3;
tx = S.f(s0) - dist*(-1i*S.df(s0)/abs(S.df(s0)));
Sunb = kernels.laplace.slp(S, 'target', tx, 'density', unb, 'side', 'i', 'closeeval', true);
Dub  = kernels.laplace.dlp(S, 'target', tx, 'density', ub,  'side', 'i', 'closeeval', true);
vt = Sunb - Dub;
pass(2) = abs(vt - u(tx)) < tol;

%% Now let's test using panels
p = 16;
S = Boundary.star(p, 'quadrature', 'panel');
S = refine(S);

% get bdry data u^-, u_n^- for test solution ...
sx  = cell2mat(S.x);      sx  = sx(:,1)  + sx(:,2)*1i;
snx = cell2mat(S.normal); snx = snx(:,1) + snx(:,2)*1i;
ub = u(sx);
unb = real(snx).*ux(sx) + imag(snx).*uy(sx);    % normal dot grad u

% ========= test far interior point via plain native quadrature eval...
tx = 0.5+0.2i; % targ: far inside pt (so plain quadr good)
Sunb = kernels.laplace.slp(S, 'target', tx, 'density', unb, 'side', 'i');
Dub  = kernels.laplace.dlp(S, 'target', tx, 'density', ub,  'side', 'i');
vt = Sunb - Dub;
pass(3) = abs(vt - u(tx)) < tol;
fprintf('GRF far int pt err = (%.3g)\n',abs(vt - u(tx)))

% ========= test close interior point
s0 = 0.3;
dist = 1e-3;
tx = S.f(s0) - dist*(-1i*S.df(s0)/abs(S.df(s0)));
Sunb = kernels.laplace.slp(S, 'target', tx, 'density', unb, 'side', 'i', 'closeeval', true);
Dub  = kernels.laplace.dlp(S, 'target', tx, 'density', ub,  'side', 'i', 'closeeval', true);
vt = Sunb - Dub;
pass(4) = abs(vt - u(tx)) < tol;
fprintf('GRF close int pt err = (%.3g)\n',abs(vt - u(tx)))

end
