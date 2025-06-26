%% FFT2 scaling

ns_fft = 2.^(4:16);
hs_fft = 1./ns_fft;
ts_fft = zeros(size(hs_fft));

for k = [1 1:length(hs_fft)]
    n = ns_fft(k);
    h = hs_fft(k);
    fprintf('n = %d, h = %g\n', n, h);
    x = rand(n);
    t = tic;
    tmp = fft2(x);
    ts_fft(k) = toc(t);
    clear x tmp
end

etas_fft = 10*hs_fft;
dof_fft = 2*ns_fft.^2;
mem_fft = dof_fft * 8 * 2 / 1e9;

%% Adaptive scaling

% Exponentially clustered Gaussians for reference solution on the raindrop
ngauss = 30;
ref = 1;
eta = 1e-8;
ymax = 0;
ymin = -0.8;

cy = logspace(log10(eta), 0, ceil(ref*ngauss));
cy = ymax - cy*(ymax-ymin);
cx = zeros(size(cy));
s = -1e2./logspace(log10(eta), 0, ceil(ref*ngauss));
s = -0.003*s.^2;
cx = reshape(cx, 1, 1, []);
cy = reshape(cy, 1, 1, []);
s  = reshape(s,  1, 1, []);
bg = @(x,y) 2*cos(3*pi*y).*x;
sol = @(x,y) sum(exp(s.*((x-cx).^2+(y-cy).^2)), 3) + bg(x,y);
syms x y
f = matlabFunction(laplacian(sol(x,y)), 'vars', [x y]);
f = @(x,y) f(x,y) + 0*x;
bc = @(x,y) sol(x,y);

n = 16;
nbulk = 400;

opts = [];
opts.debug = true;
opts.n_re  = 24;
opts.n_sem = [16 24];
opts.n_box = 16;
opts.beta = 1;
opts.stripWidth = 0.5;
opts.tol = 1e-10;
opts.boxToStripRatio = 0.5;

solve_opts = [];
solve_opts.tol = opts.tol;

alpha = 0.5;
etas = [1e0, 1e-1, 1e-2, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8];
dof       = zeros(length(etas), 1);
err       = zeros(length(etas), 1);
t_init    = zeros(length(etas), 1);
t_solve   = zeros(length(etas), 1);
t_eval    = zeros(length(etas), 1);
t_total   = zeros(length(etas), 1);
mem_init  = zeros(length(etas), 1);
mem_solve = zeros(length(etas), 1);
mem_total = zeros(length(etas), 1);
Timer.reset();

% Do it once to warm up
j = 1; k = 1;
Gamma = Boundary.raindrop2(n, quadrature='panel', eta=etas(k), alpha=alpha, tol=opts.tol);
Gamma = levelrestrict(Gamma, fac=2.1);
S = AdaptivePoissonSolver(Gamma, f, opts);
u = S.solve(bc, solve_opts);
[uu, xx, yy, inGamma, inGamma1, inStrip] = S.sampleUniform(u, nbulk);

for k = 1:length(etas)

    fprintf('k=%d\n', k);

    Gamma = Boundary.raindrop2(n, quadrature='panel', eta=etas(k), alpha=alpha, tol=opts.tol);
    Gamma = levelrestrict(Gamma, fac=2.1);

    memtic;
    Timer.tic();
    S = AdaptivePoissonSolver(Gamma, f, opts);
    t_init(k,j) = Timer.toc('Initialization');
    mem_init(k,j) = memtoc;

    dof(k,j) = numel(S.tf);

    memtic;
    Timer.tic();
    u = S.solve(bc, solve_opts);
    t_solve(k,j) = Timer.toc('Solve');
    mem_solve(k,j) = memtoc;

    t_total(k,j) = t_init(k,j) + t_solve(k,j);
    mem_total(k,j) = mem_init(k,j) + mem_solve(k,j);

    Timer.tic();
    [uu, xx, yy, inGamma, inGamma1, inStrip] = S.sampleUniform(u, nbulk);
    t_eval(k,j) = Timer.toc('Evaluation');

    uu_sol = nan(nbulk);
    uu_sol(inGamma) = sol(xx(inGamma),yy(inGamma));

    err(k,j) = max(max(abs(uu - uu_sol)));
    err(k,j)

end

etas_adap = etas;
ts_adap = t_total;
mem_adap = dof * 8 * 2 / 1e9;
file = 'figure5_2.mat';
vars = {'etas_fft', 'ts_fft', 'mem_fft', ...
        'etas_adap', 'ts_adap', 'mem_adap'};