doExport = false; %#ok<*UNRCH>

%% Make the geometry

align = 1;
amp   = 1;
delta = 1e-4;
freq  = 3;
gap   = 0.3;

% Exact solution
L = log(1/delta);
be = L*freq;
t1 = @(r) pi*(1 + log(r)/L);
mouthang = @(r) pi + 0.5*gap*amp/freq - amp/freq*sin(be*t1(r)).^2;
r = @(x,y) sqrt(x.^2 + y.^2);
real_scl = @(x,y) x.*cos(mouthang(r(x,y))) - y.*sin(mouthang(r(x,y)));
imag_scl = @(x,y) x.*sin(mouthang(r(x,y))) + y.*cos(mouthang(r(x,y)));
theta = @(x,y) atan2(imag_scl(x,y), real_scl(x,y)) - mouthang(r(x,y));
sol = @(x,y) 1/16*(x.^2+y.^2) .* (theta(x,y) + pi);
f = @(x,y) 1/16*4*(theta(x,y) + pi);
bc = @(x,y) sol(x,y);

n = 16;
geom_tol = 1e-12;
Gamma = Boundary.gobbler(n, quadrature='panel', align=align, amp=amp, ...
    delta=delta, freq=freq, gap=gap, tol=geom_tol);

% Further refinement...
Gamma = kdrefine(Gamma);
Gamma = levelrestrict(Gamma, fac=2, maxsize=0.5);

%% Generate the table

eps_ask = [1e-3 1e-6 1e-9];
ns = [4 7 10]; % Scale n with 1/eps

eps_get = zeros(size(eps_ask));
n_dof   = zeros(size(eps_ask));
mem_est = zeros(size(eps_ask));
t_setup = zeros(size(eps_ask));
t_solve = zeros(size(eps_ask));
t_homo  = zeros(size(eps_ask));
t_eval  = zeros(size(eps_ask));
t_total = zeros(size(eps_ask));
pps     = zeros(size(eps_ask));

t_evalself = zeros(size(eps_ask));
eps_get_self = zeros(size(eps_ask));
eps_get_bulk = zeros(size(eps_ask));
eps_get_strip = zeros(size(eps_ask));

% Warm up by running the first case twice
Timer.reset();
for k = [1 1:3]

    tol = eps_ask(k);
    n = ns(k);
    Gamma_k = resample(Gamma, n);

    % Initialize the solver
    opts = [];
    opts.debug = true;
    opts.n_sem = [n ceil(1.5*n)];
    opts.n_re  = ceil(1.5*n);
    opts.n_box = n;
    opts.beta = 1;
    opts.average = 5;
    opts.stripWidth = 0.5;
    opts.tol = tol;
    opts.boxToStripRatio = 0.5;

    % Initialize the solver
    Timer.tic();
    S = AdaptivePoissonSolver(Gamma_k, f, opts);
    t_setup(k) = Timer.toc('Initialize');

    % Solve
    Timer.tic();
    u = S.solve(bc);
    [~, spool] = Timer.toc('Solve');

    t_solve(k) = spool(1).toc + spool(5).toc + spool(6).toc;
    t_homo(k)  = spool(7).toc;

    % Evaluate
    nbulk = 1000;
    Timer.tic();
    dom = boundingbox(S.Gamma);
    [uu, xx, yy, inGamma, inGamma1, inStrip] = S.sampleUniform(u, nbulk, dom);
    t_eval(k) = Timer.toc('Eval');

    % Evaluate back on self
    Timer.tic();
    [u_bulk, u_strip, err_bulk, err_strip] = S.evaluatePotentials(u, sol);
    t_evalself(k) = Timer.toc('Evaluate back on self');

    % Calculate error
    uu_sol = nan(nbulk);
    uu_sol(inGamma) = sol(xx(inGamma),yy(inGamma));
    eps_get(k) = max(max(abs(uu - uu_sol)));

    eps_get_bulk(k) = err_bulk;
    eps_get_strip(k) = err_strip;
    eps_get_self(k) = max(err_bulk, err_strip);

    eps_get_self(k)

    n_dof(k) = numel(S.tf) + ...                          % Tree DOFs
               prod(opts.n_sem)*length(S.strip_dom) + ... % Strip DOFs
               2*opts.n_re*length(S.Gamma1) + ...         % Gamma' DOFs
               opts.n_re*length(S.Gamma);                 % Gamma DOFs

    mem_est(k) = 8*n_dof(k) / 1e9;

    t_total(k) = t_setup(k) + t_solve(k) + t_homo(k) + t_eval(k);

    pps(k) = n_dof(k) / t_total(k);

end

%% Save

save('table5_1.mat', 'eps_ask', 'ns', 'eps_get', 'n_dof', 'mem_est', ...
    't_setup', 't_solve', 't_homo', 't_eval', 't_total', 'pps');
