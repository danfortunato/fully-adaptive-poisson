n = 16;
Gamma = Boundary.star(n, 'quadrature', 'panel');
%Gamma = Boundary.circle(n, 'quadrature', 'panel');
D1 = kernels.laplace.dlp(Gamma);
D2 = kernels.laplace.LapDLP_nystrom(Gamma);
D3 = kernels.laplace.FMM_eval_nystrom(Gamma);
D_corr = kernels.laplace.dlp_correction(Gamma, 'i');

%% Chunkie
chnkr = Boundary.toChunker(Gamma);
fkern = @(s,t) -0.5/(2*pi) + 0*chnk.lap2d.kern(s,t,'D');
%opts = []; opts.quad = 'native';
D5 = chunkermat(chnkr, fkern);

%opts = []; opts.robust = false;
%D6 = chnk.quadadap.buildmat(chnkr, fkern, [1 1], 'log', opts);

%% Chunkie QuadJH

opdims = [1 1];
isclosed = 0;
logquad = chnk.quadggq.setuplogquad(n, opdims);
logquad.LogC = chnk.quadjh.setuplogquad(1, n, isclosed, opdims);
logquad.LogC0 = chnk.quadjh.setuplogquad([1, 0.5, 0.5], n, isclosed, opdims);
logquad.LogC1 = chnk.quadjh.setuplogquad([0.5, 0.5, 1], n, isclosed, opdims);

ilist = [];
[glnodes, glwts] = lege.exps(n);
fkern = @(s,t) chnk.lap2d.kern(s,t,'D');
D7 = chnk.quadjh.buildmat(chnkr, fkern, opdims, glwts, ilist, logquad);

%% Plotting stuff

M = 200;
box = boundingbox(Gamma);
[xx, yy] = meshgrid(linspace(box(1), box(2), M), ...
                    linspace(box(3), box(4), M));
xy = [xx(:) yy(:)];
ii = isinterior(Gamma, xx, yy);

%% BVP solve
D = D1;
I = eye(size(D));
A = D - I/2;

L = LaplaceSolver(Gamma, side='interior', bc='dirichlet', method='gmres', fmm=true);
u = L \ bc;

L = LaplaceSolver(Gamma, side='interior', bc='dirichlet');
sigma = L.op \ bc;
u = L.rep(sigma);
L.rep = @(sigma, varargin) kernels.laplace.dlp(Gamma, varargin{:}, density=sigma);

sol = @(x) real(exp(x(:,1)+x(:,2)*1i));
bc = sol(cell2mat(Gamma.x));

sigma = A \ bc;
v = nan(M);
v(ii) = kernels.laplace.dlp(Gamma, 'target',    [xx(ii) yy(ii)], ...
                                   'density',   sigma,           ...
                                   'closeeval', true,           ...
                                   'side',      'i');
v = reshape(v, size(xx));

w = nan(M);
dens = reshape(sigma, n, Gamma.np);
%opts = [];
%opts.flam = false;
fkern = @(s,t) chnk.lap2d.kern(s,t,'D');
w(ii) = chunkerkerneval(chnkr, fkern, dens, [xx(ii) yy(ii)].');
w = reshape(w, size(xx));

u = reshape(sol(xy), M, M);
u(~ii) = nan;

surf(xx, yy, log10(abs(u-v)))
shading interp
colorbar

norm(u(ii) - v(ii))