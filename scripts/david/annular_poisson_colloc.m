 %#ok<*UNRCH>
nb = 400;
ng = floor(nb/2);
M = 12;
interior = false;
chebkind = 2;

% Test problem
k = 2*pi/3;
sol_func = @(x,y) exp(sin(k*x)).*sin(k*y);
f_func = @(x,y) k^2*exp(sin(k*x)).*sin(k*y).*(cos(k*x).^2-sin(k*x)-1);

%% Construct boundary
%Gamma = Boundary.star(nb, 'wobble', 0.1, 'arms', 5);
Gamma = Boundary.circle(nb);

%% Construct grid
dom_global = [-1.5 1.5 -1.5 1.5];
[xx, yy] = meshgrid(trigpts(ng, dom_global(1:2)), ...
                    trigpts(ng, dom_global(3:4)));
xh = diff(dom_global(1:2))/ng;
yh = diff(dom_global(3:4))/ng;
h = xh*0.75;

%% Construct perturbed boundary
radial_width = M*h;
if ( interior )
    sgn = -1;
else
    sgn = 1;
end
[Gamma1, xyfun] = perturb(Gamma, sgn*radial_width);
xfun = real(xyfun);
yfun = imag(xyfun);

%% Radial grids
X = Gamma.x{1}(:,1);
Y = Gamma.x{1}(:,2);
NX = Gamma.normal{1}(:,1);
NY = Gamma.normal{1}(:,2);

% Generate the boundary fitted grid
if ( interior )
    lb = -radial_width;
    ub = 0;
else
    lb = 0;
    ub = radial_width;
end
radial_rv = chebpts(M, [lb ub], chebkind);
radial_tv = Gamma.s{1};
[radial_r, radial_t] = ndgrid(radial_rv, radial_tv);
radial_x = X.' + radial_r .* NX.';
radial_y = Y.' + radial_r .* NY.';

% Compute the approximate radius
centroid = mean(Gamma.x{1});
dxy = Gamma.x{1} - centroid;
rd = sqrt(sum(dxy.^2, 2));
approximate_radius = mean(rd);

%% Extract radial information from ebdy and construct annular solver

% Get the forces and BCs for the problem
f = f_func(radial_x, radial_y);
sol = sol_func(radial_x, radial_y);
ibc = sol_func(Gamma1.x{1}(:,1), Gamma1.x{1}(:,2));
ubc = sol_func( Gamma.x{1}(:,1),  Gamma.x{1}(:,2));
if ( ~interior )
    [ibc, ubc] = deal(ubc, ibc);
end

% Construct the preconditioner geometry
fprintf('Generating preconditioner\n')
AAG = ApproximateAnnularGeometry(nb, M, radial_width, approximate_radius, chebkind);

% Get the real geometry
if ( interior )
    sp = Gamma.speed{1};
    cur = Gamma.curvature{1};
else
    sp = Gamma1.speed{1};
    cur = Gamma1.curvature{1};
end
RAG = RealAnnularGeometry(sp, cur, M, radial_width, chebkind);

% Construct the solver
APS = AnnularPoissonSolver(RAG, AAG);

%% Solve
fprintf('Solving\n')
%u = APS.solve(f, ibc, ubc, 1e-14);
u = util.curvedPoisson(f_func, xfun, yfun, [M nb], 'periodic', ubc, ibc);

t = trigpts(nb, [0 2*pi]);
r = chebpts(M, [0 1]);
[tt,rr] = meshgrid(t,r);
u = u(tt,rr);

%% Compute error
sol = sol_func(xfun(tt,rr),yfun(tt,rr));
err_colloc = max(abs(u(:)-sol(:)));
fprintf('Error over collocation points: %g\n', err_colloc)

[tt, rr] = meshgrid(linspace(0, 2*pi, 1000), ...
                    linspace(0, 1,    100));
solf = chebfun2(sol, [0 2*pi 0 1], 'trigx');
uf = chebfun2(u, [0 2*pi 0 1], 'trigx');
ss = solf(tt,rr);
uu = uf(tt,rr);
err_uniform = max(abs(uu(:)-ss(:)));
fprintf('Error over uniform points: %g\n', err_uniform)

%% Plot the errors
clf
xf = chebfun2(radial_x, [0 2*pi 0 1], 'trigx');
yf = chebfun2(radial_y, [0 2*pi 0 1], 'trigx');
hold on
surf(xfun(tt,rr), yfun(tt,rr), log10(abs(uu-ss)))
%surf(xf(tt,rr), yf(tt,rr), uu)
plot(Gamma, 'k-', 'LineWidth', 1)
plot(Gamma1, 'k-', 'LineWidth', 1)
hold off
shading interp
axis tight equal
colorbar
shg
