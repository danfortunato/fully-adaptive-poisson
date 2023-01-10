function pass = test_laplace_bvp()
% Solve the Laplace equation on interior/exterior domains

tol = 1e-11;
pass = [];

%% First let's test using the periodic trapezoid rule
N = 200; % Number of quadrature points
S = Boundary.star(N);

% Evaluate on an M x M grid
M = 100;
L = max(abs(cat(1, S.z{:})));
[xx, yy] = meshgrid(linspace(-L, L, M));
ii = isinterior(S, xx, yy);

% Construct a solution
sol = @(x,y) real(exp(x+y*1i)); % Real part of holomorphic function
bxy = cell2mat(S.x);
bc = sol(bxy(:,1), bxy(:,2));

% Interior Laplace Dirichlet BVP solver using DLP
I = eye(N);
K = kernels.laplace.dlp(S);
sigma = (K - I/2) \ bc;
u = nan(M);
u(ii) = kernels.laplace.dlp(S, 'target',    [xx(ii) yy(ii)], ...
                               'density',   sigma,           ...
                               'closeeval', true,            ...
                               'side',      'i');
pass(1) = max(max(abs(u(ii) - sol(xx(ii),yy(ii))))) < tol;

% Interior Laplace Dirichlet BVP solver using DLP + SLP
I = eye(N);
K = kernels.laplace.dlp(S) + kernels.laplace.slp(S);
sigma = (K - I/2) \ bc;
u = nan(M);
u(ii) = kernels.laplace.dlp(S, 'target',  [xx(ii) yy(ii)], ...
                               'density', sigma, 'closeeval', true, 'side', 'i') + ...
        kernels.laplace.slp(S, 'target',  [xx(ii) yy(ii)], ...
                               'density', sigma, 'closeeval', true, 'side', 'i');
pass(2) = max(max(abs(u(ii) - sol(xx(ii),yy(ii))))) < tol;

% Exterior Laplace Dirichlet BVP solver using DLP
% Construct a solution
sol = @(x,y) real(1./(x+y*1i-0.1-0.3i)) + 2; % Point source
bxy = cell2mat(S.x);
bc = sol(bxy(:,1), bxy(:,2));
I = eye(N);
K = kernels.laplace.dlp(S, 'modified', true);
sigma = (K + I/2) \ bc;
u = nan(M);
u(~ii) = kernels.laplace.dlp(S, 'target',    [xx(~ii) yy(~ii)], ...
                                'density',   sigma,             ...
                                'modified',  true,              ...
                                'closeeval', true,              ...
                                'side',      'e');
pass(3) = max(max(abs(u(~ii) - sol(xx(~ii),yy(~ii))))) < tol;

%% Now let's test using panels
p = 16; % Number of quadrature points per panel
S = Boundary.star(p, 'quadrature', 'panel');
N  = p*S.np; % Total number of quadrature points

% Evaluate on an M x M grid
M = 100;
L = max(abs(cat(1, S.z{:})));
[xx, yy] = meshgrid(linspace(-L, L, M));
ii = isinterior(S, xx, yy);

% Construct a solution
sol = @(x,y) real(exp(x+y*1i)); % Real part of holomorphic function
bxy = cell2mat(S.x);
bc = sol(bxy(:,1), bxy(:,2));

% Interior Laplace Dirichlet BVP solver using DLP
I = eye(N);
K = kernels.laplace.dlp(S);
sigma = (K - I/2) \ bc;
u = nan(M);
u(ii) = kernels.laplace.dlp(S, 'target',    [xx(ii) yy(ii)], ...
                               'density',   sigma,           ...
                               'closeeval', true,            ...
                               'side',      'i');
pass(4) = max(max(abs(u(ii) - sol(xx(ii),yy(ii))))) < tol;

%%
% Exterior Laplace Dirichlet BVP solver using DLP
% Construct a solution
sol = @(x,y) real(1./(x+y*1i-0.1-0.3i)) + 2; % Point source
bxy = cell2mat(S.x);
bc = sol(bxy(:,1), bxy(:,2));
I = eye(N);
K = kernels.laplace.dlp(S, 'modified', true);
sigma = (K + I/2) \ bc;
u = nan(M);
u(~ii) = kernels.laplace.dlp(S, 'target',    [xx(~ii) yy(~ii)], ...
                                'density',   sigma,             ...
                                'modified',  true,              ...
                                'closeeval', true,              ...
                                'side',      'e');
pass(5) = max(max(abs(u(~ii) - sol(xx(~ii),yy(~ii))))) < tol;

end
