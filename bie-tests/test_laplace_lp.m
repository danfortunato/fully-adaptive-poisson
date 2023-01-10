function pass = test_laplace_lp()
% Test that the Laplace double-layer potential with density sigma = 1
% evaluates to -1 inside the domain and 0 outside the domain.

tol = 1e-12;
pass = [];

%% First let's test using the periodic trapezoid rule
N = 100; % Number of quadrature points
S = Boundary.star(N);
sigma = ones(N, 1);

% Evaluate on an M x M grid
M = 100;
L = max(abs(cat(1, S.z{:})));
[xx, yy] = meshgrid(linspace(-L, L, M));
xy = [xx(:) yy(:)];
ii = isinterior(S, xx, yy);

% Interior Laplace
sol = -1;
u = kernels.laplace.dlp(S, 'target',    xy,    ...
                           'density',   sigma, ...
                           'closeeval', true,  ...
                           'side',      'i');
u = reshape(u, size(xx));
pass(1) = max(max(abs(u(ii) - sol))) < tol;

% Exterior Laplace
sol = 0;
u = kernels.laplace.dlp(S, 'target',    xy,    ...
                           'density',   sigma, ...
                           'closeeval', true,  ...
                           'side',      'e');
u = reshape(u, size(xx));
pass(2) = max(max(abs(u(~ii) - sol))) < tol;
    
%% Now let's test using panels
np = 40;   % Number of panels
p  = 16;   % Number of quadrature points per panel
N  = p*np; % Total number of quadrature points
%S = Boundary.star(p, 'quadrature', 'panel', 'panels', np);
S = Boundary.star(p, 'quadrature', 'panel');
S = resample(S, 2*p);
N = S.N*S.np;
sigma = ones(N, 1);

% Evaluate on an M x M grid
M = 300;
L = max(abs(cat(1, S.z{:})));
[xx, yy] = meshgrid(linspace(-L, L, M));
xy = [xx(:) yy(:)];
ii = isinterior(S, xx, yy);

% Interior Laplace
sol = -1;
% [xk, ~, vk] = legpts(16);
% x = legpts(32);
% sigma_re = cell(S.np, 1);
% for k = 1:S.np
%     sigma_re{k} = bary(x, ones(16, 1), xk, vk);
% end
% sigma_re = cell2mat(sigma_re);
u = kernels.laplace.dlp(S, 'target',    xy,    ...
                           'density',   sigma, ...
                           'closeeval', true,  ...
                           'side',      'i');
u = reshape(u, size(xx));
pass(3) = max(max(abs(u(ii) - sol))) < 100*tol;

% Exterior Laplace
sol = 0;
u = kernels.laplace.dlp(S, 'target',    xy,    ...
                           'density',   sigma, ...
                           'closeeval', true,  ...
                           'side',      'e');
u = reshape(u, size(xx));
pass(4) = max(max(abs(u(~ii) - sol))) < 100*tol;

end
