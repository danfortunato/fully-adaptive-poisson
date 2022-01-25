%% Solve the Laplace equation on interior/exterior domains

np = 60;   % Number of panels
p  = 16;  % Number of quadrature points per panel
N  = p*np; % Total number of quadrature points
M  = 400;  % Number of plot points

%S = Boundary.star(N);
S = Boundary.star(p, 'quadrature', 'panel');
%S = Boundary.star(p, 'quadrature', 'panel', 'panels', np);
%S = Boundary.ellipse(p, 2, 1, 'quadrature', 'panel', 'panels', np);

np = S.np;
N = p*np;
I = eye(N);

L = 2;
grid = linspace(-L, L, M);
[xx, yy] = meshgrid(grid);
xy = [xx(:) yy(:)];
ii = isinterior(S, xx, yy);
plotopts = {'-k.', 'LineWidth', 2, 'MarkerSize', 15};

%% Interior Laplace Dirichlet BVP solver
% fint = @(x) real(exp(x(:,1)+x(:,2)*1i)); % Real part of holomorphic function
% 
% bc = fint(cell2mat(S.x));
% K = kernels.laplace.dlp(S);
% sigma = (K - I/2) \ bc;
% uint = nan(M);
% tic
% uint(ii) = kernels.laplace.dlp(S, 'target',    [xx(ii) yy(ii)], ...
%                                   'density',   sigma,           ...
%                                   'closeeval', false,            ...
%                                   'side',      'i');
% toc
% uint = reshape(uint, size(xx));

%% Interior DLP + SLP
fint = @(x) real(exp(x(:,1)+x(:,2)*1i)); % Real part of holomorphic function
%fint = @(x) real(exp(1i*(x(:,1)+x(:,2)*1i+1))); % Real part of holomorphic function

bc = fint(cell2mat(S.x));
A = kernels.laplace.dlp(S) - I/2;% + kernels.laplace.slp(S);
sigma = A \ bc;
uint = nan(M);
tic
uint(ii) = kernels.laplace.dlp(S, 'target',  [xx(ii) yy(ii)], ...
                                  'density', sigma, 'closeeval', true, 'side', 'i');%         + ...
           %kernels.laplace.slp(S, 'target',  [xx(ii) yy(ii)], ...
           %                       'density', sigma, 'closeeval', true, 'side', 'i');
toc
uint = reshape(uint, size(xx));

%% Exterior Laplace Dirichlet BVP solver
% fext = @(x) real(1./(x(:,1)+x(:,2)*1i-0.1-0.3i)) + 2; % Point source
% 
% bc = fext(cell2mat(S.x));
% K = kernels.laplace.dlp(S, 'modified', true);
% sigma = (K + I/2) \ bc;
% uext = nan(M);
% tic
% uext(~ii) = kernels.laplace.dlp(S, 'target',    [xx(~ii) yy(~ii)], ...
%                                  'density',   sigma,             ...
%                                  'modified',  true,              ...
%                                  'closeeval', true,              ...
%                                  'side',      'e');
% toc
% uext = reshape(uext, size(xx));

%% Exterior Laplace Dirichlet BVP solver
fext = @(x) real(1./(x(:,1)+x(:,2)*1i-0.1-0.3i)) + 2; % Point source

bc = fext(cell2mat(S.x));
A = kernels.laplace.dlp(S, 'modified', true) + I/2;
sigma = A \ bc;
uext = nan(M);
tic
uext(~ii) = kernels.laplace.dlp(S, 'target',    [xx(~ii) yy(~ii)], ...
                                   'density',   sigma,             ...
                                   'modified',  true,              ...
                                   'closeeval', true,              ...
                                   'side',      'e');
toc
uext = reshape(uext, size(xx));

%% Plot
figure(1), clf

subplot(221)
title('Interior Dirichlet solution')
hold on
imagesc(grid, grid, uint)
plot(S, plotopts{:})
hold off
colorbar
axis equal tight

subplot(222)
title('Interior Dirichlet error')
sol = reshape(fint(xy), M, M);
hold on
sol(~ii) = nan;
h = imagesc(grid, grid, log10(abs(uint-sol)));
%set(h, 'AlphaData', ~isnan(sol))
plot(S, plotopts{:})
hold off
colorbar
shading interp
axis equal tight

subplot(223)
title('Exterior Dirichlet solution')
hold on
imagesc(grid, grid, uext)
plot(S, plotopts{:})
hold off
colorbar
axis equal tight

subplot(224)
title('Exterior Dirichlet error')
sol = reshape(fext(xy), M, M);
hold on
imagesc(grid, grid, log10(abs(uext-sol)))
plot(S, plotopts{:})
hold off
colorbar
axis equal tight

colormap(gcf, 'turbo')
