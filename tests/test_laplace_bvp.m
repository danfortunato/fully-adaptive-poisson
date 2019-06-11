%% Solve the Laplace equation on interior/exterior domains

np = 15;   % Number of panels
p  = 10;   % Number of quadrature points per panel
N  = p*np; % Total number of quadrature points
M  = 300;  % Number of plot points

%S = ClosedCurve.star(p, 'quadrature', 'ptr');
S = ClosedCurve.ellipse(p, 2, 1, 'quadrature', 'panel', 'panels', np);
I = eye(N);

L = 2;
grid = linspace(-L, L, M);
[xx, yy] = meshgrid(grid);
xy = [xx(:) yy(:)];
plotopts = {'-k.', 'LineWidth', 2, 'MarkerSize', 15};

%% Interior Laplace Dirichlet BVP solver
fint = @(x) real(exp(x(:,1)+x(:,2)*1i)); % Real part of holomorphic function

bc = fint(cell2mat(S.x));
K = LaplaceKernel.dlp(S);
sigma = (K - I/2) \ bc;
uint = LaplaceKernel.dlp(S, 'target', xy, 'density', sigma);
uint = reshape(uint, size(xx));

%% Exterior Laplace Dirichlet BVP solver
fext = @(x) real(1./(x(:,1)+x(:,2)*1i-0.1-0.3i)) + 2; % Point source

bc = fext(cell2mat(S.x));
K = LaplaceKernel.dlp(S, 'modified', true);
sigma = (K + I/2) \ bc;
uext = LaplaceKernel.dlp(S, 'target', xy, 'density', sigma, 'modified', true);
uext = reshape(uext, size(xx));

%% Plot
figure(1), clf

subplot(221)
title('Interior Dirichlet solution')
hold on
imagesc(grid, grid, uint)
plot(S, plotopts{:})
hold off
colorbar, caxis([-1 2])
axis equal tight

subplot(222)
title('Interior Dirichlet error')
sol = reshape(fint(xy), M, M);
hold on
imagesc(grid, grid, log10(abs(uint-sol)))
plot(S, plotopts{:})
hold off
colorbar, caxis([-16 0])
axis equal tight

subplot(223)
title('Exterior Dirichlet solution')
hold on
imagesc(grid, grid, uext)
plot(S, plotopts{:})
hold off
colorbar, caxis([-1 2])
axis equal tight

subplot(224)
title('Exterior Dirichlet error')
sol = reshape(fext(xy), M, M);
hold on
imagesc(grid, grid, log10(abs(uext-sol)))
plot(S, plotopts{:})
hold off
colorbar, caxis([-16 0])
axis equal tight
