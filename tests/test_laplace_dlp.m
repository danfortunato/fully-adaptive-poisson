% Test that the Laplace double-layer potential with density sigma = 1
% evaluates to -1 inside the domain and 0 outside the domain.

np = 15;   % Number of panels
p  = 10;   % Number of quadrature points per panel
N  = p*np; % Total number of quadrature points
M  = 300;  % Number of plot points

S = ClosedCurve.star(p, 'quadrature', 'panel', 'panels', np);

grid = linspace(-2, 2, M);
[xx, yy] = meshgrid(grid);
xy = [xx(:) yy(:)];
u = LaplaceKernel.dlp(S, 'target', xy, 'density', ones(N,1));
u = reshape(u, size(xx));

figure(1), clf

subplot(121)
hold on
sol = -1;
imagesc(grid, grid, log10(abs(u-sol)))
plot(S, '-k.', 'LineWidth', 2, 'MarkerSize', 15)
hold off
colorbar, caxis([-16 0])
title('Interior Laplace DLP error')
axis equal tight

subplot(122)
hold on
sol = 0;
imagesc(grid, grid, log10(abs(u-sol)))
plot(S, '-k.', 'LineWidth', 2, 'MarkerSize', 15)
hold off
colorbar, caxis([-16 0])
title('Exterior Laplace DLP error')
axis equal tight
