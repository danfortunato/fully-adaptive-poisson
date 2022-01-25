% Test that the Laplace double-layer potential with density sigma = 1
% evaluates to -1 inside the domain and 0 outside the domain.

np = 40;   % Number of panels
p  = 16;   % Number of quadrature points per panel
N  = p*np; % Total number of quadrature points
M  = 800;  % Number of plot points
closeeval = true;

if ( np == 1 )
    S = Boundary.star(p);
else
    %S = Boundary.star(p, 'quadrature', 'panel', 'panels', np);
    %S = Boundary.star(p, 'quadrature', 'panel');
    %S = Boundary.wavepacket(p, 'quadrature', 'panel');
    S = Boundary.gobbler(p, 'delta', 0.1, 'align', 1, 'quadrature', 'panel');
end
%%
N = p*S.np;

L = max(abs(cat(1, S.z{:})));
grid = linspace(-L, L, M);
[xx, yy] = meshgrid(grid);
ii = isinterior(S, xx, yy);
xy = [xx(:) yy(:)];

u = kernels.laplace.dlp(S, 'target',    xy,        ...
                           'density',   ones(N,1), ...
                           'closeeval', false, ...
                           'side',      'i');
u = reshape(u, size(xx));
u(abs(u+1) < 0.9) = nan;
%u(abs(u) > 0.5) = nan;

% chnkr = makeChunker(S);
% %chnkr = chunkerfunc(@(t) starfish(t));
% kernd = @(s,t) chnk.lap2d.kern(s, t, 'd');
% dens1 = ones(chnkr.k, chnkr.nch);
% wts = weights(chnkr);
% opts = []; %opts.usesmooth = 1;
% u = chunkerkerneval(chnkr, kernd, dens1, xy.', opts);
% u = reshape(u, size(xx));

figure(1), clf

subplot(121)
hold on
sol = -1;
%imagesc(grid, grid, log10(abs(u-sol)))
%imagesc(grid, grid, u)
surf(grid, grid, u, 'EdgeColor', 'None')
plot(S, '-k.', 'LineWidth', 2, 'MarkerSize', 15)
hold off
colorbar, %caxis([-16 0])
title('Interior Laplace DLP error')
axis equal tight

% u = kernels.laplace.dlp(S, 'target',    xy,        ...
%                            'density',   ones(N,1), ...
%                            'closeeval', false, ...
%                            'side',      'e');
u = kernels.laplace.dlp(S, 'target', xy, 'density', ones(N,1));
u = reshape(u, size(xx));
%u(abs(u) < 1e-8) = nan;

subplot(122)
hold on
sol = 0;
%imagesc(grid, grid, log10(abs(u-sol)))
surf(grid, grid, u, 'EdgeColor', 'None')
plot(S, '-k.', 'LineWidth', 2, 'MarkerSize', 15)
hold off
colorbar, %caxis([-16 0])
caxis([-1 0])
title('Exterior Laplace DLP error')
axis equal tight

colormap(gcf, 'turbo')
