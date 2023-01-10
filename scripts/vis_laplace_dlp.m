% Test that the Laplace double-layer potential with density sigma = 1
% evaluates to -1 inside the domain and 0 outside the domain.

np = 40;   % Number of panels
p  = 16;   % Number of quadrature points per panel
N  = p*np; % Total number of quadrature points
M  = 400;  % Number of plot points
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
%S = C;
%S = Gamma;
%S = resample(Gamma, 20);
%S = Boundary.star(16, 'quadrature', 'panel');
%S = Boundary.star(16, quadrature='panel', arms=40);
%S = B;
%S = G;
%S = resample(S, 6);
%S = Boundary.star(16, quadrature='panel', arms=10);
%S = Gamma;

N = p*S.np;

L = max(abs(cat(1, S.z{:})));
dom = [-L L -L L]; 
%dom = [-0.35 -0.25 -0.7 -0.55];
%dom = [-0.5 0 -1 -0.5];
%dom = [ -0.2739   -0.2644   -0.5768   -0.5625 ];
%dom = [ -0.2691   -0.2678   -0.5690   -0.5670];
%dom = [-0.1732    0.1605    0.5744    0.9081];
gridx = linspace(dom(1), dom(2), M);
gridy = linspace(dom(3), dom(4), M);
[xx, yy] = meshgrid(gridx, gridy);
%ii = isinterior(resample(refine(S), 32), xx, yy);
ii = isinterior(S, xx, yy);
xy = [xx(:) yy(:)];

uu = kernels.laplace.dlp(S, target=xy(ii,:), density=ones(S.N*S.np,1), closeeval=true, side='i');
u = nan(size(xx));
u(ii) = uu;
%u = reshape(u, size(xx));
%u(abs(u+1) < 0.9) = nan;
%u(abs(u) > 0.5) = nan;
%u(~ii) = nan;

%kernels.laplace.dlp(S, target=[xi yi], density=ones(N,1), closeeval=true, side='i');

% chnkr = Boundary.toChunker(S);
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
imagesc(gridx, gridy, log10(abs(u-sol)))
%imagesc(grid, grid, u)
%surf(grid, grid, log10(abs(u-sol)), 'EdgeColor', 'None')
plot(S.zbreaks, 'k.', 'markersize', 10)
%plot(S)
%plot(S, '-k', 'LineWidth', 2, 'MarkerSize', 15)
hold off
colorbar, %caxis([-16 0])
title('Interior Laplace DLP error')
axis equal tight
axis(dom)

%caxis([-16 -11])

%%
u = kernels.laplace.dlp(S, 'target',    xy,        ...
                           'density',   ones(N,1), ...
                           'closeeval', true, ...
                           'side',      'e');
%u = kernels.laplace.dlp(S, 'target', xy, 'density', ones(N,1));
u = reshape(u, size(xx));
%u(abs(u) < 1e-8) = nan;
u(ii) = nan;

subplot(122)
hold on
sol = 0;
imagesc(gridx, gridy, log10(abs(u-sol)))
%surf(grid, grid, u, 'EdgeColor', 'None')
%plot(S, '-k.', 'LineWidth', 2, 'MarkerSize', 15)
hold off
colorbar, %caxis([-16 0])
title('Exterior Laplace DLP error')
axis equal tight

colormap(gcf, 'turbo')
