% Truncated box code
doExport = false; %#ok<*UNRCH>

% Define the ficitious domain
Bout = Boundary.ellipse(400, 0.9, 0.9, 'quadrature', 'ptr');
Bin  = Boundary.ellipse(400, 0.8, 0.8, 'quadrature', 'ptr');

% Test problem
syms x y
rng(0)
dom = [-1 1 -1 1];
sol = randnfun2(0.5, dom);
sol = sol - min2(sol);
sol = sol ./ max2(sol);
sol = 2*sol - 1;
bc = @(x,y) sol(x,y);
f = lap(sol);
f = @(x,y) f(x,y) + 0*x;
rng(0)

% Make a reference tree to start
tfref = treefun2(@(x,y) f(x,y), dom, 16);
tfref.n = 16;
cutoff = 'hard';
tol = 1e-11;

fs = 30;
lw = 2;
ticks2tex = @(x) cellfun(@(c) ['$' c '$'], x, 'UniformOutput', false);
stylefig = @(cb) chain(...
    @() set(gca, 'fontsize', fs), ...
    @() set(gca, 'linewidth', lw), ...
    @() set(gcf, 'Color', [1 1 1]), ...
    @() colormap('inferno'), ...
    @() axis('equal'), ...
    @() axis('tight'), ...
    @() axis(dom), ...
    @() axis('off'), ...
    @() shading('interp'), ...
    @() view(2), ...
    @() set(cb, 'location', 'EastOutside'), ...
    @() set(cb, 'fontsize', fs), ...
    @() set(cb, 'linewidth', 2), ...
    @() set(cb, 'TickLabelInterpreter', 'latex'), ...
    @() set(cb, 'TickLabels', ticks2tex(cb.TickLabels)), ...
    @() set(cb, 'TickDirection', 'out'), ...
    @() set(cb, 'Color', [0 0 0]), ...
    @() set(cb.Label, 'Color', [0 0 0]));

err = inf;
iter = 1;
while ( max(err(:)) > tol )

    % Compute a particular solution
    switch lower(cutoff)
        case 'hard'
            % Discontinuity occurs along box interfaces
            tf = treefun2(f, tfref);
            leafIds = leaves(tf);
            [xx, yy] = leafpts(tf);
            for k = 1:length(tf)
                if ( any(~isinterior(Bout, xx{k}, yy{k}), 'all') )
                    tf.coeffs{leafIds(k)} = zeros(tf.n);
                end
            end
        case 'soft'
            % Discontinuity occurs within a box
            tf = treefun2(@(x,y) f(x,y).*(sqrt(x.^2+y.^2) <= rout), tfref);
    end

    % Run the box code
    u_bulk = treefun2.poisson(-tf);

    % Correct the boundary data
    Bxy = cell2mat(Bin.x);
    Bx = Bxy(:,1);
    By = Bxy(:,2);

    % Compute the boundary value to correct by
    bc_corr = bc(Bx, By) - u_bulk(Bx, By);

    % Solve for double-layer density corresponding to the boundary correction
    K = kernels.laplace.dlp(Bin);
    I = eye(Bin.N * Bin.np);
    sigma = (K - I/2) \ bc_corr;
    u_bc = @(x,y) kernels.laplace.dlp(Bin, 'density', sigma, 'target', [x y], 'closeeval', true, 'side', 'i');

    % Evaluate and plot
    m = 400;
    r = 1;
    [xx, yy] = meshgrid(linspace(dom(1), dom(2), m), ...
                        linspace(dom(3), dom(4), m));
    inB = isinterior(Bin, xx, yy);
    vv = nan(m);
    vv(inB) = u_bulk(xx(inB),yy(inB)) + u_bc(xx(inB),yy(inB));
    err = abs(vv - sol(xx,yy)) ./ max(max(abs(sol(xx,yy))));
    max(err(:))

    reserr = abs(lap(u_bulk) - tf);
    normF = norm(tf, inf);
    for k = leaves(tf)
        reserr.coeffs{k} = reserr.coeffs{k} / normF;
    end

    ids = leaves(tf);
    keep = [];
    for id = ids(:).'
        if ( ~all(tf.coeffs{id} == 0, 'all') )
            keep(end+1) = id; %#ok<SAGROW>
        end
    end

    % Plot
    figure(1), clf
    plotOnly(tf, keep)
    hold on
    plot(Bout, 'k',   'linewidth', lw)
    hold off
    cb = colorbar;
    clim([-50 50])
    cb.Ticks = -50:25:50;
    stylefig(cb);
    pos = cb.Position;
    pos(3) = 0.04;
    cb.Position = pos;
    delete(cb);

    figure(2), clf
    plotOnly(log10(reserr), keep)
    hold on
    plot(Bin,  'k--', 'linewidth', lw)
    plot(Bout, 'k',   'linewidth', lw)
    hold off
    cb1 = colorbar;
    clim([-16 0])
    cb1.Ticks = -16:4:0;
    cb1.TickLabels = cellfun(@(s) ['10^{' s '}'], cb1.TickLabels, 'UniformOutput', false);
    stylefig(cb1);
    pos = cb1.Position;
    pos(3) = 0.04;
    cb1.Position = pos;

    figure(3), clf
    temp = log10(err);
    temp(isinf(temp)) = -16;
    surf(xx, yy, temp)
    hold on
    plot(Bin,  'k--', 'linewidth', lw)
    plot(Bout, 'k',   'linewidth', lw)
    hold off
    cb2 = colorbar;
    clim([-16 0])
    cb2.Ticks = -16:4:0;
    cb2.TickLabels = cellfun(@(s) ['10^{' s '}'], cb2.TickLabels, 'UniformOutput', false);
    stylefig(cb2);
    pos = cb2.Position;
    pos(3) = 0.04;
    cb2.Position = pos;

    if ( doExport )
        figure(1)
        export_fig('-r400', '-transparent', sprintf('images/truncation_only_f_%d.png', iter));
        figure(2)
        export_fig('-r400', '-transparent', sprintf('images/truncation_only_res_%d.png', iter));
        figure(3)
        export_fig('-r400', '-transparent', sprintf('images/truncation_only_err_%d.png', iter));
    end

    % Pause
    drawnow
    alignfigs
    keyboard

    % Refine
    unresolved = norm(reserr, inf, 'all') > 1e-8;
    leafs = leaves(tfref);
    tfref = refine(tfref, leafs(unresolved));

    iter = iter+1;

end
