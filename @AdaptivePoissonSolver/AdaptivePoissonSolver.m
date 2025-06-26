classdef AdaptivePoissonSolver < handle %#ok<*PROPLC>
%ADAPTIVEPOISSONSOLVER   An adaptive Poisson solver based on a partition of
% unity via function intension.

    properties

        n            % Number of nodes per panel
        n_re         % Resample original curve to use n_re per panel
        n_sem        % Number of nodes for each spectral element
        n_box        % Number of nodes for each box in tree
        m_box        % Box error is computed using m_box equispaced points

        Gamma        % Original boundary
        Gamma1       % Fictitious boundary
        Gamma_re     % Resampled boundary
        domain       % Bounding square containing Gamma

        f            % Right-hand side
        tf           % Right-hand side after function intension and boxing
        isource      % Boxes of tf which contribute to volume potential

        strip_dom    % Curvilinear mesh defining strip region
        strip_solver % Strip solver
        opts

    end

    properties ( Hidden )

        gamxy
        gamx
        gamy
        gamx_cheb
        gamy_cheb
        gam_nx
        gam_nx_cheb

        gam1xy
        gam1x
        gam1y
        gam1x_cheb
        gam1y_cheb
        gam1_nx
        gam1_nx_cheb

        gamxy_re
        gamx_re
        gamy_re

    end

    methods

        function S = AdaptivePoissonSolver(Gamma, f, opts)

            defaults = [];
            defaults.debug = true;
            defaults.n_re = 2*Gamma.N;
            defaults.n_sem = defaults.n_re;
            defaults.n_box = 16;
            defaults.m_box = 2*defaults.n_box;
            defaults.beta = 1;
            defaults.stripWidth = 0.5;
            defaults.bleed = 5;
            defaults.tol = 1e-11;
            defaults.boxToStripRatio = 0.5;
            defaults.maxBoxes = 300000;

            if ( nargin < 3 )
                opts = defaults;
            else
                opts = setDefaults(opts, defaults);
            end
            S.opts = opts;

            S.n     = Gamma.N;
            S.n_re  = opts.n_re;
            S.n_sem = opts.n_sem;
            S.n_box = opts.n_box;
            S.m_box = opts.m_box;

            S.Gamma = Gamma;

            Timer.tic();
            S.Gamma_re = resample(S.Gamma, S.n_re);
            Timer().toc('Resampling Gamma');

            Timer.tic();
            S.Gamma1 = AdaptivePoissonSolver.defineGamma1(S.Gamma_re, opts);
            Timer.toc('Defining Gamma''');

            S.domain = boundingbox(S.Gamma, 1.1);
            dx = diff(S.domain(1:2))/2; cx = sum(S.domain(1:2))/2;
            dy = diff(S.domain(3:4))/2; cy = sum(S.domain(3:4))/2;
            d = max(dx, dy);
            S.domain = [cx-d cx+d cy-d cy+d];

            S.gamxy  = cell2mat(S.Gamma.x);
            S.gamx   = S.gamxy(:,1);
            S.gamy   = S.gamxy(:,2);
            S.gam_nx = cell2mat(S.Gamma.normal);

            S.gam1xy  = cell2mat(S.Gamma1.x);
            S.gam1x   = S.gam1xy(:,1);
            S.gam1y   = S.gam1xy(:,2);
            S.gam1_nx = cell2mat(S.Gamma1.normal);

            S.gamxy_re = cell2mat(S.Gamma_re.x);
            S.gamx_re  = S.gamxy_re(:,1);
            S.gamy_re  = S.gamxy_re(:,2);

            Timer.tic();
            [S.strip_dom, S.gamx_cheb,  S.gamy_cheb,  S.gam_nx_cheb, ...
                          S.gam1x_cheb, S.gam1y_cheb, S.gam1_nx_cheb ...
            ] = buildStripGrid(S.Gamma_re, S.Gamma1, S.n_sem);
            Timer.toc('Building strip grid');

            S.f = f;
            idx = isStripResolved(S.strip_dom, S.f, opts.tol);
            if ( ~isempty(idx) )
                warning('Inhomogeneity is not resolved to requested tolerance in strip. Refining...');
                Gamma_ref = refine(Gamma, 1, idx);
                S = AdaptivePoissonSolver(Gamma_ref, f, opts);
                return
            end

            Timer.tic();
            [S.tf, S.isource] = constructTreeIntension(S, S.f, opts);
            Timer.toc('Constructing tree');

            S.strip_solver = StripSolver(S.strip_dom);

        end

    end

    methods ( Static )

        varargout = defineGamma1(varargin);

    end

end

function [strip_dom, gamx_cheb, gamy_cheb, gam_nx_cheb, gam1x_cheb, ...
    gam1y_cheb, gam1_nx_cheb] = buildStripGrid(Gamma, Gamma1, n_sem)

if ( ~isscalar(n_sem) )
    nt = n_sem(1);
    nr = n_sem(2);
else
    nt = n_sem;
    nr = n_sem;
end

n  = Gamma.N;
np = Gamma.np;
t = chebpts(nt, [0 1]);
LV2CV = legvals2chebvals(eye(n));
r = chebpts(nr);
x = chebpts(n);
RE = barymat(r, x);
xx = cell(np, 1);
yy = cell(np, 1);
gamx_cheb    = cell(np, 1);
gamy_cheb    = cell(np, 1);
gam1x_cheb   = cell(np, 1);
gam1y_cheb   = cell(np, 1);
gam_nx_cheb  = cell(np, 1);
gam1_nx_cheb = cell(np, 1);
for k = 1:np
    gamx_cheb{k}    = RE * LV2CV * real(Gamma.x{k}(:,1));
    gamy_cheb{k}    = RE * LV2CV * real(Gamma.x{k}(:,2));
    gam1x_cheb{k}   = RE * LV2CV * real(Gamma1.x{k}(:,1));
    gam1y_cheb{k}   = RE * LV2CV * real(Gamma1.x{k}(:,2));
    gam_nx_cheb{k}  = RE * LV2CV * Gamma.normal{k};
    gam1_nx_cheb{k} = RE * LV2CV * Gamma1.normal{k};
    xx{k} = t.*gam1x_cheb{k}.' + (1-t).*gamx_cheb{k}.';
    yy{k} = t.*gam1y_cheb{k}.' + (1-t).*gamy_cheb{k}.';
end
strip_dom = cell2struct([xx yy], {'x','y'}, 2);

end

function idx = isStripResolved(dom, f, tol)

[ny, nx] = size(dom(1).x);
xx = cat(3, dom.x);
yy = cat(3, dom.y);
ff = feval(f, xx, yy);
hx = max(xx, [], 1:2) - min(xx, [], 1:2);
hy = max(yy, [], 1:2) - min(yy, [], 1:2);
hh = squeeze(min(hx, hy));

eta = 0;
resolved = false(length(dom), 1);
for k = 1:length(dom)
    vals = ff(:,:,k);
    coeffs = chebtech2.vals2coeffs( chebtech2.vals2coeffs(vals).' ).';
    coeffs = coeffs(1:ny,1:nx);
    Ex = sum(abs(coeffs(end-1:end,:)), 'all') / (2*nx);
    Ey = sum(abs(coeffs(:,end-1:end)), 'all') / (2*ny);
    err_cfs = (Ex + Ey) / 2;
    h = hh(k);
    vmax = max(abs(vals(:)));
    resolved(k) = ( err_cfs * h^eta < tol * max(vmax, 1) );
end

idx = find(~resolved);

end
