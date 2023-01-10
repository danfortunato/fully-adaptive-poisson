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
            S.Gamma1 = defineGamma1(S.Gamma_re, opts);
            Timer.toc('Defining Gamma''');

            S.domain = boundingbox(S.Gamma, 1.1);
            S.domain([1 3]) = min(S.domain);
            S.domain([2 4]) = max(S.domain);

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
            Timer.tic();
            [S.tf, S.isource] = constructTreeIntension(S, S.f, opts);
            Timer.toc('Constructing tree');

            S.strip_solver = StripSolver(S.strip_dom);

        end

    end

end

function Gamma1 = defineGamma1(Gamma, opts)

defaults = [];
defaults.beta = 1;
defaults.bleed = 5;
defaults.stripWidth = 0.5;

if ( nargin < 2 )
    opts = defaults;
else
    opts = setDefaults(opts, defaults);
end

n = Gamma.N;
np = Gamma.np;
[x, y] = smoothStrip(Gamma, n, opts);
xleg = chebvals2legvals(reshape(x, n, np));
yleg = chebvals2legvals(reshape(y, n, np));
z1 = mat2cell(xleg(:) + 1i*yleg(:), repmat(n, np, 1), 1);
Gamma1 = Boundary(z1);

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
