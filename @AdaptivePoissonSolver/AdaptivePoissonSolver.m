classdef AdaptivePoissonSolver < handle %#ok<*PROPLC>
%ADAPTIVEPOISSONSOLVER   An adaptive Poisson solver based on a partition of
% unity via function intension.

    properties

        n      % Number of nodes per panel
        n_re   % Resample original curve to use n_re per panel
        n_sem  % Number of nodes for each spectral element
        n_box  % Number of nodes for each box in tree
        m_box  % Box error is computed using m_box equispaced points

        Gamma    % Original boundary
        Gamma1   % Fictitious boundary
        Gamma_re % Resampled boundary
        domain   % Bounding square containing Gamma

        f  % Right-hand side
        tf % Right-hand side after function intension and boxing

        strip_dom
        strip_solver

    end

    properties ( Hidden )

        gamx
        gamy
        gam1x
        gam1y
        gamx_re
        gamy_re
        gamx_cheb
        gamy_cheb
        gam1x_cheb
        gam1y_cheb
        gam_nx
        gam1_nx
        gam_nx_cheb
        gam1_nx_cheb

    end

    methods

        function S = AdaptivePoissonSolver(Gamma, f)

            S.n = Gamma.N;
            S.n_re  = 2*S.n;
            S.n_sem = 2*S.n;
            S.n_box = S.n;
            S.m_box = 2*S.n_box;
            %S.m_box = S.n_box+2;

            S.Gamma = Gamma;
            S.Gamma_re = resample(S.Gamma, S.n_re);
            S.Gamma1 = defineGamma1(S.Gamma_re);
            %S.Gamma1 = Boundary.ellipse(S.n_re, 0.8, 0.8, 'quadrature', 'panel', 'panels', length(S.Gamma.x));
            S.domain = boundingbox(S.Gamma, 1.2);
            S.domain([1 3]) = min(S.domain);
            S.domain([2 4]) = max(S.domain);

            gamxy    = cell2mat(S.Gamma.x);
            gam1xy   = cell2mat(S.Gamma1.x);
            gamxy_re = cell2mat(S.Gamma_re.x);
            S.gamx    = gamxy(:,1);
            S.gamy    = gamxy(:,2);
            S.gam1x   = gam1xy(:,1);
            S.gam1y   = gam1xy(:,2);
            S.gamx_re = gamxy_re(:,1);
            S.gamy_re = gamxy_re(:,2);
            S.gam_nx  = cell2mat(S.Gamma.normal);
            S.gam1_nx = cell2mat(S.Gamma1.normal);

            [S.strip_dom,    ...
             S.gamx_cheb,    ...
             S.gamy_cheb,    ...
             S.gam1x_cheb,   ...
             S.gam1y_cheb,   ...
             S.gam_nx_cheb,  ...
             S.gam1_nx_cheb] = buildStripGrid(S.Gamma_re, S.Gamma1);

            S.f = f;
            S.tf = constructTreeIntension2(S, S.f);

            S.strip_solver = StripSolver(S.strip_dom, S.f);
            build(S.strip_solver);

        end

    end

end

function Gamma1 = defineGamma1(Gamma)

n = Gamma.N;
np = Gamma.np;
%beta = 4;
beta = 6;
bleed = 10;
width = 0.4;
[x, y] = smoothStrip2(Gamma, n, beta, bleed, width);
xleg = chebvals2legvals(reshape(x, n, np));
yleg = chebvals2legvals(reshape(y, n, np));
z1 = mat2cell(xleg(:) + 1i*yleg(:), repmat(n, np, 1), 1);
Gamma1 = Boundary(z1);

% z = cell(np, 1);
% width = 0.2;
% for k = 1:np
%     nx = Gamma.normal{k};
%     nx = nx(:,1)+nx(:,2)*1i;
%     nx = nx./abs(nx);
%     z{k} = Gamma.z{k} - width*nx;
% end
% Gamma1 = Boundary(z);

end

function [strip_dom, gamx_cheb, gamy_cheb, gam1x_cheb, gam1y_cheb, ...
    gam_nx_cheb, gam1_nx_cheb] = buildStripGrid(Gamma, Gamma1)

n  = Gamma.N;
np = Gamma.np;
t = chebpts(n, [0 1]);
LV2CV = legvals2chebvals(eye(n));
xx = cell(np, 1);
yy = cell(np, 1);
gamx_cheb    = cell(np, 1);
gamy_cheb    = cell(np, 1);
gam1x_cheb   = cell(np, 1);
gam1y_cheb   = cell(np, 1);
gam_nx_cheb  = cell(np, 1);
gam1_nx_cheb = cell(np, 1);
for k = 1:np
    gamx_cheb{k}  = LV2CV * real(Gamma.x{k}(:,1));
    gamy_cheb{k}  = LV2CV * real(Gamma.x{k}(:,2));
    gam1x_cheb{k} = LV2CV * real(Gamma1.x{k}(:,1));
    gam1y_cheb{k} = LV2CV * real(Gamma1.x{k}(:,2));
    gam_nx_cheb{k} = LV2CV * Gamma.normal{k};
    gam1_nx_cheb{k} = LV2CV * Gamma1.normal{k};
    xx{k} = t.*gam1x_cheb{k}.' + (1-t).*gamx_cheb{k}.';
    yy{k} = t.*gam1y_cheb{k}.' + (1-t).*gamy_cheb{k}.';
end
strip_dom = cell2struct([xx yy], {'x','y'}, 2);

% xsem = chebpts(n_sem);
% xx = cell(Gamma_re.np,1);
% gamx_cheb_re    = cell(Gamma_re.np, 1);
% gamy_cheb_re    = cell(Gamma_re.np, 1);
% gam1x_cheb_re   = cell(Gamma_re.np, 1);
% gam1y_cheb_re   = cell(Gamma_re.np, 1);
% gam_nx_cheb_re  = cell(Gamma_re.np, 1);
% gam1_nx_cheb_re = cell(Gamma_re.np, 1);
% for k = 1:Gamma.np
%     gamx_cheb_re{k}  = bary(xsem, gamx_cheb{k});
%     gamy_cheb_re{k}  = bary(xsem, gamy_cheb{k});
%     gam1x_cheb_re{k} = bary(xsem, gam1x_cheb{k});
%     gam1y_cheb_re{k} = bary(xsem, gam1y_cheb{k});
%     gam_nx_cheb_re{k} = bary(xsem, gam_nx_cheb{k});
%     gam1_nx_cheb_re{k} = bary(xsem, gam1_nx_cheb{k});
%     xx{k} = t.*gam1x_cheb_re{k}.' + (1-t).*gamx_cheb_re{k}.';
%     yy{k} = t.*gam1y_cheb_re{k}.' + (1-t).*gamy_cheb_re{k}.';
% end

end
