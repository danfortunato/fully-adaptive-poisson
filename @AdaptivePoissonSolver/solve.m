function u = solve(S, bc, opts)
%SOLVE   Perform a partition-of-unity solve using function intension.
%   u = SOLVE(S, bc) returns a struct u containing functions such that
%
%      u = / u.bulk  + u.glue_i + u.bc  in int(S.Gamma1),
%          \ u.strip + u.glue_e + u.bc  in int(S.Gamma) \ int(S.Gamma1).
%
%   satisfies
%
%      lap(u) = S.f  in int(S.Gamma),
%           u = bc   on S.Gamma.

defaults = [];
defaults.debug = true;
defaults.useResampledCurve = false;

if ( nargin < 3 )
    opts = defaults;
else
    opts = setDefaults(opts, defaults);
end

% Pack up into a struct for easy passing later
u = struct();

Timer.tic();
u.bulk = solveBulk(S);
Timer.toc('Box code');

Timer.tic();
[u.strip, du_strip] = solveStrip(S);
Timer.toc('Strip solve');

Timer.tic();
u.glue = solveGlue(S, u.bulk, du_strip);
Timer.toc('Glue correction');

Timer.tic();
u.bc = solveBC(S, bc, u.glue.ext, opts);
Timer.toc('Boundary solve');

end

function u_bulk = solveBulk(S)
%SOLVEBULK   Solve the bulk problem.
%   u_bulk = SOLVEBULK(S) computes a function u_bulk such that
%
%      lap(u_bulk) = S.tf  in S.domain.

u_bulk = treefun2.poisson(-S.tf, S.isource);

end

function [u_strip, du_strip] = solveStrip(S)
%SOLVESTRIP   Solve the strip problem.
%   [u_strip, du_strip] = SOLVESTRIP(S) computes a function u_strip such
%   that
%
%      lap(u_strip) = S.f  in int(S.Gamma) \ int(S.Gamma1),
%           u_strip = 0    on S.Gamma + S.Gamma1,
%
%   where du_strip contains the values of the normal derivative of u on
%   S.Gamma and S.Gamma1.

[u_strip, du_strip] = S.strip_solver.solve(S.f);

end

function u_glue = solveGlue(S, u_bulk, du_strip)
%SOLVEGLUE   Solve the Neumann glue problem.
%   u_glue = SOLVEGLUE(S, u_bulk, bc_strip) computes functions u_glue.int
%   and u_glue.ext such that
%
%         lap(u_glue.int/ext) = 0                             in int/ext(S.Gamma1),
%            [u_glue.int/ext] = u_bulk - u_strip              on S.Gamma1,
%      [d(u_glue.int/ext)/dn] = d(u_bulk)/dn - d(u_strip)/dn  on S.Gamma1,
%
%   where d/dn denotes the outward normal derivative on S.Gamma1.

% Value and normal derivative of u_bulk
u_bulk_dir = feval (u_bulk, S.gam1x, S.gam1y);
[dx, dy]   = fevald(u_bulk, S.gam1x, S.gam1y);
u_bulk_neu = dx.*S.gam1_nx(:,1) + dy.*S.gam1_nx(:,2);

% Value and normal derivative of u_strip
% (Note: u_strip_dir is zero)
u_strip_neu = strip2leg(S, du_strip, 'inner');

% Dirichlet and Neumann jumps across Gamma'
dir_jump = u_bulk_dir;
neu_jump = u_bulk_neu - u_strip_neu;

% Correct jumps with layer potentials
u_glue = struct();
u_glue.int = @(x,y) kernels.laplace.dlp(S.Gamma1, 'density',  dir_jump, 'target', [x y], 'closeeval', true, 'side', 'i') + ...
                    kernels.laplace.slp(S.Gamma1, 'density', -neu_jump, 'target', [x y], 'closeeval', true, 'side', 'i');
u_glue.ext = @(x,y) kernels.laplace.dlp(S.Gamma1, 'density',  dir_jump, 'target', [x y], 'closeeval', true, 'side', 'e') + ...
                    kernels.laplace.slp(S.Gamma1, 'density', -neu_jump, 'target', [x y], 'closeeval', true, 'side', 'e');

end

function u_bc = solveBC(S, bc, bc_glue, opts)
%SOLVEBC   Solve the boundary correction problem.
%   u_bc = SOLVEBC(S, bc, bc_strip, bc_glue) computes a function u_bc such
%   that
%
%      lap(u_bc) = 0                        in S.domain \ S.Gamma,
%           u_bc = bc - bc_strip - bc_glue  on S.Gamma.

% Compute the boundary value to correct by:
bc_corr = bc(S.gamx_re, S.gamy_re) - bc_glue(S.gamx_re, S.gamy_re);

% Should we use the original boundary or the resampled boundary?
if ( opts.useResampledCurve )
    Gamma = S.Gamma_re;
else
    Gamma = S.Gamma;
    bc_corr = resample(bc_corr, S.n_re, S.n, S.Gamma.np);
end

L = LaplaceSolver(Gamma, side='interior', bc='dirichlet', method='fmm');
u_bc = L \ bc_corr;

end

%% Utilities

function bc = strip2leg(S, bc, side)
    switch lower(side)
        case 'inner'
            bc = bc(:,:,1);
        case 'outer'
            bc = bc(:,:,2);
    end
    x_re = legpts(S.n_re);
    bc = bary(x_re, bc);
    bc = bc(:);
end

function bc_re = resample(bc, m, n, np)
    x = legpts(n);
    [xk, ~, vk] = legpts(m);
    bc_re = bary(x, reshape(bc, m, np), xk, vk);
    bc_re = bc_re(:);
end
