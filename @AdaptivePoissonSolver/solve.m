function u = solve(S, bc)
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

% Pack up into a struct for easy passing later
u = struct();

u.bulk = solveBulk(S);
[u.strip, bc_strip] = solveStrip(S, u.bulk);
[u.glue_i, u.glue_e] = solveGlue(S, u.bulk, bc_strip);
u.bc = solveBC(S, bc, bc_strip, u.glue_e);
% u_i = @(x,y) u.bulk(x,y) + u.glue_i(x,y) + u.bc(x,y);
% u_e = @(x,y)   stripvals + u.glue_e(x,y) + u.bc(x,y);

end

function u_bulk = solveBulk(S)
%SOLVEBULK   Solve the bulk problem.
%   u_bulk = SOLVEBULK(S) computes a function u_bulk such that
%
%      lap(u_bulk) = S.tf  in S.domain.

u_bulk = treefun2.poisson(-S.tf);

end

function [u_strip, bc_strip] = solveStrip(S, bc)
%SOLVESTRIP   Solve the strip problem.
%   [u_strip, bc_strip] = SOLVESTRIP(S, bc) computes a function u_strip
%   such that
%
%      lap(u_strip) = S.f       in int(S.Gamma) \ int(S.Gamma1),
%           u_strip = bc_strip  on S.Gamma + S.Gamma1,
%
%   where bc_strip contains the values of bc on S.Gamma and S.Gamma1.

% This will call treefun2/feval
bc_strip = bc(S.strip_solver.patches{1}.xy(:,1), S.strip_solver.patches{1}.xy(:,2));
u_strip = S.strip_solver \ bc_strip;

end

function [u_glue_i, u_glue_e] = solveGlue(S, u_bulk, bc_strip)
%SOLVEGLUE   Solve the Neumann glue problem.
%   [u_glue_i, u_glue_e] = SOLVEGLUE(S, u_bulk, bc_strip) computes
%   functions u_glue_i and u_glue_e such that
%
%         lap(u_glue_i/e) = 0                             in int/ext(S.Gamma1),
%            [u_glue_i/e] = 0                             on S.Gamma1,
%      [d(u_glue_i/e)/dn] = d(u_strip)/dn - d(u_bulk)/dn  on S.Gamma1,
%
%   where d/dn denotes the outward normal derivative on S.Gamma1.

% Normal derivative of u_bulk
u_bulk_dn = feval(diff(u_bulk, 1, 2), S.gam1x, S.gam1y) .* S.gam1_nx(:,1) + ...
            feval(diff(u_bulk, 1, 1), S.gam1x, S.gam1y) .* S.gam1_nx(:,2);

% Normal derivative of u_strip
u_strip_dn = S.strip_solver.patches{1}.D2N * [bc_strip ; 1];
u_strip_dn = strip2leg(S, u_strip_dn, 'inner');

% Neumann jump across Gamma'
neu_jump = -(u_strip_dn + u_bulk_dn);

% Correct with single-layer potential
u_glue_i = @(x,y) kernels.laplace.slp(S.Gamma1, 'density', neu_jump, 'target', [x y], 'closeeval', true, 'side', 'i');
u_glue_e = @(x,y) kernels.laplace.slp(S.Gamma1, 'density', neu_jump, 'target', [x y], 'closeeval', true, 'side', 'e');

end

function u_bc = solveBC(S, bc, bc_strip, bc_glue)
%SOLVEBC   Solve the boundary correction problem.
%   u_bc = SOLVEBC(S, bc, bc_strip, bc_glue) computes a function u_bc such
%   that
%
%      lap(u_bc) = 0                        in S.domain \ S.Gamma,
%           u_bc = bc - bc_strip - bc_glue  on S.Gamma.

% Compute the boundary value to correct by
bc_strip = strip2leg(S, bc_strip, 'outer');
bc_corr = bc(S.gamx_re, S.gamy_re) - bc_strip - bc_glue(S.gamx_re, S.gamy_re);

% Solve for double-layer density corresponding to the boundary correction
K = kernels.laplace.dlp(S.Gamma_re);
I = eye(S.n_re * S.Gamma_re.np);
sigma = (K - I/2) \ bc_corr;
% sigma = gmres(K - I/2, bc_corr, [], 1e-14, 50);

% Correct with double-layer potential
u_bc = @(x,y) kernels.laplace.dlp(S.Gamma_re, 'density', sigma, 'target', [x y], 'closeeval', true, 'side', 'i');

end

%%% Utilities
function bc = strip2leg(S, bc, side)

outerIdx = false(size(S.strip_solver.patches{1}.xy, 1), 1);
outerIdx((1:S.n_sem-2).' + (S.n_sem-2)*(0:2:2*length(S.strip_dom)-1)) = true;
innerIdx = ~outerIdx;

switch lower(side)
    case 'inner'
        bc = bc(innerIdx);
    case 'outer'
        bc = bc(outerIdx);
end

bc = reshape(bc, S.n_sem-2, length(S.strip_dom));
% bc = util.modchebvals2legvals(bc);
bc = util.modchebvals2chebvals(bc);
x_re = legpts(S.n_re);
bc = bary(x_re, bc);
bc = bc(:);

end
