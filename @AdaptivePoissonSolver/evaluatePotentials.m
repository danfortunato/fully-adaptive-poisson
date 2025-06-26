function [u_bulk, u_strip, err_bulk, err_strip] = evaluatePotentials(S, u, sol)
%EVALUATEPOTENTIALS   Evaluate the layer potentials back on the tree and
%strip.

if ( nargin < 3 )
    sol = [];
end

%% Evaluate back on the strip

% Note: Need to use type-1 Chebyshev nodes here to avoid on-surface layer
% potential evaluation.

Timer.tic();

strip_x = cat(3, S.strip_dom.x);
strip_y = cat(3, S.strip_dom.y);
strip_x = reshape(chebvals2chebvals(strip_x, 2, 1), size(u.strip));
strip_y = reshape(chebvals2chebvals(strip_y, 2, 1), size(u.strip));

Timer.tic();
u_glue = reshape(u.glue.ext(strip_x(:),strip_y(:)), size(u.strip));
Timer.toc('Evaluate glue contribution');

Timer.tic();
u_bc = reshape(u.bc(strip_x(:),strip_y(:)), size(u.strip));
Timer.toc('Evaluate homogeneous solution');

u_strip = reshape(chebvals2chebvals(u.strip, 2, 1), size(u.strip));
u_strip = u_strip + u_glue + u_bc;

Timer.toc('Evaluate back on strip');

if ( isempty(sol) )
    sol_strip = [];
else
    sol_strip = sol(strip_x, strip_y);
    err_strip = max(abs(u_strip(:) - sol_strip(:)));
end

%% Evaluate back on the tree

Timer.tic();

[tree_x, tree_y] = leafpts(u.bulk);
n = u.bulk.n;
tree_x = reshape([tree_x{:}], n, n, []);
tree_y = reshape([tree_y{:}], n, n, []);

Timer.tic();
inGamma  = isinterior(S.Gamma,  tree_x, tree_y);
inGamma1 = isinterior(S.Gamma1, tree_x, tree_y);
inStrip = inGamma & ~inGamma1;
Timer.toc('Classify tree nodes as interior or exterior');

Timer.tic();
u_glue_int_vals = zeros(size(tree_x));
u_glue_int_vals(inGamma1) = u.glue.int(tree_x(inGamma1), tree_y(inGamma1));
% u_glue_ext_vals = zeros(size(tree_x));
% u_glue_ext_vals(inStrip) = u.glue.ext(tree_x(inStrip), tree_y(inStrip));
Timer.toc('Evaluate glue contribution');

Timer.tic();
u_bc_vals = zeros(size(tree_x));
u_bc_vals(inGamma) = u.bc(tree_x(inGamma), tree_y(inGamma));
Timer.toc('Evaluate homogeneous solution');

Timer.tic();
tree_vals = leafvals(u.bulk);
tree_vals = cat(3, tree_vals{:});
tree_vals(~inGamma) = nan;
tree_vals = tree_vals + u_glue_int_vals + u_bc_vals;
u_bulk = u.bulk;
ids = leaves(u.bulk);
for k = 1:length(ids)
    u_bulk.coeffs{ids(k)} = treefun2.vals2coeffs(tree_vals(:,:,k));
end
Timer.toc('Package into treefun2');

Timer.toc('Evaluate back on tree');

if ( isempty(sol) )
    sol_bulk = [];
else
    sol_vals = nan(size(tree_vals));
    sol_vals(inGamma) = sol(tree_x(inGamma), tree_y(inGamma));
    sol_bulk = u.bulk;
    for k = 1:length(ids)
        sol_bulk.coeffs{ids(k)} = treefun2.vals2coeffs(sol_vals(:,:,k));
    end
    err_bulk = max(abs(tree_vals(inGamma1) - sol_vals(inGamma1)));
end

end
