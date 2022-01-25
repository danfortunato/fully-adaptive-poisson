function [uu, xx, yy, inGamma] = sampleUniform(S, u, nbulk)

if ( nargin < 3 )
    nbulk = 400;
end

uu = nan(nbulk);

% Sample on a uniform grid of size nbulk x nbulk
[xx, yy] = meshgrid(linspace(S.domain(1), S.domain(2), nbulk).', ...
                    linspace(S.domain(3), S.domain(4), nbulk).');

% Partition the grid
[inGamma, inGamma1, inStrip] = S.partition(xx, yy);
xx_g = xx(inGamma); xx_g1 = xx(inGamma1); xx_s = xx(inStrip);
yy_g = yy(inGamma); yy_g1 = yy(inGamma1); yy_s = yy(inStrip);

% Sample the solution inside Gamma'
ub = u.bulk(xx,yy);
uu(inGamma1) = ub(inGamma1) + u.glue_i(xx_g1,yy_g1) + u.bc(xx_g1,yy_g1);

% Sample the solution in the strip
stripvals = zeros(size(xx_s));
[tt_s, rr_s, globalIdx] = S.assignStripPoints(xx_s, yy_s);
for k = 1:length(S.strip_dom)
    stripvals(globalIdx{k}) = util.bary2d(u.strip{k}, tt_s{k}, rr_s{k});
end
uu(inStrip) = stripvals + u.glue_e(xx_s, yy_s) + u.bc(xx_s, yy_s);

uu = real(uu);

end
