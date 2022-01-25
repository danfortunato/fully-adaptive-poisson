function fx = bary(x, fvals, xk, vk)
%BARY   Barycentric interpolation formula.
%   BARY(X, FVALS, XK, VK) uses the 2nd form barycentric formula with weights VK
%   to evaluate an interpolant of the data {XK, FVALS(:,k)} at the points X.
%   Note that XK and VK should be column vectors, and FVALS, XK, and VK should
%   have the same length.
%
%   BARY(X, FVALS) assumes XK are the 2nd-kind Chebyshev points and VK are the
%   corresponding barycentric weights.
%
%   If size(FVALS, 2) > 1 then BARY(X, FVALS) returns values in the form
%   [F_1(X), F_2(X), ...], where size(F_k(X)) = size(X) for each k.
%
%   Example:
%     x = chebpts(181);
%     f = 1./( 1 + 25*x.^2 );
%     xx = linspace(-1, 1, 1000);
%     [xx, yy] = meshgrid(xx, xx);
%     ff = bary(xx + 1i*yy, f);
%     h = surf(xx, yy, 0*xx, angle(-ff));
%     set(h, 'edgealpha', 0)
%     view(0, 90), shg
%     colormap(hsv)
%
% See also CHEBTECH.CLENSHAW.

% Parse inputs:
m = size(fvals, 2);
fx = zeros(size(x, 1), m);
for j = 1:numel(x)
    xx = vk ./ (x(j) - xk);
    fx(j,:) = (xx.'*fvals) / sum(xx);
end

end
