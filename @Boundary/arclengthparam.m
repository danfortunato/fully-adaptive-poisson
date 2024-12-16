function [s, breaks] = arclengthparam(dom)
%ARCLENGTHPARAM   Compute arc length coordinates of a given curve.
%   S = ARCLENGTHPARAM(DOM) computes arc length coordinates S of the nodes in
%   the given boundary DOM. S is a matrix such that S(j,k) is the arc length of
%   the jth node on the kth panel.
%
%   [S, BREAKS] = ARCLENGTHPARAM(DOM) also returns the arc length
%   coordinates of the breakpoints of each panel, so that the kth panel has
%   parametric endpoints BREAKS(k) and BREAKS(k+1).

% For each panel, compute the arc length of each Legendre node from the panel's
% starting endpoint:
speed = chebfun(legvals2chebvals([dom.speed{:}]));
scls = diff(dom.breaks) / 2;
f = cumsum(speed .* scls);
x = legpts(dom.N);
s = f(x);

% Shift the local panel arc lengths by the length of all previous panels:
breaks = [0 cumsum(sum([dom.w{:}]))];
s = s + breaks(1:end-1);

end
