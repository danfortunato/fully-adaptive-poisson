function foo = myblend(f1, f2, loc, dist)
% BLEND  Uses hyperbolic tangent (tanh) to smoothly blend two functions or
% smoothly combine a piecewise function into a single function handle.
% 
%   foo = BLEND(f1, f2, LOC, DIST) returns function handle foo which, when
%   evaluated, will return a smooth blend between functions provided in handles
%   f1 and f2.
% 
%   The blend is centered around LOC, and most of the smooth blend happens in
%   the region LOC +/- DIST such that F will have value ca. 0.88*f1(x) + 0.12*f2
%   at LOC - DIST and 0.12*f1(x) + 0.88*f2(x) at LOC + DIST. Default LOC = 0.
% 
%   For a pure piecewise function, set the blend distance DIST to zero. F will
%   evaluate to NaN at x = LOC. To avoid evaluating to NaN at x = LOC, use eps
%   instead of zero. Default DIST = 0.
% 
%   f1 and/or f2 may be numeric constants instead of function handles.
% 
%   Example 1: Make a single function handle for a piecewise function (no
%   smoothing).
%     foo = blend(@sin, 1, pi/2);
%     ezplot(foo, [0 pi])
% 
%   Example 2: Make a smooth function from a discontinuous piecewise function
%   and use vectorization to illustrate the effect of different blending
%   distances.
%     f1 = @(x) -x;
%     f2 = @(x)  x;
%     [dist, x] = meshgrid([0 .01 .1 .5 1], -2:.001:4);
%     foo = myblend(f1, f2, 1, dist);
%     plot(x, foo(x))
% 
%   See also tanh, piecewise.

if ( nargin < 3 ), loc = 0; end
if ( nargin < 4 ), dist = 0; end
if ( isnumeric(f1) ), f1 = @(x) f1+0*x; end
if ( isnumeric(f2) ), f2 = @(x) f2+0*x; end

persistent step
if ( isempty(step) )
    step = util.blender(30, 'pswf', [-1 1]);
    cfs = step.coeffs;
    step = @(x) chebtech.clenshaw(x, cfs);
    step = piecewise({@(x) 0*x, @(x) step(x), @(x) 1+0*x}, [-1 1]);
    %step = @(x) (erf(4*x)+1)/2;
end

fac = 4;

%foo = @(x) (1-step((x-loc)./dist/fac)).*f1(x) + step((x-loc)./dist/fac).*f2(x);
    function ifoo = ifoo(x)
        stepx = step((x-loc)./dist/fac);
        ifoo = (1-stepx).*f1(x) + stepx.*f2(x);
    end

foo = @ifoo;

end
