function C = multiscale_circle(N, varargin)
%MULTISCALE   Set up a smooth multiscale domain.

p = inputParser;
p.KeepUnmatched = true;
addRequired(p, 'N', @isfloat);
parse(p, N, varargin{:});

t = sym('t');

freq = 6*pi;
amp = 0.4;
rot = pi/2;

f = @(t) t.*sin(freq*log(t)) / pi * amp;
d1 = 1e-2; s1 = 10; damp1 = @(t) (erf(t/d1-s1)+erf(s1))./(1+erf(s1));
d2 = 0.25;  s2 = 10; damp2 = @(t) (erf((pi-t)/d2-s2)+1)./(1+erf(pi/d2-s2));
g = @(t) damp1(t).*damp2(t).*f(t);
curve = @(t) (mod(t-rot, 2*pi) <  pi) .* g(mod(t-rot, 2*pi)) + ...
             (mod(t-rot, 2*pi) >= pi) .* g(mod(2*pi-(t-rot), 2*pi));

    function r = rfun(t)
        c = curve(t);
        c(abs(t-pi/2)<eps) = 0;
        r = 1 + c;
    end

r = chebfun(@rfun, [0 2*pi], 'trig');
%r = chebfun(@rfun, [0 2*pi], 'trig', 'maxLength', 100000);
%r = chebfun(@rfun, [0 2*pi], 'maxLength', 100000);
%r = chebfun(@rfun, [0 2*pi], 'splitting', 'on');
dr = diff(r);
drr = diff(r, 2);

r = @rfun;

z   = @(t) r(t).*cos(t) + r(t).*sin(t)*1i;
dz  = @(t) (dr(t).*cos(t) - r(t).*sin(t)) + ...
           (dr(t).*sin(t) + r(t).*cos(t))*1i;
dzz = @(t) (drr(t).*cos(t) - 2*dr(t).*sin(t) - r(t).*cos(t)) + ...
           (drr(t).*sin(t) + 2*dr(t).*cos(t) - r(t).*sin(t))*1i;

levelset = @(x,y) sqrt(x.^2+y.^2) - r(atan2(y,x));

C = Boundary(N, z, dz, dzz, levelset, varargin{:});

end
