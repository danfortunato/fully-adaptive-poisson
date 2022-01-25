freq = 6*pi;
amp = 0.7;
rot = pi/2;

f = @(t) t.*sin(freq*log(t)) / pi * amp;
d1 = 1e-2; damp1 = @(t) (erf(t/d1-2)+erf(2))./(1+erf(2));
d2 = 0.5;  damp2 = @(t) (erf((pi-t)/d2-2)+erf(2))./(1+erf(2));
f = @(t) damp1(t).*damp2(t).*f(t);
curve = @(t) (mod(t-rot, 2*pi) <  pi) .* f(mod(t-rot, 2*pi)) + ...
             (mod(t-rot, 2*pi) >= pi) .* f(mod(2*pi-(t-rot), 2*pi));

n = 1e6;
t = trigpts(n, [0 2*pi]);
curve = curve(t);
curve(isnan(curve)) = 0;
r = 1 + curve;
plot(r.*cos(t), r.*sin(t))
axis equal
shg
