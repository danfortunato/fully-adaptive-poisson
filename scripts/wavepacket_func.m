function [r, d, d2] = wavepacket_func(t)

w = 1e4;

syms x
packet_f = @(x) 1/sqrt(w) * cos(pi*sqrt(w)*x) .* exp(-w*x.^2);
r = @(x) 1 + packet_f((x-pi/2)/pi);
dr = matlabFunction(diff(r(x)));
drr = matlabFunction(diff(dr(x)));

z   = @(t) r(t).*cos(t) + r(t).*sin(t)*1i;
dz  = @(t) (dr(t).*cos(t) - r(t).*sin(t)) + ...
           (dr(t).*sin(t) + r(t).*cos(t))*1i;
dzz = @(t) (drr(t).*cos(t) - 2*dr(t).*sin(t) - r(t).*cos(t)) + ...
           (drr(t).*sin(t) + 2*dr(t).*cos(t) - r(t).*sin(t))*1i;

xs = real(z(t));
ys = imag(z(t));
dxs = real(dz(t));
dys = imag(dz(t));
d2xs = real(dzz(t));
d2ys = imag(dzz(t));

r = [(xs(:)).'; (ys(:)).'];
d = [(dxs(:)).'; (dys(:)).'];
d2 = [(d2xs(:)).'; (d2ys(:)).'];

end
