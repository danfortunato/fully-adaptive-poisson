function C = squished_circle(N, r, b, rot)

if ( nargin < 2 )
    r = 1;
end
if ( nargin < 3 )
    b = 1;
end
if ( nargin < 4 )
    rot = 0;
end

t = trigpts(N, [0 2*pi]);
c = r*exp(1i*rot)*(exp(1i*t)-1i*(1-b)*sin(t).^3);
C = Boundary([real(c) imag(c)]);

end
