function C = multiscale(N, varargin)
%MULTISCALE   Set up a smooth multiscale domain.

p = inputParser;
p.KeepUnmatched = true;
addRequired(p, 'N', @isfloat);
parse(p, N, varargin{:});

e1 = 1e-4;
e2 = 0.2;
b1 = @(t) 1+erf(t/e1-5)-(1+erf(-5));
b2 = @(t) 1+erf((2*pi-t)/e2-5)-(1+erf(-5));
%y = @(t) t.*sin(10*log(t));
y = @(t) real(t.^(7i+1));
r = @(t) 0.5 + b1(t).*b2(t).*y(t)/50;

%r = @(t) 1 + 0*t;

%t = sym('t', 'real');
%assume(t>=0 & t<=2*pi)
%dr  = eval(['@(t) ' vectorize(char(diff(r(t))))  ' + 0*t']);
%drr = eval(['@(t) ' vectorize(char(diff(dr(t)))) ' + 0*t']);
%dr  = matlabFunction(diff(r(t))+0*t);
%drr = matlabFunction(diff(dr(t))+0*t);

%dr = func2str(dr);
%dr = eval(['@(t) ' dr(5:end)]);
%drr = func2str(drr);
%drr = eval(['@(t) ' drr(5:end)]);

%pref = chebfunpref();
%pref.splitting = 1;
%pref.techPrefs.maxLength = 10;
%pref.techPrefs.fixedLength = 10;
%pref.splitPrefs.splitLength = 10;
%pref.techPrefs.chebfuneps = 1e-2;
r = chebfun(r, [0 2*pi]);
dr = diff(r);
drr = diff(r,2);

% z = chebfun(@(t) r(t).*cos(t) + r(t).*sin(t)*1i, [0 2*pi], 'trig');
% dz = diff(z);
% dzz = diff(z,2);

z   = @(t) r(t).*cos(t) + r(t).*sin(t)*1i;
dz  = @(t) (dr(t).*cos(t) - r(t).*sin(t)) + ...
          (dr(t).*sin(t) + r(t).*cos(t))*1i;
dzz = @(t) (drr(t).*cos(t) - 2*dr(t).*sin(t) - r(t).*cos(t)) + ...
          (drr(t).*sin(t) + 2*dr(t).*cos(t) - r(t).*sin(t))*1i;

levelset = [];
%levelset = @(x,y) sqrt(x.^2+y.^2) - r(atan2(y,x));

C = Boundary(N, z, dz, dzz, levelset, varargin{:});

end
