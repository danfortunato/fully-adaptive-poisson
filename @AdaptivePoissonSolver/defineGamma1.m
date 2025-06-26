function Gamma1 = defineGamma1(Gamma, opts)

defaults = [];
defaults.debug = true;
defaults.beta = 1;
defaults.bleed = 5;
defaults.average = 5;
defaults.stripWidth = 0.5;

if ( nargin < 2 )
    opts = defaults;
else
    opts = setDefaults(opts, defaults);
end

n = Gamma.N;
np = Gamma.np;
[aleg, a] = arclengthparam(Gamma);
panelsize = arclength(Gamma, 1:np);
panelsize = util.windowed_average(panelsize, opts.average, 'arith');
f = [panelsize; panelsize(1)];
[x, v] = mollify_fast(a, f, opts.beta, n, aleg);

z1 = cell(np, 1);
v = [v v(:,1)];
worst = 0;
for k = 1:np
    err = abs(v(end,k) - v(1,k+1));
    worst = max(worst, err);
    % if ( err > 1e-2 )
    %     k
    %     err
    %     warning('Asymptotics not matched.');
    % end
    vk = v(:,k);
    nx = Gamma.normal{k}(:,1) + Gamma.normal{k}(:,2)*1i;
    nx = nx./abs(nx);
    z1{k} = Gamma.z{k} - opts.stripWidth * vk .* nx;
end
%dbprint(opts.debug, 'worst = %g\n', worst);
x = real(cell2mat(z1));
y = imag(cell2mat(z1));

xleg = reshape(x, n, np);
yleg = reshape(y, n, np);
z1 = mat2cell(xleg(:) + 1i*yleg(:), repmat(n, np, 1), 1);
Gamma1 = Boundary(z1);

end

%% Mollify

function [aleg, fleg, sfunsx, x] = mollify_fast(a, f, beta, n, aleg)

% We assume the interval [a(1), a(end)] is periodic.
% To handle the periodicity, we append the first panel to the end.
a = a(:).';
f = f(:).';
npanel = length(a)-1;
scl = a(end)-a(1);

nnbr = 5;
bleed = 1;
nextra = 2*nnbr;
padL = npanel-(nnbr+bleed):npanel;
padR = 2:(nnbr+bleed)+3;
aa = [a(padL)-scl  a a(padR)+scl];
ff = [f(padL)      f f(padR)];
ss = diff(ff) ./ diff(aa); % Slopes
npad = length(aa) - length(a);

aa3 = zeros(1, npanel, npad);
ff3 = zeros(1, npanel, npad);
ss3 = zeros(1, npanel, npad);
for k = 1:npanel
    aa3(:,k,:) = aa(k:k+npad-1);
    ff3(:,k,:) = ff(k:k+npad-1);
    ss3(:,k,:) = ss(k:k+npad-1);
end

x = zeros(n, npanel, npad-1);
for k = 1:npanel
    for j = 1:npad-1
        x(:,k,j) = aleg(:,k);
    end
end

% Make rounded approximations to each kink, which occurs between panels
% with endpoints [k-1, k] and [k, k+1]. The rounding scales with the local
% panel size.
kk = 1:npad-1;
scale = ff3(:,:,kk+1);
b = beta./scale; % If reparametrize(Gamma)
%b = beta;       % If not reparametrize(Gamma)
xa = x - aa3(:,:,kk+1);
rounding = (exp(-(b.*xa).^2)/sqrt(pi) + b.*xa.*(1+erf(b.*xa))) ./ (2*b);
froundx = ff3(:,:,kk+1) + ss3(:,:,kk).*xa + (ss3(:,:,kk+1)-ss3(:,:,kk)).*rounding;
linx    = ff3(:,:,kk+1) - ss3(:,:,kk+1).*xa;

% Combine the rounded approximations from each kink into functions defined
% on each panel. This will add the rounded approximations defined at the
% left and right panel endpoints using matched asymptotics.
% - Panel spans [k, k+1]
kk = 1:npad-2;
sfunsx = froundx(:,:,kk+1) + froundx(:,:,kk) + linx(:,:,kk) - 2*ff3(:,:,kk+1);

% Panel spans [k, k+1]
kk = 2:npad-3;
kl = kk-bleed;
kr = kk+bleed+1;
fsumx = froundx(:,:,kr);
for j = 0:2*bleed
    fsumx = fsumx + froundx(:,:,kl+j) + linx(:,:,kl+j) - 2*ff3(:,:,kl+j+1);
end
sfunsx(:,:,kk) = fsumx;

% Output the blended functions on the set of panel Chebyshev nodes.
step = @(x) (erf(4*x)+1)/2;
fac = 4;
% step = util.blender(30, 'pswf', [-1 1]);
% step = util.piecewise({@(x) 0*x, @(x) step(x), @(x) 1+0*x}, [-1 1]);

fleg = zeros(n, npanel);
for k = 1:npanel
    y = aleg(:,k);
    for j = 1:2*nnbr
        idx = j+k+1;
        delta = 0.5 * (aa(idx+2)-aa(idx)) / 3;
        stepx = step((y-aa(idx+1))/delta/fac);
        fleg(:,k) = (1-stepx).*fleg(:,k) + stepx.*sfunsx(:,k,j+2);
    end
end

end

