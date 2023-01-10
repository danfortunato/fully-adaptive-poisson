function Gamma1 = adaptivePerturb2(Gamma)

n = Gamma.N;
np = Gamma.np;
%beta = 4;
%beta = 6;
%beta = 20; bleed = 7; 3/8/2022
beta = 2;
%beta = 40;
%bleed = 22;
bleed = 5;
width = 0.5;
[x, y] = smoothStrip2(Gamma, beta, bleed, width);
xleg = chebvals2legvals(reshape(x, n, np));
yleg = chebvals2legvals(reshape(y, n, np));
z1 = mat2cell(xleg(:) + 1i*yleg(:), repmat(n, np, 1), 1);

Gamma1 = Boundary(z1);

end

function [x, y] = smoothStrip2(dom, beta, bleed, width)

if ( nargin < 3 ), beta  = 4;   end
if ( nargin < 4 ), bleed = 2;   end
if ( nargin < 5 ), width = 0.5; end

panelsize = arclength(dom, 1:dom.np);
panelsize = util.windowed_average(panelsize, 3, 'arith');

n = dom.N;
a = dom.breaks;
f = [panelsize; panelsize(1)];
[x, v] = mollify5(a, f, beta, n);

% plot(a, f, '-o'), hold on
% plot(cell2mat(x), cell2mat(v), '-o')
% shg

z1 = cell(dom.np, 1);
v = [v; v(1)];
CV2LV = chebvals2legvals(eye(dom.N));
LV2CV = legvals2chebvals(eye(dom.N));
worst = 0;
for k = 1:dom.np
    worst = max(worst, abs(v{k}(end) - v{k+1}(1)));
    if ( abs(v{k}(end) - v{k+1}(1)) > 1e-10 )
        abs(v{k}(end) - v{k+1}(1))
        warning('Asymptotics not matched on panel %d.', k);
    end
    vk = CV2LV * v{k};
    nx = dom.normal{k}(:,1) + dom.normal{k}(:,2)*1i;
    nx = nx./abs(nx);
    z1{k} = dom.z{k} - width * vk .* nx;
    z1{k} = LV2CV * z1{k};
end
fprintf('worst = %g\n', worst);

x = real(cell2mat(z1));
y = imag(cell2mat(z1));

end

