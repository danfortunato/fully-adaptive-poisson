function [x, y] = smoothStrip2(dom, n, beta, bleed, width)

if ( nargin == 0 )
    n = 16;
    dom = Boundary.star(n, 'quadrature', 'panel');
    dom = refine(dom, 1);
end

if ( nargin < 3 )
    beta = 4;
end

if ( nargin < 4 )
    bleed = 2;
end

if ( nargin < 5 )
    width = 0.5;
end

panelsize = sum([dom.w{:}]).';
%panelsize = [panelsize(end); panelsize];
panelsize = windowed_average(panelsize, 3, 'geom');

%panelsize = [panelsize(end); panelsize; panelsize(1)];
%panelsize = (panelsize(1:end-1) + panelsize(2:end)) / 2;
%panelsize = sqrt(panelsize(1:end-1).*panelsize(2:end));

a = dom.breaks;
%f = [panelsize(end); panelsize];
f = [panelsize; panelsize(1)];
%f = panelsize;
[x, v] = mollify(a, f, beta, bleed, n);

% plot(dom.breaks, 2*panelsize, '-o')
% axis equal
% shg

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
        k
        abs(v{k}(end) - v{k+1}(1))
        warning('Asymptotics not matched.');
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

% dom1 = Boundary(z1);
% plot(dom1)
% hold on
% plot(dom)
% axis equal
% shg

end

function avg = windowed_average(f, k, type)

if ( nargin < 2 )
    k = 2;
end

if ( nargin < 3 )
    type = 'arith';
end

if ( k == 0 || k == 1)
    avg = f;
    return
end

if ( strcmpi(type, 'geom') )
    geom = true;
    op = @times;
else
    geom = false;
    op = @plus;
end

kl = floor(k/2);
kr = floor(k/2)-mod(k+1,2);

n = numel(f);
g = f(:);
g = [g(n-kl+1:n); g(:); g(1:kr)];
avg = g(1:n);
for j = 1:kl+kr
    avg = op(avg, g((1:n)+j));
end

if ( geom )
    avg = avg.^(1/k);
else
    avg = avg / k;
end
avg = reshape(avg, size(f));

end
