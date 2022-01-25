function [x, y] = smoothStrip(dom, n, width)

% Run a test if given no arguments
if ( nargin == 0 )
    [varargout{1:nargout}] = test();
    if ( nargout > 0 )
        [x, y] = deal(varargout{:});
    end
    return
end

if ( nargin < 3 )
    width = 0.4;
end

% Perturb the panel breaks of the original domain based on the local panel
% size
z = perturbBreaks(dom, width);
%z = dom.zbreaks(1:end-1).';

% Run the 2D surface smoother to get a smooth curve
[x, y] = smoother2D([real(z) imag(z)], n, n);

end

function [x, y] = test()

n = 16;
%dom = Boundary.star(n, 'quadrature', 'panel');
%dom = Boundary.circle(n, 'quadrature', 'panel', 'panels', 20);
%dom = Boundary.ellipse(n, 5, 1, 'quadrature', 'panel');
%dom = Boundary.multiscale(n, 'quadrature', 'panel');
dom = Boundary.multiscale_circle(n, 'quadrature', 'panel');
dom = refine(dom);
[x, y] = smoothStrip(dom, n);

zz = mat2cell(x+1i*y, repmat(n,dom.np,1), 1);
rad = chebpts(n);
vol = cell(dom.np,1);
for k = 1:dom.np
    %zz{k} = legvals2chebvals(real(zz{k})) + 1i*legvals2chebvals(imag(zz{k}));
    dom.z{k} = legvals2chebvals(real(dom.z{k})) + 1i*legvals2chebvals(imag(dom.z{k}));
    vol{k} = (rad+1)/2.*dom.z{k}.' + (1-rad)/2.*zz{k}.';
    plot(vol{k}, 'k.'), hold on
end

plot(x,y,'r-','linewidth',1)
hold on
plot(dom)
axis equal
shg

end
