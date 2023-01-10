function [step, bump] = blender(w, type, dom)

if ( nargin < 1 ), w = 30;        end
if ( nargin < 2 ), type = 'pswf'; end
if ( nargin < 3 ), dom = [-1 1];  end

switch lower(type)
    case 'pswf'
        bump = pswf(0, w);
    case 'dpss'
        n = max(2*w+1, 1000);
        %x = slepian(n, w/n);
        %x = dpss(n, w/4, 1);
        x = dpss(n, w/2, 1);
        bump = chebfun(x, 'equi');
    case 'gaussian'
        bump = chebfun(@(x) exp(-w*x.^2/2));
    otherwise
        error('Unknown blending function.');
end

% Normalize
bump = bump / bump(0);

% Scale to the given domain
bump = chebfun(@(x) bump((x-dom(1))/diff(dom)*2-1), dom);

% Compute a step function
step = cumsum(bump) / sum(bump);

end

function x = slepian(N, width)

% If you have the signal processing toolbox, this is a normalized version
% of dpss(N, N*width/4, 1)

if ( nargin == 1 )
    width = 1;
end

m = 0:N-1;
H = zeros(2,N);
H(1,2:end) = m(2:end) .* (N - m(2:end)) / 2;
H(2,:) = ((N-1-2*m)/2).^2 * cos(2*pi*width/4);
A = spdiags([[H(1,2:end) 0].' H(2,:).' [0 H(1,2:end)].'], -1:1, N, N);
[x,~] = eig(full(A));
x = abs(x(:,N));
x = x/max(x);

end
