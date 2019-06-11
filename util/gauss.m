function [x,w] = gauss(N, a, b)
%GAUSS   Gauss-Legendre nodes and weights on [a,b].

if ( nargin < 3 )
    b = 1;
    if ( nargin < 2 )
        a = -1;
    end
end

beta = 0.5./sqrt(1-(2*(1:N-1)).^(-2));
T = diag(beta,1) + diag(beta,-1);
[V,D] = eig(T);
x = diag(D);
[x,i] = sort(x);
w = 2*V(1,i).'.^2;

% Map to [a,b]
x = (x+1)*(b-a)/2+a;
w = w*(b-a)/2;

end
