function A = circulant(x)
%CIRCULANT   Square circulant matrix with first row x
% Barnett 2/5/08

x = x(:);
A = toeplitz([x(1); x(end:-1:2)], x);

end
