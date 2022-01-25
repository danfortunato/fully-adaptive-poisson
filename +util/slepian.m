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
