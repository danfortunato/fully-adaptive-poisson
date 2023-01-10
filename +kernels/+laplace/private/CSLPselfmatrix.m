function A = CSLPselfmatrix(source, side)
% complex SLP Kress-split Nystrom matrix
% s = src seg, even # nodes. side = 'i'/'e' int/ext case. Barnett 10/18/13.
% only correct for zero-mean densities.

N = source.N;
z     = cell2mat(source.z);
dz    = cell2mat(source.dz);
s     = cell2mat(source.s);
speed = cell2mat(source.speed);

d = bsxfun(@minus,z,z.'); % C-# displacements mat
x = exp(1i*s); % unit circle nodes relative to which angles measured.
A = -log(d) + log(x*ones(1,N) - ones(N,1)*x.'); % NB circ angles
A(1:N+1:N^2) = log(1i*x./dz);               % complex diagonal limit
% O(N^2) hack to remove 2pi phase jumps in S (assumes enough nodes for smooth):
for i=1:numel(A)-1
    p = imag(A(i+1)-A(i)); % phase jump between pixels
    A(i+1) = A(i+1) - 2i*pi*round(p/(2*pi));
end

m = 1:N/2-1;
% cmplx Kress Rj(N/2)/4pi
if ( side == 'e')
    Rjn = ifft([0 1./m 1/N 0*m]); % imag sign dep on side
else
    Rjn = ifft([0 0*m 1/N 1./m(end:-1:1)]);
end
A = A/N + util.circulant(Rjn); % incl 1/2pi SLP prefac. drops const part in imag
A = bsxfun(@times, A, speed.');  % include speed (2pi/N weights already in)

end
