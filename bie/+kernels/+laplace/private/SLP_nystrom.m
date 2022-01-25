function A = SLP_nystrom(source, side)
%SLP_NYSTROM   Build Nystrom matrix for single-layer potential.

% Unpack things
N     = source.N * source.np;
z     = cell2mat(source.z);
dz    = cell2mat(source.dz);
s     = cell2mat(source.s);
speed = cell2mat(source.speed);
w = cell2mat(source.w);

d = bsxfun(@minus,z,z.'); % C-# displacements mat
if ( strcmpi(source.rule, 'ptr') )
    % Fill in single-layer potential using Kress quadrature
    %x = exp(1i*s); % unit circle nodes relative to which angles measured.
    %A = -log(d) + log(x*ones(1,N) - ones(N,1)*x.'); % NB circ angles
    %A(1:N+1:end) = log(1i*x./dz);               % complex diagonal limit
    A = -log(abs(d)) + util.circulant(0.5*log(4*sin(pi*(0:N-1)/N).^2));   % peri log
    A(1:N+1:end) = -log(speed);                 % diagonal limit
    m = 1:N/2-1;
    Rjn = ifft([0 1./m 2/N 1./m(end:-1:1)])/2;  % Kress Rj(N/2)/4pi
    % cmplx Kress Rj(N/2)/4pi
    %if ( side == 'e' )
    %    Rjn = ifft([0 1./m 1/N 0*m]); % imag sign dep on side
    %else
    %    Rjn = ifft([0 0*m 1/N 1./m(end:-1:1)]);
    %end
    A = A/N + util.circulant(Rjn);              % includes SLP prefac 1/2pi. Kress peri log matrix L
    A = bsxfun(@times, A, speed.');             % do speed factors (2pi/N weights already)
else
    error('SLP self-interaction is not implemented for panels.');
    %A = -(1/2/pi) * log(abs(d));
    %w = cell2mat(source.w);
    %A = A .* repmat(w.', [N 1]);
    %A(1:N+1:end) = 0;
end

end
