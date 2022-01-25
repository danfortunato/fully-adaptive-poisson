function u = DLP_eval_close_global(source, target, sigma, side)
%DLP_EVAL_CLOSE_GLOBAL   Evaluate double-layer potential at target
%points using global close evaluation.

N  = source.N;
z  = cell2mat(source.z);
s  = cell2mat(source.s);
cw = cell2mat(source.cw);

% Helsing step 1: eval bdry limits at nodes of v = complex DLP(tau)...
% (note to future fortran/C versions: this also needs to be able to handle
% complex input, real output, since the Stokes DLP feeds that in)

sfun = chebfun(sigma, [0 2*pi], 'trig');
dsigma = feval(diff(sfun), s); % numerical deriv of density

% BLAS3 sped-up version of Eqn (4.2) in [lsc2d], beats socks off bsxfun...
Y = 1 ./ bsxfun(@minus, z, z.');
Y(1:N+1:end) = 0;              % j.ne.i Cauchy mat
Y = Y .* (cw*ones(1,N));       % include complex wei over 1st index
Y(1:N+1:end) = -sum(Y).';      % set diag to: -sum_{j.ne.i} w_j/(y_j-y_i)
vb = Y.'*sigma*(1/(-2i*pi));   % v @ nodes, size N, w/ prefactor
vb = vb - (1/(1i*N))*dsigma;   % diagonal sigma' correction, last term in (4.2).
if ( side == 'i' )
    vb = vb - sigma; % JR's add for v^-, cancel for v^+ (Eqn 4.3)
end

% Helsing step 2: compensated close-evaluation of v, followed by take Re:
v = Cauchy_close_global(source, target, vb, side); % does Sec. 3 of [lsc2d]
u = real(v);

end
