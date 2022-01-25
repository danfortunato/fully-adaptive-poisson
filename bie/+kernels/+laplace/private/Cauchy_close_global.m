function [vc, vcp, vcpp] = Cauchy_close_global(source, target, vb, side)
%CAUCHY_CLOSE_GLOBAL   Globally compensated barycentric int/ext Cauchy integral
%
% This is a spectrally-accurate close-evaluation scheme for Cauchy integrals.
%  It returns approximate values (and possibly first derivatives) of a function
%  either holomorphic inside of, or holomorphic and decaying outside of, a
%  closed curve, given a set of its values on nodes of a smooth global
%  quadrature rule for the curve (such as the periodic trapezoid rule).
%  This is done by approximating the Cauchy integral
%
%       v(x) =  +- (1/(2i.pi)) integral_Gamma v(y) / (x-y) dy,
%
%  where Gamma is the curve, the sign is + (for x interior) or - (exterior),
%  using special barycentric-type formulae which are accurate arbitrarily close
%  to the curve.
%
%  By default, for the value these formulae are (23) for interior and (27) for
%  exterior, from [hel08], before taking the real part.  The interior case
%  is originally due to [ioak]. For the derivative, the Schneider-Werner formula
%  (Prop 11 in [sw86]; see [berrut]) is used to get v' at the nodes, then since
%  v' is also holomorphic, it is evaluated using the same scheme as v.  This
%  "interpolate the derivative" suggestion of Trefethen (personal communication,
%  2014) contrasts [lsc2d] which "differentes the interpolant". The former gives
%  around 15 digits for values v, 14 digits for interior derivatives v', but
%  only 13 digits for exterior v'. (The other [lsc2d] scheme gives 14 digits
%  in the last case; see options below).  The paper [lsc2d] has key background,
%  and is helpful to understand [hel08] and [sw86].  This code replaces code
%  referred to in [lsc2d].
%
%  The routine can (when vb is empty) instead return the full M-by-N dense
%  matrices mapping v at the nodes to values (and derivatives) at targets.
%
% Basic use:  v = Cauchy_close_global(x,s,vb,side)
%             [v, vp] = Cauchy_close_global(x,s,vb,side)
%             [v, vp, vpp] = Cauchy_close_global(x,s,vb,side)
%
% Inputs:
%  source = CloseCurve object
%  target = row or col vec of M target points in complex plane
%  vb = col vec (or stack of such) of N boundary values of holomorphic function
%       v. If empty, causes outputs to be the dense matrix/matrices.
%  side = 'i' or 'e' specifies if all targets interior or exterior to curve.
%
% Outputs:
%  v  = col vec (or stack of such) approximating the homolorphic function v
%       at the M targets
%  vp = col vec (or stack of such) approximating the complex first derivative
%       v' at the M targets
%  vpp = col vec (or stack of such) approximating the complex 2nd derivative
%       v'' at the M targets
%
% Without input arguments, a self-test is done outputting errors at various
%  distances from the curve for a fixed N. (Needs setupquad.m)
%
% Notes:
% 1) For accuracy, the smooth quadrature must be accurate for the boundary data
%  vb, and it must come from a holomorphic function (be in the right Hardy
%  space).
% 2) For the exterior case, v must vanish at infinity.
% 3) The algorithm is O(NM) in time. In order to vectorize in both sources and
%  targets, it is also O(NM) in memory - using loops this could of course be
%  reduced to O(N+M). If RAM is a limitation, targets should be blocked into
%  reasonable numbers and a separate call done for each).
%
% If vb is empty, the outputs v and vp are instead the dense evaluation
%  matrices, and the time cost is O(N^2M) --- this should be rewritten.
%
% References:
%
%  [sw86]  C. Schneider and W. Werner, Some new aspects of rational
%          interpolation, Math. Comp., 47 (1986), pp. 285–299
%
%  [ioak]  N. I. Ioakimidis, K. E. Papadakis, and E. A. Perdios, Numerical
%          evaluation of analytic functions by Cauchy’s theorem, BIT Numer.
%          Math., 31 (1991), pp. 276–285
%
%  [berrut] J.-P. Berrut and L. N. Trefethen, Barycentric Lagrange
%          interpolation, SIAM Review, 46 (2004), pp. 501-517
%
%  [hel08] J. Helsing and R. Ojala, On the evaluation of layer potentials close
%          to their sources, J. Comput. Phys., 227 (2008), pp. 2899–292
%
%  [lsc2d] Spectrally-accurate quadratures for evaluation of layer potentials
%          close to the boundary for the 2D Stokes and Laplace equations,
%          A. H. Barnett, B. Wu, and S. Veerapaneni, SIAM J. Sci. Comput.,
%          37(4), B519-B542 (2015)   https://arxiv.org/abs/1410.2187
%
% Todo: * allow mixed interior/exterior targets, and/or auto-detect this.
% * O(N) faster matrix filling version!
% * Think about if interface should be t.x.
% Note in/output format changed to col vecs, 6/27/16

% (c) Alex Barnett, June 2016, based on code from 10/22/13. Blocked 8/2/16
% 2nd deriv added, Oct. 2018, Jun Wang, Flatiron Inst.
% Barnett sped up by at least 20x, via all O(N^2.Nc) via GEMM, tidied, 5/4/19
% Adapted by Dan Fortunato, 2019

M = size(target,1);
N = source.N;
Nc = size(vb,2); % # input col vecs
x  = cell2mat(source.x);
nx = cell2mat(source.normal);
cw = cell2mat(source.cw);

% Complexify
if ( size(x,2) == 2 ), x = x(:,1) + x(:,2)*1i; end
if ( size(nx,2) == 2 ), nx = nx(:,1) + nx(:,2)*1i; end
if ( size(cw,2) == 2 ), cw = cw(:,1) + cw(:,2)*1i; end
if ( size(target,2) == 2 ), target = target(:,1) + target(:,2)*1i; end
 
% Note: only multi-vector version now (faster even for Nc=1 single-vec).
% (Note: is non-optimal as method for matrix filling when case vb=sparse)
% Do bary interp for value outputs:
% Precompute weights in O(NM)... note sum along 1-axis faster than 2-axis...
comp = repmat(cw, [1 M]) ./ (repmat(x,[1 M]) - repmat(target(:).',[N 1]));
% mult input vec version (transp of Wu/Marple): comp size N*M, I0 size M*Nc
I0 = blockedinterp(vb,comp);    % local func, directly below
J0 = sum(comp).';  % size N*1, Ioakimidis notation
if side=='e', J0 = J0-2i*pi; end                      % Helsing exterior form
vc = I0./(J0*ones(1,Nc));                 % bary form (multi-vec), size M*Nc
[jj, ii] = ind2sub(size(comp),find(~isfinite(comp)));  % node-targ coincidences
for l=1:numel(jj), vc(ii(l),:) = vb(jj(l),:); end     % replace each hit w/ corresp vb

if nargout>1   % 1st deriv also wanted... Trefethen idea first get v' @ nodes
  Y = 1 ./ bsxfun(@minus,x,x.'); Y(1:N+1:N^2)=0;  % j.ne.i Cauchy mat
  Y = Y .* (cw*ones(1,N));      % include complex wei over 1st index
  Y(1:N+1:N^2) = -sum(Y).';     % set diag to: -sum_{j.ne.i} w_j/(y_j-y_i)
  vbp = Y.'*vb;                 % v' @ nodes, size N*Nc
  if side=='e', vbp = vbp + 2i*pi*vb; end    % S-W variant derived 6/12/16...
  vbp = ((-1./cw)*ones(1,Nc)).*vbp;
  % now again do bary interp of v' using its value vbp at nodes...
  I0 = blockedinterp(vbp,comp);
  J0 = sum(comp).';
  if side=='e', J0 = J0-2i*pi; end                    % Helsing exterior form
  vcp = I0./(J0*ones(1,Nc));                          % bary form
  for l=1:numel(jj), vcp(ii(l),:) = vbp(jj(l),:); end % replace hits w/ vbp
end

if nargout>2   % 2nd deriv, mult-col version; we use vcp and Y from 1st deriv
  vbpp = Y.'*vbp;                 % v'' @ nodes, size N*Nc
  if side=='e', vbpp = vbpp + 2i*pi*vbp; end  % S-W variant derived 6/12/16
  vbpp = ((-1./cw)*ones(1,Nc)).*vbpp;
  % now again do bary interp of v'' 
  I0 = blockedinterp(vbpp,comp);
  % J0 computed above in the nargout>1 case
  vcpp = I0./(J0*ones(1,Nc));                          % bary form
  for l=1:numel(jj), vcpp(ii(l),:) = vbpp(jj(l),:); end % replace hits w/ vbp
end

end

function I0 = blockedinterp(vb,comp)   % ....................................
% perform barycentric interpolation using precomputed comp wei mat, used in
% multi-density vec version above. Output: I0 (size M*Nc).
% Barnett 5/4/19
I0 = (vb.'*comp).';    % unbelievable we didn't notice this GEMM earlier :)
end
