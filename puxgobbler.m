function z = puxgobbler(t,p)
% PUXGOBBLER  Smooth parameterized boundary with exterior near-intersections
%
% z = puxgobbler(t,p) returns boundary points of multiscale domain, with
%  fine features accumulating at the origin in space, and in t.
% Inputs:
%  t: parameter in [0,2pi), or array of such.
%  p: optional parameter struct with optional fields (defaults used otherwise):
%     p.tol - allowed jump discontinuity relative error size ("noise floor")
%     p.delta - smallest spatial scale, say 1e-12 to 1e-1
%     p.align : +1 for everywhere-close-touching-teeth, -1 for canine-only close
%     p.a - amplitude (aspect ratio of teeth), around 1 is good
%     p.f - freq (per exp growth), around 4
%     p.g - gap, rel. to ampl, eg << 1
%     p.ws - blend width at small end, in radians of osc. Around 10 is good
%     p.wb - blend width at large "
% Outputs z, same shape as t, giving the points in complex plane whose
%  parameters are t.
%
%  Notes:
%  1) coords only, no t-derivs.
%  2) outer part is swept over fast in t.
%
% Without arguments, does self-test

% Barnett 12/13/21
if nargin==0, test_puxgobbler; return; end
if nargin<2, p=[]; end
if ~isfield(p,'tol'), p.tol = 1e-12; end
if ~isfield(p,'delta'), p.delta = 1e-3; end
if ~isfield(p,'align'), p.align = 1; end
if ~isfield(p,'a'), p.a = 1; end
if ~isfield(p,'f'), p.f = 3; end
if ~isfield(p,'g'), p.g = 0.1; end
if ~isfield(p,'ws'), p.ws = 10; end
if ~isfield(p,'wb'), p.wb = 10; end

L = log(1/p.delta);   % horiz width before exponentiating
be = L*p.f;   % freq wrt t, so that p.f is wrt Re y
             % note 1/4-period shift in break-pt for Re part here: (see periblender test)
t1 = @(t) mod(t+pi/2,2*pi)-pi/2;  % no mod-breaks in nbhd of [0,pi]
f1 = @(t) L*(-1+t1(t)/pi) + 1i*p.a/p.f*sin(be*t1(t)).^2;
t2 = @(t) mod(t-pi/2,2*pi)+pi/2;  % no mod-breaks in nbhd of [pi,2pi]
f2 = @(t) L*(1-t2(t)/pi) + 1i*(2*pi - p.g*p.a/p.f) + (p.align*1i*p.a/p.f)*sin(be*(2*pi-t2(t))).^2;
% Set up blending intervals...
as = {-p.ws/be, pi-p.wb/be};
bs = {pi+p.wb/be, p.ws/be};
y = periblender(t,{f1,f2}, as, bs, p.tol);   % do POU
z = exp(y);         % boom!

%%%%%
function test_puxgobbler
n=1e4; t = (1:n)*(2*pi/n);
figure(1); clf; figure(2); clf;
for align = [-1,1]     % two shapes
  p = []; p.align = align;
  z = puxgobbler(t,p);
  figure(1);
  subplot(2,2,2+align); plot(log(z),'.'); axis tight; title('log z')
  subplot(2,2,3+align); plot(z,'.'); axis equal tight
  title(sprintf('z, puxgobbler (align=%d)',align))
  figure(2); subplot(2,1,1+(1+align)/2);
  ff = abs(fft(z));  % noise floor should be tol
  semilogy(ff(1:n/2+1),'.'); axis tight; title('Fourier coeffs, series in t')
end
