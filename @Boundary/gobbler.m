function C = gobbler(N, varargin)
%GOBBLER   Set up a gobbler.

parser = inputParser;
parser.KeepUnmatched = true;
addRequired(parser, 'N', @isfloat);
addParameter(parser, 'delta', 1e-3, @isfloat);
addParameter(parser, 'align', 1,    @isfloat);
addParameter(parser, 'amp',   1,    @isfloat);
addParameter(parser, 'freq',  3,    @isfloat);
addParameter(parser, 'gap',   0.1,  @isfloat);
parse(parser, N, varargin{:});

p = struct();
p.delta = parser.Results.delta;
p.align = parser.Results.align;
p.a     = parser.Results.amp;
p.f     = parser.Results.freq;
p.g     = parser.Results.gap;
p.tol = 1e-12;

z = chebfun(@(t) puxgobbler(t,p), [0 2*pi], 'trig');
dz = diff(z);
dzz = diff(dz);

levelset = [];
C = Boundary(N, z, dz, dzz, levelset, varargin{:});

end

%% puxgobbler

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

end

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

end

%% periblender

function y = periblender(t,fs,as,bs,tol)
% PERIBLENDER  Blend a set of 1D functions periodically using POU via erfs
%
% y = periblender(t,fs,as,bs,tol) evaluates blended function y = f(t), returning
%  array same shape as t the list of arguments, interpreted 2pi-periodically.
%  fs is a cell array of function handles. The function fs{p} is used over the
%  (2pi-periodically-wrapped) interval [as{p}, bs{p}), which must overlap
%  the corresponding interval for p-1 and p+1. Thus the endpoints must
%  alternate a2, b1, a3, b2, a4, b3, ... cyclically when taken mod 2pi, and
%  cyclically with respect to index in {1,2,..,numel(fs)}. This is not yet
%  tested.
%  tol sets the blending function approximate relative jump discontinuity.
%  The periodicity is 2pi, also assumed for the funcs; they are eval
%  in [0,2pi) or one or two subintervals of it.
%
% Without argument, does self-test (see end of code for demo)
%
% Self-test needs: VLINE

% Barnett 12/13/21
if nargin==0, test_periblender; return; end

erftail = sqrt(log(0.1/tol));        % how far to go into tail of erf
blend = @(t) (1+erf((t-0.5)*(2*erftail)))/2;   % maps [0,1] to about [tol,1-tol]
% figure; t=linspace(0,1,1e3); plot(t,blend(t),'-');   % test

siz = size(t);
t = mod(t(:),2*pi);   % col vec, wrap it
y = 0*t;
n = numel(fs);
for p=1:n    % loop over partitions
  prevp = mod(p-2,n)+1;
  nextp = mod(p,n)+1;
  % the blending part...
  a = mod(as{p},2*pi);        % blend fs{p} with fs{prevp} over [a,b)...
  b = mod(bs{prevp},2*pi);
  %fprintf('blend a,b = %.3g,%.3g\n',a,b)
  if b>a     % no wrap-around of blend interval
    jj = t>=a & t<b;
    psi = blend((t(jj)-a)/(b-a));
    y(jj) = psi.*fs{p}(t(jj)) + (1-psi).*fs{prevp}(t(jj));    % do blend
  else           % interval wraps around at 2pi
    jj = t>=a;
    psi = blend((t(jj)-a)/(2*pi+b-a));   % part just before wrapping, [a,2pi)
    y(jj) = psi.*fs{p}(t(jj)) + (1-psi).*fs{prevp}(t(jj));    % do blend
    jj = t<b;
    psi = blend((t(jj)-a+2*pi)/(2*pi+b-a));   % after wrapping, [0,b)
    y(jj) = psi.*fs{p}(t(jj)) + (1-psi).*fs{prevp}(t(jj));    % do blend
  end
  % the pure func part...
  a = mod(bs{prevp},2*pi);
  b = mod(as{nextp},2*pi);        % use fs{p} over [a,b)... (note names swap)
  %fprintf('f%d  a,b = %.3g,%.3g\n',p,a,b)
  if b>a     % no wrap-around of interval
    jj = t>=a & t<b;
  else       % wrapping
    jj = t>=a | t<b;
  end
  y(jj) = fs{p}(t(jj));
end
y = reshape(y,siz);

end

function test_periblender

fs = {@(t) sin(t), @(t) 0*t + 2};
as = {1,2.5};   % start pts of valid region for each func
bs = {3, 1.7};  % end pts  "

fs = {@(t) -1+t/pi, @(t) 1-t/pi};  % use case for puxgobbler.. tricky
fs = {@(t) -1+mod(t/pi+.5,2)-.5, @(t) 1-mod(t/pi-.5,2)-.5};  % ...fixed :)
as = {-1, pi-.5}; bs = {pi+.5, 1};

n=1e3; tt = (2*pi/n)*(1:n);
tol = 1e-12;
figure;

%sh = 0.1;   % delta shift
%for i=1:1    % test periodic shifting
%for p=1:numel(as), as{p}=mod(as{p}+sh,2*pi); bs{p}=mod(bs{p}+sh,2*pi); end

y = periblender(tt,fs,as,bs,tol);
plot(tt,[fs{1}(tt); fs{2}(tt)], '-'); hold on;
plot(tt,y,'-','linewidth',3); hold off;
legend('f_1','f_2','blended');
vline(vertcat(as{:}),'b:');
vline(vertcat(bs{:}),'r:');
drawnow; %pause(1/60);  % matlab's plot really slow anyway :(

end
