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

%%%%%%
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
