a = 1;
b = 3;
%x = [-1-a -1 1 1+b];
sj = [2.^(-20:0), 2:6];
wj = 0*sj; wj(2:end-1) = (sj(3:end)-sj(1:end-2))/2; wj(end) = 1;
w = wj;
w = [1 1 1 1 2 2 4 8 8 4 8 4 2 2 1 1];
x = cumsum(w);
x = rescale(x, -1, 1);
w = diff(x);
%%w = [(w(1)+w(end))/2 (w(1:end-1)+w(2:end))/2 (w(1)+w(end))/2];
w = [geomean([w(1) w(end)]) geomean([w(1:end-1); w(2:end)]) geomean([w(1) w(end)])];

funs = cell(numel(x)-1,1);
for k = 1:numel(x)-1
    funs{k} = @(t) (w(k+1)-w(k))/(x(k+1)-x(k)).*(t-x(k)) + w(k);
end
chebfun(funs, x);

n = 1000;
[step, bump] = util.makeMask(n, 10);

% minx = -1-a;
% maxx = 1+b;
% bump1 = chebfun(@(x) bump(x+2), [x(1) x(2)]);
% step1 = chebfun(@(x) step(x+2), [x(1) x(2)]);
% 
% bump2 = chebfun(@(x) bump((x+1)/2), [x(2) x(3)]);
% step2 = chebfun(@(x) 1-step((x+1)/2), [x(2) x(3)]);
% 
% bump3 = chebfun(@(x) bump((x+1)/2), [x(2) x(3)]);
% step3 = chebfun(@(x) step((x+1)/2), [x(2) x(3)]);
% 
% bump4 = chebfun(@(x) bump((x-1)/3), [x(3) x(4)]);
% step4 = chebfun(@(x) 1-step((x-1)/3), [x(3) x(4)]);

%pp = spline(x,w);
pp = csape(x, w, 'periodic');
%xx = linspace(x(1), x(end), 1000);
xx = logspace(log10(min(sj)/2), log10(max(sj)),1e4);
ww = ppval(pp,xx);

figure(1)
plot(x,w,'o-')
%loglog(x,w,'bo-')
hold on
%plot(xx,ww,'-')
loglog(xx,ww,'r-')
% plot(w(2)*step1,'b')
% plot(w(2)*step2,'b')
% plot(w(3)*step3,'b')
% plot(w(3)*step4,'b')
% xx = linspace(x(2), x(3), 100);
% plot(xx, w(2)*step2(xx) + w(3)*step3(xx))
hold off

%% RBF interpolation
dom = [0 10];
x = [2.^(-20:0), 2:6];
w = 0*x; w(2:end-1) = (x(3:end)-x(1:end-2))/2;
%xx = logspace(log10(min(x)/2), log10(max(x)),1e4);
%xx = logspace(-10, log10(dom(2)), 1e4);
xx = linspace(dom(1), dom(2), 1e4);

figure(1), clf
op = rbfcreate(x, w, 'RBFFunction', 'multiquadric', 'RBFConstant', 0.5+0*w, 'RBFSmooth', w);
%op = rbfcreate(x, w, 'RBFFunction', 'multiquadric', 'RBFConstant', 0.1+w, 'RBFSmooth', 0.1+0*w);
%op = rbfcreate(x, w, 'RBFFunction', 'gaussian', 'RBFConstant', 0.1+0.5*w, 'RBFSmooth', 0.1*w);
%op = rbfcreate(x, w, 'RBFFunction', 'invmultiquadric', 'RBFConstant', 0.5+0*w, 'RBFSmooth', w);
%op = rbfcreate(x, w, 'RBFFunction', 'thinplate', 'RBFConstant', 0.5+0*w, 'RBFSmooth', w);
%op = rbfcreate(x, w, 'RBFFunction', 'gaussian', 'RBFConstant', 1e-6+w, 'RBFSmooth', 0*w);
rbfcheck(op);
ww = rbfinterp(xx, op);

subplot(211), plot(xx,ww,'-'), hold on
plot([x; x],[0*w; w], 'k.-', 'MarkerSize', 20);
%set(gca,'ylim',[0 2]);
subplot(212), loglog(xx,ww,'-'), hold on
plot(x, w, 'k.', 'MarkerSize', 20);
%set(gca,'ylim',[1e-7 1]);

fun = chebfun(@(x) rbfinterp(x, op), [x(1), x(end)]);
%fun = chebfun(@(x) rbfinterp(x, op), dom, 'trig');
figure(2)
length(fun)
plotcoeffs(fun)
alignfigs

%% Multiscale mollifer
lambda = 1;
%D = diff(x);
%c = 0.5*(x(1:end-1)+x(2:end));
D = w;
c = x;
s = D/lambda;
s0 = 0.1;

dom = [-1 1];
periodize = @(x) mod(x+dom(2), dom(2)-dom(1)) + dom(1);

xx = linspace(-1,1,1000).';
gg = @(x) exp(-0.5*(periodize(x-c)./s).^2);
sigma = @(x) sum(s.*gg(x), 2) ./ sum(gg(x), 2);

plot(xx, sigma(xx)), hold on
plot(c,D,'o'), hold off
shg

%% Gaussian convolution
figure(3)
%sigma = @(w) 0.1;
gh = @(t, w) 1./(sigma(w)*sqrt(2*pi)) .* exp(-0.5*(t./sigma(w)).^2);
%sigma = 0.05;
%g = chebfun(@(t) 1/(sigma*sqrt(2*pi)) * exp(-0.5*(t/sigma).^2), 'trig');
funs = cell(numel(x)-1,1);
for k = 1:numel(x)-1
    funs{k} = @(t) (w(k+1)-w(k))/(x(k+1)-x(k)).*(t-x(k)) + w(k);
end
f = chebfun(funs, x);
h = mycircconv(f, gh);
plot(x,w,'o')
hold on
plot(xx,h(xx),'-')
%xx = chebpts(200);
%plot(xx,h,'-')
hold off

alignfigs

%% Sum of Gaussians
%xx = linspace(x(1), x(end), 1000).';
xx = chebpts(500);
ff = 0*xx;
dd = 0*xx;
close all
hold on

dom = [-1 1];
periodize = @(x) mod(x+dom(2), dom(2)-dom(1)) + dom(1);

%phi = @(t) (1-t.^2).^10;
%phi = @(t) phi(t).*(t.^2<1);
phi = @(t,s) exp(-0.5*(t./s).^2);
%phi = @(t) exp(-1./(1-t.^2));

% Make it periodic
%phi = chebfun(phi, 'trig');

for k = 1:numel(w)
    s = w(k)/2;
    a = w(k);
    %d = exp(-0.5*((xx-x(k))./s).^2);
    p = chebfun(@(t) phi(t,s), 'trig');
    d = p(xx-x(k));
    %plot(xx, a*d)
    ff = ff + a*d;
    dd = dd + d;
end

plot(x,w,'o')
plot(xx,ff./dd,'-')
hold off
shg
