% Build an adaptive strip
%#ok<*UNRCH>

panels = true;

if ( panels )
    p = 8;
    dom = Boundary.star(p, 'quadrature', 'panel');
    %dom = Boundary.hand(p, 'quadrature', 'panel');
    %dom = Boundary.circle(p, 'quadrature', 'panel');
else
    p = 100;
    dom = Boundary.star(p);
end

%%
x = dom.breaks;
scl = min(x(x>0));
w = diff(x);
scl = min(w);
w = [(w(1)+w(end))/2 (w(1:end-1)+w(2:end))/2 (w(1)+w(end))/2];

const = 0.1;
smooth = 0.1;
op = rbfcreate(x, w, 'RBFFunction', 'multiquadric', 'RBFConstant', 1+0*const, 'RBFSmooth', smooth);
rbfcheck(op);
xx = linspace(x(1), x(end), 10000);
ww = rbfinterp(xx, op);
plot(x,w,'o')
hold on
plot(xx,ww,'-')
hold off
%ylim([0 4*scl])
shg

fun = chebfun(@(x) rbfinterp(x, op), [x(1), x(end)]);
length(fun)

%%
dom_strip = [0 2*pi 0 1];
dom1 = perturb(dom, @(t) -fun(t), dom_strip);

%%
close all
plot(dom)
hold on
plot(dom1)
hold off
shg

%% Define the strip region
width = -0.2;
dom_strip = [0 2*pi 0 1];
curv = chebfun(dom.curvature, [0 2*pi], 'trig');
dff = chebfun(dom.dff, [0 2*pi], 'trig');
%[Gamma1, xyfun] = perturb(Gamma, @(t) width-0.01.*abs(Gamma.dff(t)), dom_strip);
[dom1, xyfun] = perturb(dom, @(t) width-0.1.*(max(abs(curv))-abs(curv(t))), dom_strip);
%[Gamma1, xyfun] = perturb(Gamma, @(t) width-0.01.*(max(abs(dff))-abs(Gamma.dff(t))), dom_strip);