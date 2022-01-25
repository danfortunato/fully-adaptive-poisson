function pass = test_periodic_poisson()

tol = 1e-16;
pass = [];

n = 30;
u = chebfun2(@(x,y) -cos(x).*exp(sin(x)).*sin(y), [-pi pi -2*pi 2*pi], 'trig');
f = lap(u);
v = util.periodicPoisson(f, n);
pass(1) = ( norm(u - v) < tol );

end
