function pass = test_periodic_poisson()

tol = 1e2 * eps;
pass = [];

n = 30;
u = chebfun2(@(x,y) -cos(x).*exp(sin(x)).*sin(y), [-pi pi -pi pi], 'trig');
f = lap(u);
v = periodicPoisson(f, n);
pass(end+1) = ( norm(u - v) < tol );

end
