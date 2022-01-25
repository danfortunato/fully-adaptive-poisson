% Test divide by 1

N = 30;
S = Boundary.circle(N);
f = @(z) exp(z);
z = S.x{1}(:,1) + S.x{1}(:,2)*1i;
z0 = z(1)-0.1;

dz = S.dz(S.s{1});
dz = dz(:,1) + dz(:,2)*1i;
w = 2*pi/N*dz;
%w = S.w;
nx = S.normal{1}(:,1) + S.normal{1}(:,2)*1i;

% Find the density for this test case
%sigma = ones(1,N);
D = kernels.laplace.dlp(S);
I = eye(N);
rhs = real(f(z));
sigma = (D - I/2) \ rhs;
%sigma = sigma.';

% Evaluate the solution
sol = f(z0);
%u = -1/(2*pi) * sum( w.*real(conj(nx).*(z-z0))./abs(z-z0).^2 .* sigma );
u = -1/(2*pi*1i) * sum(w.*sigma./(z-z0));
cauchy = abs( u - sol )
% u = -1 * sum( w.*real(conj(nx).*(z-z0)).*sigma./abs(z-z0).^2 ) / ...
%          sum( 1/(1i) * 2*pi/N*dz./(z-z0) );
% u = -1 * sum( w.*real(conj(nx).*(z-z0))./abs(z-z0).^2 .* sigma ) / ...
%          sum( w.*real(conj(nx).*(z-z0))./abs(z-z0).^2 );
u = -sum(w.*sigma./(z-z0)) ./ sum(w./(z-z0));
interp = abs( u - sol )

%% Cauchy's integral formula

N = 100;
S = Boundary.star(N);
f = @(z) exp(z);
z = S.x{1}(:,1) + S.x{1}(:,2)*1i;
z0 = z(1) - 0.1;

dz = S.dz(S.s{1});
dz = dz(:,1) + dz(:,2)*1i;
w = 2*pi/N*dz;

sol = f(z0);
cauchy = abs( 1/(2*pi*1i) * sum(w.*f(z)./(z-z0))     - sol )
trick  = abs( sum(w.*f(z)./(z-z0)) ./ sum(w./(z-z0)) - sol )

%%
N = 10;
S = Boundary.star(N);
f = @(z) exp(z-0.9);
z0 = 0.9;
z0 = -0.31 - 0.95i;
z = exp(2*pi*1i*(1:N)/N);
w = 2*pi/N;

sol = f(z0);
cauchy = abs( 1/(2*pi*1i)*sum( w.*1i.*z.*f(z)./(z-z0) ) - sol )
interp = abs( 1/(2*pi*1i)*sum( w.*1i*z.*f(z)./(z-z0)) ./ (1/(2*pi*1i)*sum(w.*1i*z./(z-z0))) - sol )
%cauchy = abs( mean(z.*f(z)./(z-z0)) - sol )
%interp = abs( mean(z.*f(z)./(z-z0))./mean(z./(z-z0)) - sol )
