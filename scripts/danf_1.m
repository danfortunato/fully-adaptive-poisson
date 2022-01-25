% panel-aware strip width blending function. Barnett + Fortunato   10/28/20
%phi = @(t) exp(1./(1-1./t.^2));
p = 10; phi = @(t) (1-t.^2).^p;   % p even, order of the zero at +-1
phi = @(t) phi(t).*(t.^2<1);   % truncate it to [-1,1]
t = linspace(-1.5,1.5,1e3);
figure(1); clf; plot(t,phi(t),'-'); title('\phi')

sj = [2.^(-20:0), 2:6];  % some breakpts in params: dyadic then uniform
% just make local width wj for all sj except endpoints...
wj = 0*sj; wj(2:end-1) = (sj(3:end)-sj(1:end-2))/2;   % local widths - maybe not needed?
ss = logspace(log10(min(sj)/2), log10(max(sj)),1e4);  % param plot vals
relwid = 1.5;   % not too large otherwise can touch unlimited panels

% log xform version - failed due to rel err getting exp'ed...
%ll=0*ss; for j=2:numel(sj)-1, ll = ll + log(fj(j))*phi(fac*(ss-sj(j))/fj(j)); end, ff = exp(ll);

% plain version
%ff =0*ss; for j=2:numel(sj)-1, ff = ff + wj(j)*phi((ss-sj(j))/wj(j)/relwid); end

ff = spline(sj,wj,ss);  % cubic spline, decent but only C^1

figure(2); clf;
subplot(2,1,1); plot(ss,ff,'-'); hold on;
plot([sj;sj],[0*wj;wj],'k.-','markersize',20);   % wj stick plots
subplot(2,1,2); loglog(ss,ff,'-'); hold on; plot(sj,wj,'k.','markersize',20);
set(gca,'ylim',[1e-7 1]);

