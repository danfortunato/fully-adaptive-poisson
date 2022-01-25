function [xcheb, fcheb] = mollify(a, f, beta, bleed, n)

if ( nargin == 0 ), test(); return; end
if ( nargin < 3 ), beta = 4;  end
if ( nargin < 4 ), bleed = 2; end
if ( nargin < 5 ), n = 16;    end

npanel = length(a)-1;
scl = a(end)-a(1);
a = a(:).';
f = f(:).';
padL = npanel-bleed+1:npanel;
padR = 2:bleed+3;
a = [a(padL)-scl  a a(padR)+scl];
f = [f(padL)      f f(padR)];

nextra = 2*bleed+1;
fround = cell(npanel+nextra, 1);
lin    = cell(npanel+nextra, 1);
s = diff(f) ./ diff(a); % Slopes
for k = 1:npanel+nextra
    j = k+1;
    %scale = mean([a(j)-a(j-1) a(j+1)-a(j)]);
    scale = f(j);
    b = beta/scale;
    srelu = @(x) (exp(-(b*x).^2)/sqrt(pi) + b*x.*(1+erf(b*x))) / (2*b);
    fround{k} = @(x) f(j) + s(j-1)*(x-a(j)) + (s(j)-s(j-1))*srelu(x-a(j));
    lin{k}    = @(x) f(j) - s(j)*(x-a(j));
end

xcheb = cell(npanel, 1);
fcheb = cell(npanel, 1);
for k = 1:npanel
    % This panel has endpoints k and k+1
    kl = k;
    kr = k+2*bleed+1;
    fsum = sprintf('@(x) fround{%d}(x)', kr);
    for j = kl:kr-1
        fsum = sprintf('%s + fround{%d}(x) + lin{%d}(x) - 2*f(%d)', fsum, j, j, j+1);
    end
    fsum = eval(fsum);
    xcheb{k} = chebpts(n, [a(k+bleed) a(k+bleed+1)]);
    fcheb{k} = fsum(xcheb{k});
end

end

function test()

%a = [1 2 3 3.5 4.5 6.5 7.5];
sj = [2.^(-20:0), 2:8 10 11:14];  % some breakpts in params: dyadic then uniform
% just make width-func f for all sj except endpoints
wj = 0*sj; wj(2:end-1) = (sj(3:end)-sj(1:end-2))/2;   % local widths - maybe not needed?
ss = logspace(log10(min(sj)/2), log10(max(sj)),1e4);  % param grid
a = sj;

w = diff(a);
f = [(w(1)+w(end))/2 (w(1:end-1)+w(2:end))/2 (w(1)+w(end))/2];

n = 16;
%beta = 1;
%bleed = 2;
beta = 4;
bleed = 2;
[x, v] = mollify(a, f, beta, bleed, n);

x = x(bleed+2:end-bleed-1);
v = v(bleed+2:end-bleed-1);
worst = 0;
for k = 1:length(v)-1
    worst = max(worst, abs(v{k}(end) - v{k+1}(1)));
    if ( abs(v{k}(end) - v{k+1}(1)) > 1e-10 )
        k
        abs(v{k}(end) - v{k+1}(1))
        warning('Asymptotics not matched.');
    end
end
worst
x = cell2mat(x);
v = cell2mat(v);

plot(a, f, '-o'), hold on
plot(x, v, '-o')
shg

figure(2); clf;
subplot(2,1,1); plot(x,v,'-'); hold on;
plot([sj;sj],[0*wj;wj],'k.-','markersize',20);   % wj stick plots
subplot(2,1,2); loglog(x,v,'-'); hold on; plot(sj,wj,'k.','markersize',20);
set(gca,'ylim',[1e-7 1]);

end
