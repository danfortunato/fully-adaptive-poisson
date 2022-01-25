beta = 2;
x = linspace(-1, 1, 101);

relu = @(x) max(0,x);
srelu = restrict(cumsum(chebfun(@(x) 1+erf(beta*x), [-100 1])), [-1 1]);
splus = chebfun(@(x) softplus(beta*x)/beta);

plot(x, relu(x)), hold on
plot(x, srelu(x))
plot(x, splus(x))
axis equal
shg

%%
close all

% w = [1 2 2 1];
% a = [0 cumsum(w)];
% a = rescale(a, -1, 1);
% w = diff(a);
% f = [(w(1)+w(end))/2 (w(1:end-1)+w(2:end))/2 (w(1)+w(end))/2];

a = 2.^(0:-1:-5); a = [a; a]; a = a(:).';
a = [a flip(a)];
a = [0 cumsum(a)];
w = diff(a);
f = [(w(1)+w(end))/2 (w(1:end-1)+w(2:end))/2 (w(1)+w(end))/2];

% sj = [2.^(-4:0), 2:6];
% wj = 0*sj; wj(2:end-1) = (sj(3:end)-sj(1:end-2))/2;
% a = sj;
% f = wj;

% a = [-2 -1 0 1 1.5 2 3];
% f = [2 1 2 1 1.5 1];
% a = [-3 -1 0 1 1.5 2 3];
% f = [2 1.5 1 .75 .5 .75];
a = [1 2 3 3.5 4.5 6.5 7.5];
%f = [1 1 0.5 1 2 1];
%f = [f f(1)];
w = diff(a);
f = [(w(1)+w(end))/2 (w(1:end-1)+w(2:end))/2 (w(1)+w(end))/2];

plot(a, f), hold on

dom = [a(1) a(end)];
x = linspace(dom(1), dom(2), 1000);
nvert  = length(f)-1;
npanel = nvert-1;

beta = 1;
relu = @(x) max(0, x);
%srelu = restrict(cumsum(chebfun(@(x) (1+erf(beta*x))/2, [-100 10])), [-10 10]);
%srelu = chebfun(@(x) softplus(beta*x)/beta, [-10 10]); % Softmax

fun    = cell(nvert, 1);
fround = cell(nvert, 1);
lin    = cell(nvert, 1);
s = diff(f) ./ diff(a); % Slopes
srelu = cell(nvert, 1);
for k = 2:nvert
    %beta = 4*f(k-1);
    beta = 5;
    srelu{k} = restrict(cumsum(chebfun(@(x) (1+erf(beta*x))/2, [-100 10])), [-10 10]);
    fun{k}    = @(x) f(k) + s(k-1)*(x-a(k)) + (s(k)-s(k-1))*relu(x-a(k));
    fround{k} = @(x) f(k) + s(k-1)*(x-a(k)) + (s(k)-s(k-1))*srelu{k}(x-a(k));
    lin{k}    = @(x) f(k) - s(k)*(x-a(k));
end
% The first and last ones wrap around
scl = diff(dom);
fun{1}    = @(x) f(1) + s(end)*(x-a(1)) + (s(1)-s(end))*relu(x-a(1));
fround{1} = @(x) f(1) + s(end)*(x-a(1)) + (s(1)-s(end))*srelu{end}(x-a(1));
lin{1}    = @(x) f(1) - s(1)*(x-a(1));

fsum12 = @(x) fround{1}(x) + fround{2}(x)     + lin{1}(x) - 2*f(1);
fsum23 = @(x) fround{2}(x) + fround{3}(x)     + lin{2}(x) - 2*f(2);
fsum34 = @(x) fround{3}(x) + fround{4}(x)     + lin{3}(x) - 2*f(3);
fsum45 = @(x) fround{4}(x) + fround{5}(x)     + lin{4}(x) - 2*f(4);
fsum56 = @(x) fround{5}(x) + fround{6}(x)     + lin{5}(x) - 2*f(5);
fsum61 = @(x) fround{6}(x) + fround{1}(x-scl) + lin{6}(x) - 2*f(6);

fsum123 = @(x) fround{1}(x) + fround{2}(x) + fround{3}(x) + ...
               lin{1}(x) + lin{2}(x) - 2*(f(1)+f(2));
fsum234 = @(x) fround{2}(x) + fround{3}(x) + fround{4}(x) + ...
               lin{2}(x) + lin{3}(x) - 2*(f(2)+f(3));
fsum345 = @(x) fround{3}(x) + fround{4}(x) + fround{5}(x) + ...
               lin{3}(x) + lin{4}(x) - 2*(f(3)+f(4));
fsum456 = @(x) fround{4}(x) + fround{5}(x) + fround{6}(x) + ...
               lin{4}(x) + lin{5}(x) - 2*(f(4)+f(5));
fsum561 = @(x) fround{5}(x) + fround{6}(x) + fround{1}(x-scl) + ...
               lin{5}(x) + lin{6}(x) - 2*(f(5)+f(6));
fsum612 = @(x) fround{6}(x) + fround{1}(x-scl) + fround{2}(x-scl) + ...
               lin{6}(x) + lin{1}(x-scl) - 2*(f(6)+f(1));

fsum1234 = @(x) fround{1}(x) + fround{2}(x) + fround{3}(x) + fround{4}(x) + ...
                lin{1}(x) + lin{2}(x) + lin{3}(x) - 2*(f(1)+f(2)+f(3));
fsum2345 = @(x) fround{2}(x) + fround{3}(x) + fround{4}(x) + fround{5}(x) + ...
                lin{2}(x) + lin{3}(x) + lin{4}(x) - 2*(f(2)+f(3)+f(4));
fsum3456 = @(x) fround{3}(x) + fround{4}(x) + fround{5}(x) + fround{6}(x) + ...
                lin{3}(x) + lin{4}(x) + lin{5}(x) - 2*(f(3)+f(4)+f(5));

fsum12345 = @(x) fround{1}(x) + fround{2}(x) + fround{3}(x) + fround{4}(x) + fround{5}(x) + ...
                 lin{1}(x) + lin{2}(x) + lin{3}(x) + lin{4}(x) - 2*(f(1)+f(2)+f(3)+f(4));

fsum123456 = @(x) fround{1}(x) + fround{2}(x) + fround{3}(x) + fround{4}(x) + fround{5}(x) + fround{6}(x) + ...
                  lin{1}(x) + lin{2}(x) + lin{3}(x) + lin{4}(x) + lin{5}(x) - 2*(f(1)+f(2)+f(3)+f(4)+f(5));

fsum1234561 = @(x) fround{1}(x) + fround{2}(x) + fround{3}(x) + fround{4}(x) + fround{5}(x) + fround{6}(x) + fround{1}(x-scl) + ...
                   lin{1}(x) + lin{2}(x) + lin{3}(x) + lin{4}(x) + lin{5}(x) + lin{6}(x) - 2*(f(1)+f(2)+f(3)+f(4)+f(5)+f(6));

fsumall = @(x) fround{1}(x) + fround{2}(x) + fround{3}(x) + fround{4}(x) + fround{5}(x) + fround{6}(x) + ...
               fround{1}(x-scl) + fround{2}(x-scl) + fround{3}(x-scl) + fround{4}(x-scl) + fround{5}(x-scl) + fround{6}(x-scl) + ...
               lin{1}(x) + lin{2}(x) + lin{3}(x) + lin{4}(x) + lin{5}(x) + lin{6}(x) + ...
               lin{1}(x-scl) + lin{2}(x-scl) + lin{3}(x-scl) + lin{4}(x-scl) + lin{5}(x-scl) + ...
               - 2*(f(1)+f(2)+f(3)+f(4)+f(5)+f(6)) + ...
               - 2*(f(1)+f(2)+f(3)+f(4)+f(5));

% fsum = fsum456;
% plot(x, fsum(x), 'r'), plot(x, fsum(x+scl), 'r')
plot(x, fround{4}(x))
plot(x, fround{5}(x))
% plot(x, fsum12(x))
% plot(x, fsum23(x))
%plot(x, fsum23(x))
% plot(x, fsum45(x))
% plot(x, fsum56(x))
% plot(x, fsum61(x))
% plot(x, fsum123(x))
% plot(x, fsum234(x))
% plot(x, fsum345(x))
% plot(x, fsum456(x))
% plot(x, fsum561(x)), plot(x, fsum561(x+scl))
% plot(x, fsum612(x)), plot(x, fsum612(x+scl))
% plot(x, fsum1234(x))
%plot(x, fsum2345(x))
%plot(x, fsum3456(x))
%plot(x, fsum12345(x))
%plot(x, fsum123456(x))
% plot(x, fsum1234561(x))
%plot(x, fsumall(x))
axis equal tight
hold on
plot(a, f, 'ro')

%g = chebfun(fsum23, dom, 'trig');
%g = chebfun(fsumall, dom);

figure(2)
%err = abs(fround{3}(x) - fsum23(x));
err = @(x) abs(fsum3456(x) - fsum123456(x));
% err = abs(fsum123456(x) - fsumall(x));
subplot(211), plot(x, err(x))
subplot(212), semilogy(x, err(x)), xlim(dom)
hold on
plot(a, err(a), 'ro')

h = findobj([figure(1) figure(2)], 'Type', 'Line');
set(h, 'LineWidth', 2)
figure(1)
alignfigs
