function adaptivePerturb(S)
%ADAPTIVEPERTURB   Construct a Boundary by perturbing another.

plot(S)
hold on

np = S.np;
len = arcLength(S, 1:np);

% Wrap around
len = [len; len(1)];
avglen = (len(1:np)+len(2:np+1))/2;
nx = [S.normal; S.normal{1}];

r = 0.4;
x = [real(S.zbreaks).' imag(S.zbreaks).'];

for k = 1:np
    % Get the normal vector at the interface between panels k and k+1
    nxL = util.feval(nx{k},    1); % The left panel evaluated on the right side
    nxR = util.feval(nx{k+1}, -1); % The right panel evaluated on the left side
    % Take the mean and normalize
    nxk = (nxL+nxR)/2;
    nxk = nxk ./ norm(nxk);
    x(k+1,:) = x(k+1,:) - r * avglen(k) * nxk;
end
x(1,:) = x(np+1,:);

xx = trigpts(np, [0 2*pi]);
yy = x(1:np,1)+x(1:np,2)*1i;
f = chebfun(x(1:np,1)+x(1:np,2)*1i, [0 2*pi], 'trig', 2*np);
%pp = pchip(xx,yy);
%f = chebfun(avglen, [0 2*pi], 'trig');
plot(f, 'r', 'LineWidth', 1)

%[psa,nt] = smthpoly(x(1:np,:), 0.1, 10, 2048);
%plot(psa(1:nt,1), psa(1:nt,2), '-', 'LineWidth', 2)

% xx = [xx; 2*pi];
% yy = [yy; yy(1,:)];
% cs = spaps(xx, yy,  0.05);
% xxx = linspace(0, 2*pi, 1000);
% vvv = fnval(cs, xxx);
% plot(vvv)

plot(x(:,1)+x(:,2)*1i, 'bo', 'LineWidth', 1)
hold off

end
