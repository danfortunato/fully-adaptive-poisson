figure(1)
clf
set(gcf, 'Position', [1 485 1792 532]);

n = 16;
%dom = Boundary.wavepacket(n, 'quadrature', 'panel');
%dom = refine(dom);
%dom = Boundary.star(n, 'quadrature', 'panel');
dom = Boundary.circle(n, 'quadrature', 'panel');
dom = refine(dom, 1);

% Make a smooth strip
%[x, y] = smoothStrip(dom, n);
% zz = cell2mat(dom.z);
% normals = cell2mat(dom.normal);
% th = 0;
% R = [cos(th) -sin(th);
%      sin(th)  cos(th)];
% normals = (R*normals.').';
% cn = normals(:,1) + 1i*normals(:,2);
% zz1 = zz - 0.3*cn;
% x = real(zz1);
% y = imag(zz1);
beta = 4;
bleed = 10;
width = 0.75;
[x, y] = smoothStrip2(dom, n, beta, bleed, width);
xleg = chebvals2legvals(reshape(x, n, dom.np)); xleg = xleg(:);
yleg = chebvals2legvals(reshape(y, n, dom.np)); yleg = yleg(:);
z1 = mat2cell(xleg + 1i*yleg, repmat(n, dom.np, 1), 1);
Gamma1 = Boundary(z1);

% Build the strip grid
% cx = mat2cell(x, repmat(n,dom.np,1), 1);
% cy = mat2cell(y, repmat(n,dom.np,1), 1);
% t = chebpts(n, [0 1]);
% xx = cell(dom.np, 1);
% yy = cell(dom.np, 1);
% for k = 1:dom.np
%     dom.x{k}(:,1) = legvals2chebvals(real(dom.x{k}(:,1)));
%     dom.x{k}(:,2) = legvals2chebvals(real(dom.x{k}(:,2)));
%     %cx{k} = legvals2chebvals(cx{k});
%     %cy{k} = legvals2chebvals(cy{k});
%     xx{k} = (1-t).*cx{k}.' + t.*dom.x{k}(:,1).';
%     yy{k} = (1-t).*cy{k}.' + t.*dom.x{k}(:,2).';
%     %xx{k} = (t).*cx{k}.' + (1-t).*dom.x{k}(:,1).';
%     %yy{k} = (t).*cy{k}.' + (1-t).*dom.x{k}(:,2).';
% end
% dom = cell2struct([xx yy], {'x','y'}, 2);
% Build the strip grid
xsem = chebpts(n);
t = chebpts(n, [0 1]);
xx = cell(dom.np,1);
yy = cell(dom.np,1);
gamx_cheb  = cell(dom.np, 1);
gamy_cheb  = cell(dom.np, 1);
gam1x_cheb = cell(dom.np, 1);
gam1y_cheb = cell(dom.np, 1);
gam_nx_cheb = cell(dom.np, 1);
gam1_nx_cheb = cell(dom.np, 1);
LV2CV = legvals2chebvals(eye(n));
for k = 1:dom.np
    gamx_cheb{k}  = LV2CV * real(dom.x{k}(:,1));
    gamy_cheb{k}  = LV2CV * real(dom.x{k}(:,2));
    gam1x_cheb{k} = LV2CV * real(Gamma1.x{k}(:,1));
    gam1y_cheb{k} = LV2CV * real(Gamma1.x{k}(:,2));
    gam_nx_cheb{k} = LV2CV * dom.normal{k};
    gam1_nx_cheb{k} = LV2CV * Gamma1.normal{k};
    xx{k} = t.*gam1x_cheb{k}.' + (1-t).*gamx_cheb{k}.';
    yy{k} = t.*gam1y_cheb{k}.' + (1-t).*gamy_cheb{k}.';
end
dom = cell2struct([xx yy], {'x','y'}, 2);
K = length(dom);

fs = 16; % Font size

subplot(131)
for k = 1:K
    surf(dom(k).x, dom(k).y, 0*dom(k).x, 'FaceAlpha', 0), hold on
end
hold off, view(2), axis equal tight
title('Mesh')
ax = gca;
pos = ax.Position;
ax.FontSize = fs;
drawnow, shg

% Test against a random smooth function
rng(0)
rfun = randnfun2(1, [-2 2 -2 2]);
rhs = lap(rfun);
%bc = @(x,y) rfun(x,y);

tic
S = StripSolver(dom, rhs);
build(S)
bc = rfun(S.patches{1}.xy(:,1), S.patches{1}.xy(:,2));
u = S \ bc;
toc

subplot(132)
for k = 1:K
    surf(dom(k).x, dom(k).y, u{k}), hold on
end
hold off, shading interp, view(2), colorbar
axis equal tight
title('Solution')
ax = gca;
ax.Position(3:4) = pos(3:4);
ax.FontSize = fs;
drawnow, shg

subplot(133)
err = 0;
for k = 1:K
    rr = rfun(dom(k).x, dom(k).y);
    surf(dom(k).x, dom(k).y, log10(abs(rr - u{k}))), hold on
    err = max(err, max(max(abs(rr-u{k}))));
end
err
hold off, shading interp, view(2), colorbar
axis equal tight
title('Error')
ax = gca;
ax.Position(3:4) = pos(3:4);
ax.FontSize = fs;
drawnow, shg
