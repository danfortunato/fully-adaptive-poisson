n = 16;
dom = Boundary.star(n, 'quadrature', 'panel');
[x, y] = smoothStrip(dom, n);

cx = mat2cell(x, repmat(n,dom.np,1), 1);
cy = mat2cell(y, repmat(n,dom.np,1), 1);

figure(1)
plot(dom)
hold on

t = chebpts(n, [0 1]);
xx = cell(dom.np,1);
yy = cell(dom.np,1);
for k = 1:dom.np
    xx{k} = (1-t).*cx{k}.' + t.*dom.x{k}(:,1).';
    yy{k} = (1-t).*cy{k}.' + t.*dom.x{k}(:,2).';
    scatter(xx{k}(:),yy{k}(:),'k')
end
shg