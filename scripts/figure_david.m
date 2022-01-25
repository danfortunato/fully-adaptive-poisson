export = true;

figure(1)
plot(Gamma)
axis equal tight
ax = gca;
set(ax, 'fontsize', 18)
ax.XAxis.LineWidth = 1;
ax.YAxis.LineWidth = 1;
box on
shg

if ( export )
    export_fig -r400 gamma.png
end


figure(2)
plot(Gamma)
axis equal
scl = 2e-3;
axis([-scl scl 1-scl 1+scl])
ax = gca;
set(ax, 'fontsize', 18)
ax.XAxis.LineWidth = 1;
ax.YAxis.LineWidth = 1;
box on
shg

if ( export )
    export_fig -r400 gamma_zoom.png
end

figure(3)
plot(Gamma)
hold on
x = Gamma1.x;
for k = 1:Gamma1.np
    plot(x{k}(:,1), x{k}(:,2), 'k--', 'LineWidth', 1);
end
plot(Gamma1.zbreaks, 'k.', 'MarkerSize', 15)
axis equal tight
ax = gca;
set(ax, 'fontsize', 18)
ax.XAxis.LineWidth = 1;
ax.YAxis.LineWidth = 1;
box on
shg

if ( export )
    export_fig -r400 gamma1.png
end

figure(4)
plot(Gamma)
hold on
x = Gamma1.x;
for k = 1:Gamma1.np
    plot(x{k}(:,1), x{k}(:,2), 'k--', 'LineWidth', 1);
end
plot(Gamma1.zbreaks, 'k.', 'MarkerSize', 15)
axis equal
scl = 2e-3;
axis([-scl scl 1-scl 1+scl])
ax = gca;
set(ax, 'fontsize', 18)
ax.XAxis.LineWidth = 1;
ax.YAxis.LineWidth = 1;
box on
shg

if ( export )
    export_fig -r400 gamma1_zoom.png
end

figure(5)
plot(Gamma)
hold on
x = Gamma1.x;
for k = 1:Gamma1.np
    plot(x{k}(:,1), x{k}(:,2), 'k--', 'LineWidth', 1);
end
plot(Gamma1.zbreaks, 'k.', 'MarkerSize', 15)
axis equal
scl = 4e-4;
axis([-scl scl 1-scl 1+scl])
ax = gca;
set(ax, 'fontsize', 18)
ax.XAxis.LineWidth = 1;
ax.YAxis.LineWidth = 1;
box on
shg

if ( export )
    export_fig -r400 gamma1_zoom2.png
end

alignfigs

%%
figure(6)
plot_f = tf;
nplotpts = 400;
[x, y] = meshgrid(linspace(plot_f.boxes(1).domain(1), plot_f.boxes(1).domain(2), nplotpts), ...
                  linspace(plot_f.boxes(1).domain(3), plot_f.boxes(1).domain(4), nplotpts));
plot_u = feval(plot_f, x, y);
surface(x, y, 0*plot_u, plot_u, 'EdgeAlpha', 1);
shading interp
view(2)
ax = gca;
axis equal tight
set(ax, 'fontsize', 18)
ax.XAxis.LineWidth = 1;
ax.YAxis.LineWidth = 1;
box on

hold on
plot(Gamma)
hold on
x = Gamma1.x;
for k = 1:Gamma1.np
    plot(x{k}(:,1), x{k}(:,2), 'k--', 'LineWidth', 1);
    hold on
end
hold on
plot(Gamma1.zbreaks, 'k.', 'MarkerSize', 15)
colorbar
axis equal tight


if ( export )
    export_fig -r400 rhs_nogrid.png
end

%%
figure(7)
plot_f = tf;
nplotpts = 400;
[x, y] = meshgrid(linspace(plot_f.boxes(1).domain(1), plot_f.boxes(1).domain(2), nplotpts), ...
                  linspace(plot_f.boxes(1).domain(3), plot_f.boxes(1).domain(4), nplotpts));
plot_u = feval(plot_f, x, y);
surface(x, y, 0*plot_u, plot_u, 'EdgeAlpha', 1);
hold on
plot(u_bulk)
shading interp
view(2)
ax = gca;
axis equal tight
set(ax, 'fontsize', 18)
ax.XAxis.LineWidth = 1;
ax.YAxis.LineWidth = 1;
box on

hold on
plot(Gamma)
hold on
x = Gamma1.x;
for k = 1:Gamma1.np
    plot(x{k}(:,1), x{k}(:,2), 'k--', 'LineWidth', 1);
    hold on
end
hold on
plot(Gamma1.zbreaks, 'k.', 'MarkerSize', 15)
colorbar
axis equal tight

if ( export )
    export_fig -r400 rhs_grid.png
end

%%
figure(8)
plot_f = tf;
nplotpts = 400;
scl = 2e-3;
[x, y] = meshgrid(linspace(-scl, scl, nplotpts), ...
                  linspace(1-scl, 1+scl, nplotpts));
plot_u = feval(plot_f, x, y);
surface(x, y, 0*plot_u, plot_u, 'EdgeAlpha', 1);
axis equal
shading interp
view(2)
ax = gca;
set(ax, 'fontsize', 18)
ax.XAxis.LineWidth = 1;
ax.YAxis.LineWidth = 1;
box on

hold on
plot(Gamma)
hold on
x = Gamma1.x;
vals = leafvals(tf);
for k = 1:Gamma1.np
    plot(x{k}(:,1), x{k}(:,2), 'k--', 'LineWidth', 1);
    hold on
end
hold on
plot(Gamma1.zbreaks, 'k.', 'MarkerSize', 15)
colorbar
axis([-scl scl 1-scl 1+scl])
colorbar

if ( export )
    export_fig -r400 rhs_nogrid_zoom.png
end

%%
figure(9)
plot_f = tf;
nplotpts = 400;
scl = 2e-3;
[x, y] = meshgrid(linspace(-scl, scl, nplotpts), ...
                  linspace(1-scl, 1+scl, nplotpts));
plot_u = feval(plot_f, x, y);
surface(x, y, 0*plot_u, plot_u, 'EdgeAlpha', 1);
hold on
plot(u_bulk)
axis equal
shading interp
view(2)
ax = gca;
set(ax, 'fontsize', 18)
ax.XAxis.LineWidth = 1;
ax.YAxis.LineWidth = 1;
box on

hold on
plot(Gamma)
hold on
x = Gamma1.x;
for k = 1:Gamma1.np
    plot(x{k}(:,1), x{k}(:,2), 'k--', 'LineWidth', 1);
    hold on
end
hold on
plot(Gamma1.zbreaks, 'k.', 'MarkerSize', 15)
colorbar
axis([-scl scl 1-scl 1+scl])
colorbar

if ( export )
    export_fig -r400 rhs_grid_zoom.png
end

%%
figure(10)
plot_f = tf;
boxes = leaves(plot_f);
nplotpts = 400;
scl = 4e-4;
[x, y] = meshgrid(linspace(-scl, scl, nplotpts), ...
                  linspace(1-scl, 1+scl, nplotpts));
plot_u = feval(plot_f, x, y);
h = surface(x, y, 0*plot_u, plot_u, 'EdgeAlpha', 1);
hold on
plot(u_bulk)
axis equal
shading interp
view(2)
ax = gca;
set(ax, 'fontsize', 18)
ax.XAxis.LineWidth = 1;
ax.YAxis.LineWidth = 1;
box on

hold on
plot(Gamma)
hold on
x = Gamma1.x;
for k = 1:Gamma1.np
    plot(x{k}(:,1), x{k}(:,2), 'k--', 'LineWidth', 1);
    hold on
end
hold on
plot(Gamma1.zbreaks, 'k.', 'MarkerSize', 15)
axis([-scl scl 1-scl 1+scl])
colorbar

if ( export )
    export_fig -r400 rhs_grid_zoom2.png
end

%%
figure(11)
plot_f = tf;
boxes = leaves(plot_f);
nplotpts = 400;
scl = 4e-4;
[x, y] = meshgrid(linspace(-scl, scl, nplotpts), ...
                  linspace(1-scl, 1+scl, nplotpts));
plot_u = feval(plot_f, x, y);
h = surface(x, y, 0*plot_u, plot_u, 'EdgeAlpha', 1);
hold on
axis equal
shading interp
view(2)
ax = gca;
set(ax, 'fontsize', 18)
ax.XAxis.LineWidth = 1;
ax.YAxis.LineWidth = 1;
box on

hold on
plot(Gamma)
hold on
x = Gamma1.x;
for k = 1:Gamma1.np
    plot(x{k}(:,1), x{k}(:,2), 'k--', 'LineWidth', 1);
    hold on
end
hold on
plot(Gamma1.zbreaks, 'k.', 'MarkerSize', 15)
axis([-scl scl 1-scl 1+scl])
colorbar

% if ( export )
%     export_fig -r400 rhs_nogrid_zoom2.png
% end

%%
figure(12)
K = length(dom);
for k = 1:K
    surf(dom(k).x, dom(k).y, 0*dom(k).x, 'FaceAlpha', 0), hold on
end

hold on
plot(Gamma)
hold on
plot(Gamma1)

axis equal tight

view(2)
ax = gca;
set(ax, 'fontsize', 18)
ax.XAxis.LineWidth = 1;
ax.YAxis.LineWidth = 1;
box on

if ( export )
    export_fig -r400 sem.png
end

%%
figure(13)
K = length(dom);
for k = 1:K
    surf(dom(k).x, dom(k).y, 0*dom(k).x, 'FaceAlpha', 0), hold on
end

hold on
plot(Gamma)
hold on
plot(Gamma1)

axis equal tight

view(2)
ax = gca;
set(ax, 'fontsize', 18)
ax.XAxis.LineWidth = 1;
ax.YAxis.LineWidth = 1;
box on

scl = 5e-5;
axis([-scl scl 1.00028-scl 1.00028+scl])
if ( export )
    export_fig -r400 sem_zoom.png
end

%%

alignfigs