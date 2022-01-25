load('multiscale.mat', 'f', 'breaks')

% Plot the true smooth curve, f:
tt = linspace(0, 1, 100).';
tt = (breaks(2:end)-breaks(1:end-1)).*tt + breaks(1:end-1);
plot(f(tt), 'k-', 'LineWidth', 1), hold on

% Plot the panel break points:
xy = f(breaks).';
xy = [real(xy) imag(xy)];
plot(xy(:,1), xy(:,2), 'k.', 'MarkerSize', 15)

% Run the 2D surface smoother on the panel break points to get a smooth
% curve:
n = 8;
[x_sm, y_sm] = find_smooth_curve_2D(xy, n, n, []);

plot(x_sm, y_sm, 'r-', 'LineWidth', 1)
axis equal
shg
