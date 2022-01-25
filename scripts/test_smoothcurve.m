%% Make a panelized curve

p = 16; % Points per panel
dom = Boundary.star(p, 'quadrature', 'panel');
%dom = Boundary.multiscale(p, 'quadrature', 'panel');
%dom = newdom;

% These vertices define a piecewise linear polygon that approximates the
% original smooth curve.
v = [real(dom.zbreaks).' imag(dom.zbreaks).'];
v = v(1:end-1,:);

width = 0.4;

%% Define parameters for corner rounding

z = perturbBreaks(dom, width);
v = [real(z) imag(z)];
v = v(1:end-1,:);

nsmooth = 10; % Points per smoothed corner
k = 2;        % Exponent in smoothing kernel

% Smoothing widths (one per panel junction)
hl = vecnorm(v - circshift(v,  1), 2, 2) / 4; % Left neighbor widths
hr = vecnorm(v - circshift(v, -1), 2, 2) / 4; % Right neighbor widths
h = (hl + hr) / 2;                            % Average
h = h / 1.5;                                  % Scale by a constant

%% Run the smoother

figure(1), clf, hold on
samp = smoothcurve(v, h, k, nsmooth);
% plot(samp(:,1), samp(:,2), '-', 'LineWidth', 2), hold on
%plot(dom)

axis equal
hold off
shg
