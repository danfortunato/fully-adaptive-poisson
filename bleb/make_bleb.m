% Get 2D curve from blebbing image.

im_orig  = imread('bleb_orig.png');
im_trace = imread('bleb_trace.png');
BW = logical(im_trace(:,:,1));
B = bwboundaries(rot90(BW, -1), 'noholes');

% Spectral approximation
x = B{1}(:,end:-1:1);
n = size(x, 1);
xhat = ifft(x);

% Apply a filter to roll off high frequencies
kfilter = cos((0:n-1)'/n*pi).^2;
xhat = xhat .* kfilter;

cfs = fftshift(xhat);
xfun = real(chebfun(cfs(:,1), 'coeffs', 'trig'));
yfun = real(chebfun(cfs(:,2), 'coeffs', 'trig'));
zfun = xfun + 1i*yfun;

% Convolve with Gaussian to smooth
s = 1/200;
g = chebfun(@(x) exp(-0.5*x.^2/s.^2)/(s*sqrt(2*pi)), 'trig');
xfun_conv = circconv(xfun, g);
yfun_conv = circconv(yfun, g);
zfun_conv = xfun_conv + 1i*yfun_conv;

plot(zfun_conv)
axis equal
axis([100 2300 200 2200])
xticks([])
yticks([])
shg

save('bleb_z.mat', 'zfun_conv')
