function samp = smoothcurve(v, h, k, n)
%SMOOTHCURVE   Smooth a curve.
%   SAMP = SMOOTHCURVE(V, H, K, M) returns a smoothed curve defined by the
%   sample points SAMP, given the vertices of the curve V and smoothing
%   parameters H and K. The number of sample points used for each
%   smoothed vertex is 2*N.
%
%   The method is based on the paper:
%
%     C. L. Epstein and M. O'Neil, "Smoothed corners and scattered waves",
%     SIAM J. Sci. Comput., 38 (2016), pp. A2665-A2698,
%     https://doi.org/10.1137/15M1028248.

    nv = length(v);

    % We now set up the local geometries for the edges
    for j = 1:nv
        j1 = 1 + mod(j, nv);
        j0 = 1 + mod(j-2, nv);
        X(j,1:2) = v(j1,1:2) - v(j,1:2);
        X(j,1:2) = X(j,1:2) / sqrt(X(j,1:2)*X(j,1:2)');
        Y(j,1:2) = v(j0,1:2) - v(j,1:2);
        Y(j,1:2) = Y(j,1:2) / sqrt(Y(j,1:2)*Y(j,1:2)');
    end

    [as0, j0, j1] = smoothabs(0.5, k, 2*n);
    s0 = (0:4*n-1)/(2*n)-1;
    x = zeros(length(j0:j1+1), nv);
    y = zeros(length(j0:j1+1), nv);
    for i = 1:nv
        as = 2 * h(i) * as0(j0:j1+1);
        s  = 2 * h(i) *  s0(j0:j1+1);
        x(:,i) = as + s;
        y(:,i) = as - s;
    end
    
    % The number of points on the smoothed edges.
    npt = j1 - j0 + 2;
    
    % Loop over vertices
    samp = [];
    for i = 1:nv
        start = length(samp) + 1;

        % We get samples of the function used to smooth the vertices
        % We use the part with indices from j0 to j1+1

        % The number of sample points on all the curved parts.
        npt1 = nv*npt; 
        nx2 = ceil(log2(npt1));
        nflp = nx2;

        samp(end+(1:npt),1:2) = ones(npt,1)*v(i,1:2) + x(1:npt,i)*X(i,1:2) + y(1:npt,i)*Y(i,1:2);

        % The last point on the current edge  
        Z1 = v(i,1:2) + x(npt,i)*X(i,1:2) + y(npt,i)*Y(i,1:2); 

        % Next panel index (in cyclic order)
        next = 1 + mod(i, nv); 

        % The first point on the next edge
        Z2 = v(next,1:2) + x(1,next)*X(next,1:2) + y(1,next)*Y(next,1:2);
            
        % Test code:
        Zdif = Z2-Z1;
        nZ = sqrt(Zdif*Zdif');

        % Construct the straight segment joining the two curved ones
%         psamp(end+(1:nflp),1:2) = ...
%            ones(nflp,1)*Z1(1:2)+((1:nflp)'/(1+nflp))*Zdif(1:2);
        %lin = ones(nflp,1)*Z1(1:2)+((1:nflp)'/(1+nflp))*Zdif(1:2);
        lin = [Z2; Z1];
        
        figure(1)
        plot(samp(start:end,1), samp(start:end,2), 'r-', 'LineWidth', 2)
        hold on
        plot(lin(:,1), lin(:,2), 'k-', 'LineWidth', 2)
    end
    samp(end+1,1:2) = samp(1,1:2);

end
  
function [samps, j0, j1] = smoothabs(h, k, n)
%SMOOTHABS   Smoothed absolute value.
%   [SAMP, J0, J1] = SMOOTHABS(H, K, N) uses the FFT to return 2*N samples
%   SAMP of the convolution between |x| and (1-(x/h)^2)^k*chi_{[-1,1]}(x).
%   Note that the minimum value of the smoothed function occurs at index
%   N+1. We also return the indices of the endpoints of the curved section,
%   J0 and J1, which coincide with the intersections of the smoothed curve
%   with the faces.

% First compute the Fourier coefficients of |x| exactly
ah(1)       = pi^2;                                 % Zero frequency
ah(2:n)     = -2*(1-(-1).^(1:(n-1)))./(1:(n-1)).^2; % Positive frequency
ah(n+1:2*n) = -2*(1-(-1).^(n:-1:1))./(n:-1:1).^2;   % Negative frequency

% Compute samples of the (unscaled) smoothing kernel
j0 = floor(n*(1-h)+1);
j1 = floor(n*(1+h)+1);

% jj = j0:(j1-1);
jj = chebpts(n, [j0 j1-1]);

% Take account of the characteristic function of the interval [-h, h]
fk(1:j0) = 0;
fk(j0+1:j1) = (1-(jj./(n*h)-1/h).^2).^k;
fk(j1+1:2*n) = 0;
ftfk = fft(fk);

% Multiply the FT of |x| by the FT of the filter
filftfk = ah.*ftfk; 

% Invert to get the spatial samples
samps = n * real(ifft(filftfk)) / (pi^2 * ftfk(1));

end
