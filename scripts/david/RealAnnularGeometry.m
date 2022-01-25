classdef RealAnnularGeometry
% Approximate Annular Geometry for solving PDE in annular regions
%     n: number of discrete points in tangential direction
%     M: number of chebyshev modes in radial direction
%     width: width of radial region
%     approx_r: approximate radius of annulus

    properties ( Access = public )

        psi0
        psi1
        psi2
        inv_psi0
        inv_psi1
        inv_psi2
        DR_psi2
        ipsi_DR_ipsi_DT_psi2
        ipsi_DT_ipsi_DR_psi2

    end

    methods

        function self = RealAnnularGeometry(speed, curvature, M, width, chebkind)
            n = size(curvature,1);
            k = fftfreq(n, 1/n);
            dt_curvature = real(ifft(fft(curvature.').*1i.*k));
            rv0 = chebpts(M,   [-width 0], chebkind);
            rv1 = chebpts(M-1, [-width 0], chebkind);
            rv2 = chebpts(M-2, [-width 0], chebkind);
            self.psi0 = speed.' .* (1 + rv0.*curvature.');
            self.psi1 = speed.' .* (1 + rv1.*curvature.');
            self.psi2 = speed.' .* (1 + rv2.*curvature.');
            self.inv_psi0 = 1 ./ self.psi0;
            self.inv_psi1 = 1 ./ self.psi1;
            self.inv_psi2 = 1 ./ self.psi2;
            self.DR_psi2 = speed.' .* curvature.' .* ones(size(rv2));
            denom2 = speed.' .* (1 + rv2.*curvature.').^3;
            idenom2 = 1 ./ denom2;
            % these are what i think it should be? need to check computation
            %self.ipsi_DR_ipsi_DT_psi2 = (curvature.' - dt_curvature) .* idenom2;
            %self.ipsi_DT_ipsi_DR_psi2 = -dt_curvature .* idenom2;
            % these are what work...
            self.ipsi_DR_ipsi_DT_psi2 = dt_curvature .* idenom2;
            self.ipsi_DT_ipsi_DR_psi2 = dt_curvature .* idenom2;

        end

    end

end

function out = fftfreq(n, d)

if ( nargin == 1 )
    d = 1;
end

if ( mod(n,1) ~= 0)
    error('n should be an integer.');
end

val = 1/(n*d);
results = zeros(1,n);
k = floor((n-1)/2)+1;
results(1:k) = 0:k-1;
results(k+1:n) = -floor(n/2):-1;
out = results * val;

end