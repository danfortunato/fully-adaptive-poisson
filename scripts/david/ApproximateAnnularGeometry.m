classdef ApproximateAnnularGeometry
% Approximate Annular Geometry for solving PDE in annular regions
%     n: number of discrete points in tangential direction
%     M: number of chebyshev modes in radial direction
%     width: width of radial region
%     approx_r: approximate radius of annulus

    properties ( Access = public )

        n
        M
        radius
        width
        radial_h
        tangent_h
        ns
        n2
        k
        ks
        iks
        rv0
        rv1
        rv2
        ratio
        approx_psi0
        approx_psi1
        approx_psi2
        approx_inv_psi0
        approx_inv_psi1
        approx_inv_psi2
        CO

    end

    methods

        function self = ApproximateAnnularGeometry(n, M, width, approx_r, chebkind)

            self.n = n;
            self.M = M;
            self.radius = approx_r;
            self.width = width;
            self.radial_h = width/M;
            self.tangent_h = 2*pi/n;
            self.ns = n-1;
            self.n2 = floor(n/2);
            self.k = fftfreq(n, 1/n);
            self.ks = [ravel(self.k(1:self.n2)); ravel(self.k(self.n2+2:end))];
            self.iks = 1i*self.ks;

            % r grids
            self.rv0 = chebpts(M,   [-width 0], chebkind);
            self.rv1 = chebpts(M-1, [-width 0], chebkind);
            self.rv2 = chebpts(M-2, [-width 0], chebkind);
            self.ratio = -width/2;

            % Coordinate transformations
            self.approx_psi0 = self.radius + self.rv0;
            self.approx_psi1 = self.radius + self.rv1;
            self.approx_psi2 = self.radius + self.rv2;
            self.approx_inv_psi0 = 1 ./ self.approx_psi0;
            self.approx_inv_psi1 = 1 ./ self.approx_psi1;
            self.approx_inv_psi2 = 1 ./ self.approx_psi2;

            % Chebyshev Operators
            self.CO = ChebyshevOperators(M, self.ratio, chebkind);

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