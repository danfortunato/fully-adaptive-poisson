classdef ChebyshevOperators
% Construct Chebyshev operators to be used in annular solvers
% 
% Inputs:
%     M (int),     number of modes in chebyshev grid
%     rat (float), ratio giving width of annulus to [-1,1]

    properties ( Access = public )

        M
        V0
        V1
        V2
        VI0
        VI1
        VI2
        D00
        D01
        D12
        ibc_dirichlet
        obc_dirichlet
        ibc_neumann
        obc_neumann
        R01
        R12
        R02
        P10
        
    end

    methods

        function self = ChebyshevOperators(M, rat, chebkind)

            self.M = M;
            [x0, ~, v0, t0] = chebpts(M,   chebkind); x0 = flip(x0);
            [x1, ~, v1, t1] = chebpts(M-1, chebkind); x1 = flip(x1);
            [x2, ~, v2, t2] = chebpts(M-2, chebkind); x2 = flip(x2);

            % Vandermonde and inverse Vandermonde matrices
            self.V0 = chebvander(x0, M-1);
            self.V1 = chebvander(x1, M-2);
            self.V2 = chebvander(x2, M-3);
            self.VI0 = inv(self.V0);
            self.VI1 = inv(self.V1);
            self.VI2 = inv(self.V2);

            % Differentiation matrices
            %DC01 = chebder(eye(M))   / rat;
            %DC12 = chebder(eye(M-1)) / rat;
            %DC00 = [DC01; zeros(1,M)];
            %self.D00 = self.V0 * (DC00 * self.VI0);
            %self.D01 = self.V1 * (DC01 * self.VI0);
            %self.D12 = self.V2 * (DC12 * self.VI1);

            dom = rat * [-1 1];
            self.D00 = diffmat(M,         dom, ['chebkind' num2str(chebkind)]);
            self.D01 = diffmat([M-1 M],   dom, ['chebkind' num2str(chebkind)]);
            self.D12 = diffmat([M-2 M-1], dom, ['chebkind' num2str(chebkind)]);

            % Boundary condition operators
            %self.ibc_dirichlet = chebvander( 1, M-1) * self.VI0;
            %self.obc_dirichlet = chebvander(-1, M-1) * self.VI0;
            self.ibc_dirichlet = barymat(-1, x0, v0, pi, t0);
            self.obc_dirichlet = barymat( 1, x0, v0, 0,  t0);
            self.ibc_neumann = self.ibc_dirichlet * self.D00;
            self.obc_neumann = self.obc_dirichlet * self.D00;

            % Rank reduction operators
            %self.R01 = self.V1 * (speye(M-1, M)   * self.VI0);
            %self.R12 = self.V2 * (speye(M-2, M-1) * self.VI1);
            %self.R02 = self.R12 * self.R01;
            self.R01 = barymat(x1, x0, v0, t1, t0, 1);
            self.R12 = barymat(x2, x1, v1, t2, t1, 1);
            self.R02 = barymat(x2, x0, v0, t2, t0, 1);

            % Get poof operator from M-1 --> M
            %self.P10 = self.V0 * (speye(M, M-1) * self.VI1);
            self.P10 = barymat(x0, x1, v1, t0, t1, 1);

        end

    end

end

function out = chebvander(x, n)
    out = feval(chebpoly(0:n), x);
end

function out = chebder(c)
    [n, m] = size(c);
    out = zeros(n-1, m);                        % Initialize vector {c_r}
    w = repmat(2*(1:n-1)', 1, m);
    v = w.*c(2:end,:);                          % Temporal vector
    out(n-1:-2:1,:) = cumsum(v(n-1:-2:1,:), 1); % Compute c_{n-2}, c_{n-4}, ...
    out(n-2:-2:1,:) = cumsum(v(n-2:-2:1,:), 1); % Compute c_{n-3}, c_{n-5}, ...
    out(1,:) = .5*out(1,:);                     % Adjust the value for c_0
end

