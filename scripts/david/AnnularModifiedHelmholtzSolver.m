classdef AnnularModifiedHelmholtzSolver
% Spectrally accurate Modified Helmholtz solver on annular domain
%
%     Solves (k^2-L)u = f in the annulus described by the Annular Geometry AG
%     Subject to the Robin boundary condition:
%     ia*u(ri) + ib*u_r(ri) = ig (boundary condition at the inner radius)
%     oa*u(ro) + ob*u_r(ro) = og (boundary condition at the outer radius)
% 
%     On instantionation, a preconditioner is formed with ia, ib, ua, ub
%         defining the boundary conditions
%     These can be changed at solvetime, but preconditioning may not work so well

    properties ( Access = public )

        RAG
        AAG
        ia
        ib
        oa
        ob
        k
        M
        ns
        n
        NB
        small_shape
        shape
        iterations_last_call
        KLUS

    end

    methods ( Access = public )

        function self = AnnularModifiedHelmholtzSolver(RAG, AAG, k, ia, ib, oa, ob)

            self.RAG = RAG;
            self.AAG = AAG;
            self.k = k;
            
            if ( nargin == 3 )
                ia = 1; ib = 0;
                oa = 1; ob = 0;
            end

            self.ia = ia; self.ib = ib;
            self.oa = oa; self.ob = ob;
            self.M = AAG.M;
            self.ns = AAG.ns;
            self.n = AAG.n;
            self.NB = self.M * self.ns;
            self.small_shape = [self.M, self.ns];
            self.shape = [self.M, self.n];
            
            CO = AAG.CO;
            apsi1 =  AAG.approx_psi1;
            aipsi1 = AAG.approx_inv_psi1;
            aipsi2 = AAG.approx_inv_psi2;
            ks = AAG.ks;
            D01 = CO.D01;
            D12 = CO.D12;
            R01 = CO.R01;
            R12 = CO.R12;
            R02 = CO.R02;
            ibcd = CO.ibc_dirichlet;
            ibcn = CO.ibc_neumann;
            obcd = CO.obc_dirichlet;
            obcn = CO.obc_neumann;
            ns = self.ns;
            M = self.M;

            self.KLUS = cell(ns,1);
            for i = 1:ns
                K = zeros(M);
                LL = aipsi2.*(D12*(apsi1.*D01)) - ks(i)^2.*(R12*(aipsi1.*R01));
                K(1:M-2,:) = self.k^2*R02 - LL;
                K(M-1,:) = self.ia*ibcd + self.ib*ibcn;
                K(M,:)   = self.oa*obcd + self.ob*obcn;
                self.KLUS{i} = decomposition(K, 'lu');
            end

        end

        function out = solve(self, f, ig, og, tol)
            R02 = self.AAG.CO.R02;
            ff = [R02*f; ig.'; og.']; ffh = mfft(ff); ffh = ffh(:);
            
            % Iterative solve:
            %[x, ~, ~, ~, resvec] = right_gmres(@(x) self.apply(x), ffh, [], tol, 100, @(x) self.preconditioner(x));
            %[x, ~, ~, ~, resvec] = gmres(@(x) self.apply(x), ffh, [], tol, 100, @(x) self.preconditioner(x));
            %fprintf('GMRES took %i iterations.\n', numel(resvec))
            
            % Direct solve:
            A = self.make_matrix();
            %ffh = self.preconditioner(ffh);
            x = A \ ffh;

            out = real( mifft(reshape(x, self.small_shape)) );
        end

    end
    
    methods ( Access = public )

        function y = preconditioner(self, fh)
            fh = reshape(fh, self.small_shape);
            fo = zeros(self.small_shape);
            for i = 1:self.ns
                fo(:,i) = self.KLUS{i} \ fh(:,i);
            end
            y = fo(:);
        end

        function y = apply(self, uh)
            CO = self.AAG.CO;
            ibcd = CO.ibc_dirichlet; ibcn = CO.ibc_neumann;
            obcd = CO.obc_dirichlet; obcn = CO.obc_neumann;
            uh = reshape(uh, self.small_shape);
            luh = scalar_laplacian(CO, self.AAG, self.RAG, uh);
            fuh = self.k^2*CO.R02*uh - luh;
            ibc = (self.ia*ibcd + self.ib*ibcn)*uh;
            obc = (self.oa*obcd + self.ob*obcn)*uh;
            y = [fuh; ibc; obc];
            y = y(:);
        end
        
        function A = make_matrix(self)
            N = self.M * self.ns;
            A = zeros(N);
            x = zeros(N,1);
            for i = 1:N
                x(i) = 1;
                A(:,i) = self.apply(x);
                %A(:,i) = self.preconditioner(A(:,i));
                x(i) = 0;
            end
        end

    end

end

function out = mfft(f)
    [m,n] = size(f);
    n2 = floor(n/2);
    fh = fft(f, [], 2);
    out = zeros(m,n-1);
    out(:,1:n2) = fh(:,1:n2);
    out(:,n2+1:end) = fh(:,n2+2:end);
end

function out = mifft(fh)
    [m,ns] = size(fh);
    n = ns+1;
    n2 = floor(n/2);
    out = zeros(m,n);
    out(:,1:n2) = fh(:,1:n2);
    out(:,n2+2:end) = fh(:,n2+1:end);
    out = ifft(out, [], 2);
end

function out = fourier_multiply(fh, m)
    out = mfft(m.*mifft(fh));
end

function luh = scalar_laplacian(CO, AAG, RAG, uh)
    uh_t  = CO.R01 * (uh .* AAG.iks.');
    uh_tt = CO.R12 * (fourier_multiply(uh_t, RAG.inv_psi1) .* AAG.iks.');
    uh_rr = CO.D12 * (fourier_multiply(CO.D01*uh, RAG.psi1));
    luh = fourier_multiply(uh_rr+uh_tt, RAG.inv_psi2);
end