classdef LaplaceSolver
%LAPLACESOLVER

    properties
        op
        rep
        side
        bc
        method
        n
    end

    methods

        function L = LaplaceSolver(Gamma, opts)

            arguments
                Gamma Boundary
                opts.side   (1,:) char {mustBeMember(opts.side, {'interior', 'exterior'})} = 'interior'
                opts.bc     (1,:) char {mustBeMember(opts.bc, {'dirichlet', 'neumann'})}   = 'dirichlet'
                opts.method (1,:) char {mustBeMember(opts.method, {'direct', 'fmm'})}      = 'fmm'
            end

            L.side   = opts.side;
            L.bc     = opts.bc;
            L.method = opts.method;

            if ( strcmpi(L.side, 'interior') && strcmpi(L.bc, 'dirichlet') )
                c = -1/2;
                modified = false;
                op  = 'dlp';
                rep = 'dlp';
            elseif ( strcmpi(L.side, 'exterior') && strcmpi(L.bc, 'dirichlet') )
                c = 1/2;
                modified = true;
                op  = 'dlp';
                rep = 'dlp';
            elseif ( strcmpi(L.side, 'interior') && strcmpi(L.bc, 'neumann') )
                c = 1/2;
                modified = true;
                op  = 'splp';
                rep = 'slp';
            elseif ( strcmpi(L.side, 'exterior') && strcmpi(L.bc, 'neumann') )
                c = -1/2;
                modified = true;
                op  = 'splp';
                rep = 'slp';
            end

            side = L.side(1);
            L.n = numel(Gamma);

            % Generate near corrections as a sparse matrix
            K_corr = kernels.laplace.([op '_correction'])(Gamma, side);

            if ( strcmpi(L.method, 'direct') )
                I = eye(L.n);
                K_smooth = kernels.laplace.(op)(Gamma, modified=modified);
                L.op = c*I + K_corr + K_smooth;
            else
                K_smooth = @(sigma) kernels.laplace.([op '_apply'])(Gamma, sigma, modified=modified);
                L.op = @(sigma) c*sigma + K_corr*sigma + K_smooth(sigma);
            end

            L.rep = @(sigma) @(x,y) kernels.laplace.(rep)(Gamma, ...
                side=side, ...
                density=sigma, ...
                modified=modified, ...
                target=[x(:) y(:)], ...
                closeeval=true ...
            );

        end

    end

end
