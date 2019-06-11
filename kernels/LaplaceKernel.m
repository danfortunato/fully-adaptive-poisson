classdef LaplaceKernel < Kernel
%LAPLACEKERNEL   Represent the Laplace kernel.

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% STATIC METHODS
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    methods ( Access = public, Static = true )

        function u = dlp(source, varargin)
        %DLP   Evaluate the Laplace double-layer potential.

            % Parse optional inputs
            p = inputParser;
            addRequired(p,  'source', @(s) isa(s,'ClosedCurve'));
            addParameter(p, 'target',   [],    @isnumeric);
            addParameter(p, 'density',  [],    @isnumeric);
            addParameter(p, 'modified', false, @islogical);
            parse(p, source, varargin{:});
            target   = p.Results.target;
            density  = p.Results.density;
            modified = p.Results.modified;

            dipole = false;
            if ( isempty(target) )
                target = cell2mat(source.x);
                dipole = true;
            end
            dipole = dipole || isequal(cell2mat(source.x), target);

            % Unpack things
            M     = size(target,1);
            N     = source.N * source.np;
            x     = cell2mat(source.x);
            w     = cell2mat(source.w);
            nx    = cell2mat(source.normal);
            kappa = cell2mat(source.curvature);

            % Fill in double-layer potential

            u = zeros(M,N);
            for j = 1:N
                r = target - x(j,:);
                u(:,j) = r./sum(r.^2,2)*nx(j,:).'/(2*pi)*w(j);
            end

            % Limiting value when source == target (dipole kernel)
            if ( dipole )
                u(1:N+1:end) = -kappa/(4*pi).*w;
            end

            % Correction for modified DLP
            if ( modified )
                u = u + w.';
            end

            % Include the density
            if ( ~isempty(density) )
                u = u * density;
            end

        end

        function u = slp(source, varargin)
        %SLP   Evaluate the Laplace single-layer potential.
        
            % TODO: Needs singular quadrature.

            % Parse optional inputs
            p = inputParser;
            addRequired(p,  'source', @(s) isa(s,'ClosedCurve'));
            addParameter(p, 'target',  [], @isnumeric);
            addParameter(p, 'density', [], @isnumeric);
            parse(p, source, varargin{:});
            target  = p.Results.target;
            density = p.Results.density;

            dipole = false;
            if ( isempty(target) )
                target = cell2mat(source.x);
                dipole = true;
            end
            dipole = dipole || isequal(cell2mat(source.x), target);

            % Unpack things
            M     = size(target,1);
            N     = source.N * source.np;
            x     = cell2mat(source.x);
            w     = cell2mat(source.w);
            nx    = cell2mat(source.normal);
            kappa = cell2mat(source.curvature);

            % Fill in single-layer potential
            u = zeros(M,N);
            for j = 1:N
                r = target - x(j,:);
                u(:,j) = log(1./sqrt(sum(r.^2,2)))/(2*pi)*w(j);
            end

            if ( ~isempty(density) )
                u = u * density;
            end

        end

    end

end
