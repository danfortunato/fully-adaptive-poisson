classdef ClosedCurve
%CLOSEDCURVE   Represent a closed curve in 2D through parametrization.

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% CLASS PROPERTIES:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    properties ( Access = public )

        z   % Parametrization of curve on [0, 2pi)
        dz  % Parametrization of 1st derivative
        dzz % Parametrization of 2nd derivative

        N   % Number of nodes
        np  % Number of nodes per panel
        s   % Quadrature nodes (in reference space)
        x   % Quadrature nodes (in real space)
        w   % Quadrature weights
        cw  % Quadrature weights (vector-valued)

        velocity  % Parametrization velocity at quadrature nodes
        speed     % Parametrization speed at quadrature nodes
        normal    % Normal vectors at quadrature nodes
        curvature % Curvature at quadrature nodes

    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% CLASS CONSTRUCTOR
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    methods

        function obj = ClosedCurve(N, z, dz, dzz, varargin)
        %CLOSEDCURVE   Class constructor for the @ClosedCurve class.
            obj.z   = z;
            obj.dz  = dz;
            obj.dzz = dzz;

            p = inputParser;
            p.KeepUnmatched = true;
            checkquad = @(s) any(strcmp(s, {'ptr', 'panel'}));
            checkpanels = @(s) strcmp(s, 'auto') || isfloat(s);
            addRequired(p, 'N', @isfloat);
            addRequired(p, 'z',   @(f) isa(f, 'function_handle'));
            addRequired(p, 'dz',  @(f) isa(f, 'function_handle'));
            addRequired(p, 'dzz', @(f) isa(f, 'function_handle'));
            addParameter(p, 'quadrature', 'ptr', checkquad);
            addParameter(p, 'panels', 5, checkpanels);
            parse(p, N, z, dz, dzz, varargin{:});

            % Quadrature rule
            [s,w] = quadrature(N, p.Results.quadrature, p.Results.panels);

            np = size(s,1);
            obj.x         = cell(np,1);
            obj.w         = cell(np,1);
            obj.velocity  = cell(np,1);
            obj.speed     = cell(np,1);
            obj.normal    = cell(np,1);
            obj.curvature = cell(np,1);
            for k = 1:np
                ds  = dz(s{k});
                dss = dzz(s{k});
                obj.velocity{k}  = ds;
                obj.speed{k}     = sqrt(sum(ds.^2,2));
                obj.normal{k}    = [ ds(:,2), -ds(:,1) ] ./ obj.speed{k};
                obj.curvature{k} = -dot(dss, obj.normal{k}, 2) ./ obj.speed{k}.^2;
                obj.x{k} = z(s{k});
                obj.w{k}  = w{k} .* obj.speed{k};
                obj.cw{k} = w{k} .* obj.velocity{k};
            end

            obj.N  = N;
            obj.np = np;
            obj.s  = s;

        end
        
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% METHODS
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    methods ( Access = public, Static = false )

        function varargout = plot(S, varargin)
        %PLOT   Plot a ClosedCurve.
        %   PLOT(S) plots the ClosedCurve S.
        %
        %   PLOT(S, C), where C is a single character string chosen from
        %   the list 'r', 'g', 'b', 'c', 'm', 'y', 'w', 'k', or an RGB row
        %   vector triple, [r g b], fills the domain with the constant
        %   specified color.
        %
        %   PLOT(S, C, PROP1, VAL1, PROP2, VAL2) allows adjusting of the
        %   plot in any way allowed by the built-in PLOT() method.
        %
        %   H = PLOT(S, ...) returns a figure handle of the form returned
        %   by H = PLOT(...), where PLOT() is the built-in MATLAB method.

            % Choose a style:
            if ( nargin > 1 )
                style = varargin{1};
                varargin(1) = [];
            else
                style = '-o';
            end

            holdState = ishold();

            %x = [S.x; S.x(1,:)]; % Wrap around
            x = S.x;
            x{end} = [x{end}; x{1}(1,:)];
            for k = 1:S.np
                plot(x{k}(:,1), x{k}(:,2), style, varargin{:})
            end
            axis equal

            % Return hold state:
            if ( ~holdState), hold off, end
            % Assign output if required:
            if ( nargout > 0 ), varargout = {h}; end
            
        end

    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% STATIC METHODS
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    methods ( Access = public, Static = true )

        function obj = circle(N, varargin)
        %CIRCLE   Set up a circle.

            obj = ClosedCurve.ellipse(N, 1, 1, varargin{:});

        end

        function obj = ellipse(N, a, b, varargin)
        %ELLIPSE   Set up an ellipse.

            if ( nargin == 1 )
                a = 1;
                b = 1;
            end

            z   = @(t) [  a.*cos(t),  b.*sin(t) ];
            dz  = @(t) [ -a.*sin(t),  b.*cos(t) ];
            dzz = @(t) [ -a.*cos(t), -b.*sin(t) ];

            obj = ClosedCurve(N, z, dz, dzz, varargin{:});

        end

        function obj = star(N, varargin)
        %SMOOTHSTAR   Set up a smooth star.

            p = inputParser;
            p.KeepUnmatched = true;
            addRequired(p, 'N', @isfloat);
            addParameter(p, 'arms', 5, @isfloat);
            addParameter(p, 'wobble', 0.3, @isfloat);
            parse(p, N, varargin{:});
            k = p.Results.arms;
            a = p.Results.wobble;

            r   = @(t)  1 + a*cos(k*t);
            dr  = @(t)   -k*a*sin(k*t);
            drr = @(t) -k^2*a*cos(k*t);

            z   = @(t) [ r(t).*cos(t), ...
                         r(t).*sin(t) ];
            dz  = @(t) [ dr(t).*cos(t) - r(t).*sin(t), ...
                         dr(t).*sin(t) + r(t).*cos(t) ];
            dzz = @(t) [ drr(t).*cos(t) - 2*dr(t).*sin(t) - r(t).*cos(t), ...
                         drr(t).*sin(t) + 2*dr(t).*cos(t) - r(t).*sin(t) ];

            obj = ClosedCurve(N, z, dz, dzz, varargin{:});

        end

    end

end
