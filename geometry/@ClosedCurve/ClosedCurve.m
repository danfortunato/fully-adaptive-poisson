classdef ClosedCurve
%CLOSEDCURVE   Represent a closed curve in 2D through parametrization.

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% CLASS PROPERTIES:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    properties ( Access = public )

        % Boundary info

        z        % Parametrization of curve on [0, 2pi)
        dz       % Parametrization of 1st derivative
        dzz      % Parametrization of 2nd derivative
        levelset % Level set defining boundary
        polygon  % Polygon of boundary

        % Discretization info

        rule      % Type of quadrature
        N         % Number of nodes per panel
        np        % Number of panels
        s         % Quadrature nodes (in reference space)
        x         % Quadrature nodes (in real space)
        w         % Quadrature weights
        cx        % Complex quadrature nodes
        cw        % Complex quadrature weights
        velocity  % Parametrization velocity at quadrature nodes
        speed     % Parametrization speed at quadrature nodes
        normal    % Normal vectors at quadrature nodes
        curvature % Curvature at quadrature nodes

    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% CLASS CONSTRUCTOR
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    methods

        function obj = ClosedCurve(varargin)
        %CLOSEDCURVE   Class constructor for a CLOSEDCURVE.

            if ( nargin == 1 || nargin == 2 )

                % Call is CLOSEDCURVE(x)
                x = varargin{1};

                pref = chebfunpref();
                if ( iscell(x) )
                    obj.x = x;
                    obj.N = size(x{1},1);
                    obj.np = size(x,1);
                    obj.rule = 'panel';
                    pref.tech = @chebtech2;
                else
                    obj.x = {x};
                    obj.N = size(x,1);
                    obj.np = 1;
                    obj.rule = 'ptr';
                    pref.tech = @trigtech;
                end

                [s, w, breaks] = quadrature(obj.N, obj.rule, obj.np);

                obj.s = s;
                obj.w         = cell(obj.np,1);
                obj.cx        = cell(obj.np,1);
                obj.cw        = cell(obj.np,1);
                obj.velocity  = cell(obj.np,1);
                obj.speed     = cell(obj.np,1);
                obj.normal    = cell(obj.np,1);
                obj.curvature = cell(obj.np,1);
                for k = 1:obj.np
                    obj.cx{k} = obj.x{k}(:,1)+obj.x{k}(:,2)*1i;
                end
                vals = obj.cx;
                if ( strcmpi(obj.rule, 'panel') )
                    vals = cellfun(@legvals2chebvals, vals, 'UniformOutput', false);
                end
                xfun   = chebfun(vals, breaks, pref);
                dxfun  = diff(xfun);
                dxxfun = diff(dxfun);
                for k = 1:obj.np
                    ds  = feval(dxfun,  s{k}); ds  = [ real(ds)  imag(ds)  ];
                    dss = feval(dxxfun, s{k}); dss = [ real(dss) imag(dss) ];
                    obj.velocity{k}  = ds;
                    obj.speed{k}     = sqrt(sum(ds.^2,2));
                    obj.normal{k}    = [ ds(:,2), -ds(:,1) ] ./ obj.speed{k};
                    obj.curvature{k} = -dot(dss, obj.normal{k}, 2) ./ obj.speed{k}.^2;
                    obj.w{k}  = w{k} .* obj.speed{k};
                    obj.cw{k} = w{k} .* obj.velocity{k};
                end

            elseif ( nargin > 4 )

                % Call is CLOSEDCURVE(N, z, dz, dzz, levelset, ...)

                p = inputParser;
                p.KeepUnmatched = true;
                checkquad = @(s) any(strcmp(s, {'ptr', 'panel'}));
                checkpanels = @(s) strcmp(s, 'auto') || isfloat(s);
                addRequired(p, 'N', @isfloat);
                addRequired(p, 'z',   @(f) isa(f, 'function_handle'));
                addRequired(p, 'dz',  @(f) isa(f, 'function_handle'));
                addRequired(p, 'dzz', @(f) isa(f, 'function_handle'));
                addRequired(p, 'levelset', @(f) isa(f, 'function_handle'));
                addParameter(p, 'quadrature', 'ptr', checkquad);
                addParameter(p, 'panels', 5, checkpanels);
                parse(p, varargin{:});

                obj.N = p.Results.N;
                obj.z = p.Results.z;
                obj.dz = p.Results.dz;
                obj.dzz = p.Results.dzz;
                obj.levelset = p.Results.levelset;
                obj.rule = p.Results.quadrature;

                % Quadrature rule
                [s,w] = quadrature(obj.N, obj.rule, p.Results.panels);

                np = size(s,1);
                obj.x         = cell(np,1);
                obj.w         = cell(np,1);
                obj.cx        = cell(np,1);
                obj.cw        = cell(np,1);
                obj.velocity  = cell(np,1);
                obj.speed     = cell(np,1);
                obj.normal    = cell(np,1);
                obj.curvature = cell(np,1);
                for k = 1:np
                    ds  = obj.dz(s{k});
                    dss = obj.dzz(s{k});
                    obj.velocity{k}  = ds;
                    obj.speed{k}     = sqrt(sum(ds.^2,2));
                    obj.normal{k}    = [ ds(:,2), -ds(:,1) ] ./ obj.speed{k};
                    obj.curvature{k} = -dot(dss, obj.normal{k}, 2) ./ obj.speed{k}.^2;
                    obj.x{k} = obj.z(s{k});
                    obj.w{k}  = w{k} .* obj.speed{k};
                    obj.cx{k} = obj.x{k}(:,1)+obj.x{k}(:,2)*1i;
                    obj.cw{k} = w{k} .* obj.velocity{k};
                end

                obj.np = np;
                obj.s  = s;

            end

            nodes = cell2mat(obj.x);
            obj.polygon = polyshape(nodes(:,1), nodes(:,2));

        end

    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% METHODS
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    methods ( Access = public, Static = false )

        function ii = isinterior(S, x, y)
        %ISINTERIOR   Find interior points.

            if ( ~all(size(x) == size(y)) )
                error('Number of X and Y points must match.');
            end

            if ( ~isempty(S.levelset) )
                ii = S.levelset(x,y) < 10*eps;
            else
                warning('Level set not found. Defaulting to polygon.');
                ii = isinterior(S.polygon, x(:), y(:));
                ii = reshape(ii, size(x));
            end

        end

        function obj = perturb(S, width)
        %PERTURB   Construct a ClosedCurve by perturbing another.

            newxy = cell(S.np,1);
            for k = 1:S.np
                newxy{k} = S.x{k} - width .* S.normal{k};
            end

            if ( strcmpi(S.rule, 'ptr') )
                newxy = newxy{1};
            end

            obj = ClosedCurve(newxy);

        end

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
            hold on
            for k = 1:S.np
                plot(x{k}(:,1), x{k}(:,2), style, varargin{:})
            end
            hold off
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

            levelset = @(x,y) sqrt((x/a).^2+(y/b).^2) - 1;

            obj = ClosedCurve(N, z, dz, dzz, levelset, varargin{:});

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

            levelset = @(x,y) sqrt(x.^2+y.^2) - r(atan2(y,x));

            obj = ClosedCurve(N, z, dz, dzz, levelset, varargin{:});

        end

    end

end
