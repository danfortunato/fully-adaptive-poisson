classdef Boundary
%BOUNDARY   Represent a closed curve in 2D through parametrization.

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% CLASS PROPERTIES:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    properties ( Access = public )

        % Boundary info

        f        % Parametrization of curve on [0, 2pi)
        df       % Parametrization of 1st derivative
        dff      % Parametrization of 2nd derivative
        bend
        levelset % Level set defining boundary
        polygon  % Polygon of boundary
        polynodes
        polyedges

        % Discretization info

        rule      % Type of quadrature
        N         % Number of nodes per panel
        np        % Number of panels
        breaks    % Locations of panel breaks
        zbreaks   % Locations of complex panel breaks
        x         % Quadrature nodes        (vector form)
        dx        % 1st derivative at nodes (vector form)
        dxx       % 2nd derivative at nodes (vector form)
        z         % Quadrature nodes        (complex form)
        dz        % 1st derivative at nodes (complex form)
        dzz       % 2nd derivative at nodes (complex form)
        s         % Quadrature nodes in reference space
        w         % Quadrature weights
        cw        % Complex quadrature weights
        speed     % Parametrization speed at quadrature nodes
        accel     % Parametrization acceleration at quadrature nodes
        normal    % Normal vectors at quadrature nodes
        curvature % Curvature at quadrature nodes

    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% CLASS CONSTRUCTOR
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    methods

        function obj = Boundary(varargin)
        %BOUNDARY   Class constructor for a BOUNDARY.
        
            if ( nargin == 0 )
                return
            end

            if ( nargin == 1 || nargin == 2 )

                if ( nargin == 1 )
                    % Call is BOUNDARY(x)
                    x = varargin{1};
                    quadinfo = [];
                else
                    % Call is BOUNDARY(x, quadinfo)
                    x = varargin{1};
                    quadinfo = varargin{2};
                end

                pref = chebfunpref();
                if ( iscell(x) )
                    if ( size(x{1},2) == 1 )
                        x = cellfun(@(c) [real(c) imag(c)], x, 'UniformOutput', false);
                    end
                    obj.x = x;
                    obj.N = size(x{1},1);
                    obj.np = size(x,1);
                    obj.rule = 'panel';
                    pref.tech = @chebtech2;
                else
                    if ( size(x,2) == 1 )
                        x = [real(x) imag(x)];
                    end
                    obj.x = {x};
                    obj.N = size(x,1);
                    obj.np = 1;
                    obj.rule = 'ptr';
                    pref.tech = @trigtech;
                end

                if ( ~isempty(quadinfo) )
                    % We were given quadrature info
                    s = quadinfo.s;
                    w = quadinfo.w;
                    obj.breaks = quadinfo.breaks;
                else
                    % Compute quadrature info
                    [s, w, obj.breaks] = util.quadrature(obj.N, obj.rule, obj.np);
                end

                obj.s = s;
                obj.z              = cell(obj.np,1);
                obj.dz             = cell(obj.np,1);
                obj.dzz            = cell(obj.np,1);
                obj.dx             = cell(obj.np,1);
                obj.dxx            = cell(obj.np,1);
                obj.w              = cell(obj.np,1);
                obj.cw             = cell(obj.np,1);
                obj.speed          = cell(obj.np,1);
                obj.accel          = cell(obj.np,1);
                obj.normal         = cell(obj.np,1);
                obj.curvature      = cell(obj.np,1);

                for k = 1:obj.np
                    obj.z{k} = obj.x{k}(:,1)+obj.x{k}(:,2)*1i;
                end

%                 vals = obj.z;
%                 if ( strcmpi(obj.rule, 'panel') )
%                     vals = cellfun(@legvals2chebvals, vals, 'UniformOutput', false);
%                     %vals = cellfun(@(x) chebvals2chebvals(x, 1, 2), vals, 'UniformOutput', false);
%                 end
%                 obj.f   = chebfun(vals, obj.breaks, pref);
%                 obj.df  = diff(obj.f);
%                 obj.dff = diff(obj.df);

                [~, ~, D] = util.gauss(obj.N);
                for k = 1:obj.np
                    scl = 2 / (obj.breaks(k+1) - obj.breaks(k));
                    obj.dz{k}  = scl * (D * obj.z{k});
                    obj.dzz{k} = scl * (D * obj.dz{k});
                    %obj.dz{k}  = obj.df(s{k});
                    %obj.dzz{k} = obj.dff(s{k});
                    obj.dx{k}  = [ real(obj.dz{k})  imag(obj.dz{k})  ];
                    obj.dxx{k} = [ real(obj.dzz{k}) imag(obj.dzz{k}) ];
                    obj.speed{k} = sqrt(sum(obj.dx{k}.^2, 2));
                    obj.accel{k} = real(obj.dz{k} ./ obj.speed{k} .* conj(obj.dzz{k}));
                    obj.normal{k} = [ obj.dx{k}(:,2), -obj.dx{k}(:,1) ] ./ obj.speed{k};
                    obj.curvature{k} = -dot(obj.dxx{k}, obj.normal{k}, 2) ./ obj.speed{k}.^2;
                    obj.w{k}  = w{k} .* obj.speed{k};
                    obj.cw{k} = w{k} .* obj.dz{k};
                end
                
                [xleg, ~, vleg] = legpts(obj.N);
                obj.zbreaks = bary(-1, [obj.z{:}], xleg, vleg);
                obj.zbreaks(end+1) = obj.zbreaks(1);
                
                % Reparametrize by arc length
%                 alens = arclength(obj, 1:obj.np).';
%                 obj.breaks = 2*pi * [0 cumsum(alens)./sum(alens)];
%                 obj.f   = chebfun(vals, obj.breaks, pref);
%                 obj.df  = diff(obj.f);
%                 obj.dff = diff(obj.df);
%                 for k = 1:obj.np
%                     obj.dz{k}  = obj.df(s{k});
%                     obj.dzz{k} = obj.dff(s{k});
%                     obj.dx{k}  = [ real(obj.dz{k})  imag(obj.dz{k})  ];
%                     obj.dxx{k} = [ real(obj.dzz{k}) imag(obj.dzz{k}) ];
%                     obj.speed{k}     = sqrt(sum(obj.dx{k}.^2,2));
%                     obj.normal{k}    = [ obj.dx{k}(:,2), -obj.dx{k}(:,1) ] ./ obj.speed{k};
%                     obj.curvature{k} = -dot(obj.dxx{k}, obj.normal{k}, 2) ./ obj.speed{k}.^2;
%                     obj.w{k}  = w{k} .* obj.speed{k};
%                     obj.cw{k} = w{k} .* obj.dz{k};
%                 end

            elseif ( nargin > 4 )

                % Call is BOUNDARY(N, z, dz, dzz, levelset, ...)

                p = inputParser;
                p.KeepUnmatched = true;
                checkquad = @(s) any(strcmp(s, {'ptr', 'panel'}));
                checkpanels = @(s) strcmp(s, 'auto') || isfloat(s);
                addRequired(p, 'N', @isfloat);
                addRequired(p, 'f',   @(f) isa(f, 'function_handle') || isa(f, 'chebfun'));
                addRequired(p, 'df',  @(f) isa(f, 'function_handle') || isa(f, 'chebfun'));
                addRequired(p, 'dff', @(f) isa(f, 'function_handle') || isa(f, 'chebfun'));
                addRequired(p, 'levelset', @(f) isa(f, 'function_handle') || isa(f, 'chebfun2') || isempty(f));
                addParameter(p, 'quadrature', 'ptr', checkquad);
                addParameter(p, 'panels', 'auto', checkpanels);
                parse(p, varargin{:});

                obj.N = p.Results.N;
                obj.f = p.Results.f;
                obj.df = p.Results.df;
                obj.dff = p.Results.dff;
                obj.levelset = p.Results.levelset;
                obj.rule = p.Results.quadrature;

                % Quadrature rule
%                 zp = stdpan.D * f{1}(a+L/2*(1+stdpan.nodes));
%                 zpp = stdpan.D^2 * f{1}(a+L/2*(1+stdpan.nodes));
%                 ben = sum(stdpan.weights .* imag(zpp./zp).^2 ./ abs(zp)) / (2*pi); % bending energy (normalized by 2pi)
%                 if ben > 5/(log10(1/tol))  % bw - temporary rule restricting bending energy
%                     split = inhalf;
%                 end
                obj.bend = @(t) abs(imag(obj.dff(t)./obj.df(t))).^2 ./ abs(obj.df(t));
                [s, w, obj.breaks] = util.quadrature(obj.N, obj.rule, p.Results.panels, ...
                    {obj.f, @(t) abs(obj.df(t)), @(t) obj.df(t)./abs(obj.df(t)), obj.bend});
                %[s, w, obj.breaks] = util.quadrature(obj.N, obj.rule, p.Results.panels, ...
                %  {obj.f, @(t) abs(obj.df(t)), obj.dff});
                %[s, w, obj.breaks] = util.quadrature(obj.N, obj.rule, p.Results.panels, {obj.f});

                np = size(s,1);
                obj.x         = cell(np,1);
                obj.dx        = cell(np, 1);
                obj.dxx       = cell(np, 1);
                obj.z         = cell(np, 1);
                obj.dz        = cell(np, 1);
                obj.dzz       = cell(np, 1);
                obj.w         = cell(np, 1);
                obj.cw        = cell(np, 1);
                obj.speed     = cell(np, 1);
                obj.accel     = cell(np, 1);
                obj.normal    = cell(np, 1);
                obj.curvature = cell(np, 1);
                for k = 1:np
                    obj.z{k}   = obj.f(s{k});
                    obj.dz{k}  = obj.df(s{k});
                    obj.dzz{k} = obj.dff(s{k});
                    obj.x{k}   = [ real(obj.z{k})   imag(obj.z{k})   ];
                    obj.dx{k}  = [ real(obj.dz{k})  imag(obj.dz{k})  ];
                    obj.dxx{k} = [ real(obj.dzz{k}) imag(obj.dzz{k}) ];
                    obj.speed{k} = sqrt(sum(obj.dx{k}.^2,2));
                    obj.accel{k} = real(obj.dz{k} ./ obj.speed{k} .* conj(obj.dzz{k}));
                    obj.normal{k}    = [ obj.dx{k}(:,2), -obj.dx{k}(:,1) ] ./ obj.speed{k};
                    obj.curvature{k} = -dot(obj.dxx{k}, obj.normal{k}, 2) ./ obj.speed{k}.^2;
                    obj.w{k}  = w{k} .* obj.speed{k};
                    obj.cw{k} = w{k} .* obj.dz{k};
                end

                obj.np = np;
                obj.s  = s;
                obj.zbreaks = obj.f(obj.breaks(:)).';

            end

%             nodes = cell2mat(obj.x);
%             warnstate = warning('off', 'MATLAB:polyshape:repairedBySimplify');
%             obj.polygon = polyshape(nodes(:,1), nodes(:,2));
%             warning(warnstate);

            nodes = cell2mat(obj.x);
            npts = size(nodes, 1);
            obj.polynodes = nodes;
            obj.polyedges = [1:npts ; 2:npts 1].';
            
%             linnodes = cell(obj.np, 1);
%             [x, ~, w] = legpts(obj.N);
%             y = trigpts(4*obj.N);
%             B = barymat(y, x, w);
%             for k = 1:obj.np
%                 linnodes{k} = B * obj.x{k};
%             end
%             linnodes = cell2mat(linnodes);
%             npts = size(linnodes, 1);
%             obj.polynodes = linnodes;
%             obj.polyedges = [1:npts ; 2:npts 1].';

        end

    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% STATIC METHODS
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    methods ( Static )

        C = circle(N, varargin);
        C = squished_circle(N, r, b, rot);
        C = squircle(N, v, a, b, varargin);
        C = ellipse(N, a, b, varargin);
        C = star(N, varargin);
        C = hand(N, varargin);
        C = multiscale(N, varargin);
        C = multiscale_circle(N, varargin);
        C = wavepacket(N, varargin);
        C = c_shape(N, varargin);
        C = gobbler(N, varargin);
        C = spikey(N, varargin);
        C = saw(N, varargin);
        C = spiral(N, varargin);
        C = larrycup(N, varargin);
        C = bleb(N, varargin);

        C = fromChunker(chnkr);
        chnkr = toChunker(C);

    end

end
