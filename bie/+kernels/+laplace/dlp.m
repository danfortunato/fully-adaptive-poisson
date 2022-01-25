function u = dlp(source, varargin)
%DLP   Evaluate the Laplace double-layer potential.

% Parse optional inputs
p = inputParser;
addRequired(p,  'source', @(s) isa(s,'Boundary'));
addParameter(p, 'target',    [],    @isnumeric);
addParameter(p, 'density',   [],    @isnumeric);
addParameter(p, 'modified',  false, @islogical);
addParameter(p, 'closeeval', false, @islogical);
addParameter(p, 'side',      [],    @(s) any(strcmp(s,{'i','e'})) || isempty(s));
parse(p, source, varargin{:});
target    = p.Results.target;
density   = p.Results.density;
modified  = p.Results.modified;
closeeval = p.Results.closeeval;
side      = p.Results.side;

if ( ~isreal(target) )
    target = [real(target(:)) imag(target(:))];
end

nystrom = false;
if ( isempty(target) )
    nystrom = true;
end
nystrom = nystrom || isequal(cell2mat(source.x), target);

if ( nystrom )
    % Build Nystrom matrix directly
    u = DLP_nystrom(source);
    % Correction for modified DLP
    if ( modified )
        w = cell2mat(source.w);
        u = u + w.';
    end
elseif ( ~isempty(density) )
    % Evaluate density at targets
    if ( closeeval )
        % Use close evaluation
        if ( isempty(side) )
            % If side was not given, try to determine it
            side = 'e';
            if ( isinterior(source, target(1,1), target(1,2)) )
                side = 'i';
            end
        end
        if ( strcmpi(source.rule, 'ptr') )
            u = DLP_eval_close_global(source, target, density, side);
        elseif ( strcmpi(source.rule, 'panel') )
            u = DLP_eval_close_panel(source, target, density, side);
        end
    else
        % Use FMM
        u = FMM_eval(source, target, 'dipstr', density);
    end
    % Correction for modified DLP
    if ( modified )
        w = cell2mat(source.w);
        u = u + w.' * density;
    end
end

% if ( closeeval && ~isempty(density) )
% 
%     if ( strcmpi(source.rule, 'ptr') )
% 
%         if ( isempty(side) )
%             % If side was not given, try to determine it
%             side = 'e';
%             if ( isinterior(source, target(1,1), target(1,2)) )
%                 side = 'i';
%             end
%         end
%         % Use global close evaluation based on Helsing
%         u = Laplace_DLP_close_global(source, target, density, side);
% 
%     elseif ( strcmpi(source.rule, 'panel') )
% 
%         % Use panel close evaluation from Ludvig
%         u = Laplace_DLP_close_panel(source, target, density, side);
% 
%     end
% 
%     % Correction for modified DLP
%     if ( modified )
%         u = u + w.' * density;
%     end
% 
% else
% end

end