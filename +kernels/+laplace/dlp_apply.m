function u = dlp_apply(source, density, varargin)
%DLP_APPLY   Apply the Laplace double-layer potential.

% Parse optional inputs
p = inputParser;
addRequired(p,  'source', @(s) isa(s,'Boundary'));
addRequired(p,  'density',   @isnumeric);
addParameter(p, 'modified',  false, @islogical);
addParameter(p, 'closeeval', false, @islogical);
addParameter(p, 'side',      'i',    @(s) any(strcmp(s,{'i','e'})));
parse(p, source, density, varargin{:});
modified  = p.Results.modified;
closeeval = p.Results.closeeval;
side      = p.Results.side;

% Evaluate density at targets
if ( closeeval )
    % Use close evaluation
    if ( strcmpi(source.rule, 'ptr') )
        u = DLP_apply_close_global(source, density, side);
    elseif ( strcmpi(source.rule, 'panel') )
        u = DLP_apply_close_panel(source, density, side);
    end
else
    % Use FMM
    target = cell2mat(source.x);
    u = FMM_eval(source, target, 'dipstr', density);
end

% Correct self
w     = cell2mat(source.w);
kappa = cell2mat(source.curvature);
self_corr = -kappa/(4*pi).*w;
u = u + self_corr .* density;

% Correction for modified DLP
if ( modified )
    w = cell2mat(source.w);
    u = u + w.' * density;
end

end