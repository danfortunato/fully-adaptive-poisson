function u = slp(source, varargin)
%SLP   Evaluate the Laplace single-layer potential.

% TODO: Needs singular quadrature.

% Parse optional inputs
p = inputParser;
addRequired(p,  'source', @(s) isa(s,'Boundary'));
addParameter(p, 'target',    [],    @isnumeric);
addParameter(p, 'density',   [],    @isnumeric);
addParameter(p, 'closeeval', false, @islogical);
addParameter(p, 'side',      [],    @(s) any(strcmp(s,{'i','e'})) || isempty(s));
parse(p, source, varargin{:});
target    = p.Results.target;
density   = p.Results.density;
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
    if ( isempty(side) ), side = 'i'; end
    u = SLP_nystrom(source, side);
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
            u = SLP_eval_close_global(source, target, density, side);
        elseif ( strcmpi(source.rule, 'panel') )
            u = SLP_eval_close_panel(source, target, density, side);
            %error('SLP panel-based close evaluation is not implemented.');
        end
    else
%                     x = cell2mat(source.x);
%                     w = cell2mat(source.w);
%                     t = target(:,1)+target(:,2)*1i;
%                     x = x(:,1)+x(:,2)*1i;
%                     d = bsxfun(@minus,t,x.'); % C-# displacements mat
%                     u = bsxfun(@times, log(abs(d)), -(1/2/pi)*w(:)');  % prefactor & wei
%                     u = u * density;
        % Use FMM
        u = FMM_eval(source, target, 'charge', density);
    end
end

end