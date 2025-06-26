function varargout = plot(S, varargin)
%PLOT   Plot a Boundary.
%   PLOT(S) plots the Boundary S.
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
    style = 'k-';
end

lw = 2;

parser = inputParser;
parser.KeepUnmatched = true;
addParameter(parser, 'endpoints', true, @islogical);
parse(parser, varargin{:});
keys = fieldnames(parser.Unmatched);
vals = struct2cell(parser.Unmatched);
v = [keys(:).' ; vals(:).'];
varargin = v(:).';

showEndpoints = parser.Results.endpoints;

holdState = ishold();

% Plot a smooth curve, if we can
% if ( isa(S.f, 'function_handle') )
%     tt = linspace(0, 2*pi, 1000);
%     plot(S.f(tt), 'k-', 'LineWidth', 1);
% else
    x = S.x;
    if ( strcmpi(S.rule, 'ptr') )
        % Wrap around
        x{end} = [x{end}; x{1}(1,:)];
    end
    hold on
    for k = 1:S.np
        x{k} = [real(S.zbreaks(k)) imag(S.zbreaks(k)) ; x{k} ; real(S.zbreaks(k+1)) imag(S.zbreaks(k+1))];
        plot(x{k}(:,1), x{k}(:,2), style, 'LineWidth', lw, 'MarkerSize', 15, varargin{:});
    end
    if ( strcmpi(S.rule, 'panel') && showEndpoints )
        % Plot the panel endpoints
        plot(S.zbreaks, 'k.', 'MarkerSize', 15, varargin{:})
    end
    hold off
% end

axis equal

% Return hold state:
if ( ~holdState ), hold off, end
% Assign output if required:
if ( nargout > 0 ), varargout = {h}; end

end
