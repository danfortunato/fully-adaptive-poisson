function varargout = mldivide(varargin)
%MLDIVIDE   Solve a strip problem.
%   MLDIVIDE is a shorthand for SOLVE().
%
%   See also SOLVE.

[varargout{1:nargout}] = solve(varargin{:});

end
