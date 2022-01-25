function varargout = mldivide(varargin)
%MLDIVIDE   Perform a partition-of-unity solve using function intension.
%   MLDIVIDE is a shorthand for SOLVE().
%
%   See also SOLVE.

[varargout{1:nargout}] = solve(varargin{:});

end
