function initialize(S, varargin)
%INITIALIZE   Initialize and solve local subproblems.
%   INITIALIZE(S, RHS) will initialize the STRIPSOLVER object S on each of
%   the subpatches of S.domain using the righthand side RHS.
%
%   INITIALIZE(S) assumes the problem is homogeneous (i.e., RHS = 0).
%
%   The full sequence for solving a problem using a STRIPSOLVER object S
%   is:
%
%      initialize(S, RHS)
%      build(S)   % (optional)
%      sol = S\bc % or sol = solve(S, bc)
%
%   See also BUILD, SOLVE.

tic
% Initialize all leaf patches:
S.patches = StripSolver.Leaf.initialize(S.domain, varargin{:});
%S.patches = StripSolver.Leaf.parinitialize(S.domain, varargin{:});
toc

end
