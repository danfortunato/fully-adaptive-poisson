function [S, D2N, Aii] = constructOperators(x, y, rhs, B, nthreads)

if ( nargin < 5 )
    nthreads = feature('numcores');
end

[n, ~, numPatches] = size(x);
numIntPts = (n-2)*(n-2) + 4;
numBdyPts = 4*(n-2);
D = diffmat(n);
S   = zeros(n^2,       numBdyPts+1, numPatches);
D2N = zeros(numBdyPts, numBdyPts+1, numPatches);
Aii = zeros(numIntPts, numIntPts,   numPatches);
mex_id_ = 'constructOperators32(i int, i int, i double[], i double[], i double[], i double[], i double[], io double[], io double[], io double[])';
[S, D2N, Aii] = gateway(mex_id_, numPatches, nthreads, x, y, rhs, D, B, S, D2N, Aii);
S   = reshape(S,   [n^2       numBdyPts+1 numPatches]);
D2N = reshape(D2N, [numBdyPts numBdyPts+1 numPatches]);
Aii = reshape(Aii, [numIntPts numIntPts   numPatches]);

end
