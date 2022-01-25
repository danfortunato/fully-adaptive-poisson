function [inGamma, inGamma1, inStrip] = partition(S, x, y)
%PARTITION   Partition a set of points.
%   [INGAMMA, INGAMMA1, INSTRIP] = PARTITION(S, X, Y) returns logical
%   arrays of size SIZE(X) which partition the arrays of points X and Y
%   according the partition-of-unity regions defined by S. The k-th entry
%   of each array is true if the point (X(k), Y(k)) is inside S.Gamma,
%   inside S.Gamma1, or in the strip between the two, respectively. The
%   arrays X and Y must have the same size.

if ( ~all(size(x) == size(y)) )
    error('X and Y arrays must have the same size.'); 
end

inGamma  = isinterior(S.Gamma_re,  x, y);
inGamma1 = isinterior(S.Gamma1, x, y);
inStrip  = inGamma & ~inGamma1;

end
