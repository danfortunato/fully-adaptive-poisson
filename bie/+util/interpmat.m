function R = interpmat(s, t)
%INTERPMAT   Interpolation matrix between node sets.
%   INTERPMAT(S, T) returns a matrix which interpolates values from source
%   nodes S to target nodes T.

N = numel(s); s = s(:);
M = numel(t); t = t(:);
S = ones(N);
T = ones(M,N);
for j = 2:N
    S(:,j) = S(:,j-1).*s; % polyval matrix at source nodes
    T(:,j) = T(:,j-1).*t; % polyval matrix at target nodes
end
R = (S' \ T')'; % Backwards stable

end
