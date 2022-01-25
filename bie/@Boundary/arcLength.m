function out = arcLength(S, kk)
%ARCLENGTH   Arc length of a Boundary.
%   ARCLENGTH(S) computes the arc length of the Boundary S.
%
%   ARCLENGTH(S, KK) computes the arc length of the panels of S with
%   indices given by KK.

if ( nargin == 1 )
    kk = 1:S.np;
end

if ( strcmpi(S.rule, 'ptr') )
    out = sum(S.w{1});
else
    out = zeros(numel(kk),1);
    idx = 1;
    for k = kk(:).'
        out(idx) = sum(S.w{k});
        idx = idx+1;
    end
end

if ( nargin == 1 )
    out = sum(out);
end

end
