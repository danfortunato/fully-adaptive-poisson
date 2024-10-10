function [uu, xx, yy, inGamma, inGamma1, inStrip] = sampleUniform(S, u, nbulk, dom)

if ( nargin < 3 )
    nbulk = 400;
end

if ( nargin < 4 )
    dom = S.domain;
end

% Sample on a uniform grid of size nbulk x nbulk
[xx, yy] = meshgrid(linspace(dom(1), dom(2), nbulk), ...
                    linspace(dom(3), dom(4), nbulk));

[uu, inGamma, inGamma1, inStrip] = S.sample(u, xx, yy);

end
