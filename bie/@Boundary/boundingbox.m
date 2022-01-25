function dom = boundingbox(S, scale)
%BOUNDINGBOX   Compute the bounding box of a Boundary.

if ( nargin == 1 )
    scale = 1;
end

xy = cell2mat(S.x);
x = xy(:,1);
y = xy(:,2);
dom = scale * [min(x) max(x) min(y) max(y)];

end
