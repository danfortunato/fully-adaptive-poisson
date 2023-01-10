function dom = boundingbox(S, scale)
%BOUNDINGBOX   Compute the bounding box of a Boundary.

if ( nargin == 1 )
    scale = 1;
end

xy = cell2mat(S.x);
x = xy(:,1);
y = xy(:,2);
dom = [min(x) max(x) min(y) max(y)];

% Shift to origin
cx = mean(dom(1:2));
cy = mean(dom(3:4));
dom = dom - [cx cx cy cy];

% Scale
dom = scale * dom;

% Shift back
dom = dom + [cx cx cy cy];

end
