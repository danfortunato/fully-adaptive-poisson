function chebvals = modchebvals2chebvals(chebvals)

% Interpolation operator for corner values:
n = size(chebvals,1)+2;
X = chebpts(n);
Xii = X(2:(n-1));
B = [Xii-1, -1-Xii].';
B(:,1:2:end) = -B(:,1:2:end);
if ( mod(n-1, 2) )
    B(2,:) = -B(2,:);
end

corners = B*chebvals;
chebvals = [corners(1,:); chebvals; corners(end,:)];

end
