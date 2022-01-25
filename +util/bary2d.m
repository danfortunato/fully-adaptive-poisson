function vals = bary2d(chebvals, x, y)

vals = 0*x;
C = chebtech2.vals2coeffs( chebtech2.vals2coeffs( chebvals ).' ).';
Cy = chebtech2.clenshaw(y(:), C).';
for k = 1:numel(x)
    vals(k) = chebtech2.clenshaw(x(k), Cy(:,k));
end

end
