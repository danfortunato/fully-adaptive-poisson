C = Boundary.star(7, 'wobble', 0.2, 'quadrature', 'panel');
plot(C)
hold on
for k = 1:C.np
    x = C.x{k};
    nx = C.normal{k};
    %sqrt(nx(:,1).^2 + nx(:,2).^2)
    quiver(x(:,1), x(:,2), nx(:,1), nx(:,2), 2)
end
hold off
shg

%%
C = Boundary.star(7, 'wobble', 0.2, 'quadrature', 'panel');
plot(C)
hold on
x = C.zbreaks;
for k = 1:C.np-1
    % Get the normal vector at the interface between panels k and k+1
    nxL = util.feval(C.normal{k},    1); % The left panel evaluated on the right side
    nxR = util.feval(C.normal{k+1}, -1); % The right panel evaluated on the left side
    % Take the mean and normalize
    nx = (nxL+nxR)/2;
    nx = nx ./ norm(nx);
    quiver(real(x(k+1)), imag(x(k+1)), nx(1), nx(2))
end
shg
