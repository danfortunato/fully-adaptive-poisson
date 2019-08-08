function [tt, rr] = cart2curv(Gamma, xfun, yfun, xx, yy)
%CART2CURV   Transform Cartesian to curved coordinates.

tol = eps;
maxiter = 5; % Sometimes Newton's method stalls... why?

xd = xx(:) - Gamma.x{1}(:,1).';
yd = yy(:) - Gamma.x{1}(:,2).';
D = sqrt(xd.^2 + yd.^2);
[~, guess_ind] = min(D,[],2);
tt = Gamma.s{1}(guess_ind);
xdg = xx - Gamma.x{1}(guess_ind,1);
ydg = yy - Gamma.x{1}(guess_ind,2);
rr = sqrt(xdg.^2 + ydg.^2);
xo = xfun(tt, rr);
yo = yfun(tt, rr);
remx = xo - xx;
remy = yo - yy;
rem = max(abs(remx.^2 + remy.^2));

iter = 0;
while ( rem > tol && iter < maxiter )

    % Compute the Jacobian
    dxdt = feval(diff(xfun,1,2), tt, rr);
    dxdr = feval(diff(xfun,1,1), tt, rr);
    dydt = feval(diff(yfun,1,2), tt, rr);
    dydr = feval(diff(yfun,1,1), tt, rr);
    J = dxdr.*dydt - dxdt.*dydr;
    % Compute d[t,r]/d[x,y] as the inverse of the Jacobian
    dtdx = -dydr ./ J;
    dtdy =  dxdr ./ J;
    drdx =  dydt ./ J;
    drdy = -dxdt ./ J;

    % Line search
    delt = -[ dtdx.*remx + dtdy.*remy, drdx.*remx + drdy.*remy ];
    line_factor = 1;
    while ( true )
        t_new = tt + line_factor*delt(:,1);
        r_new = rr + line_factor*delt(:,2);
        xo = xfun(t_new, r_new);
        yo = yfun(t_new, r_new);
        remx = xo - xx;
        remy = yo - yy;
        rem_new = max(abs(remx.^2 + remy.^2));
        if ( rem_new < (1-0.5*line_factor)*rem || line_factor < 1e-4 )
            tt = t_new;
            rr = r_new;
            rem = rem_new;
            break
        end
        line_factor = 0.5*line_factor;
    end
    iter = iter+1;
end

tt(tt<0) = tt(tt<0) + 2*pi;

end
