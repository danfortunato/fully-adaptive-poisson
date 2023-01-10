function [tt, rr, converged] = cart2curv_element(xfun, yfun, xx, yy, init, varargin)

tol = 1e-30;
maxiter = 10; % Sometimes Newton's method stalls... why?

% Initial guess is the center of the cell
%tt = 0*xx;
%rr = 0*yy;
tt = init(:,1) + 0*xx;
rr = init(:,2) + 0*xx;
xo = xfun(tt,rr);
yo = yfun(tt,rr);
remx = xo - xx;
remy = yo - yy;
%rem = max(abs(remx.^2 + remy.^2), [], 'all');
rem = abs(remx.^2 + remy.^2);

if ( nargin == 8 )
    [dxfun_dt, dxfun_dr, dyfun_dt, dyfun_dr] = deal(varargin{:});
else
    dxfun_dt = diff(xfun, 1, 2);
    dxfun_dr = diff(xfun, 1, 1);
    dyfun_dt = diff(yfun, 1, 2);
    dyfun_dr = diff(yfun, 1, 1);
end

t_new = tt;
r_new = rr;
%t_new = zeros(size(xx));
%r_new = zeros(size(xx));
rem_new = zeros(size(xx));

iter = 0;
while ( any(rem > tol) && iter < maxiter )

    % Compute the Jacobian
    dxdt = feval(dxfun_dt, tt, rr);
    dxdr = feval(dxfun_dr, tt, rr);
    dydt = feval(dyfun_dt, tt, rr);
    dydr = feval(dyfun_dr, tt, rr);
    J = dxdr.*dydt - dxdt.*dydr;

    % Compute d[t,r]/d[x,y] as the inverse of the Jacobian
    dtdx = -dydr ./ J;
    dtdy =  dxdr ./ J;
    drdx =  dydt ./ J;
    drdy = -dxdt ./ J;
    
%     tt = tt - (dtdx.*xfun(tt,rr) + dtdy.*yfun(tt,rr));
%     rr = rr - (drdx.*xfun(tt,rr) + drdy.*yfun(tt,rr));
%     xo = xfun(tt,rr);
%     yo = yfun(tt,rr);
%     remx = xo - xx;
%     remy = yo - yy;
%     rem = abs(remx.^2 + remy.^2);

    % Line search
    delt_t = dtdx.*remx + dtdy.*remy;
    delt_r = drdx.*remx + drdy.*remy;
    line_factor = ones(size(xx));
    refined = false(size(xx));
    refined(rem <= tol) = true;
    while ( true )
        t_new(~refined) = tt(~refined) - line_factor(~refined).*delt_t(~refined);
        r_new(~refined) = rr(~refined) - line_factor(~refined).*delt_r(~refined);
        
        
%         for k = 1:numel(t_new)
%             while t_new(k) < -1 || t_new(k) > 1
%                 if t_new(k) < -1
%                     t_new(k) = -1 - (t_new(k) + 1);
%                     t_new(k)
%                 elseif t_new(k) > 1
%                     t_new(k) = 1 - (t_new(k) - 1);
%                     t_new(k)
%                 end
%             end
%         end
        
        xo(~refined) = xfun(t_new(~refined), r_new(~refined));
        yo(~refined) = yfun(t_new(~refined), r_new(~refined));
        remx(~refined) = xo(~refined) - xx(~refined);
        remy(~refined) = yo(~refined) - yy(~refined);

        rem_new(~refined) = abs(remx(~refined).^2 + remy(~refined).^2);
        refined(~refined) = (rem_new(~refined) < (1-0.5*line_factor(~refined)).*rem(~refined)) ...
                          | (line_factor(~refined) < 1e-4);
        tt(refined) = t_new(refined);
        rr(refined) = r_new(refined);
        rem(refined) = rem_new(refined);
        if ( all(refined) )
            break
        end

%         rem_new = max(abs(remx.^2 + remy.^2), [], 'all');
%         if ( rem_new < (1-0.5*line_factor)*rem || line_factor < 1e-4 )
%             tt = t_new;
%             rr = r_new;
%             rem = rem_new;
%             break
%         end

        line_factor(~refined) = 0.5*line_factor(~refined);
    end
    iter = iter+1;
end

converged = rem <= tol;
tt(~converged) = nan;
rr(~converged) = nan;

if ( iter == maxiter )
    warning('Maximum number of Newton iterations reached.');
end

end
