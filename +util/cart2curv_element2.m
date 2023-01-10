function [tt, rr, converged] = cart2curv_element2(xfun, yfun, xx, yy, init, varargin)

if ( nargin == 8 )
    [dxfun_dt, dxfun_dr, dyfun_dt, dyfun_dr] = deal(varargin{:});
else
    dxfun_dt = diff(xfun, 1, 2);
    dxfun_dr = diff(xfun, 1, 1);
    dyfun_dt = diff(yfun, 1, 2);
    dyfun_dr = diff(yfun, 1, 1);
end

solver_type = 'NewtonElementProject';
ref_tol = 1e-14;
phys_rtol = 1e-10;
maxiter = 10;
%init = [0 0];
tt = init(:,1) + 0*xx;
rr = init(:,2) + 0*xx;
flag = zeros(size(xx));

for k = 1:numel(xx)

    hit_bdr = false;
    prev_hit_bdr = false;
    phys_tol = phys_rtol * max(abs([xx(k) yy(k)]));
    %phys_tol = phys_rtol;

    it = 0;
    while ( true )
        x = xfun(tt(k),rr(k));
        y = yfun(tt(k),rr(k));
        x = xx(k) - x;
        y = yy(k) - y;

        err_phys = max(max(abs(x), [], 'all'), max(abs(y), [], 'all'));

        if ( err_phys < phys_tol )
            if ( ~strcmpi(solver_type, 'Newton') )
                flag(k) = 1;
            else
                error('unimplemented')
                %return Geometry::CheckPoint(geom, ip, ip_tol) ? Inside : Outside;
            end
            break
        end
    
        if ( hit_bdr )
            if ( prev_hit_bdr || it == maxiter )
                dt = tt(k) - prev_tt;
                dr = rr(k) - prev_rr;
                real_dt_norm = max(max(abs(dt), [], 'all'), max(abs(dr), [], 'all'));
                %real_dt_norm
                if ( prev_hit_bdr && real_dt_norm < ref_tol )
                    warning('Newton: *** stuck on boundary!');
                    flag(k) = -1;
                    break
                end
            end
        end
    
        if ( it == maxiter )
            warning('maxiter reached');
            break
        end
    
        if ( ~strcmpi(solver_type, 'Newton') )
            prev_tt = tt(k);
            prev_rr = rr(k);
            prev_hit_bdr = hit_bdr;
        end

        % Compute the Jacobian
        dxdt = feval(dxfun_dt, tt(k), rr(k));
        dxdr = feval(dxfun_dr, tt(k), rr(k));
        dydt = feval(dyfun_dt, tt(k), rr(k));
        dydr = feval(dyfun_dr, tt(k), rr(k));
        J = dxdr.*dydt - dxdt.*dydr;

        % Compute d[t,r]/d[x,y] as the inverse of the Jacobian
        dtdx = -dydr ./ J;
        dtdy =  dxdr ./ J;
        drdx =  dydt ./ J;
        drdy = -dxdt ./ J;

        dt = dtdx.*x + dtdy.*y;
        dr = drdx.*x + drdy.*y;
        tt(k) = tt(k) + dt;
        rr(k) = rr(k) + dr;
        it = it + 1;

        % Perform projection based on solver_type:
        switch ( solver_type )
            case 'Newton'
                break
            case 'NewtonSegmentProject'
                error('unimplemented')
                %hit_bdr = !Geometry::ProjectPoint(geom, prev_xip, xip);
            case 'NewtonElementProject'
                [tt(k), rr(k), inside] = projectPointElement(tt(k), rr(k));
                hit_bdr = ~inside;
            otherwise
                error("Invalid solver type.");
        end

        % Check for convergence in reference coordinates:
        dt_norm = max(max(abs(dt), [], 'all'), max(abs(dr), [], 'all'));
        if ( dt_norm < ref_tol )
            if ( ~strcmpi(solver_type, 'Newton') )
                flag(k) = 1;
            else
                error('unimplemented')
                %flag = CheckPoint(geom, ip, ip_tol) ? Inside : Outside;
            end
        end
    end
    
    %it
end

converged = (flag == 1);

end

function [t, r, inside] = projectPointElement(t, r)

if (t < -1)
    in_x = false;
    t = -1;
elseif (t > 1)
    in_x = false;
    t = 1;
else
    in_x = true;
end

if (r < -1)
    in_y = false;
    r = -1;
elseif (r > 1)
    in_y = false;
    r = 1;
else
    in_y = true;
end

inside = in_x && in_y;

end

function [t, r, inside] = projectPointSegment(t_beg, r_beg, t, r)

lend = [t, r, 1-t, 1-r];
lbeg = [t_beg, r_beg, 1-t_beg, 1-r_beg];
s = 1;
out = false;
for i = 1:4
    lbeg(i) = max(lbeg(i), 0); % remove round-off
    if (lend(i) < 0)
        out = true;
        s = min(s, lbeg(i)/(lbeg(i)-lend(i)));
    end
end

if (out)
    t = s*lend(1) + (1-s)*lbeg(1);
    r = s*lend(2) + (1-s)*lbeg(2);
end
inside = ~out;
end
