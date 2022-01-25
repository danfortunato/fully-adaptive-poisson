function [u, f] = pfreegauss(d, x0, y0)

    function f = ffunc(x, y)
        f = exp(-((x-x0).^2 + (y-y0).^2)/d);
    end

    function u = ufunc(x, y)
        gamma = 0.57721566490153286061;
        r = (x-x0).^2 + (y-y0).^2;
        if ( r/d < 1e-7 )
            u = -(d/4)*(-gamma + log(d) + r/d);
        else
            v = expint(r/2);
            u = -(d/4)*(v + log(r));
        end
        u = -u;
    end

u = @ufunc;
f = @ffunc;

end
