function x = levelrestrict(x, N, df)

fac = 2;
kind = 'arclength';
%kind = '';
npmax = 2000;
maxiter = 1000;
[xs, ws] = legpts(N);

np = length(x)-1;

ab = zeros(2,npmax);
ab(1,1:np) = x(1:end-1);
ab(2,1:np) = x(2:end);

adjs = zeros(2,npmax);
adjs(1,1:np) = [np 1:np-1];
adjs(2,1:np) = [2:np 1];

iter = 0;
done = false;
while ( ~done && iter < maxiter )

    done = true;
    npold = np;
    for k = 1:npold
        kl = adjs(1,k);
        kr = adjs(2,k);

        a = ab(1,k);
        b = ab(2,k);

        if ( strcmpi(kind, 'arclength') )
            self  = arclength(df, a, b, xs, ws);
            left  = arclength(df, ab(1,kl), ab(2,kl), xs, ws);
            right = arclength(df, ab(1,kr), ab(2,kr), xs, ws);
        else
            self  = b-a;
            left  = ab(2,kl)-ab(1,kl);
            right = ab(2,kr)-ab(1,kr);
        end

        if ( self > fac*left || self > fac*right )

            if ( np == npmax )
                error('Too many panels.');
            end

            done = false;

            knew = np+1;
            adjs(2,k) = knew;
            adjs(:,knew) = [k kr];
            adjs(1,kr) = knew;

            ab(:,k)    = [a (a+b)/2];
            ab(:,knew) = [(a+b)/2 b];

            np = np+1;
        end

    end

end

%ab = ab(:,1:np);
%adjs = adjs(:,1:np);
x = [sort(ab(1,1:np),2), x(end)];

end

function len = arclength(df, a, b, xs, ws)
    ts = a + (b-a)*(xs+1)/2;
    %dsdt = sqrt(sum(df(ts).^2,1));
    dsdt = df(ts);
    len = dot(dsdt, ws) * (b-a)/2;
end
