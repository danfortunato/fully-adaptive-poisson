function out = fftfreq(n, d)

    if ( nargin == 1 )
        d = 1;
    end

    if ( mod(n,1) ~= 0)
        error('n should be an integer.');
    end

    val = 1/(n*d);
    results = zeros(n,1);
    k = floor((n-1)/2)+1;
    results(1:k) = 0:k-1;
    results(k+1:n) = -floor(n/2):-1;
    out = results * val;

end
