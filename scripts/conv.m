function ff = conv(f, x, xx)
%CONV   Circular convolution.
%   H = CONV(F, X, XX) produces the convolution of F:
%                     - 
%                    /
%           H(x) =   |    F(t) G(x-t) dt,  x in [0, 2*pi]
%                    /
%                   -

sigma = 0.1;
g = chebfun(@(t) 1/(sigma*sqrt(2*pi)) * exp(-0.5*(t/sigma).^2), 'trig', [-pi pi]);

end
