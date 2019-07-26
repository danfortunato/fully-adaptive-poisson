function [x,w,breaks] = quadrature(N, rule, np)
%QUADRATURE   Quadrature rule on [0,2pi].

if ( nargin < 3 )
    np = 1;
    if ( nargin < 2 )
        rule = 'ptr';
    end
end

switch rule
    case 'ptr'
        x{1} = 2*pi/N*(0:N-1).';
        w{1} = 2*pi/N*ones(N,1);
        breaks = [0 2*pi];
    case 'panel'
        if ( strcmp(np, 'auto') )
            % Adaptive panels
            % f =  ...
            % t = split(f, 0, 2*pi, N);
        else
            % Uniform panels
            breaks = 2*pi*(0:np)/np;
        end
        [x01,w01] = gauss(N, 0, 1);
        x = cell(np,1);
        w = cell(np,1);
        for k = 1:np
            x{k} = (breaks(k+1)-breaks(k))*x01 + breaks(k);
            w{k} = (breaks(k+1)-breaks(k))*w01;
        end
    otherwise
        error('Unknown quadrature rule.');
end
