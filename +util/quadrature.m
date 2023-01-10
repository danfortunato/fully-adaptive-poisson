function [x, w, breaks] = quadrature(N, rule, np, resfun)
%QUADRATURE   Quadrature rule on [0, 2pi].

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
            %breaks = util.split(resfun, 0, 2*pi, N);
            if ( isa(resfun{1}, 'chebfun') )
                dom = domain(resfun{1});
                a = dom(1);
                b = dom(end);
            else
                a = 0;
                b = 2*pi;
            end
            breaks = util.split(resfun, a, b, N);
            %breaks = util.levelrestrict(breaks, N, resfun{2});
        else
            % Uniform panels
            breaks = 2*pi*(0:np)/np;
        end
        [x01, w01] = util.gauss(N, 0, 1);
%         [x01, w01] = chebpts(N, [0 1]); w01 = w01.';
%         [x01, w01] = chebpts(N, [0 1], 1); w01 = w01.';
        np = numel(breaks)-1;
        x = cell(np,1);
        w = cell(np,1);
        for k = 1:np
            x{k} = (breaks(k+1)-breaks(k))*x01 + breaks(k);
            w{k} = (breaks(k+1)-breaks(k))*w01;
        end
    otherwise
        error('Unknown quadrature rule.');
end
