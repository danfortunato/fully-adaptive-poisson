function handle = getElemCoordinates_(n)
    if (n < 4 || n > 32)
        error(['Template has not been instantiated for n = ' num2str(n)]);
    end
    handle = eval(['@getElemCoordinates' num2str(n)]);
end

% function varargout = getElemCoordinates(n, varargin)
%     persistent handle n
%     if ( isempty(handle) )
%         handle = eval(['@getElemCoordinates' num2str(n)]);
%     end
%     [varargout{1:nargout}] = handle(varargin{:});
% end

function [t, r, found] = getElemCoordinates12(x, y, gamx, gamy, gam1x, gam1y, nx, ny, xcheb, vcheb, nthreads)

npts = size(x, 1);
t = zeros(npts, 1);
r = zeros(npts, 1);
found = zeros(npts, 1);
mex_id_ = 'getElemCoordinates12(i int, i int, i double[], i double[], i double[], i double[], i double[], i double[], i double[], i double[], i double[], i double[], io double[], io double[], io int[])';
[t, r, found] = util.gateway(mex_id_, npts, nthreads, x, y, gamx, gamy, gam1x, gam1y, nx, ny, xcheb, vcheb, t, r, found);
found = logical(found);

end

function [t, r, found] = getElemCoordinates16(x, y, gamx, gamy, gam1x, gam1y, nx, ny, xcheb, vcheb, nthreads)

npts = size(x, 1);
t = zeros(npts, 1);
r = zeros(npts, 1);
found = zeros(npts, 1);
mex_id_ = 'getElemCoordinates16(i int, i int, i double[], i double[], i double[], i double[], i double[], i double[], i double[], i double[], i double[], i double[], io double[], io double[], io int[])';
[t, r, found] = util.gateway(mex_id_, npts, nthreads, x, y, gamx, gamy, gam1x, gam1y, nx, ny, xcheb, vcheb, t, r, found);
found = logical(found);

end

function [t, r, found] = getElemCoordinates24(x, y, gamx, gamy, gam1x, gam1y, nx, ny, xcheb, vcheb, nthreads)

npts = size(x, 1);
t = zeros(npts, 1);
r = zeros(npts, 1);
found = zeros(npts, 1);
mex_id_ = 'getElemCoordinates24(i int, i int, i double[], i double[], i double[], i double[], i double[], i double[], i double[], i double[], i double[], i double[], io double[], io double[], io int[])';
[t, r, found] = util.gateway(mex_id_, npts, nthreads, x, y, gamx, gamy, gam1x, gam1y, nx, ny, xcheb, vcheb, t, r, found);
found = logical(found);

end

function [t, r, found] = getElemCoordinates30(x, y, gamx, gamy, gam1x, gam1y, nx, ny, xcheb, vcheb, nthreads)

npts = size(x, 1);
t = zeros(npts, 1);
r = zeros(npts, 1);
found = zeros(npts, 1);
mex_id_ = 'getElemCoordinates30(i int, i int, i double[], i double[], i double[], i double[], i double[], i double[], i double[], i double[], i double[], i double[], io double[], io double[], io int[])';
[t, r, found] = util.gateway(mex_id_, npts, nthreads, x, y, gamx, gamy, gam1x, gam1y, nx, ny, xcheb, vcheb, t, r, found);
found = logical(found);

end

function [t, r, found] = getElemCoordinates32(x, y, gamx, gamy, gam1x, gam1y, nx, ny, xcheb, vcheb, nthreads)

npts = size(x, 1);
t = zeros(npts, 1);
r = zeros(npts, 1);
found = zeros(npts, 1);
mex_id_ = 'getElemCoordinates32(i int, i int, i double[], i double[], i double[], i double[], i double[], i double[], i double[], i double[], i double[], i double[], io double[], io double[], io int[])';
[t, r, found] = util.gateway(mex_id_, npts, nthreads, x, y, gamx, gamy, gam1x, gam1y, nx, ny, xcheb, vcheb, t, r, found);
found = logical(found);

end
