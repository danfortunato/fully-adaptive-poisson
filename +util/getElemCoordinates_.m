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

function [t, r, found] = getElemCoordinates4(x, y, gamx, gamy, gam1x, gam1y, nx, ny, xcheb, vcheb, nthreads, tol)

    npts = size(x, 1);
    t = zeros(npts, 1);
    r = zeros(npts, 1);
    found = zeros(npts, 1);
    mex_id_ = 'getElemCoordinates4(i int, i int, i double, i double[], i double[], i double[], i double[], i double[], i double[], i double[], i double[], i double[], i double[], io double[], io double[], io int[])';
[t, r, found] = util.getElemCoordinates_mex(mex_id_, npts, nthreads, tol, x, y, gamx, gamy, gam1x, gam1y, nx, ny, xcheb, vcheb, t, r, found);
    found = logical(found);
    
end

function [t, r, found] = getElemCoordinates6(x, y, gamx, gamy, gam1x, gam1y, nx, ny, xcheb, vcheb, nthreads, tol)

    npts = size(x, 1);
    t = zeros(npts, 1);
    r = zeros(npts, 1);
    found = zeros(npts, 1);
    mex_id_ = 'getElemCoordinates6(i int, i int, i double, i double[], i double[], i double[], i double[], i double[], i double[], i double[], i double[], i double[], i double[], io double[], io double[], io int[])';
[t, r, found] = util.getElemCoordinates_mex(mex_id_, npts, nthreads, tol, x, y, gamx, gamy, gam1x, gam1y, nx, ny, xcheb, vcheb, t, r, found);
    found = logical(found);
    
end

function [t, r, found] = getElemCoordinates8(x, y, gamx, gamy, gam1x, gam1y, nx, ny, xcheb, vcheb, nthreads, tol)

    npts = size(x, 1);
    t = zeros(npts, 1);
    r = zeros(npts, 1);
    found = zeros(npts, 1);
    mex_id_ = 'getElemCoordinates8(i int, i int, i double, i double[], i double[], i double[], i double[], i double[], i double[], i double[], i double[], i double[], i double[], io double[], io double[], io int[])';
[t, r, found] = util.getElemCoordinates_mex(mex_id_, npts, nthreads, tol, x, y, gamx, gamy, gam1x, gam1y, nx, ny, xcheb, vcheb, t, r, found);
    found = logical(found);
    
end

function [t, r, found] = getElemCoordinates9(x, y, gamx, gamy, gam1x, gam1y, nx, ny, xcheb, vcheb, nthreads, tol)

    npts = size(x, 1);
    t = zeros(npts, 1);
    r = zeros(npts, 1);
    found = zeros(npts, 1);
    mex_id_ = 'getElemCoordinates9(i int, i int, i double, i double[], i double[], i double[], i double[], i double[], i double[], i double[], i double[], i double[], i double[], io double[], io double[], io int[])';
[t, r, found] = util.getElemCoordinates_mex(mex_id_, npts, nthreads, tol, x, y, gamx, gamy, gam1x, gam1y, nx, ny, xcheb, vcheb, t, r, found);
    found = logical(found);
    
end

function [t, r, found] = getElemCoordinates10(x, y, gamx, gamy, gam1x, gam1y, nx, ny, xcheb, vcheb, nthreads, tol)

    npts = size(x, 1);
    t = zeros(npts, 1);
    r = zeros(npts, 1);
    found = zeros(npts, 1);
    mex_id_ = 'getElemCoordinates10(i int, i int, i double, i double[], i double[], i double[], i double[], i double[], i double[], i double[], i double[], i double[], i double[], io double[], io double[], io int[])';
[t, r, found] = util.getElemCoordinates_mex(mex_id_, npts, nthreads, tol, x, y, gamx, gamy, gam1x, gam1y, nx, ny, xcheb, vcheb, t, r, found);
    found = logical(found);
    
end

function [t, r, found] = getElemCoordinates11(x, y, gamx, gamy, gam1x, gam1y, nx, ny, xcheb, vcheb, nthreads, tol)

    npts = size(x, 1);
    t = zeros(npts, 1);
    r = zeros(npts, 1);
    found = zeros(npts, 1);
    mex_id_ = 'getElemCoordinates11(i int, i int, i double, i double[], i double[], i double[], i double[], i double[], i double[], i double[], i double[], i double[], i double[], io double[], io double[], io int[])';
[t, r, found] = util.getElemCoordinates_mex(mex_id_, npts, nthreads, tol, x, y, gamx, gamy, gam1x, gam1y, nx, ny, xcheb, vcheb, t, r, found);
    found = logical(found);
    
end

function [t, r, found] = getElemCoordinates12(x, y, gamx, gamy, gam1x, gam1y, nx, ny, xcheb, vcheb, nthreads, tol)

    npts = size(x, 1);
    t = zeros(npts, 1);
    r = zeros(npts, 1);
    found = zeros(npts, 1);
    mex_id_ = 'getElemCoordinates12(i int, i int, i double, i double[], i double[], i double[], i double[], i double[], i double[], i double[], i double[], i double[], i double[], io double[], io double[], io int[])';
[t, r, found] = util.getElemCoordinates_mex(mex_id_, npts, nthreads, tol, x, y, gamx, gamy, gam1x, gam1y, nx, ny, xcheb, vcheb, t, r, found);
    found = logical(found);

end

function [t, r, found] = getElemCoordinates13(x, y, gamx, gamy, gam1x, gam1y, nx, ny, xcheb, vcheb, nthreads, tol)

    npts = size(x, 1);
    t = zeros(npts, 1);
    r = zeros(npts, 1);
    found = zeros(npts, 1);
    mex_id_ = 'getElemCoordinates13(i int, i int, i double, i double[], i double[], i double[], i double[], i double[], i double[], i double[], i double[], i double[], i double[], io double[], io double[], io int[])';
[t, r, found] = util.getElemCoordinates_mex(mex_id_, npts, nthreads, tol, x, y, gamx, gamy, gam1x, gam1y, nx, ny, xcheb, vcheb, t, r, found);
    found = logical(found);

end

function [t, r, found] = getElemCoordinates14(x, y, gamx, gamy, gam1x, gam1y, nx, ny, xcheb, vcheb, nthreads, tol)

    npts = size(x, 1);
    t = zeros(npts, 1);
    r = zeros(npts, 1);
    found = zeros(npts, 1);
    mex_id_ = 'getElemCoordinates14(i int, i int, i double, i double[], i double[], i double[], i double[], i double[], i double[], i double[], i double[], i double[], i double[], io double[], io double[], io int[])';
[t, r, found] = util.getElemCoordinates_mex(mex_id_, npts, nthreads, tol, x, y, gamx, gamy, gam1x, gam1y, nx, ny, xcheb, vcheb, t, r, found);
    found = logical(found);
    
end

function [t, r, found] = getElemCoordinates15(x, y, gamx, gamy, gam1x, gam1y, nx, ny, xcheb, vcheb, nthreads, tol)

    npts = size(x, 1);
    t = zeros(npts, 1);
    r = zeros(npts, 1);
    found = zeros(npts, 1);
    mex_id_ = 'getElemCoordinates15(i int, i int, i double, i double[], i double[], i double[], i double[], i double[], i double[], i double[], i double[], i double[], i double[], io double[], io double[], io int[])';
[t, r, found] = util.getElemCoordinates_mex(mex_id_, npts, nthreads, tol, x, y, gamx, gamy, gam1x, gam1y, nx, ny, xcheb, vcheb, t, r, found);
    found = logical(found);
    
end

function [t, r, found] = getElemCoordinates16(x, y, gamx, gamy, gam1x, gam1y, nx, ny, xcheb, vcheb, nthreads, tol)

    npts = size(x, 1);
    t = zeros(npts, 1);
    r = zeros(npts, 1);
    found = zeros(npts, 1);
    mex_id_ = 'getElemCoordinates16(i int, i int, i double, i double[], i double[], i double[], i double[], i double[], i double[], i double[], i double[], i double[], i double[], io double[], io double[], io int[])';
[t, r, found] = util.getElemCoordinates_mex(mex_id_, npts, nthreads, tol, x, y, gamx, gamy, gam1x, gam1y, nx, ny, xcheb, vcheb, t, r, found);
    found = logical(found);

end

function [t, r, found] = getElemCoordinates18(x, y, gamx, gamy, gam1x, gam1y, nx, ny, xcheb, vcheb, nthreads, tol)

    npts = size(x, 1);
    t = zeros(npts, 1);
    r = zeros(npts, 1);
    found = zeros(npts, 1);
    mex_id_ = 'getElemCoordinates18(i int, i int, i double, i double[], i double[], i double[], i double[], i double[], i double[], i double[], i double[], i double[], i double[], io double[], io double[], io int[])';
[t, r, found] = util.getElemCoordinates_mex(mex_id_, npts, nthreads, tol, x, y, gamx, gamy, gam1x, gam1y, nx, ny, xcheb, vcheb, t, r, found);
    found = logical(found);
    
    end

function [t, r, found] = getElemCoordinates20(x, y, gamx, gamy, gam1x, gam1y, nx, ny, xcheb, vcheb, nthreads, tol)

    npts = size(x, 1);
    t = zeros(npts, 1);
    r = zeros(npts, 1);
    found = zeros(npts, 1);
    mex_id_ = 'getElemCoordinates20(i int, i int, i double, i double[], i double[], i double[], i double[], i double[], i double[], i double[], i double[], i double[], i double[], io double[], io double[], io int[])';
[t, r, found] = util.getElemCoordinates_mex(mex_id_, npts, nthreads, tol, x, y, gamx, gamy, gam1x, gam1y, nx, ny, xcheb, vcheb, t, r, found);
    found = logical(found);

end

function [t, r, found] = getElemCoordinates21(x, y, gamx, gamy, gam1x, gam1y, nx, ny, xcheb, vcheb, nthreads, tol)

    npts = size(x, 1);
    t = zeros(npts, 1);
    r = zeros(npts, 1);
    found = zeros(npts, 1);
    mex_id_ = 'getElemCoordinates21(i int, i int, i double, i double[], i double[], i double[], i double[], i double[], i double[], i double[], i double[], i double[], i double[], io double[], io double[], io int[])';
[t, r, found] = util.getElemCoordinates_mex(mex_id_, npts, nthreads, tol, x, y, gamx, gamy, gam1x, gam1y, nx, ny, xcheb, vcheb, t, r, found);
    found = logical(found);
    
    end

function [t, r, found] = getElemCoordinates24(x, y, gamx, gamy, gam1x, gam1y, nx, ny, xcheb, vcheb, nthreads, tol)

    npts = size(x, 1);
    t = zeros(npts, 1);
    r = zeros(npts, 1);
    found = zeros(npts, 1);
    mex_id_ = 'getElemCoordinates24(i int, i int, i double, i double[], i double[], i double[], i double[], i double[], i double[], i double[], i double[], i double[], i double[], io double[], io double[], io int[])';
[t, r, found] = util.getElemCoordinates_mex(mex_id_, npts, nthreads, tol, x, y, gamx, gamy, gam1x, gam1y, nx, ny, xcheb, vcheb, t, r, found);
    found = logical(found);

end

function [t, r, found] = getElemCoordinates30(x, y, gamx, gamy, gam1x, gam1y, nx, ny, xcheb, vcheb, nthreads, tol)

    npts = size(x, 1);
    t = zeros(npts, 1);
    r = zeros(npts, 1);
    found = zeros(npts, 1);
    mex_id_ = 'getElemCoordinates30(i int, i int, i double, i double[], i double[], i double[], i double[], i double[], i double[], i double[], i double[], i double[], i double[], io double[], io double[], io int[])';
[t, r, found] = util.getElemCoordinates_mex(mex_id_, npts, nthreads, tol, x, y, gamx, gamy, gam1x, gam1y, nx, ny, xcheb, vcheb, t, r, found);
    found = logical(found);

end

function [t, r, found] = getElemCoordinates32(x, y, gamx, gamy, gam1x, gam1y, nx, ny, xcheb, vcheb, nthreads, tol)

    npts = size(x, 1);
    t = zeros(npts, 1);
    r = zeros(npts, 1);
    found = zeros(npts, 1);
    mex_id_ = 'getElemCoordinates32(i int, i int, i double, i double[], i double[], i double[], i double[], i double[], i double[], i double[], i double[], i double[], i double[], io double[], io double[], io int[])';
[t, r, found] = util.getElemCoordinates_mex(mex_id_, npts, nthreads, tol, x, y, gamx, gamy, gam1x, gam1y, nx, ny, xcheb, vcheb, t, r, found);
    found = logical(found);

end
