function out = feval(vals, x)

N = size(vals,1);
pts = util.gauss(N);
R = util.interpmat(pts, x);
out = R * vals;

end
