cheboppref.setDefaults('discretization', 'coeffs')

p = 6;
L = chebop(@(u) diff(u,2*(p+1)));
diffs = '';
for k = 1:p
    diffs = [diffs ' ; diff(u,' num2str(k) ')']; %#ok<AGROW>
end
L.lbc = eval(['@(u) [u'   diffs ']']);
L.rbc = L.lbc;
%L.rbc = eval(['@(u) [u-1' diffs ']']);
x = chebfun('x');
rhs = exp(-10*x.^2);
% rhs = 0;
rhs = 1;
v = L \ rhs;
v = v / max(abs(minandmax(v)));
plot(v)
shg

vals = zeros(p+1, 2);
for k = 1:p+1
    vals(k,:) = feval(diff(v,k-1), [-1 1]);
end
vals