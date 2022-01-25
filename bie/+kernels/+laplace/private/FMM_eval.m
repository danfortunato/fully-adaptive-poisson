function u = FMM_eval(source, target, varargin)
%FMM_eval   Evaluate Laplace layer potentials at target points
%using the fast multipole method.

p = inputParser;
addRequired(p, 'source', @(s) isa(s,'Boundary'));
addRequired(p, 'target', @isnumeric);
addParameter(p, 'charge',  [],  @isnumeric);
addParameter(p, 'dipstr',  [],  @isnumeric);
addParameter(p, 'tol',     eps, @isscalar);
parse(p, source, target, varargin{:});
charge  = p.Results.charge;
dipstr  = p.Results.dipstr;
tol     = p.Results.tol;

% 1 => tolerance =.5d-3  =>  3 digits
% 2 => tolerance =.5d-6  =>  6 digits
% 3 => tolerance =.5d-9  =>  9 digits
% 4 => tolerance =.5d-12 => 12 digits
% 5 => tolerance =.5d-15 => 15 digits
if tol > .5d-3
    iprec = 1;
elseif tol > .5d-6
    iprec = 2;
elseif tol > .5d-9
    iprec = 3;
elseif tol > .5d-12
    iprec = 4;
else 
    iprec = 5;
end

nsource = source.N * source.np;
x = cell2mat(source.x).';
y = target.';
w = cell2mat(source.w);
ifcharge = ~isempty(charge);
if ( ifcharge )
    charge = charge .* w;
end
ifdipole = ~isempty(dipstr);
if ( ifdipole )
    dipstr = dipstr .* w;
end
dipvec = cell2mat(source.normal).';
ntarget = size(target,1);
ifpot  = 0;
ifgrad = 0;
ifhess = 0;
ifpottarg  = 1;
ifgradtarg = 0;
ifhesstarg = 0;
U = rfmm2dpart( iprec, nsource, x, ifcharge, charge, ifdipole, dipstr, dipvec, ifpot, ifgrad, ifhess, ...
                       ntarget, y, ifpottarg, ifgradtarg, ifhesstarg );
u = -reshape(U.pottarg, ntarget, 1) / (2*pi);

end
