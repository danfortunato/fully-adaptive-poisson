function c = merge(a, b)
%MERGELR   Merge two patch objects.
%   C = MERGE(A, B) returns a patch C formed by merging the two patches A
%   and B. Typically A and B will be adjacent and any common edges will be
%   eliminated by enforcing continuity and continuity of the derivative
%   across the boundary.

% Parse inputs:
if ( nargin == 0 )
    c = [];
    return
elseif ( nargin == 1 )
    c = a;
    return
end

% Compute the indices of intersecting points in a and b.
[i1, i2, s1, s2, dom, edges] = intersect(a, b);

% Extract D2N maps:
D2Na = a.D2N; D2Nb = b.D2N;

% Compute new solution operator:
A = -( D2Na(s1,s1) + D2Nb(s2,s2) );
z = [ D2Na(s1,i1), D2Nb(s2,i2), D2Na(s1,end) + D2Nb(s2,end) ];
%                              |----------- rhs -----------|

% Store the decomposition for reuse in updateRHS().
dA = decomposition(A, 'lu');
S = dA \ z;

if ( isempty(s1) )
    keyboard
end

% Compute new D2N maps:
sparse_tol = 1e-14;
if ( issparse(D2Na) && issparse(D2Nb) && issparse(S) )

    nz = sum(abs(S) >= sparse_tol, 'all');
    if ( nz / numel(S) < 0.75 )
        S(abs(S) < sparse_tol) = 0;
        S = sparse(S);
    end

    n1 = numel(i1);
    n2 = numel(i2);
    D2Na_i1 = D2Na(i1,i1);
    D2Nb_i2 = D2Nb(i2,i2);
    nzmax = nnz(D2Na_i1) + nnz(D2Nb_i2) + n1 + n2;
    % TODO: Use block diag or somesuch to speed up indexing
    D2N = sparse(n1+n2, n1+n2+1, nzmax);
    D2N(1:n1,1:n1) = D2Na_i1;
    D2N(n1+(1:n2),n1+(1:n2)) = D2Nb_i2;
    D2N(1:n1,end) = D2Na(i1,end);
    D2N(n1+(1:n2),end) = D2Nb(i2,end);
    D2N = D2N + [ D2Na(i1,s1) ; D2Nb(i2,s2) ] * S;
else    
    Z12 = zeros(numel(i1), numel(i2));
    %                                 |--- rhs ----|
    D2N = [ D2Na(i1,i1),  Z12,         D2Na(i1,end) ;
            Z12.',        D2Nb(i2,i2), D2Nb(i2,end) ] ...
        + [ D2Na(i1,s1) ; D2Nb(i2,s2) ] * S;

    nz = sum(abs(D2N) >= sparse_tol, 'all');
    if ( nz / numel(D2N) < 0.75 )
        D2N(abs(D2N) < sparse_tol) = 0;
        D2N = sparse(D2N);
        S(abs(S) < sparse_tol) = 0;
        S = sparse(S);
    end
end

% if ( isa(a, 'StripSolver.Parent') && ...
%      isa(a.child1, 'StripSolver.Parent') && ...
%      isa(a.child1.child1, 'StripSolver.Parent') )
%     keyboard
% end

% Construct the new patch:
xy = [a.xy(i1,:) ; b.xy(i2,:)];
c = StripSolver.Parent(dom, S, D2N, dA, edges, xy, a, b, {i1, s1}, {i2, s2});

end
