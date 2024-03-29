function pass = test_boundary()

tol = 1e-8;
pass = [];

%% Construct a periodic boundary from given nodes
N = 100;
A = Boundary.star(N);
B = Boundary(A.x{1});
pass(end+1) = ( norm(A.x{1}         - B.x{1})         < tol ) && ...
              ( norm(A.dx{1}        - B.dx{1})        < tol ) && ...
              ( norm(A.dxx{1}       - B.dxx{1})       < tol ) && ...
              ( norm(A.z{1}         - B.z{1})         < tol ) && ...
              ( norm(A.dz{1}        - B.dz{1})        < tol ) && ...
              ( norm(A.dzz{1}       - B.dzz{1})       < tol ) && ...
              ( norm(A.s{1}         - B.s{1})         < tol ) && ...
              ( norm(A.w{1}         - B.w{1})         < tol ) && ...
              ( norm(A.cw{1}        - B.cw{1})        < tol ) && ...
              ( norm(A.speed{1}     - B.speed{1})     < tol ) && ...
              ( norm(A.normal{1}    - B.normal{1})    < tol ) && ...
              ( norm(A.curvature{1} - B.curvature{1}) < tol );

%% Construct a panelized boundary from given nodes
N  = 10;
np = 70;
A = Boundary.star(N, 'quadrature', 'panel', 'panels', np);
B = Boundary(A.x);
pass(end+1) = ( norm(cell2mat(A.x)         - cell2mat(B.x))         < tol ) && ...
              ( norm(cell2mat(A.dx)        - cell2mat(B.dx))        < tol ) && ...
              ( norm(cell2mat(A.dxx)       - cell2mat(B.dxx))       < tol ) && ...
              ( norm(cell2mat(A.z)         - cell2mat(B.z))         < tol ) && ...
              ( norm(cell2mat(A.dz)        - cell2mat(B.dz))        < tol ) && ...
              ( norm(cell2mat(A.dzz)       - cell2mat(B.dzz))       < tol ) && ...
              ( norm(cell2mat(A.s)         - cell2mat(B.s))         < tol ) && ...
              ( norm(cell2mat(A.w)         - cell2mat(B.w))         < tol ) && ...
              ( norm(cell2mat(A.cw)        - cell2mat(B.cw))        < tol ) && ...
              ( norm(cell2mat(A.speed)     - cell2mat(B.speed))     < tol ) && ...
              ( norm(cell2mat(A.normal)    - cell2mat(B.normal))    < tol ) && ...
              ( norm(cell2mat(A.curvature) - cell2mat(B.curvature)) < tol );

end
