classdef AnnularPoissonSolver < AnnularModifiedHelmholtzSolver
% Spectrally accurate Poisson solver on annular domain
% 
%     Solves Lu = f in the annulus described by the Annular Geometry AG
%     Subject to the Robin boundary condition:
%     ia*u(ri) + ib*u_r(ri) = ig (boundary condition at the inner radius)
%     oa*u(ro) + ob*u_r(ro) = og (boundary condition at the outer radius)
% 
%     On instantionation, a preconditioner is formed with ia, ib, ua, ub
%         defining the boundary conditions
%     These can be changed at solvetime, but preconditioning may not work so well

    methods ( Access = public )

        function self = AnnularPoissonSolver(RAG, AAG, ia, ib, oa, ob)
            
            if ( nargin == 2 )
                ia = 1; ib = 0;
                oa = 1; ob = 0;
            end
            
            self = self@AnnularModifiedHelmholtzSolver(RAG, AAG, 0, ia, ib, oa, ob);

        end

        function out = solve(self, f, ig, og, tol)
            out = solve@AnnularModifiedHelmholtzSolver(self, -f, ig, og, tol);
        end

    end

end
