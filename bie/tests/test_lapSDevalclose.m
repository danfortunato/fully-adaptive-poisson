function pass = test_lapSDevalclose
% Solve 2D Dirichlet Laplace BVPs to test global close-evaluation formulae
% of Barnett which generalize Helsing/Ioakimidis.
% Currently DLP and SLP (lp='s','d'), int and ext cases (side='i','e').
% We use Ch.6 of LIE Kress book for formulation and uniqueness of the 4 BVPs,
% except we only consider u_infty=0 for ext Dir BVP (and don't modify DLP here).
% For ext Neu BVP only the zero-mean SLP is tested for close-eval.
% Needs: lapDevalclose.m, lapSevalclose.m, quadr.m, and whatever they needs.
% Alex Barnett for Shravan + Wu 10/18/13-11/4/13, based on 4/20/13.

tol = 1e-15;
pass = [];
N = 200;
S = Boundary.star(N);

for lp = ['s' 'd']
    for side = ['i' 'e']

        if ( side == 'i' )
            % Interior
            x = 0.2 + 0.3i;
        else
            % Exterior (note has O(1/|z|) decay so not general Dirichlet BVP solution)
            x = 1.2 + 1i;
        end

        % Initial test given density, if value and derivative agree in the
        % far field where standard trapezoid rule is accurate. Can test if
        % non-zero mean SLP is good:
        t = cell2mat(S.s);
        sigma = exp(sin(4*t+1)); % Generic smooth density with non-zero mean

        if ( lp == 'd' )
            A = DLPmatrix(S, x);
            u = kernels.laplace.dlp(S, 'target',    x,     ...
                                       'density',   sigma, ...
                                       'closeeval', true,  ...
                                       'side',      side);
        else
            A = SLPmatrix(S, x);
            u = kernels.laplace.slp(S, 'target',    x,     ...
                                       'density',   sigma, ...
                                       'closeeval', true,  ...
                                       'side',      side);
        end
        pass(end+1) = abs(u - A*sigma) < tol; %#ok<AGROW>
    end
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% end main %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function A = DLPmatrix(source, t) % double-layer kernel matrix & targ n-deriv
% t = target seg (x,nx cols), s = src seg, a = optional translation of src seg
% No jump included on self-interaction (ie principal-value integral)
N = source.N;
z = cell2mat(source.z);
w = cell2mat(source.w);
nx = cell2mat(source.normal); nx = nx(:,1)+nx(:,2)*1i;
cur = cell2mat(source.curvature);
M = numel(t);
d = repmat(t, [1 N]) - repmat(z.', [M 1]);    % C-# displacements mat
ny = repmat(nx.', [M 1]);      % identical rows given by src normals
A = (1/2/pi) * real(ny./d);      % complex form of dipole
if numel(z)==numel(t) && max(abs(z-t))<1e-14
    n = size(A,1);
    A(1:n+1:n^2) = -cur/4/pi;  % self? diagonal term for Laplace
end
A = A .* repmat(w(:)', [M 1]);

end

function A = SLPmatrix(source, t) % single-layer kernel matrix & targ n-deriv
% t = target seg (x,nx cols), s = src seg, a = optional translation of src seg
% No jump included on self-interaction of derivative (ie PV integral).
N = source.N;
z = cell2mat(source.z);
w = cell2mat(source.w);
M = numel(t);
d = repmat(t, [1 N]) - repmat(z.', [M 1]);    % C-# displacements mat
A = -(1/2/pi) * log(abs(d)) .* repmat(w(:)', [M 1]);  % infty for self-int
end
