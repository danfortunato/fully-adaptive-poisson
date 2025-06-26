function func = piecewise(funs, ends)
%PIECEWISE   Piecewise function with vectorization.
%   F = PIECEWISE({F_1, ..., F_(n+1)}}, [X_1, ..., X_n]) returns a function
%   handle F defined as the piecewise function
%
%             / F_1(x),      for x < X_1,
%      F(x) = | F_k(x),      for X_k < x < X_{k+1} and k = 1, ..., n-1
%             \ F_(n+1)(x),  for x > X_n,
%
%   but constructed in a way so that F_k(x) is computed only if x satisfies
%   the piecewise criteria. This avoids the so-called "NaN problem" that
%   occurs when constructing piecewise functions in MATLAB using logical
%   indexing, but still allows the evaluation of F_k(x) to be vectorized.
%   Each function F_k should be given as a vectorized function handle.

%ends = [-inf ends inf];

    function out = ipiecewise(x)
        out = zeros(size(x));
        for k = 1:length(ends)-1
            idx = ( x >= ends(k) & x <= ends(k+1) );
            out(idx) = funs{k}(x(idx));
        end
    end

func = @ipiecewise;

end
