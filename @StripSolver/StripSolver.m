classdef StripSolver
%STRIPSOLVER

    properties
        dom
    end

    methods

        function L = StripSolver(dom)
            assert(isstruct(dom), 'Invalid domain.');
            L.dom = dom;
        end

    end

end
