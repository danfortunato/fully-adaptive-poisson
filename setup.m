function setup

fprintf('Initializing path directories...')
root = fileparts(mfilename('fullpath'));
addpath(root)
addpath([root '/scripts'])

% External dependencies
addpath([root '/external/treefun'])
addpath([root '/external/inpoly'])
addpath([root '/external/kdtree/toolbox'])

if ( isempty(which('chebfun')) )
    error('Chebfun not found. Please install chebfun.');
end

if ( isempty(which('fmm2d')) )
    error('FMMLIB2D not found. Please install FMMLIB2D.');
end

fprintf(' done.\n');

end
