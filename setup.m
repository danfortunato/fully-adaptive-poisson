function setup

warnstate = warning('off', 'MATLAB:dispatcher:nameConflict');

fprintf('Initializing path directories...')
root = fileparts(mfilename('fullpath'));
addpath(root)
addpath([root '/bie'])
addpath([root '/tests'])
addpath([root '/scripts'])

% External dependencies
addpath([root '/external/treefun'])
addpath([root '/external/inpoly'])
addpath([root '/external/kdtree/toolbox'])
addpath([root '/external/BIE2D/panels'])
addpath([root '/external/BIE2D/kernels'])
addpath([root '/external/BIE2D/utils'])
addpath([root '/external/linequad/matlab/bin'])
addpath([root '/external/linequad/matlab/external/legtools'])
addpath(genpath([root '/external/linequad/matlab/src']))

if ( isempty(which('chebfun')) )
    error('Chebfun not found. Please install chebfun.');
end

if ( isempty(which('fmm2d')) )
    error('FMMLIB2D not found. Please install FMMLIB2D.');
end

fprintf(' done.\n');

warning(warnstate)

end
