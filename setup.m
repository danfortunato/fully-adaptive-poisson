function setup

warnstate = warning('off', 'MATLAB:dispatcher:nameConflict');

fprintf('Initializing path directories...')
root = pwd;
addpath(root)
addpath([root '/bie'])
addpath([root '/sem'])
addpath([root '/tests'])
addpath([root '/scripts'])

% External dependencies
addpath([root '/external/inpoly'])
addpath([root '/external/kdtree/toolbox'])
addpath([root '/external/BIE2D/panels'])
addpath([root '/external/BIE2D/kernels'])
addpath([root '/external/BIE2D/utils'])
addpath('~/Research/Barnett/treefun')
addpath('~/Research/fmmlib2d/matlab')

% Add linequad for close panel evaluation
linequadroot = '~/Research/Barnett/linequad/matlab';
addpath([linequadroot '/bin'])
addpath(genpath([linequadroot '/src']))
addpath([linequadroot '/external/legtools'])

fprintf(' done.\n');

warning(warnstate)

end
