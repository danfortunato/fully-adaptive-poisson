function setup

warnstate = warning('off', 'MATLAB:dispatcher:nameConflict');

fprintf('Initializing path directories...')
root = pwd;
addpath(root)
addpath([root, '/bie'])
addpath([root, '/sem'])
addpath([root, '/tests'])
addpath([root, '/scripts'])
addpath('~/Research/Barnett/treefun')
addpath('~/Research/fmmlib2d/matlab')
addpath('~/Research/kdtree/toolbox')
addpath('~/Research/Barnett/inpoly')

% Add linequad for close panel evaluation
linequadroot = '~/Research/Barnett/linequad/matlab';
addpath([linequadroot '/bin'])
addpath(genpath([linequadroot '/src']))
addpath([linequadroot '/external/legtools'])

fprintf(' done.\n');

warning(warnstate)

end
