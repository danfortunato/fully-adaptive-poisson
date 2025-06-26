function setup

fprintf('Initializing path directories...')
root = fileparts(mfilename('fullpath'));
addpath(root)
addpath([root '/scripts'])

% Add external dependencies
external = [root '/external'];
addpath([external '/chebfun'])
addpath([external '/treefun'])
addpath([external '/treefun/external/boxcode2d/matlab'])
addpath([external '/inpoly'])
addpath([external '/kdtree/toolbox'])
addpath([external '/fmmlib2d/matlab'])
addpath([external '/alignfigs'])
addpath([external '/export_fig'])

fprintf(' done.\n');

end
