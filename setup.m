fprintf('Initializing path directories...')
root = fileparts([pwd, filesep]);
addpath(root)
addpath([root, '/geometry']);
addpath([root, '/kernels']);
addpath([root, '/tests']);
addpath([root, '/util']);
clear root
fprintf(' done.\n');
