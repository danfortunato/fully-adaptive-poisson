function mem = memtoc(start)
%MEMTOC   Read the memory timer.

if ( nargin == 0 )
    file = '/tmp/memtic';
    fid = fopen('/tmp/memtic_pid');
    pid = cell2mat(textscan(fid, '%f'));
    fclose(fid);
else
    file = start.file;
    pid = start.pid;
end

fid = fopen(file);
mems = textscan(fid, '%f');
mems = cell2mat(mems);
fclose(fid);

unix(sprintf('kill %d', pid));
unix(sprintf('rm -f %s /tmp/memtic_pid', file));

mem = max(mems) - mems(1);
mem = mem / 1.024e6;

end
