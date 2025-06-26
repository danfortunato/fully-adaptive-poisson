function varargout = memtic
%MEMTIC   Start a memory timer.

if ( nargout == 0 )
    file = '/tmp/memtic';
else
    file = sprintf('/tmp/memtic_%x', randi([0 intmax('uint32')], 'uint32'));
end

if ( isfile(file) )
    unix(['rm -f ' file]);
end

pid = feature('getpid');
dt = 1;

cmd = sprintf(['rm -f file; ' ...
               'while true; do ' ...
                 'ps -o rss -p %d | tail -1 >> %s; ' ...
                 'sleep %.1f; ' ...
               'done'], pid, file, dt);
[~, out] = unix([cmd ' & jobs -p']);
args = textscan(out, '%s %s %d');
cmdpid = args{3};
unix(sprintf('rm -f memtic_pid; echo %d > /tmp/memtic_pid', cmdpid));

if ( nargout == 1 )
    varargout{1} = struct('file', file, 'pid', cmdpid);
end

end
