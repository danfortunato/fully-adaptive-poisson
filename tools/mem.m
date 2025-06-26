function varargout = mem()
%MEM   Get the current memory being used by MATLAB.

pid = feature('getpid');
[~, realmem] = unix(sprintf('ps -o rss -p %d | tail -1', pid));
[~, macmem]  = unix(sprintf('top -stats mem -pid %d -l 1 | tail -1', pid));

realmem = sanitize([realmem 'K']);
macmem  = sanitize(macmem);

if ( nargout == 0 )
    fprintf('Mac memory:  %.2fGB\n', macmem);
    fprintf('Real memory: %.2fGB\n', realmem);
elseif ( nargout == 1 )
    varargout{1} = macmem;
else
    varargout{1} = macmem;
    varargout{2} = realmem;
end

end

function mem = sanitize(mem)

mem = strip(mem);
if ( any(strcmp(mem(end), {'+', '-'})) )
    mem = mem(1:end-1);
end

% Convert to GB
switch upper(mem(end))
    case 'G', mem = str2double(mem(1:end-1));           % Already in GB
    case 'M', mem = str2double(mem(1:end-1)) / 1.024e3; % MB to GB
    case 'K', mem = str2double(mem(1:end-1)) / 1.024e6; % KB to GB
    otherwise mem = str2double(mem)          / 1.024e9; % Bytes to GB
end

end
