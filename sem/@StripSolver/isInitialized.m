function out = isInitialized(S)
%ISINITIALIZED   Check to see if a STRIPSOLVER has been initialized.

out = ~isempty(S.patches) && size(S.patches{1}.S, 2) > 0;

end
