function opts = setDefaults(opts, defaults)

fs = fieldnames(defaults);
for field = fs(:).'
    if ( ~isfield(opts, field{1}) )
        opts.(field{1}) = defaults.(field{1});
    end
end

end
