function xClean = cleanDiffs(x)
    % Make x a vector
    x = reshape(x, numel(x), 1);
    
    % Remove NaNs and Infs
    xClean = x(~isnan(x) & ~isinf(x));
end