function [madClean, stdadClean, taxaNamesClean] = clean_MAD_STDAD_andUpdateTaxaNames(mad, stdad, taxaNames)
    % Make x a vector
    mad = reshape(mad, numel(mad), 1);
    
    % Remove NaNs and Infs
    madClean = mad(~isnan(mad) & ~isinf(mad));
    stdadClean = stdad(~isnan(mad) & ~isinf(mad));

    taxaNamesClean = taxaNames(~isnan(mad) & ~isinf(mad));

end