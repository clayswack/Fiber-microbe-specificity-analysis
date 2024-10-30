function rbNull = makeNullDataset(rb)

% Objectice:
% Randomly shuffle the log changes of relative abundance one row at a time
% This preserves the original distribution of changes to taxa for that
% sample, but assigns those changes to random taxa

rbNull = NaN(size(rb));

for i = 1:height(rb)

    logRBnow = rb(i,:); % logRB vector for one sample

    % Find indices of Inf, -Inf or NaN
    infIxs = isinf(logRBnow); % indices of infs
    nanIxs = isnan(logRBnow); % indices of nans

    % Find indices of zeros
    zeroPositions = find(logRBnow == 0); % positions of zeros (but not a logical array)
    zeroIxs = false(1, length(logRBnow)); % preallocate a logical array
    zeroIxs(zeroPositions) = 1; % weirdest line of code in my life

    % Find indices of everything else (the data we want to shuffle)
    badIxs = or(infIxs, nanIxs); % indices of nans or infs
    badIxs = or(badIxs, zeroIxs); % indices of nans, infs, or zeros
    keepIxs = find(~badIxs); % indices of good values
    logRBnowKeep = logRBnow(keepIxs); % good values (to shuffle)

    % Shuffle them
    logRBnowShuffled = logRBnowKeep(randperm(length(logRBnowKeep))); % good values shuffled
    rbNull(i,keepIxs) = logRBnowShuffled; % add the shuffled values to the null dataset

end

end