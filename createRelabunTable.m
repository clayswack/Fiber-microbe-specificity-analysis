function [rbTable, meta] = createRelabunTable(raw_reads_table, meta_raw, minReadsForSample, minReadsForTaxa, minRelabunForTaxa)
%% Objective: start with raw reads table
% Filter out samples with low total reads
% Then calculate relative abundance

%% FILTER 1: Remove samples with total reads < minReadsForSample

    k = 1;
    for i = 1:height(raw_reads_table)
    
        % Get total reads for the sample (one row per sample)
        sampleTotalReads(i) = sum(raw_reads_table(i, :));
    
        % Now find indices of samples with few reads
        if sampleTotalReads(i) < minReadsForSample
            exclusionIxs(k) = i;
            k = k + 1;
        end
    
    end
    
    % Create copy of raw reads table and metadata
    rr_filt1 = raw_reads_table;
    meta_filt1 = meta_raw;
    
    % Now remove samples with with total reads < minReadsForSample
    if exist('exclusionIxs', 'var')
        rr_filt1(exclusionIxs,:) = [];
        meta_filt1(exclusionIxs,:) = [];
        N_samples_removed = numel(exclusionIxs)
    end
    
    
    %% FILTER 2: Remove taxa where reads < minReadsForTaxa
    
    rr_filt2 = rr_filt1;
    rr_filt2(rr_filt2 < minReadsForTaxa) = 0;
    
    
    %% First pass calculate relative abundance
    
    rb_filt2 = zeros(size(rr_filt2));
    for i = 1:height(rb_filt2)
    
        % Get total reads for the sample (one row per sample)
        sampleTotalReads = sum(rr_filt2(i, :));
    
        % Now calculate relative abundance
        rb_filt2(i,:) = rr_filt2(i, :)/sampleTotalReads;
    
    end
    
    
    %% FILTER 3: Remove taxa where relabun < minRelabunForTaxa
    
    rr_filt3 = rr_filt2;
    rr_filt3(rb_filt2 < minRelabunForTaxa) = 0;
    
    
    %% Second pass calculate relative abundance (after removing rare taxa)
    
    rb_filt3 = zeros(size(rr_filt3));
    for i = 1:height(rb_filt3)
    
        % Get total reads for the sample (one row per sample)
        sampleTotalReads = sum(rr_filt3(i, :));
    
        % Now calculate relative abundance
        rb_filt3(i,:) = rr_filt3(i, :)/sampleTotalReads;
    
    end
    
    
    %% Optional: sum check to make sure relabun is calculated properly
    
    % for i = 1:height(rb_filt2)
    %     sumCheckrbfilt2(i) = nansum(rb_filt2(i,:));
    %     sumCheckrbfilt3(i) = nansum(rb_filt3(i,:));
    % end
    % figure
    % plot(sumCheckrbfilt2) % should be = 1 +/- rounding error
    % hold on
    % plot(sumCheckrbfilt3)
    % leg = legend('filt2', 'filt3');

    %% Assign output
    rbTable = rb_filt3;
    meta = meta_filt1;
end