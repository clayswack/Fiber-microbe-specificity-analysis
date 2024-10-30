%% cantu_jungles_hamaker_pipeline5.m


%% Objective:

% Calculate delta, fold, and log change of relative abundance
% Calculate alpha diversity of the changes
% Calculate the mean among donors (MAD), and coefficient of variance among
% donors (COVAD) for the delta, fold, and log changes

% Pipeline 5 has ONLY calculations involving relative abundance
% Pipeline 5 has NO subtracting the blank

%% Setup

clear
close all
clc
tic


%% User inputs

rawDataName = 'cantu_jungles_hamaker_raw'
saveFilePrefix = 'cantu_jungles_hamaker_';
N_metadata_cols = 3; % number of colums on the left side with metadata

minReadsForSample = 0; % remove samples with fewer total raw reads than this
minReadsForTaxa = 0; % remove taxa with < this number of reads
minRelabunForTaxa = 0; % old default = 0.001; % (0.005 = 0.5%) remove taxa with < this relative abundance
minDonorsToCalcSTDAD = 3; % min number of donors to calculate the stdev of a taxa between donors
absFlag = false; % express delta, fold, and log change as absolute vals


%% Get raw data

opts = detectImportOptions(strcat('D:\Specificity\Raw data\', rawDataName, '.csv'));
varTypes = opts.VariableTypes;
varTypes(1:N_metadata_cols) = repmat({'categorical'},1,N_metadata_cols);
varTypes(N_metadata_cols+1:end) = repmat({'double'},1,length(varTypes) - N_metadata_cols);
opts.VariableTypes = varTypes;
rr_raw = readtable(strcat('D:\Specificity\Raw data\', rawDataName, '.csv'), opts); % raw reads table including metadata

% Change some variable types
meta_raw = rr_raw(:,1:N_metadata_cols); % metadata
rr_raw_double = rr_raw{:, N_metadata_cols+1:end}; % raw reads without metadata

% Remove blanks we don't need them
blankRowIxs = find(rr_raw.fiber_type == 'Blank');
rr_raw(blankRowIxs,:) = [];



%% Calculate relative abundance from the raw reads table

% Also, filter out samples that had total reads < minReadsForSample
[rbTableFilt, metaFilt] = createRelabunTable(rr_raw_double, meta_raw,  minReadsForSample, minReadsForTaxa, minRelabunForTaxa);


% Make output table holding filtered relative abundance

% Now convert relative abundance to table with appropriate taxa names
taxaNames = rr_raw.Properties.VariableNames;
taxaNames = taxaNames(N_metadata_cols+1:end);
rbt = array2table(rbTableFilt, "VariableNames", taxaNames);

% Matlab modifies the taxa names to make valid table variable names
% Here, store the original, full names
fullTaxaNames = rr_raw.Properties.VariableDescriptions(N_metadata_cols+1:end);

% Append the metadata
% Filter 1 is the only step that gets rid of whole samples (if sampleTotalReads < minReadsForSample),
% so it is the only step that can change the metadata
rbTable = horzcat(metaFilt, rbt);
rbDouble = rbTable{:, N_metadata_cols+1:end};

% % Optional sum check relabun again: should be = 1 +/- rounding error
% for i = 1:height(rbTable)
%     sumCheck(i) = nansum(rbTable{i, N_metadata_cols+1:end});
% end
% figure
% plot(sumCheck)

%% Separate relabun table into final and initial samples

% RB Final
rbFinal = rbTable(rbTable.time == "after" & rbTable.fiber_type ~= 'Blank', :);
rbFinalDouble = rbFinal{:, N_metadata_cols+1:end};
finalMeta = rbFinal(:, 1:N_metadata_cols);

% RB Initial
rbInitial = rbTable(rbTable.time == "before",:);
rbInitialDouble = rbInitial{:, N_metadata_cols+1:end};
initialMeta = rbInitial(:, 1:N_metadata_cols);


%% Average out technical replicates from the same donor

% Get unique fibers and donors from the final (after treatment) samples
unqFibers = unique(finalMeta.fiber_type)
nFibers = numel(unqFibers);
unqDonors = categorical(unique(double(finalMeta.donor)))
nDonors = numel(unqDonors);

[rbFinal, rbInitial] = avgTechReps(rbFinalDouble, rbInitialDouble, initialMeta, finalMeta, taxaNames, unqFibers, nFibers, unqDonors, nDonors);


%% Calculate delta, fold, and log change in relabun

% The metadata of these is different
rbFinalNoTechReps = rbFinal{:, N_metadata_cols+1:end};
rbInitialNoTechReps = rbInitial{:, N_metadata_cols:end};

nTaxa = length(taxaNames)
nSamples = height(rbFinal) % N samples (after the treatment)

% Preallocate
deltaRB = zeros(size(rbFinalNoTechReps));
foldRB = zeros(size(rbFinalNoTechReps));

% Run outer loop once per sample (row) in the final (after fermentation) results table
for i = 1:nSamples
    donorNow = rbFinal.donor(i); % donor for this row (one row per sample)
    
    % Run inner loop once for each taxa in the ith sample
    for k = 1:nTaxa

        % Get all changes to this this taxa across all donors 
        donor_init_reps_rb = rbInitialNoTechReps(rbInitial.donor == donorNow, k);
        avg_rb_init_donor = mean(donor_init_reps_rb);

        % Run some tests, these should never print anything
        if isempty(avg_rb_init_donor)
            disp('empty')
        elseif isnan(avg_rb_init_donor)
            disp('NaN')
        end

        % If the donor's initial RB for that taxa was zero but the final
        % was nonzero, we can compute delta RB but not fold or log RB
        if avg_rb_init_donor == 0 && rbFinalNoTechReps(i,k) ~= 0
            deltaRB(i,k) = rbFinalNoTechReps(i,k) - avg_rb_init_donor;
            foldRB(i,k) = NaN;

        % If the donor's initial RB for that taxa was zero and the final
        % was also zero, we cannot compute anything
        elseif avg_rb_init_donor == 0 && rbFinalNoTechReps(i,k) == 0
            deltaRB(i,k) = NaN;
            foldRB(i,k) = NaN;

        % If the initial and final RB are nonzero, we can compute
        % everything
        else
            deltaRB(i,k) = rbFinalNoTechReps(i,k) - avg_rb_init_donor;
            foldRB(i,k) = rbFinalNoTechReps(i,k)/avg_rb_init_donor;
        end

    end

end

logRB = log10(foldRB);


%% Create null dataset (absence of fiber/microbe alignment)

% Randomly shuffle the log changes of relative abundance one row at a time
% This preserves the original distribution of changes to taxa for that
% sample, but assigns those changes to random taxa

logRBNull = makeNullDataset(logRB);


%% Optional: save delta, fold, and log change in relabun tables as spreadsheets

    %
% Convert deltaRB, foldRB, and logRB to tables with appropriate taxa names
samplesMeta = table;
samplesMeta.donor = rbFinal.donor;
samplesMeta.fiber_type = rbFinal.fiber_type;
samplesMeta.time = rbFinal.time;
deltaRBtable = horzcat(samplesMeta, array2table(deltaRB,"VariableNames",taxaNames));
foldRBtable = horzcat(samplesMeta, array2table(foldRB,"VariableNames",taxaNames));
logRBtable = horzcat(samplesMeta, array2table(logRB,"VariableNames",taxaNames));
logRBtableNull = horzcat(samplesMeta, array2table(logRBNull,"VariableNames",taxaNames));

filename = strcat(saveFilePrefix, 'delta_relabun.csv');
writetable(deltaRBtable,filename);

filename = strcat(saveFilePrefix, 'fold_relabun.csv');
writetable(foldRBtable,filename);

filename = strcat(saveFilePrefix, 'log_relabun.csv');
writetable(logRBtable,filename);

filename = strcat(saveFilePrefix, 'log_relabun_null.csv');
writetable(logRBtableNull,filename);


%% Calculate alpha diversity of each change vector
   
% Preallocate
shan_rawRB = zeros(height(rbFinal),1);
shan_deltaRB = zeros(height(rbFinal),1);
shan_foldRB = zeros(height(rbFinal),1);
shan_logRB = zeros(height(rbFinal),1);

rich_rawRB = zeros(height(rbFinal),1);
rich_deltaRB = zeros(height(rbFinal),1);
rich_foldRB = zeros(height(rbFinal),1);
rich_logRB = zeros(height(rbFinal),1);

simp_rawRB = zeros(height(rbFinal),1);
simp_deltaRB = zeros(height(rbFinal),1);
simp_foldRB = zeros(height(rbFinal),1);
simp_logRB = zeros(height(rbFinal),1);

even_rawRB = zeros(height(rbFinal),1);
even_deltaRB = zeros(height(rbFinal),1);
even_foldRB = zeros(height(rbFinal),1);
even_logRB = zeros(height(rbFinal),1);

for i = 1:height(rbFinal)

    % Raw final sample relabun (NOT the change in relabun)
    wholeRow = rbFinal{i,N_metadata_cols+1:end};
    rowClean = wholeRow(~isnan(wholeRow) & ~isinf(wholeRow) & wholeRow~=0);
    Nelements = numel(rowClean);
    probs = histcounts(rowClean, Nelements, 'Normalization','probability');
    shan_rawRB(i) = calc_shannon(probs(probs~=0));
    rich_rawRB(i) = Nelements;
    simp_rawRB(i) = 1/sum(rowClean.^2);
    even_rawRB(i) = shan_rawRB(i)/log(Nelements);

    % Delta relabun
    wholeRow = deltaRB(i,:);
    rowClean = wholeRow(~isnan(wholeRow) & ~isinf(wholeRow));
    Nelements = numel(rowClean);
    probs = histcounts(rowClean, Nelements, 'Normalization','probability');
    shan_deltaRB(i) = calc_shannon(probs(probs~=0));
    rich_deltaRB(i) = Nelements;
    simp_deltaRB(i) = 1/sum(rowClean.^2);
    even_deltaRB(i) = shan_deltaRB(i)/log(Nelements);

    % Fold relabun
    wholeRow = foldRB(i,:);
    rowClean = wholeRow(~isnan(wholeRow) & ~isinf(wholeRow));
    Nelements = numel(rowClean);
    probs = histcounts(rowClean, Nelements, 'Normalization','probability');
    shan_foldRB(i) = calc_shannon(probs(probs~=0));
    rich_foldRB(i) = Nelements;
    simp_foldRB(i) = 1/sum(rowClean.^2);
    even_foldRB(i) = shan_foldRB(i)/log(Nelements);

    % Log relabun
    wholeRow = logRB(i,:);
    rowClean = wholeRow(~isnan(wholeRow) & ~isinf(wholeRow));
    Nelements = numel(rowClean);
    if Nelements~= 0
        probs = histcounts(rowClean, Nelements, 'Normalization','probability');
        shan_logRB(i) = calc_shannon(probs(probs~=0));
        rich_logRB(i) = Nelements;
        simp_logRB(i) = 1/sum(rowClean.^2);
        even_logRB(i) = shan_logRB(i)/log(Nelements);
    else
        shan_logRB(i) = NaN;
        rich_logRB(i) = NaN;
        simp_logRB(i) = NaN;
        even_logRB(i) = NaN;
    end

end

% Make a few plots to check for anything weird
figure
subplot(2,2,1)
plot(shan_rawRB)
hold on
plot(shan_deltaRB)
plot(shan_foldRB)
plot(shan_logRB)
title('Shannon')

subplot(2,2,2)
plot(rich_rawRB)
hold on
plot(rich_deltaRB)
plot(rich_foldRB)
plot(rich_logRB)
title('Richness')

subplot(2,2,3)
plot(simp_rawRB)
hold on
plot(simp_deltaRB)
plot(simp_foldRB)
plot(simp_logRB)
title('Simpson')

subplot(2,2,4)
plot(even_rawRB)
hold on
plot(even_deltaRB)
plot(even_foldRB)
plot(even_logRB)
title('Evenness')

alphaDivVals = [shan_rawRB, shan_deltaRB, shan_foldRB, shan_logRB,...
                rich_rawRB, rich_deltaRB, rich_foldRB, rich_logRB,...
                simp_rawRB, simp_deltaRB, simp_foldRB, simp_logRB,...
                even_rawRB, even_deltaRB, even_foldRB, even_logRB];

colNames = {'shan_rawRB', 'shan_deltaRB', 'shan_foldRB', 'shan_logRB',...
            'rich_rawRB', 'rich_deltaRB', 'rich_foldRB', 'rich_logRB',...
            'simp_rawRB', 'simp_deltaRB', 'simp_foldRB', 'simp_logRB',...
            'even_rawRB', 'even_deltaRB', 'even_foldRB', 'even_logRB'};

alphaDiv = horzcat(samplesMeta, array2table(alphaDivVals,"VariableNames",colNames));

%  save alpha diversity output table
filename = strcat(saveFilePrefix, 'alphaDiv.csv');
writetable(alphaDiv,filename)


%% Find the variance between donors for each taxa on a certain fiber
    
    %
% Set up loop

% Preallocate
MAD_deltaRB = NaN(nTaxa, nFibers);
MAD_foldRB = NaN(nTaxa, nFibers);
MAD_logRB = NaN(nTaxa, nFibers);
MAD_logRB_null = NaN(nTaxa, nFibers);

COVAD_deltaRB = NaN(nTaxa, nFibers);
COVAD_foldRB = NaN(nTaxa, nFibers);
COVAD_logRB = NaN(nTaxa, nFibers);
COVAD_logRB_null = NaN(nTaxa, nFibers);

STDAD_deltaRB = NaN(nTaxa, nFibers);
STDAD_foldRB = NaN(nTaxa, nFibers);
STDAD_logRB = NaN(nTaxa, nFibers);
STDAD_logRB_null = NaN(nTaxa, nFibers);

% Go one taxa at at time
for i = 1:nTaxa

    % Get the delta, fold, and log change in relabun values for taxa_i
    taxaNow_allFibers_allDonors_deltaRB = deltaRB(:,i);
    taxaNow_allFibers_allDonors_foldRB = foldRB(:,i);
    taxaNow_allFibers_allDonors_logRB = logRB(:,i);
    taxaNow_allFibers_allDonors_logRB_null = logRBNull(:,i);
    
    % Now loop through the individual fibers
    for k = 1:nFibers

        fiberNow = unqFibers(k);

        % Get Delta/Fold/Log relabun values for this taxa and fiber, all donors
        taxaNow_fiberNow_allDonors_deltaRB = taxaNow_allFibers_allDonors_deltaRB(rbFinal.fiber_type == fiberNow);
        taxaNow_fiberNow_allDonors_foldRB = taxaNow_allFibers_allDonors_foldRB(rbFinal.fiber_type == fiberNow);
        taxaNow_fiberNow_allDonors_logRB = taxaNow_allFibers_allDonors_logRB(rbFinal.fiber_type == fiberNow);
        taxaNow_fiberNow_allDonors_logRB_null = taxaNow_allFibers_allDonors_logRB_null(rbFinal.fiber_type == fiberNow);

        % Get the average delta/fold/log change in relabun for this taxa and fiber across all donors that had the taxa at time zero
        % MAD = Mean Among Donors
        MAD_deltaRB(i,k) = mean(taxaNow_fiberNow_allDonors_deltaRB, 'omitnan');
        MAD_foldRB(i,k) = mean(taxaNow_fiberNow_allDonors_foldRB, 'omitnan');
        MAD_logRB(i,k) = mean(taxaNow_fiberNow_allDonors_logRB, 'omitnan');
        MAD_logRB_null(i,k) = mean(taxaNow_fiberNow_allDonors_logRB_null, 'omitnan');

        % Number of donors that had this taxa in samples with the current fiber
        nSamples = numel(taxaNow_fiberNow_allDonors_logRB(~isnan(taxaNow_fiberNow_allDonors_logRB)));

        % If we have at least three or so samples that have this taxa with the
        % current fiber, then take the std/cov of the change
        zzReal(i,k) = nSamples;
        if nSamples >= minDonorsToCalcSTDAD
            % COVAD = Coefficient Of Variation Among Donors
            COVAD_deltaRB(i,k) = std_noNaN_noInf(taxaNow_fiberNow_allDonors_deltaRB)/abs(MAD_deltaRB(i,k));
            COVAD_foldRB(i,k) = std_noNaN_noInf(taxaNow_fiberNow_allDonors_foldRB)/abs(MAD_foldRB(i,k));
            COVAD_logRB(i,k) = std_noNaN_noInf(taxaNow_fiberNow_allDonors_logRB)/abs(MAD_logRB(i,k));

            % STDAD = STAndard Deviation Among Donors
            STDAD_deltaRB(i,k) = std_noNaN_noInf(taxaNow_fiberNow_allDonors_deltaRB);
            STDAD_foldRB(i,k) = std_noNaN_noInf(taxaNow_fiberNow_allDonors_foldRB);
            STDAD_logRB(i,k) = std_noNaN_noInf(taxaNow_fiberNow_allDonors_logRB);
        else
            COVAD_deltaRB(i,k) = NaN;
            COVAD_foldRB(i,k) = NaN;
            COVAD_logRB(i,k) = NaN;

            STDAD_deltaRB(i,k) = NaN;
            STDAD_foldRB(i,k) = NaN;
            STDAD_logRB(i,k) = NaN;
        end

        % Calculate COVAD and STDAD for the null model
        nSamplesNull = numel(taxaNow_fiberNow_allDonors_logRB_null(~isnan(taxaNow_fiberNow_allDonors_logRB_null)));
        zzNull(i,k) = nSamplesNull;
        if nSamplesNull >= minDonorsToCalcSTDAD
            % COVAD = Coefficient Of Variation Among Donors
            COVAD_logRB_null(i,k) = std_noNaN_noInf(taxaNow_fiberNow_allDonors_logRB_null)/abs(MAD_logRB_null(i,k));

            % STDAD = STAndard Deviation Among Donors
            STDAD_logRB_null(i,k) = std_noNaN_noInf(taxaNow_fiberNow_allDonors_logRB_null);
        else
            COVAD_logRB_null(i,k) = NaN;
            STDAD_logRB_null(i,k) = NaN;
        end


    end
end


%% Optional: save MAD and COVAD of each metric as spreadsheets
    
    %
% Save the MAD (mean among donors) and 
% COVAD (coefficient of variation among donors
% These all share common metadata

save(strcat(saveFilePrefix,'meta.mat'), 'finalMeta') % all final (after fermentation) all technical reps
save(strcat(saveFilePrefix,'metaNoTechReps.mat'), 'samplesMeta') % all samples with tech reps averaged out
save(strcat(saveFilePrefix,'taxaNames.mat'), 'taxaNames') % valid table variable names
save(strcat(saveFilePrefix,'fullTaxaNames.mat'), 'fullTaxaNames') % full taxa names

%%% MAD
save(strcat(saveFilePrefix,'MAD_deltaRB.mat'),'MAD_deltaRB')
save(strcat(saveFilePrefix,'MAD_foldRB.mat'),'MAD_foldRB')
save(strcat(saveFilePrefix,'MAD_logRB.mat'),'MAD_logRB')

%%% COVAD
save(strcat(saveFilePrefix,'COVAD_deltaRB.mat'),'COVAD_deltaRB')
save(strcat(saveFilePrefix,'COVAD_foldRB.mat'),'COVAD_foldRB')
save(strcat(saveFilePrefix,'COVAD_logRB.mat'),'COVAD_logRB')

%%% STDAD
save(strcat(saveFilePrefix,'STDAD_deltaRB.mat'),'STDAD_deltaRB')
save(strcat(saveFilePrefix,'STDAD_foldRB.mat'),'STDAD_foldRB')
save(strcat(saveFilePrefix,'STDAD_logRB.mat'),'STDAD_logRB')

%%% NUll MODEL
save(strcat(saveFilePrefix,'MAD_logRB_null.mat'),'MAD_logRB_null')
save(strcat(saveFilePrefix,'COVAD_logRB_null.mat'),'COVAD_logRB_null')
save(strcat(saveFilePrefix,'STDAD_logRB_null.mat'),'STDAD_logRB_null')
%}


%% Make plots of STDAD vs MAD (real and null)

minLogChange = 0; % throw out any log changes in relabun lower than this

markerStyle = {'o', 's', 'd', '^', 'x', '*', 'p', '+', 'v', '<', '>',...
            'o', 's', 'd', '^', 'x', '*', 'p', '+', 'v', '<', '>'};
ms = 6;
colors = lines;

%%% Real data
[MAD_logRB, STDAD_logRB] = process_MAD_and_COVAD(MAD_logRB, STDAD_logRB, minLogChange, absFlag);

figure

subplot(1,2,1)
hold on
box on

% Make the plot
for i = 1:width(MAD_logRB)

    scatter(MAD_logRB(:,i), STDAD_logRB(:,i), 8, "filled", markerStyle{i})

end

% Combine points from all fibers
STDAD_all = reshape(STDAD_logRB, numel(STDAD_logRB), 1);
STDAD_all = STDAD_all(~isnan(STDAD_all));

MAD_all = reshape(MAD_logRB, numel(MAD_logRB), 1);
MAD_all = MAD_all(~isnan(STDAD_all));

% Linear regression
[mdl, xModel, yModel, pValTxt, R2txt] = linearReg(MAD_all, STDAD_all);

% Plot linear regression
p1 = plot(xModel, yModel, 'r-', 'LineWidth', 1);

legend(p1, strcat('Linear regression: ', R2txt, pValTxt));

ylabel('STD among donors')
xlabel('Log change relative abundance')
title('Real data')



%%% Null model
[MAD_logRB_null, STDAD_logRB_null] = process_MAD_and_COVAD(MAD_logRB_null, STDAD_logRB_null, minLogChange, absFlag);

subplot(1,2,2)
hold on
box on

% Make the plot
for i = 1:width(MAD_logRB)

    scatter(MAD_logRB_null(:,i), STDAD_logRB_null(:,i), 8, "filled", markerStyle{i})

end

% Combine points from all fibers
STDAD_all_null = reshape(STDAD_logRB_null, numel(STDAD_logRB_null), 1);
STDAD_all_null = STDAD_all_null(~isnan(STDAD_all_null));

MAD_all_null = reshape(MAD_logRB_null, numel(MAD_logRB_null), 1);
MAD_all_null = MAD_all_null(~isnan(STDAD_all_null));

% Linear regression
[mdlNull, xModelNull, yModelNull, pValTxt, R2txt] = linearReg(MAD_all_null, STDAD_all_null);

% Plot linear regression
p2 = plot(xModelNull, yModelNull, 'r-', 'LineWidth', 1);

legend(p2, strcat('Linear regression: ', R2txt, pValTxt));

ylabel('STD among donors')
xlabel('Log change relative abundance')
title('Null model (randomized data)')

set(gcf, 'Position', [365	429.800000000000	1144	420.000000000000])

% Get number of points on each graph
nReal = numel(MAD_all(~isnan(MAD_all)))
nNull = numel(MAD_all_null(~isnan(MAD_all_null)))

%}

%}
elapsed_time = round(toc,3) / 60












