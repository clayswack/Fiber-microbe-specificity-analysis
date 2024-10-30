%% de_campos_costa_hamaker_change_dist_plots.m

%% Objective:

% Make 3 violin plots of the distribution of changes to taxa:
% deltaRB 
% foldRB
% logRB

% Make 3 bar graphs all based on average log change between donors:
% N distinct taxa that increased by > 1 log
% N distinct taxa increased by > 1.5 log
% N distinct taxa  increased by > 2 log

% Make 1 big bar graph:
% Swowing the top five taxa ranked by average log change between donors
% Make a separate cluster for each fiber in the experiment
% Make sure to show the standard dev between donors as an errorbar

% Make 2 scatterplots
% COVAD vs MAD
% STDAD vs MAD

%% Setup

clear
clc
close all


%% User inputs
sgTitleTxt = 'Cantu-Jungles and Hamaker 2021'
fibers = {'FOS', 'Glucan', 'Pectin', 'RS'};
N_metadata_cols = 4;
barYmax = 10; % upper y lim for bar graphs of "N taxa changing > thresh"
bigBarGraphLims = [0, 2.5]; % limits for the big bar graph


%% Constants
colors = lines;
ms = 1; % markersize for all points
outlierSize = 10; % markersize for points with logChangeRB > thresh
threshs = [0.5, 0.75, 1] % log changes above this are plotted separately and the taxa id queried
NTopTaxa = 5; % grab the taxa names of the top NtopTaxa from each fiber


%% Get data

% Metadata and taxa names
load cantu_jungles_hamaker_meta.mat % finalMeta
% load cantu_jungles_hamaker_taxaNames.mat % taxaNames as valid table variable names
load cantu_jungles_hamaker_fullTaxaNames.mat % full taxaNames

% For VIOLIN plots
% DeltaRB, foldRB, logRB full tables (no averaging)
deltaRB = readtable('cantu_jungles_hamaker_delta_relabun.csv');
deltaRB.fiber_type = categorical(deltaRB.fiber_type);

foldRB = readtable('cantu_jungles_hamaker_fold_relabun.csv');
foldRB.fiber_type = categorical(foldRB.fiber_type);

logRB = readtable('cantu_jungles_hamaker_log_relabun.csv');
logRB.fiber_type = categorical(logRB.fiber_type);

% For BAR plots and SCATTER plots
% Average, stdev, and coefficient of variance for log increase between donors
load cantu_jungles_hamaker_MAD_logRB.mat % MAD_logRB
load cantu_jungles_hamaker_COVAD_logRB.mat % COVAD_logRB
load cantu_jungles_hamaker_STDAD_logRB.mat % STDAD_logRB

%% Set up plots

t = tiledlayout(2,3,'TileSpacing','Compact','Padding','Compact');

%% VIOLIN PLOTS

%%% Delta RB
% Preallocate
violinX = categorical;
violinY = [];
fibersCat = categorical(fibers);

% Put together data for violin plot
for f = 1:length(fibers)

    yNow = cleanDiffs(deltaRB{deltaRB.fiber_type == fibersCat(f), N_metadata_cols+1:end});
    violinY = vertcat(violinY, yNow);
    violinX = vertcat(violinX, repmat(fibersCat(f), length(yNow), 1));

end

% Make violin plot
nexttile
vs = violinplot(violinY, violinX, 'MarkerSize', ms, 'ShowBox', false, 'ShowMedian', false, 'ShowWhiskers', false, 'EdgeColor', [0, 0, 0], 'ViolinAlpha', 0);

% Mark top 5 changing taxa with large "X"
for f = 1:length(fibers)
    yNow = cleanDiffs(deltaRB{deltaRB.fiber_type == fibersCat(f), N_metadata_cols+1:end});
    plot(f*ones(5,1), maxk(yNow,5), 'x', 'MarkerEdgeColor', colors(f,:), 'MarkerFaceColor', colors(f,:), 'MarkerSize', outlierSize)
end

title('Delta relative abundance')


%%% Fold RB
% Preallocate
violinX = categorical;
violinY = [];
fibersCat = categorical(fibers);

% Put together data for violin plot
for f = 1:length(fibers)

    yNow = cleanDiffs(foldRB{foldRB.fiber_type == fibersCat(f), N_metadata_cols+1:end});
    violinY = vertcat(violinY, yNow);
    violinX = vertcat(violinX, repmat(fibersCat(f), length(yNow), 1));

end

% Make violin plot
nexttile
vs = violinplot(violinY, violinX, 'MarkerSize', ms, 'ShowBox', false, 'ShowMedian', false, 'ShowWhiskers', false, 'EdgeColor', [0, 0, 0], 'ViolinAlpha', 0);

% Mark top 5 changing taxa with large "X"
for f = 1:length(fibers)
    yNow = cleanDiffs(foldRB{foldRB.fiber_type == fibersCat(f), N_metadata_cols+1:end});
    plot(f*ones(5,1), maxk(yNow,5), 'x', 'MarkerEdgeColor', colors(f,:), 'MarkerFaceColor', colors(f,:), 'MarkerSize', outlierSize)
end
title('Fold change relative abundance')


%%% Log RB
% Preallocate
violinX = categorical;
violinY = [];
fibersCat = categorical(fibers);

% Put together data for violin plot
for f = 1:length(fibers)

    yNow = cleanDiffs(logRB{logRB.fiber_type == fibersCat(f), N_metadata_cols+1:end});
    violinY = vertcat(violinY, yNow);
    violinX = vertcat(violinX, repmat(fibersCat(f), length(yNow), 1));

end

% Make violin plot
nexttile
vs = violinplot(violinY, violinX, 'MarkerSize', ms, 'ShowBox', false, 'ShowMedian', false, 'ShowWhiskers', false, 'EdgeColor', [0, 0, 0], 'ViolinAlpha', 0);

% Mark top 5 changing taxa with large "X"
for f = 1:length(fibers)
    yNow = cleanDiffs(logRB{logRB.fiber_type == fibersCat(f), N_metadata_cols+1:end});
    plot(f*ones(5,1), maxk(yNow,5), 'x', 'MarkerEdgeColor', colors(f,:), 'MarkerFaceColor', colors(f,:), 'MarkerSize', outlierSize)
end

hold off
title('Log change relative abundance')

%% BAR PLOTS

% Remove NaN, Inf, and -Inf from the log changes lists
fos_MAD_logRB = cleanDiffs(MAD_logRB(:,1));
glucan_MAD_logRB = cleanDiffs(MAD_logRB(:,2));
pectin_MAD_logRB = cleanDiffs(MAD_logRB(:,3));
rs_MAD_logRB = cleanDiffs(MAD_logRB(:,4));


%%% N changes > threshs(1)
N = [];
for f = 1:length(fibersCat)
    yNow = cleanDiffs(MAD_logRB(:,f));
    N(f) = length(yNow(yNow>threshs(1)));
end

nexttile
b = bar(N);

set(gca, 'xticklabel', fibers, 'Fontsize', 10)
rotateXLabels( gca(), -45)

b.FaceColor = 'flat';

for f = 1:length(fibersCat)
    b(1).CData(f,:) = colors(f,:);
end

title(strcat('N distinct taxa increasing by > ', num2str(threshs(1)), '{ }', 'log'))
ylim([0, barYmax])

%%% N changes > threshs(2)
N = [];
for f = 1:length(fibersCat)
    yNow = cleanDiffs(MAD_logRB(:,f));
    N(f) = length(yNow(yNow>threshs(2)));
end

nexttile
b = bar(N);

set(gca, 'xticklabel', fibers, 'Fontsize', 10)
rotateXLabels( gca(), -45)

b.FaceColor = 'flat';
for f = 1:length(fibersCat)
    b(1).CData(f,:) = colors(f,:);
end

title(strcat('N distinct taxa increasing by > ', num2str(threshs(2)), '{ }', 'log'))
ylim([0, barYmax])

%%% N changes > threshs(3)
N = [];
for f = 1:length(fibersCat)
    yNow = cleanDiffs(MAD_logRB(:,f));
    N(f) = length(yNow(yNow>threshs(3)));
end

nexttile
b = bar(N);

set(gca, 'xticklabel', fibers, 'Fontsize', 10)
rotateXLabels( gca(), -45)

b.FaceColor = 'flat';
for f = 1:length(fibersCat)
    b(1).CData(f,:) = colors(f,:);
end

title(strcat('N distinct taxa increasing by > ', num2str(threshs(3)), '{ }', 'log'))
ylim([0, barYmax])

sgtitle(sgTitleTxt, 'FontSize', 12)
set(gcf, 'Position', [378.600000000000	278.600000000000	1208.80000000000	666.400000000000])

%% BIG BAR GRAPH OF TOP 5 INCREASING TAXA WITH NAMES
    
% Restructure data for bar graph with arbitrary number of fibers
% fibers = categorical(fibers);

tix = 1;
for f = 1:length(fibersCat)

    % Start by removing Nan, Inf, and -Inf from the data and remove the
    % corresponding taxa names
    [cleanMAD, cleanSTDAD, updatedTaxaNames] = clean_MAD_STDAD_andUpdateTaxaNames(MAD_logRB(:,f), STDAD_logRB(:,f), fullTaxaNames);

    % Next, sort the changes from high to low
    [sortedChanges, ixs] = sort(cleanMAD, 'descend');

    % Get the taxa names
    sortedTaxaNames = updatedTaxaNames(ixs);
    taxaNamesLabels(tix:tix+NTopTaxa-1) = sortedTaxaNames(1:5);
    tix = tix + NTopTaxa;

    % Get the standard dev among donors (STDAD)
    sortedSTDAD = cleanSTDAD(ixs);

    % Now, get the top 5 changes
    MADForBar(f,:) = sortedChanges(1:NTopTaxa);
    STDADForBar(f,:) = sortedSTDAD(1:5);
end

figure

b2 = barwitherr(STDADForBar, MADForBar)

% Make all bars from the same fiber have the same color
for f = 1:length(fibersCat)

    for j = 1:NTopTaxa
    
        b2(j).FaceColor = 'flat';
        b2(j).CData(f,:) = colors(f,:);
    
    end
end

% Now get X coordinate of the bars
for j = 1:NTopTaxa

    % Get bar X ticks
    endpoints(j,:) = b2(j).XEndPoints;

end

ticks = reshape(endpoints, numel(endpoints), 1);

% Create Tick Labels
% xLab = zLogOut.taxaName;
xLab = taxaNamesLabels;

% Set individual ticks
set(gca, 'XTick', ticks, 'XTickLabel', xLab, 'xticklabelrotation', 90)

% Set x and y tick mark properties
xAX = get(gca,'XAxis');
set(xAX,'FontSize', 6); %X tick fontsize

% % yAX = get(gca,'YAxis');
% % set(yAX,'FontSize', 12); %Y tick fontsize

ylim(bigBarGraphLims)

sgtitle(strcat(sgTitleTxt, '{: }', 'Top 5 distinct taxa by log change, for each fiber'), 'FontSize', 12)
set(gcf, 'Position', [541.800000000000	85	646.400000000000	924])


%}





