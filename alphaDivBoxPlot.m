function alphaDivBoxPlot(g,y,titleTxt)

% alphaDivBoxPlot.m
% Graph alpha diversity indices using a boxplot
% See: whistlerPoster_study1_shanbar @ D:\Specificity\Pipeline2\Study1

unqFibers = unique(g);
Nfibers = length(unqFibers);

% Make boxplot
positions = linspace(1, (1+Nfibers*0.5), Nfibers);
boxplot(y, g, 'positions', positions, 'Symbol', '', 'Whisker', inf) %'Notch', 'on', 

% Change whisker from dashed to solid line
h=findobj('LineStyle','--');
set(h, 'LineStyle','-');

% Now add all points
hold on
spread = 0.04; % 0=no spread:0.5=random spread within box bounds

% Add all points
for i = 1:Nfibers

    yPts = y(g == unqFibers(i));

    x = rand(size(yPts))*spread -(spread/2) + positions(i);

    plot(x, yPts, 'ko', 'linewidth', 1, 'MarkerSize', 3)%, 'MarkerFaceColor', 'k')
    clear yPts
    
end

% Redraw medians so they are on top
boxes = findobj(gca,'Tag','Box');
medians = findobj(gca,'tag','Median');

for j = 1:length(boxes)
   line(get(medians(j),'XData'), get(medians(j),'YData'), 'Color', 'r','LineWidth', 2);
end

% Plot style
set(gca,'linewidth', 1) %change box line width

title(titleTxt, 'FontSize', 10)

% yticks(0:30:240)
% ylim([0 250])

% Set x and y tick mark properties
xAX = get(gca,'XAxis');
yAX = get(gca,'YAxis');
set(xAX,'FontSize', 11); %X tick fontsize
set(yAX,'FontSize', 11); %Y tick fontsize


end







