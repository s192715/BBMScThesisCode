%This script finds the lat vs height distribution of Es events.
%ALSO, finds height vs. local time relative Es occurrence
%Displays as relative occurrence i.e. bins Es events into lat vs height and
%divides each bin by total number of Es detected.


%for LocalTime test
% LocalTimeIndex = find(TotalEventLocalTime >= 16 & TotalEventLocalTime <= 22);
% TotalEventLat = TotalEventLat(LocalTimeIndex);
% TotalEventHeight = TotalEventHeight(LocalTimeIndex);
% TotalEventLocalTime = TotalEventLocalTime(LocalTimeIndex);

%load('/zhome/e8/9/144512/Desktop/MScThesis/MatlabScripts/SavedMatlabVariables/ForGlobalMaps/2011seasons/2011winter.mat')
year = '2007-2019';
season = 'fall';
%lat = '-5° to 5°';

%-------------------------------------------------------
%to only use measurements in a specified latitude band:
%IndexToUse = find(TotalEventLat > 30 & TotalEventLat < 40);
%IndexToUse = find(TotalEventLat > -40 & TotalEventLat < -30);
% IndexToUse = find(TotalEventLat > -5 & TotalEventLat < 5);
% 
% TotalEventLat = TotalEventLat(IndexToUse);
% TotalEventLocalTime = TotalEventLocalTime(IndexToUse);
% TotalEventHeight = TotalEventHeight(IndexToUse);
%-------------------------------------------------------



%In script GlobalMaps.m the grid is made according to lat/lon or lat/time of RO
%measurements and values TotalS4_present were split into the grid bins,
%now the grid is made according to lat/height of Es detections and no data needs to be split in bins.
%latitude vs height
minlat = -90;     %adjust these if needed. They determine where the bins start
minheight = 1; 
gridspacing_lat = 1; %degree
gridspacing_time = 1; %one km
latidx = 1 + floor((TotalEventLat - minlat) ./ gridspacing_lat);
heightdx = 1 + floor((TotalEventHeight - minheight) ./ gridspacing_time);

A = accumarray( [latidx(:), heightdx(:)], 1); %number of detections within each bin

%IF I want to smoothen data
% K = (1/9)*ones(3); %defining kernel
% A = conv2(A,K,'same');
% [row column] = find(A==max(max(A)));
%A(row,column) = 0;

xticks = 1:20:length(A(1,:)); %deciding tick placements
xticks(2:end) = xticks(2:end) -1;
%yticks = -length(A(:,1))/2:length(A(:,1))/6:length(A(:,1))/2; %for 5 degree lat resolution
%yticks(1) = -length(A(:,1))/2 +1;
yticks = -(length(A(:,1))+1)/2:(length(A(:,1))+1)/6:(length(A(:,1))+1)/2; %for 1 deg lat resolution

figure(1); imagesc(xticks, yticks, (A/max(max(A)))) %detection in each bin divided by total detections
%imagesc(-179:1:180,-89:1:90,A_S4_perc)
h = colorbar;
%ylabel(h, '%')
colormap(jet)
%caxis([0 0.5]) %for 5 deg lat resolution
%caxis([0 0.07]) %for 1 deg lat resolution
xlabel('Height [km]')
ylabel('Latitude [°]')

ylabels = {'-90' '-60' '-30' '0' '30' '60' '90'};  %time labels
set(gca, 'YTick', yticks, 'YTickLabel', ylabels);
set(gca,'YDir','normal')
%xlim([20 130]) %20 because Es events were cut below 20

title(['Normalized Es occurrence, ' year ' ' season '.'])


% %saving figure
% baseFileName = sprintf([num2str(year) num2str(season) '.png'],1);
% % Specify some particular, specific folder:
% fullFileName = fullfile(['/zhome/e8/9/144512/Desktop/MScThesis/Figures/EsLatvsHeight/seasons/', num2str(year)], baseFileName);  
% saveas(figure(1),fullFileName);
% close()


%% Height vs LocalTime


minheight = 1;  
mintime = 0;
gridspacing_height = 1;
gridspacing_time = 0.25; %15min
heightdx = 1 + floor((TotalEventHeight - minheight) ./ gridspacing_height);
timedx = 1 + floor((TotalEventLocalTime - mintime) ./ gridspacing_time);

A = accumarray( [heightdx(:), timedx(:)], 1); %number of detections within each bin


%IF I want to smoothen data
% K = (1/9)*ones(3); %defining kernel
% A = conv2(A,K,'same');
% [row column] = find(A==max(max(A)));
%A(row,column) = 0;

yticks = 1:20:length(A(:,1)); %deciding tick placements
yticks(2:end) = yticks(2:end) -1;

xticks = 0:16:96; %for 15 min res
xlabels = 0:4:24; 

Norm = 667;
%Norm = max(max(A));
figure(103); imagesc(xticks, yticks, (A/Norm)) %detections in each bin divided by total detections
%imagesc(-179:1:180,-89:1:90,A_S4_perc)
h = colorbar;
%ylabel(h, '%')
colormap(jet)
caxis([0 0.25])
xlabel('Local Time [hr]')
ylabel('Height [km]')

%ylabels = {'-90' '-60' '-30' '0' '30' '60' '90'};  %time labels
set(gca, 'YTick', yticks, 'YTickLabel', yticks);
set(gca,'YDir','normal')
set(gca, 'XTick', xticks, 'XTickLabel', xlabels);
ylim([20 120]) %20 because Es events were cut below 20

%title({['Relative Es occurrence, ' year ' ' season '.'],['Latitudes: ' lat '.']})
title(['Normalized Es occurrence, ' year ' ' season '.'])

%% LatVsLocalTime Relative Es Occurrence

minlat = -90;     %adjust these if needed. They determine where the bins start
mintime = 0;
gridspacing_lat = 1;
gridspacing_time = 0.25; %15 min
latidx = 1 + floor((TotalEventLat - minlat) ./ gridspacing_lat);
timedx = 1 + floor((TotalEventLocalTime - mintime) ./ gridspacing_time);

A = accumarray( [latidx(:), timedx(:)], 1); %number of measurements within each bin

%IF I want to smoothen data
% K = (1/9)*ones(3); %defining kernel
% A = conv2(A,K,'same');
% [row column] = find(A==max(max(A)));

figure(3); imagesc(A./length(TotalEventHeight))
h = colorbar;
ylabel(h, '%')
colormap(jet)
%caxis([0 50])
xlabel('Local Time [hr]')
ylabel('Latitude [°]')

yticks = 0:length(A(:,1))/6:length(A(:,1)); %adjust as appropriate, positive integers only
yticks(1) = 1;
ylabels = {'-90' '-60' '-30' '0' '30' '60' '90'};  %time labels
set(gca, 'YTick', yticks, 'YTickLabel', ylabels);
set(gca,'YDir','normal')
%xticks = 0:4:24 %for 1hr res
xticks = 0:16:96; %for 15 min res
xlabels = 0:4:24; %adjust as appropriate, positive integers only
%xlabels = {'0' '4' '-60' '0' '60' '120' '180'};  %time labels
set(gca, 'XTick', xticks, 'XTickLabel', xlabels);
%caxis([0 0.1])

title(['Relative Es occurrence, ' year ' ' season '.'])


