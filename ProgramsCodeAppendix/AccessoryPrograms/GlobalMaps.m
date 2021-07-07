%This script loads lat,lon,Time, of each RO measurement along with if S4 and
%Es present in measurement. Plots global maps of S4 and Es occurrence in %,
%in lat vs lon and lat vs TimeOfDay.

%for LocalTime test
% LocalTimeIndex = find(TotalROLocalTime >= 16 & TotalROLocalTime <= 22);
% TotalEs_present = TotalEs_present(LocalTimeIndex);
% TotalS4_present = TotalS4_present(LocalTimeIndex);
% TotalRO_lat = TotalRO_lat(LocalTimeIndex);
% TotalRO_lon = TotalRO_lon(LocalTimeIndex);
% TotalROLocalTime = TotalROLocalTime(LocalTimeIndex);
% TotalRO_MaxS4 = TotalRO_MaxS4(LocalTimeIndex);

% load('/zhome/e8/9/144512/Desktop/MScThesis/MatlabScripts/SavedMatlabVariables/ForGlobalMaps/2009.mat')
% TotalRO_lat = TotalRO_lat_2009;
% TotalRO_lon = TotalRO_lon_2009;
% TotalRO_Time = TotalRO_Time_2009;
% TotalEs_present = TotalEs_present_2009;
% TotalS4_present = TotalS4_present_2009;
%load('/zhome/e8/9/144512/Desktop/MScThesis/MatlabScripts/SavedMatlabVariables/ForGlobalMaps/2010.mat')
%year = 2010; %which year we're plotting, used in titles in the plots.
season = 'Summer';
year = '2007-2019';
%load('/zhome/e8/9/144512/Desktop/MScThesis/MatlabScripts/SavedMatlabVariables/ForGlobalMaps/2011seasons/2011winter.mat')
%year = 2012;
%season = 'winter';

%load('/zhome/e8/9/144512/Desktop/MScThesis/MatlabScripts/SavedMatlabVariables/ForGlobalMaps/2013seasons/2013fall.mat')

%converting Time to LocalTime
% TotalROTimeUTC  = datevec(datetime(TotalRO_Time,'ConvertFrom','juliandate'));
% TotalROTimeUTC_HourOfDay = 24*datenum(hours(TotalROTimeUTC(:,4)) + minutes(TotalROTimeUTC(:,5)) + seconds(TotalROTimeUTC(:,6)));
% 
% [zd] = timezone(TotalRO_lon);
% TotalROLocalTime = (TotalROTimeUTC_HourOfDay - zd)';
% for i = 1:length(TotalROLocalTime)
%     if TotalROLocalTime(i) < 0
%         TotalROLocalTime(i) = 24 + TotalROLocalTime(i);
%     elseif TotalROLocalTime(i) >= 24
%         TotalROLocalTime(i) = TotalROLocalTime(i) - 24;
%     end
% end



minlat = -90;     %adjust these if needed. They determine where the bins start
minlong = -180;
gridspacing = 1; %spatial resolution in degrees
latidx = 1 + floor((TotalRO_lat - minlat) ./ gridspacing);
longidx = 1 + floor((TotalRO_lon - minlong) ./ gridspacing);

A = accumarray( [latidx(:), longidx(:)], 1); %number of measurements within each bin
A_S4 = accumarray( [latidx(:), longidx(:)], TotalS4_present'); %number of S4 detections within bin
A_Es = accumarray( [latidx(:), longidx(:)], TotalEs_present'); %number of Es detections within bin

A_S4_perc = 100*(A_S4./A); %some values become nan because there was no A was zero (no measurements)
A_Es_perc = 100*(A_Es./A);

%If I want to smooth the data:
% K = (1/9)*ones(3); %defining kernel
% A_S4_perc = conv2(A_S4_perc,K,'same');
% A_Es_perc = conv2(A_Es_perc,K,'same');

%A_S4_perc(isnan(A_S4_perc)) = 0; %sets all nan values to zero
% 
% xticks = -length(A(1,:))/2:length(A(1,:))/6:length(A(1,:))/2; %deciding tick placements
% xticks(1) = -length(A(1,:))/2 +1;
% yticks = -length(A(:,1))/2:length(A(:,1))/6:length(A(:,1))/2;
% yticks(1) = -length(A(:,1))/2 +1;
% 
% figure(1); imagesc(xticks, yticks, A_S4_perc)
% %imagesc(-179:1:180,-89:1:90,A_S4_perc)
% h = colorbar;
% ylabel(h, '%')
% colormap(jet)
% caxis([0 100])
% xlabel('Longitude [°]')
% ylabel('Latitude [°]')
% 
% ylabels = {'-90' '-60' '-30' '0' '30' '60' '90'};  %time labels
% set(gca, 'YTick', yticks, 'YTickLabel', ylabels);
% set(gca,'YDir','normal')
% xlabels = {'-180' '-120' '-60' '0' '60' '120' '180'};  %time labels
% set(gca, 'XTick', xticks, 'XTickLabel', xlabels);
% 
% title(['Global S4 scintillation occurrence ' num2str(year)])

%plotting magnetic equator, only works correctly for spatial resolution 1°x1°
%the magnetic equator here is defined as the lat and lon where the magnetic
%dip angle at 60 km height is zero.
load('/zhome/e8/9/144512/Desktop/MScThesis/MatlabScripts/SavedMatlabVariables/MagneticEquator.mat')
% hold on;
% plot(lon,lat(LatIndex),'w-')


%In case I want to plot NaN values different color than zero.
% imAlpha=ones(size(A_S4_perc));
% imAlpha(isnan(A_S4_perc))=0;
% imagesc(A_S4_perc,'AlphaData',imAlpha);
% set(gca,'color',[0.2 0.2 0.2]);
% colorbar
% colormap(jet)



% figure(2); imagesc(A_Es_perc)
% h = colorbar;
% ylabel(h, '%')
% colormap(jet)
% caxis([0 100])
% xlabel('Longitude [°]')
% ylabel('Latitude [°]')
% 
% yticks = 0:length(A(:,1))/6:length(A(:,1)); %adjust as appropriate, positive integers only
% yticks(1) = 1;
% ylabels = {'-90' '-60' '-30' '0' '30' '60' '90'};  %time labels
% set(gca, 'YTick', yticks, 'YTickLabel', ylabels);
% set(gca,'YDir','normal')
% xticks = 0:length(A(1,:))/6:length(A(1,:)); %adjust as appropriate, positive integers only
% xticks(1) = 1;
% xlabels = {'-180' '-120' '-60' '0' '60' '120' '180'};  %time labels
% set(gca, 'XTick', xticks, 'XTickLabel', xlabels);
% 
% title(['Global Es occurrence ' num2str(year)])

%--------------------------------------------------------------------------
%% %latitude vs local time
minlat = -90;     %adjust these if needed. They determine where the bins start
mintime = 0;
gridspacing_lat = 1;
gridspacing_time = 0.25; %one hour
latidx = 1 + floor((TotalRO_lat - minlat) ./ gridspacing_lat);
timedx = 1 + floor((TotalROLocalTime - mintime) ./ gridspacing_time);

A2 = accumarray( [latidx(:), timedx(:)], 1); %number of measurements within each bin
A2_S4 = accumarray( [latidx(:), timedx(:)], TotalS4_present'); %number of S4 detections within bin
A2_Es = accumarray( [latidx(:), timedx(:)], TotalEs_present'); %number of Es detections within bin

A2_S4_perc = 100*(A2_S4./A2);
A2_Es_perc = 100*(A2_Es./A2);


%IF I want to smoothen data
% K = (1/9)*ones(3); %defining kernel
% A2_S4_perc = conv2(A2_S4_perc,K,'same');
% A2_Es_perc = conv2(A2_Es_perc,K,'same');

figure(3); imagesc(A2_S4_perc)
h = colorbar;
ylabel(h, '%')
colormap(jet)
caxis([0 100])
xlabel('Local Time [hr]')
ylabel('Latitude [°]')

yticks = 0:length(A2(:,1))/6:length(A2(:,1)); %adjust as appropriate, positive integers only
yticks(1) = 1;
ylabels = {'-90' '-60' '-30' '0' '30' '60' '90'};  %time labels
set(gca, 'YTick', yticks, 'YTickLabel', ylabels);
set(gca,'YDir','normal')
%xticks = 0:4:24 %for 1hr res
xticks = 0:16:96; %for 15 min res
xlabels = 0:4:24; %adjust as appropriate, positive integers only
%xlabels = {'0' '4' '-60' '0' '60' '120' '180'};  %time labels
set(gca, 'XTick', xticks, 'XTickLabel', xlabels);

title(['Global S4 scintillation occurrence '  season ' ' num2str(year)])


figure(4); imagesc(A2_Es_perc)
h = colorbar;
ylabel(h, '%')
colormap(jet)
caxis([0 50])
xlabel('Local Time [hr]')
ylabel('Latitude [°]')

yticks = 0:length(A2(:,1))/6:length(A2(:,1)); %adjust as appropriate, positive integers only
yticks(1) = 1;
ylabels = {'-90' '-60' '-30' '0' '30' '60' '90'};  %time labels
set(gca, 'YTick', yticks, 'YTickLabel', ylabels);
set(gca,'YDir','normal')
%xticks = 0:4:24 %for 1hr res
xticks = 0:16:96; %for 15 min res
xlabels = 0:4:24; %adjust as appropriate, positive integers only
%xlabels = {'0' '4' '-60' '0' '60' '120' '180'};  %time labels
set(gca, 'XTick', xticks, 'XTickLabel', xlabels);

title(['Global Es occurrence '  season ' ' num2str(year)])


%Plotting for S4 max amplitude
A2_S4max = accumarray( [latidx(:), timedx(:)], TotalRO_MaxS4'); %sum of MaxS4 amplitude within each bin
A2_S4max_mean = (A2_S4max./A2); %A2 is number of measurement within each bin

%If I want to smoothen the data
%A2_S4max_mean = conv2(A2_S4max_mean,K,'same');

figure(30); imagesc(A2_S4max_mean)
h = colorbar;
ylabel(h, 'S4')
colormap(jet)
caxis([0 1])
xlabel('Local Time [hr]')
ylabel('Latitude [°]')

yticks = 0:length(A2(:,1))/6:length(A2(:,1)); %adjust as appropriate, positive integers only
yticks(1) = 1;
ylabels = {'-90' '-60' '-30' '0' '30' '60' '90'};  %time labels
set(gca, 'YTick', yticks, 'YTickLabel', ylabels);
set(gca,'YDir','normal')
%xticks = 0:4:24 %for 1hr res
xticks = 0:16:96; %for 15 min res
xlabels = 0:4:24; %adjust as appropriate, positive integers only
%xlabels = {'0' '4' '-60' '0' '60' '120' '180'};  %time labels
set(gca, 'XTick', xticks, 'XTickLabel', xlabels);

title(['Global S4 amplitude ' season ' ' num2str(year)])

%close ALL
%% 

%--------------------------------------------------------------------------
%Show on actual global maps


%--------------------------------------------------------------------------
%TEST SECTION, trying to plot NaN values in different color than white
% Backup = A_S4_perc;
% 
% A_S4_perc = Backup;
% 
% imAlpha=zeros(size(A_S4_perc));
% imAlpha(isnan(A_S4_perc)) = 1; %imAlpha are the indexes are NaN
% A_S4_perc(isnan(A_S4_perc)) = -100; %plots NaN values in gray color
% 
% Background = 0.01*ones(size(A_S4_perc'));
% Background(1,1) = 0;
% Background(1,2) = 1;
% 
% 
% theta = [-179:180];
% phi = [-89:90];
% %obs % 360x180 2d array of data
% [lat_grid,lon_grid] = meshgrid(phi,theta); % should both be 360x180 to match data
% load coast % loads lat and long variables that define the coastline
% figure(11); ax1 = worldmap('World'); % also try axesm as it gives more options
% ax2 = worldmap('World')
% b = pcolorm(lat_grid,lon_grid,Background,'axes',ax1);
% 
% 
% c = colormap(ax1,gray);
% 
% c2 = colormap(ax2,jet);
% b2 = pcolorm(lat_grid,lon_grid,A_S4_perc');
% 
% 
% coastlines = plotm(lat,long,'w-') %draw the coastlines
% load('/zhome/e8/9/144512/Desktop/MScThesis/MatlabScripts/SavedMatlabVariables/MagneticEquator.mat')
% MagnEq = plotm(lat(LatIndex),lon,'w-') %plotting the magnetic equator
% colormap(ax2,jet)
% h = colorbar;
% ylabel(h, '%')
% caxis([0 100])
% title({['Global S4 scintillations ' num2str(year) ' ' season '.'],['Total Occurrence: ' num2str((100*sum(TotalS4_present)/length(TotalS4_present)),3) '%.']})

%133018 measurements in winter 2011
%12602 total NaN values (1°x1° bin where no measurements were recoreded)

%246464 measurements in winter 2009
%5400 total NaN values

%179361 mælingar summer 2011
%8631 total NaN values
%--------------------------------------------------------------------------


theta = [-179:180];
phi = [-89:90];
%obs % 360x180 2d array of data
[lat_grid,lon_grid] = meshgrid(phi,theta); % should both be 360x180 to match data
load coast % loads lat and long variables that define the coastline
figure(11); worldmap('World'); % also try axesm as it gives more options
pcolorm(lat_grid,lon_grid,A_S4_perc');
plotm(lat,long,'w-') %draw the coastlines
load('/zhome/e8/9/144512/Desktop/MScThesis/MatlabScripts/SavedMatlabVariables/MagneticEquator.mat')
plotm(lat(LatIndex),lon,'w-') %plotting the magnetic equator
h = colorbar;
ylabel(h, '%')
colormap(jet)
caxis([0 100])
%title({['Global S4 scintillations ' num2str(year) '.'],['Total Occurrence: ' num2str((100*sum(TotalS4_present)/length(TotalS4_present)),3) '%.']})
title({['Global S4 scintillations ' num2str(year) ' ' season '.'],['Total Occurrence: ' num2str((100*sum(TotalS4_present)/length(TotalS4_present)),3) '%.']})
%set(gca,'Visible','on','color',[0 0 0])

% %saving figure
% baseFileName = sprintf(['S4' num2str(year) num2str(season) '.png'],1);
% % Specify some particular, specific folder:
% fullFileName = fullfile(['/zhome/e8/9/144512/Desktop/MScThesis/Figures/GlobalMaps/', num2str(year), '/seasons/'], baseFileName);  
% saveas(figure(10),fullFileName);

figure(12); worldmap('World') % also try axesm as it gives more options
hold on;
pcolorm(lat_grid,lon_grid,A_Es_perc')
load coast
plotm(lat,long,'w') %draw the coastlines
load('/zhome/e8/9/144512/Desktop/MScThesis/MatlabScripts/SavedMatlabVariables/MagneticEquator.mat')
plotm(lat(LatIndex),lon,'w-') %plotting the magnetic equator
h = colorbar;
ylabel(h, '%')
colormap(jet)
caxis([0 50])
%title({['Global Es scintillations ' num2str(year) '.'],['Total Occurrence: ' num2str((100*sum(TotalEs_present)/length(TotalEs_present)),3) '%.']})
title({['Global Es scintillations ' num2str(year) ' ' season '.'],['Total Occurrence: ' num2str((100*sum(TotalEs_present)/length(TotalEs_present)),3) '%.']})


% %saving figure
% baseFileName = sprintf(['Es' num2str(year) num2str(season) '.png'],1);
% % Specify some particular, specific folder:
% fullFileName = fullfile(['/zhome/e8/9/144512/Desktop/MScThesis/Figures/GlobalMaps/', num2str(year), '/seasons/'], baseFileName);  
% print(figure(12),fullFileName);

%%  Making global maps of Scintillation amplitude
%The maps only takes into account measurements with maxS4 > 0.3
%(because if majority is below 0.3, we don't see much difference in mean global amplitude between years,
%better to have only take into account scintillations above 0.3 so we
%better see how the amplitude differs with solar activity between the years).

minlat = -90;     %adjust these if needed. They determine where the bins start
minlong = -180;
gridspacing = 1; %spatial resolution in degrees
latidx = 1 + floor((TotalRO_lat - minlat) ./ gridspacing);
longidx = 1 + floor((TotalRO_lon - minlong) ./ gridspacing);

Index_Above03 = TotalRO_MaxS4 > 0; %index of measurements that have S4max above 0.3

A = accumarray( [latidx(Index_Above03)', longidx(Index_Above03)'], 1); %number of measurements above 0.3 within each bin
AmpS4 = accumarray( [latidx(Index_Above03)', longidx(Index_Above03)'], TotalRO_MaxS4(Index_Above03)'); %sum of TotalRO_MaxS4 > 0.3 within each bin

AmpS4mean = AmpS4./A; %the mean of the max S4 amplitude within each bin


theta = [-179:180];
phi = [-89:90];
%obs % 360x180 2d array of data
[lat_grid,lon_grid] = meshgrid(phi,theta); % should both be 360x180 to match data
load coast % loads lat and long variables that define the coastline
figure(13); worldmap('World'); % also try axesm as it gives more options
pcolorm(lat_grid,lon_grid,AmpS4mean');
plotm(lat,long,'w-') %draw the coastlines
load('/zhome/e8/9/144512/Desktop/MScThesis/MatlabScripts/SavedMatlabVariables/MagneticEquator.mat')
plotm(lat(LatIndex),lon,'w-') %plotting the magnetic equator
h = colorbar;
ylabel(h, 'S4')
colormap(jet)
caxis([0 1])
title({['Global S4 amplitude ' num2str(year) ' ' season '.'],['Global mean amplitude: ' num2str((sum(TotalRO_MaxS4(Index_Above03))/length(TotalRO_MaxS4(Index_Above03))),4) '.']})
%The S4 amplitude is the observed maximum amplitude of a profile at h_t > 20km.











