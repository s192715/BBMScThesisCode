%COMBINING MONTHLY DATA
%Edit MyDataDir, path to directory containing the calculated variables for each month. 
%Outputs cell matrix, for each variable there is a cell element. 
%ALSO includes some plotting routines for MonthlyTimeSeries i.e. to investigate function of
%solar activity

MyDataDir = dir("/zhome/e8/9/144512/Desktop/MScThesis/MatlabScripts/SavedMatlabVariables/MonthsMainRunAll/CombineMonths");
MyDataDirString = "/zhome/e8/9/144512/Desktop/MScThesis/MatlabScripts/SavedMatlabVariables/MonthsMainRunAll/CombineMonths";

CombinedData = cell(1,22); %allocating the nubmer of variables that are combined
CombinedFileStampCell = cell(1,1);

%Following variables are stored in cells in CombinedData:
% TotalDelta_H
% TotalEventHeight
% TotalEventLat
% TotalEventLon
% TotalEventTime
% TotalEventAbsTime
% n  ,a number that tracks from which file an event occured
% TotalEventAzi
% TotalEventLocalTime
% 
% For RO measurement statistics:
% TotalEs_present
% TotalNumberOfEvents ,shows number of Es for each RO measurement
% TotalS4_present
% TotalBoth_present
% TotalNeither_present
% TotalOnlyEs_present 
% TotalOnlyS4_present

% TotalRO_lat
% TotalRO_lon
% TotalRO_Time
% TotalRO_MeanS4
% TotalRO_MaxS4    
% TotalRO_HeightOfS4Max

TotalRO_MonthNumber = [];
TotalEs_MonthNumber = [];
for k = 3:length(MyDataDir) %for each directory in MyDataDir (skips first 2 elements because they are just . and ..)

    Month = dir(MyDataDirString + "/" + MyDataDir(k).('name'));
    MonthString = MyDataDirString + "/" + MyDataDir(k).('name');
    load(MonthString)
    
    
    %Combines all the data into the variables, just like the program has
    %been run for more months at a time.
    CombinedData{1,1} = [CombinedData{1,1} TotalDelta_H];
    CombinedData{1,2} = [CombinedData{1,2} TotalEventHeight];
    CombinedData{1,3} = [CombinedData{1,3} TotalEventLat];
    CombinedData{1,4} = [CombinedData{1,4} TotalEventLon];
    CombinedData{1,5} = [CombinedData{1,5} TotalEventTime];
    CombinedData{1,6} = [CombinedData{1,6} TotalEventAbsTime];
    CombinedData{1,7} = [CombinedData{1,7} n];
    CombinedData{1,8} = [CombinedData{1,8} TotalEventAzi];
    CombinedData{1,21} = [CombinedData{1,21} TotalEventLocalTime];
    
    CombinedData{1,9} = [CombinedData{1,9} TotalEs_present];
    CombinedData{1,10} = [CombinedData{1,10} TotalNumberOfEvents];
    CombinedData{1,11} = [CombinedData{1,11} TotalS4_present];
    CombinedData{1,12} = [CombinedData{1,12} TotalBoth_present];
    CombinedData{1,13} = [CombinedData{1,13} TotalNeither_present];
    CombinedData{1,14} = [CombinedData{1,14} TotalOnlyEs_present];
    CombinedData{1,15} = [CombinedData{1,15} TotalOnlyS4_present];
    
    CombinedData{1,16} = [CombinedData{1,16} TotalRO_lat];
    CombinedData{1,17} = [CombinedData{1,17} TotalRO_lon];
    CombinedData{1,18} = [CombinedData{1,18} TotalRO_Time];
    CombinedData{1,19} = [CombinedData{1,19} TotalRO_MeanS4];
    CombinedData{1,20} = [CombinedData{1,20} TotalRO_MaxS4];
    CombinedData{1,22} = [CombinedData{1,22} TotalRO_HeightOfS4Max];
    
    CombinedFileStampCell{1,1} = [CombinedFileStampCell{1,1} FilestampCell];
    %Here, I can make variables if I want monthly statistics (i.e. take
    %some averages of each month).
    
    TotalRO_MonthNumber = [TotalRO_MonthNumber (k-2).*ones(1,length(TotalRO_lat))]; %each measurement has assigned a month number.
    
    TotalEs_MonthNumber = [TotalEs_MonthNumber (k-2).*ones(1,length(TotalEventLat))]; %each measurement has assigned a month number.
    
    EsMeanThickness(k-2) = mean(TotalDelta_H);
    EsMeanHeight(k-2) = mean(TotalEventHeight);
    EsMeanLocalTime(k-2) = mean(TotalEventLocalTime);
    EsMeanLatitude(k-2) = mean(TotalEventLat);
    EsMeanLongitude(k-2) = mean(TotalEventLon);
    
    
    %extracting for a certain local time
%     LocalTimeMorningIndex = find(TotalEventLocalTime < 12 & TotalEventLocalTime > 8);
%     LocalTimeEveningIndex = find(TotalEventLocalTime < 23 & TotalEventLocalTime > 19);
    
    
     %to select occurrence peaks with local time and height limits------------------
%     LocalTimeMorningIndex1Time = find(TotalEventLocalTime < 13 & TotalEventLocalTime > 8);
%     LocalTimeMorningIndex2Time = find(TotalEventLocalTime < 19 & TotalEventLocalTime > 13);
%     
%     LocalTimeMorningIndex1Height = find(TotalEventHeight > 80);
%     LocalTimeMorningIndex2Height = find(TotalEventHeight> 80 & TotalEventHeight < 90);
%     
%     LocalTimeMorningIndex1 = intersect(LocalTimeMorningIndex1Time, LocalTimeMorningIndex1Height);
%     LocalTimeMorningIndex2 = intersect(LocalTimeMorningIndex2Time, LocalTimeMorningIndex2Height);
%     LocalTimeMorningIndex = [LocalTimeMorningIndex1 LocalTimeMorningIndex2]; 
%        
%     LocalTimeEveningIndexTime = find(TotalEventLocalTime < 23 & TotalEventLocalTime > 16);
%     LocalTimeEveningIndexHeight = find(TotalEventHeight > 95);
%     LocalTimeEveningIndex = intersect(LocalTimeEveningIndexTime,LocalTimeEveningIndexHeight);

    
    %Just with height limits
%     LocalTimeMorningIndex = find(TotalEventHeight > 80);
%     LocalTimeEveningIndex = find(TotalEventHeight < 80);
    %-------------------------
    
    %To compare morning and evening below 80 km    
    LocalTimeMorningIndexTime = find(TotalEventLocalTime < 14 & TotalEventLocalTime > 8);
    LocalTimeEveningIndexTime = find(TotalEventLocalTime < 23 & TotalEventLocalTime > 16);
    LocalTimeMorningIndexHeight = find(TotalEventHeight < 80);
    LocalTimeEveningIndexHeight = find(TotalEventHeight < 80);
    
    LocalTimeMorningIndex = intersect(LocalTimeMorningIndexTime,LocalTimeMorningIndexHeight);
    LocalTimeEveningIndex = intersect(LocalTimeEveningIndexTime,LocalTimeEveningIndexHeight);
    %----------------------------------------
    
    
    EsMeanThicknessMorning(k-2) = mean(TotalDelta_H(LocalTimeMorningIndex));
    EsMeanHeightMorning(k-2) = mean(TotalEventHeight(LocalTimeMorningIndex));
    EsMeanLocalTimeMorning(k-2) = mean(TotalEventLocalTime(LocalTimeMorningIndex));
    EsMeanLatitudeMorning(k-2) = mean(TotalEventLat(LocalTimeMorningIndex));
    EsMeanLongitudeMorning(k-2) = mean(TotalEventLon(LocalTimeMorningIndex));
    
    EsMeanThicknessEvening(k-2) = mean(TotalDelta_H(LocalTimeEveningIndex));
    EsMeanHeightEvening(k-2) = mean(TotalEventHeight(LocalTimeEveningIndex));
    EsMeanLocalTimeEvening(k-2) = mean(TotalEventLocalTime(LocalTimeEveningIndex));
    EsMeanLatitudeEvening(k-2) = mean(TotalEventLat(LocalTimeEveningIndex));
    EsMeanLongitudeEvening(k-2) = mean(TotalEventLon(LocalTimeEveningIndex));
    
end

TotalDelta_H = CombinedData{1,1};
TotalEventHeight = CombinedData{1,2};
TotalEventLat = CombinedData{1,3};
TotalEventLon = CombinedData{1,4};
TotalEventTime = CombinedData{1,5};
TotalEventAbsTime = CombinedData{1,6};
n = CombinedData{1,7};
TotalEventAzi = CombinedData{1,8};
TotalEventLocalTime = CombinedData{1,21};


TotalEs_present = CombinedData{1,9};
TotalNumberOfEvents = CombinedData{1,10};
TotalS4_present = CombinedData{1,11};
TotalBoth_present = CombinedData{1,12};
TotalNeither_present = CombinedData{1,13};
TotalOnlyEs_present = CombinedData{1,14};
TotalOnlyS4_present = CombinedData{1,15};

TotalRO_lat = CombinedData{1,16};
TotalRO_lon = CombinedData{1,17};
TotalRO_Time = CombinedData{1,18};
TotalRO_MeanS4 = CombinedData{1,19};
TotalRO_MaxS4  = CombinedData{1,20};
TotalRO_HeightOfS4Max = CombinedData{1,22};

%converting Time to LocalTime
TotalROTimeUTC  = datevec(datetime(TotalRO_Time,'ConvertFrom','juliandate'));
TotalROTimeUTC_HourOfDay = 24*datenum(hours(TotalROTimeUTC(:,4)) + minutes(TotalROTimeUTC(:,5)) + seconds(TotalROTimeUTC(:,6)));

[zd] = timezone(TotalRO_lon);
TotalROLocalTime = (TotalROTimeUTC_HourOfDay - zd)';
for i = 1:length(TotalROLocalTime)
    if TotalROLocalTime(i) < 0
        TotalROLocalTime(i) = 24 + TotalROLocalTime(i);
    elseif TotalROLocalTime(i) >= 24
        TotalROLocalTime(i) = TotalROLocalTime(i) - 24;
    end
end
clear TotalROTimeUTC
clear TotalROTimeUTC_HourOfDay
clear zd


%-------------------------------------------
%First, sunspot number plot
%data from http://www.sidc.be/silso/newdataset

%SNdata = readmatrix("/zhome/e8/9/144512/Desktop/MScThesis/DataMisc/SunspotNumber2007_2013.txt");
SNdata = readmatrix("/zhome/e8/9/144512/Desktop/MScThesis/DataMisc/SunspotNumber2007_2019.txt");

DecimalYear = SNdata(:,3);
SNnumber = SNdata(:,4);
SNerror = SNdata(:,5);
%figure(10); errorbar(DecimalYear,SNnumber,SNerror)
%xlabel('year')
%ylabel('SN number')
%title('Monthly average SN number')
%hold on;
RM = movmean(SNnumber,12); %12 month running mean
%plot(DecimalYear,RM,'r-')
%grid on;
%legend('SN data', '12 month running mean')
%xline(2007,'handlevisibility','off'); xline(2008,'handlevisibility','off'); xline(2009,'handlevisibility','off'); xline(2010,'handlevisibility','off');
%xline(2011,'handlevisibility','off'); xline(2012,'handlevisibility','off'); xline(2013,'handlevisibility','off'); xline(2014,'handlevisibility','off'); 
%xline(2015,'handlevisibility','off'); xline(2016,'handlevisibility','off'); xline(2017,'handlevisibility','off'); xline(2018,'handlevisibility','off'); xline(2019,'handlevisibility','off'); 
%-------------------------------------------------------

%% 



minlat = -90;     %adjust these if needed. They determine where the bins start
mintime = 1;
gridspacing_lat = 1;
gridspacing_month = 1; %one month
latidx = 1 + floor((TotalRO_lat - minlat) ./ gridspacing_lat);
timedx = 1 + floor((TotalRO_MonthNumber- mintime) ./ gridspacing_month);

A2 = accumarray( [latidx(:), timedx(:)], 1); %number of measurements within each bin
A2_S4 = accumarray( [latidx(:), timedx(:)], TotalS4_present'); %number of S4 detections within bin
A2_Es = accumarray( [latidx(:), timedx(:)], TotalEs_present'); %number of Es detections within bin

A2_S4_perc = 100*(A2_S4./A2);
A2_Es_perc = 100*(A2_Es./A2);

%IF I want to smoothen data
K = (1/9)*ones(3); %defining kernel
A2_S4_perc = conv2(A2_S4_perc,K,'same');
A2_Es_perc = conv2(A2_Es_perc,K,'same');

figure(1); imagesc(A2_S4_perc) %imagesc(Zsmooth1)
h = colorbar;
ylabel(h, '%')
colormap(jet)
caxis([0 100])
xlabel('Year')
ylabel('Latitude [°]')

yticks = 0:length(A2(:,1))/6:length(A2(:,1)); %adjust as appropriate, positive integers only
yticks(1) = 1;
ylabels = {'-90' '-60' '-30' '0' '30' '60' '90'};  %time labels
set(gca, 'YTick', yticks, 'YTickLabel', ylabels);
set(gca,'YDir','normal')
%xticks = [1 4:4:84]; %adjust as appropriate, positive integers only
xticks = 1:12:156;
xlabels = 2007:2019; 
set(gca, 'XTick', xticks, 'XTickLabel', xlabels);
xtickangle(45)
title(['S4 occurrence vs. latitude 2007-2019'])

%adding scaled monthly SNnumber which denotes solar activity
hold on;
plot(1:length(RM),RM*1.5,'w-','Linewidth',1.5)




figure(2); imagesc(A2_Es_perc)
h = colorbar;
ylabel(h, '%')
colormap(jet)
%caxis([0 100])
caxis([0 50])
xlabel('Year')
ylabel('Latitude [°]')

yticks = 0:length(A2(:,1))/6:length(A2(:,1)); %adjust as appropriate, positive integers only
yticks(1) = 1;
ylabels = {'-90' '-60' '-30' '0' '30' '60' '90'};  %time labels
set(gca, 'YTick', yticks, 'YTickLabel', ylabels);
set(gca,'YDir','normal')
%xticks = [1 4:4:84]; %adjust as appropriate, positive integers only
xticks = 1:12:156;
xlabels = 2007:2019; 
set(gca, 'XTick', xticks, 'XTickLabel', xlabels);
xtickangle(45)
title(['Es occurrence vs. latitude 2007-2019'])

%adding scaled monthly SNnumber which denotes solar activity
hold on;
plot(1:length(RM),RM*1.5,'w-','Linewidth',1.5)


Index_Above03 = TotalRO_MaxS4 > 0.3; %index of measurements that have S4max above 0.3

A = accumarray( [latidx(Index_Above03)', timedx(Index_Above03)'], 1); %number of measurements above 0.3 within each bin
AmpS4 = accumarray( [latidx(Index_Above03)', timedx(Index_Above03)'], TotalRO_MaxS4(Index_Above03)'); %sum of TotalRO_MaxS4 > 0.3 within each bin

AmpS4mean = AmpS4./A; %the mean of the max S4 amplitude within each bin


%if i want to smoothen the data
AmpS4mean = conv2(AmpS4mean,K,'same');

figure(3); imagesc(AmpS4mean) %imagesc(Zsmooth1)
h = colorbar;
ylabel(h, 'S4')
colormap(jet)
caxis([0.3 1.3])
xlabel('Year')
ylabel('Latitude [°]')

yticks = 0:length(A2(:,1))/6:length(A2(:,1)); %adjust as appropriate, positive integers only
yticks(1) = 1;
ylabels = {'-90' '-60' '-30' '0' '30' '60' '90'};  %time labels
set(gca, 'YTick', yticks, 'YTickLabel', ylabels);
set(gca,'YDir','normal')
%xticks = [1 4:4:84]; %adjust as appropriate, positive integers only
xticks = 1:12:156;
xlabels = 2007:2019; 
set(gca, 'XTick', xticks, 'XTickLabel', xlabels);
xtickangle(45)
title(['S4 amplitude vs. latitude 2007-2019'])

%adding scaled monthly SNnumber which denotes solar activity
hold on;
plot(1:length(RM),RM*1.5,'w-','Linewidth',1.5)




%S4 map with all measurements included, not just those with S4 above 0.3
A = accumarray( [latidx(:), timedx(:)], 1); %number of measurements within each bin
AmpS4 = accumarray( [latidx(:), timedx(:)], TotalRO_MaxS4'); %sum of TotalRO_MaxS4 within each bin

AmpS4mean = AmpS4./A; %the mean of the max S4 amplitude within each bin

%if i want to smoothen the data:
AmpS4mean = conv2(AmpS4mean,K,'same');

figure(4); imagesc(AmpS4mean)
h = colorbar;
ylabel(h, 'S4')
colormap(jet)
caxis([0 1])
xlabel('Year')
ylabel('Latitude [°]')

yticks = 0:length(A2(:,1))/6:length(A2(:,1)); %adjust as appropriate, positive integers only
yticks(1) = 1;
ylabels = {'-90' '-60' '-30' '0' '30' '60' '90'};  %time labels
set(gca, 'YTick', yticks, 'YTickLabel', ylabels);
set(gca,'YDir','normal')
%xticks = [1 4:4:84]; %adjust as appropriate, positive integers only
xticks = 1:12:156;
xlabels = 2007:2019; 
set(gca, 'XTick', xticks, 'XTickLabel', xlabels);
xtickangle(45)
title(['S4 amplitude vs. latitude 2007-2019'])

%adding scaled monthly SNnumber which denotes solar activity
hold on;
plot(1:length(RM),RM*1.5,'w-','Linewidth',1.5)


%% plotting at a certain local times only


%variables use this index to only include local times between 19 and 23
%LocalTimeIndex = find(TotalROLocalTime >= 0 & TotalROLocalTime <= 24);
%LocalTimeIndex = find(TotalROLocalTime >= 19 & TotalROLocalTime <= 23);
LocalTimeIndex = find(TotalROLocalTime >= 8 & TotalROLocalTime <= 12);
LT = '8-12';

minlat = -90;     %adjust these if needed. They determine where the bins start
mintime = 1;
gridspacing_lat = 1;
gridspacing_month = 1; %one month
latidx = 1 + floor((TotalRO_lat(LocalTimeIndex) - minlat) ./ gridspacing_lat);
timedx = 1 + floor((TotalRO_MonthNumber(LocalTimeIndex) - mintime) ./ gridspacing_month);

%S4 map with all measurements included, not just those with S4 above 0.3
A = accumarray( [latidx(:), timedx(:)], 1); %number of measurements within each bin
AmpS4 = accumarray( [latidx(:), timedx(:)], TotalRO_MaxS4(LocalTimeIndex)'); %sum of TotalRO_MaxS4 within each bin

AmpS4mean = AmpS4./A; %the mean of the max S4 amplitude within each bin

% if i want to smoothen the data:
K = (1/9)*ones(3); %defining kernel
AmpS4mean = conv2(AmpS4mean,K,'same');

figure(5); imagesc(AmpS4mean)
h = colorbar;
ylabel(h, 'S4')
colormap(jet)
caxis([0 1])
xlabel('Year')
ylabel('Latitude [°]')

yticks = 0:length(A2(:,1))/6:length(A2(:,1)); %adjust as appropriate, positive integers only
yticks(1) = 1;
ylabels = {'-90' '-60' '-30' '0' '30' '60' '90'};  %time labels
set(gca, 'YTick', yticks, 'YTickLabel', ylabels);
set(gca,'YDir','normal')
%xticks = [1 4:4:84]; %adjust as appropriate, positive integers only
xticks = 1:12:156;
xlabels = 2007:2019; 
set(gca, 'XTick', xticks, 'XTickLabel', xlabels);
xtickangle(45)
title(['S4 amplitude vs. latitude 2007-2019. LT: ' LT])

%adding scaled monthly SNnumber which denotes solar activity
hold on;
plot(1:length(RM),RM*1.5,'w-','Linewidth',1.5)


%-------------For S4 and Es occurrece-------------


A2 = accumarray( [latidx(:), timedx(:)], 1); %number of measurements within each bin
A2_S4 = accumarray( [latidx(:), timedx(:)], TotalS4_present(LocalTimeIndex)'); %number of S4 detections within bin
A2_Es = accumarray( [latidx(:), timedx(:)], TotalEs_present(LocalTimeIndex)'); %number of Es detections within bin

A2_S4_perc = 100*(A2_S4./A2);
A2_Es_perc = 100*(A2_Es./A2);

%IF I want to smoothen data
K = (1/9)*ones(3); %defining kernel
A2_S4_perc = conv2(A2_S4_perc,K,'same');
A2_Es_perc = conv2(A2_Es_perc,K,'same');

figure(6); imagesc(A2_S4_perc) %imagesc(Zsmooth1)
h = colorbar;
ylabel(h, '%')
colormap(jet)
caxis([0 100])
%caxis([-20 20])
xlabel('Year')
ylabel('Latitude [°]')

yticks = 0:length(A2(:,1))/6:length(A2(:,1)); %adjust as appropriate, positive integers only
yticks(1) = 1;
ylabels = {'-90' '-60' '-30' '0' '30' '60' '90'};  %time labels
set(gca, 'YTick', yticks, 'YTickLabel', ylabels);
set(gca,'YDir','normal')
%xticks = [1 4:4:84]; %adjust as appropriate, positive integers only
xticks = 1:12:156;
xlabels = 2007:2019; 
set(gca, 'XTick', xticks, 'XTickLabel', xlabels);
xtickangle(45)
title(['S4 occurrence vs. latitude 2007-2019. LT: ' LT])
%title(['S4 occurrence difference: Morning - evening.'])

%adding scaled monthly SNnumber which denotes solar activity
hold on;
plot(1:length(RM),RM*1.5,'w-','Linewidth',1.5)


figure(7); imagesc(A2_Es_perc)
h = colorbar;
ylabel(h, '%')
colormap(jet)
%caxis([0 100])
caxis([0 50])
%caxis([-8 8])
xlabel('Year')
ylabel('Latitude [°]')

yticks = 0:length(A2(:,1))/6:length(A2(:,1)); %adjust as appropriate, positive integers only
yticks(1) = 1;
ylabels = {'-90' '-60' '-30' '0' '30' '60' '90'};  %time labels
set(gca, 'YTick', yticks, 'YTickLabel', ylabels);
set(gca,'YDir','normal')
%xticks = [1 4:4:84]; %adjust as appropriate, positive integers only
xticks = 1:12:156;
xlabels = 2007:2019; 
set(gca, 'XTick', xticks, 'XTickLabel', xlabels);
xtickangle(45)
title(['Es occurrence vs. latitude 2007-2019. LT: ' LT])
%title(['Es occurrence difference: Morning - evening.'])

%adding scaled monthly SNnumber which denotes solar activity
hold on;
plot(1:length(RM),RM*1.5,'w-','Linewidth',1.5)

%% Plotting COSMIC data coverage in LocalTime and Time

figure(10); plot(TotalRO_Time,TotalRO_lat,'k.','MarkerSize', 0.01)
xlabel('Year')
ylabel('Latitude [°]')
%title('COSMIC data coverage')
xticks = TotalRO_Time(1):(TotalRO_Time(end)-TotalRO_Time(1))/13:TotalRO_Time(end);
xlabels = 2007:2020;
set(gca, 'XTick', xticks, 'XTickLabel', xlabels);
xtickangle(45)


%%Good for showing that the max S4 is less in spring and fall and peaks in
%%summer and winter (i.e. when summer is on either hemisphere).
% figure(11); plot(TotalRO_Time,TotalRO_MaxS4,'k.','MarkerSize', 0.1)
% xlabel('Year')
% ylabel('Max S4')
% title('COSMIC data coverage')
% xticks = TotalRO_Time(1):(TotalRO_Time(end)-TotalRO_Time(1))/13:TotalRO_Time(end);
% xlabels = 2007:2020;
% set(gca, 'XTick', xticks, 'XTickLabel', xlabels);


%% Es mean thickness and mean height after months

%EsMeanThickness = movmean(EsMeanThickness,3); %3 month running mean
%EsMeanHeight = movmean(EsMeanHeight,3); %Good to look at the 3 month running mean (is it changes with seasons, thus 3 months), shows more clearly the effects.


figure(20);
subplot(2,1,1); plot(unique(TotalRO_MonthNumber),EsMeanThickness,'b.-')
xlabel('Year')
ylabel('Es thickness [km]')
xticks = 0:12:156
xlabels = 2007:2020;
set(gca, 'XTick', xticks, 'XTickLabel', xlabels);
xtickangle(45)
ylim([1.8 2.2])
hold on;
SummerStart = 5:12:156;
SummerStop = 8:12:156;
WinterStart = 11:12:154; %excluding first winter in 2007
WinterStop = 14:12:156;
for i = 1:length(SummerStart)
    patch([SummerStart(i) SummerStart(i) SummerStop(i) SummerStop(i)],[1.4 2.2 2.2 1.4],'y','EdgeColor','y','FaceAlpha',0.25,'EdgeAlpha',0.25)
end
for j = 1:length(WinterStart)
    patch([WinterStart(j) WinterStart(j) WinterStop(j) WinterStop(j)],[1.4 2.2 2.2 1.4],'c','EdgeColor','c','FaceAlpha',0.25,'EdgeAlpha',0.25)
end
patch([0 0 2 2],[1.41 2.19 2.19 1.41],'c','EdgeColor','c','FaceAlpha',0.25,'EdgeAlpha',0.25)


subplot(2,1,2); plot(unique(TotalRO_MonthNumber),EsMeanHeight,'b.-')
xlabel('Year')
ylabel('Es altitude [km]')
xticks = 0:12:156;
xlabels = 2007:2020;
set(gca, 'XTick', xticks, 'XTickLabel', xlabels);
xtickangle(45)
hold on;
%adding scaled monthly SNnumber (12 month running mean) which denotes solar activity
hold on;
plot(1:length(RM),RM*0.1 + mean(EsMeanHeight)-4,'r-')

%adding patches
ylim([94.62,107])
suptitle('Monthly means of Es thickness and altitude')

SummerStart = 5:12:156;
SummerStop = 8:12:156;
WinterStart = 11:12:154; %excluding first winter in 2007
WinterStop = 14:12:156;
for i = 1:length(SummerStart)
    patch([SummerStart(i) SummerStart(i) SummerStop(i) SummerStop(i)],[94.63 106.99 106.99 94.63],'y','EdgeColor','y','FaceAlpha',0.25,'EdgeAlpha',0.25)
end
for j = 1:length(WinterStart)
    patch([WinterStart(j) WinterStart(j) WinterStop(j) WinterStop(j)],[94.63 106.99 106.99 94.63],'c','EdgeColor','c','FaceAlpha',0.25,'EdgeAlpha',0.25)
end
patch([0 0 2 2],[94.63 106.99 106.99 94.63],'c','EdgeColor','c','FaceAlpha',0.25,'EdgeAlpha',0.25)


%Fluctuations are due to the seasons, 

%Thickness does not seem to be affected with solar activity. The low Es
%thickness in September 2019 is due to sparse data collected, where there were only measurements
%from one day in September 2019.

%-Also, yellow patches mark summers and blue patches mark winter, altitude seems to be lower during summer and winter.


%% Thickness and height plots monthly means, at morning and evening local times.

% %---------If I want to smooth the data--------
% %doing a 3 month running mean because altitude and thickness varies with
% %seasons.
% EsMeanHeightMorning = movmean(EsMeanHeightMorning,3);
% EsMeanHeightEvening = movmean(EsMeanHeightEvening,3);
% EsMeanThicknessMorning = movmean(EsMeanThicknessMorning,3);
% EsMeanThicknessEvening = movmean(EsMeanThicknessEvening,3);
% %-----------------------------



figure(30);
%subplot(2,1,1);
plot(unique(TotalRO_MonthNumber),EsMeanHeightMorning,'b.-')
xlabel('Year')
ylabel('Es altitude [km]')
xticks = 0:12:156;
xlabels = 2007:2020;
set(gca, 'XTick', xticks, 'XTickLabel', xlabels);
xtickangle(45)
hold on;
%adding scaled monthly SNnumber (12 month running mean) which denotes solar activity
plot(1:length(RM),RM*0.1 + mean(EsMeanHeightMorning(1:50))-6,'r-','Handlevisibility','Off')

%adding patches
%ylim([96 108])
ylim([88 108])
%ylim([50 75])

SummerStart = 5:12:156;
SummerStop = 8:12:156;
WinterStart = 11:12:154; %excluding first winter in 2007
WinterStop = 14:12:156;
for i = 1:length(SummerStart)
    patch([SummerStart(i) SummerStart(i) SummerStop(i) SummerStop(i)],[88 108 108 88],'y','EdgeColor','y','FaceAlpha',0.25,'EdgeAlpha',0.25,'Handlevisibility','Off')
end
for j = 1:length(WinterStart)
    patch([WinterStart(j) WinterStart(j) WinterStop(j) WinterStop(j)],[88 108 108 88],'c','EdgeColor','c','FaceAlpha',0.25,'EdgeAlpha',0.25,'Handlevisibility','Off')
end
patch([0 0 2 2],[88 108 108 88],'c','EdgeColor','c','FaceAlpha',0.25,'EdgeAlpha',0.25,'Handlevisibility','Off')
plot(unique(TotalRO_MonthNumber),EsMeanHeightEvening,'k.-')
%legend('Morning 8-12','Evening 19-23')
legend('Morning','Evening')


 %to split the plots in subplot format:
% subplot(2,1,2); plot(unique(TotalRO_MonthNumber),EsMeanHeightEvening,'b.-')
% xlabel('Year')
% ylabel('Es altitude [km]')
% xticks = 0:12:156;
% xlabels = 2007:2020;
% set(gca, 'XTick', xticks, 'XTickLabel', xlabels);
% xtickangle(45)
% hold on;
% %adding scaled monthly SNnumber (12 month running mean) which denotes solar activity
% plot(1:length(RM),RM*0.15 + mean(EsMeanHeightMorning(1:50))-11,'r-')
% 
% %adding patches
% %ylim([94.62,107])
% ylim([85 108])
% suptitle('Monthly means of Es altitude at local times 8-12 (above) and 19-23 (below)')
% 
% SummerStart = 5:12:156;
% SummerStop = 8:12:156;
% WinterStart = 11:12:154; %excluding first winter in 2007
% WinterStop = 14:12:156;
% for i = 1:length(SummerStart)
%     patch([SummerStart(i) SummerStart(i) SummerStop(i) SummerStop(i)],[85 108 108 85],'y','EdgeColor','y','FaceAlpha',0.25,'EdgeAlpha',0.25)
% end
% for j = 1:length(WinterStart)
%     patch([WinterStart(j) WinterStart(j) WinterStop(j) WinterStop(j)],[85 108 108 85],'c','EdgeColor','c','FaceAlpha',0.25,'EdgeAlpha',0.25)
% end
% patch([0 0 2 2],[85 108 108 85],'c','EdgeColor','c','FaceAlpha',0.25,'EdgeAlpha',0.25)
% 
% 
% 
% 
% %-----------------------thickness----------------
figure(31);
%subplot(2,1,1); 
plot(unique(TotalRO_MonthNumber),EsMeanThicknessMorning,'b.-')
xlabel('Year')
ylabel('Es thickness [km]')
xticks = 0:12:156;
xlabels = 2007:2020;
set(gca, 'XTick', xticks, 'XTickLabel', xlabels);
xtickangle(45)
hold on;
%adding scaled monthly SNnumber (12 month running mean) which denotes solar activity
plot(1:length(RM),RM*0.005 + mean(EsMeanThicknessMorning(1:50)) - 0.1,'r-','Handlevisibility','Off')
plot(unique(TotalRO_MonthNumber),EsMeanThicknessEvening,'k.-')
%legend('Morning 8-12','Evening 19-23')
legend('Morning','Evening')

%adding patches
ylim([1.6 2.4])
%ylim([1.4 2.8])

SummerStart = 5:12:156;
SummerStop = 8:12:156;
WinterStart = 11:12:154; %excluding first winter in 2007
WinterStop = 14:12:156;
for i = 1:length(SummerStart)
    patch([SummerStart(i) SummerStart(i) SummerStop(i) SummerStop(i)],[1.6 2.4 2.4 1.6],'y','EdgeColor','y','FaceAlpha',0.25,'EdgeAlpha',0.25,'Handlevisibility','Off')
end
for j = 1:length(WinterStart)
    patch([WinterStart(j) WinterStart(j) WinterStop(j) WinterStop(j)],[1.6 2.4 2.4 1.6],'c','EdgeColor','c','FaceAlpha',0.25,'EdgeAlpha',0.25,'Handlevisibility','Off')
end
patch([0 0 2 2],[1.6 2.4 2.4 1.6],'c','EdgeColor','c','FaceAlpha',0.25,'EdgeAlpha',0.25,'Handlevisibility','Off')

 %To split the plot in subplot format:
% subplot(2,1,2); plot(unique(TotalRO_MonthNumber),EsMeanThicknessEvening,'b.-')
% xlabel('Year')
% ylabel('Es thickness [km]')
% xticks = 0:12:156;
% xlabels = 2007:2020;
% set(gca, 'XTick', xticks, 'XTickLabel', xlabels);
% xtickangle(45)
% hold on;
% %adding scaled monthly SNnumber (12 month running mean) which denotes solar activity
% plot(1:length(RM),RM*0.005 + mean(EsMeanThicknessMorning(1:50)) - 0.1,'r-')
% 
% %adding patches
% %ylim([94.62,107])
% ylim([1.6 2.4])
% suptitle('Monthly means of Es thickness at local times 8-12 (above) and 19-23 (below)')
% 
% SummerStart = 5:12:156;
% SummerStop = 8:12:156;
% WinterStart = 11:12:154; %excluding first winter in 2007
% WinterStop = 14:12:156;
% for i = 1:length(SummerStart)
%     patch([SummerStart(i) SummerStart(i) SummerStop(i) SummerStop(i)],[1.6 2.4 2.4 1.6],'y','EdgeColor','y','FaceAlpha',0.25,'EdgeAlpha',0.25)
% end
% for j = 1:length(WinterStart)
%     patch([WinterStart(j) WinterStart(j) WinterStop(j) WinterStop(j)],[1.6 2.4 2.4 1.6],'c','EdgeColor','c','FaceAlpha',0.25,'EdgeAlpha',0.25)
% end
% patch([0 0 2 2],[1.6 2.4 2.4 1.6],'c','EdgeColor','c','FaceAlpha',0.25,'EdgeAlpha',0.25)




%% Es monthly mean lat,lon and local time as function of time i.e. solar activity.

figure(25);
subplot(3,1,1); plot(unique(TotalRO_MonthNumber),EsMeanLatitude,'b.-')
xlabel('Year')
ylabel('Es Latitude [°]')
xticks = 0:12:156;
xlabels = 2007:2020;
set(gca, 'XTick', xticks, 'XTickLabel', xlabels);
xtickangle(45)
ylim([-40 40])
hold on;
SummerStart = 5:12:156;
SummerStop = 8:12:156;
WinterStart = 11:12:154; %excluding first winter in 2007
WinterStop = 14:12:156;
for i = 1:length(SummerStart)
    patch([SummerStart(i) SummerStart(i) SummerStop(i) SummerStop(i)],[-40 40 40 -40],'y','EdgeColor','y','FaceAlpha',0.25,'EdgeAlpha',0.25)
end
for j = 1:length(WinterStart)
    patch([WinterStart(j) WinterStart(j) WinterStop(j) WinterStop(j)],[-40 40 40 -40],'c','EdgeColor','c','FaceAlpha',0.25,'EdgeAlpha',0.25)
end
patch([0 0 2 2],[-40 40 40 -40],'c','EdgeColor','c','FaceAlpha',0.25,'EdgeAlpha',0.25)


subplot(3,1,2); plot(unique(TotalRO_MonthNumber),EsMeanLongitude,'b.-')
xlabel('Year')
ylabel('Es longitude [°]')
xticks = 0:12:156;
xlabels = 2007:2020;
set(gca, 'XTick', xticks, 'XTickLabel', xlabels);
xtickangle(45)
hold on;
%adding scaled monthly SNnumber (12 month running mean) which denotes solar activity
% hold on;
% plot(1:length(RM),RM*0.1 + mean(EsMeanHeight)-4,'r-')

%adding patches
ylim([-40 40])


SummerStart = 5:12:156;
SummerStop = 8:12:156;
WinterStart = 11:12:154; %excluding first winter in 2007
WinterStop = 14:12:156;
for i = 1:length(SummerStart)
    patch([SummerStart(i) SummerStart(i) SummerStop(i) SummerStop(i)],[-40 40 40 -40],'y','EdgeColor','y','FaceAlpha',0.25,'EdgeAlpha',0.25)
end
for j = 1:length(WinterStart)
    patch([WinterStart(j) WinterStart(j) WinterStop(j) WinterStop(j)],[-40 40 40 -40],'c','EdgeColor','c','FaceAlpha',0.25,'EdgeAlpha',0.25)
end
patch([0 0 2 2],[-40 40 40 -40],'c','EdgeColor','c','FaceAlpha',0.25,'EdgeAlpha',0.25)


subplot(3,1,3); plot(unique(TotalRO_MonthNumber),EsMeanLocalTime,'b.-')
xlabel('Year')
ylabel('Es local time [hours]')
xticks = 0:12:156;
xlabels = 2007:2020;
set(gca, 'XTick', xticks, 'XTickLabel', xlabels);
xtickangle(45)
ylim([11 19])
hold on;
SummerStart = 5:12:156;
SummerStop = 8:12:156;
WinterStart = 11:12:154; %excluding first winter in 2007
WinterStop = 14:12:156;
for i = 1:length(SummerStart)
    patch([SummerStart(i) SummerStart(i) SummerStop(i) SummerStop(i)],[11 19 19 11],'y','EdgeColor','y','FaceAlpha',0.25,'EdgeAlpha',0.25)
end
for j = 1:length(WinterStart)
    patch([WinterStart(j) WinterStart(j) WinterStop(j) WinterStop(j)],[11 19 19 11],'c','EdgeColor','c','FaceAlpha',0.25,'EdgeAlpha',0.25)
end
patch([0 0 2 2],[11 19 19 11],'c','EdgeColor','c','FaceAlpha',0.25,'EdgeAlpha',0.25)
suptitle('Monthly means of Es latitude, longitude, and local time')


%finding means of all winters and all summers
% vecW = [EsMeanThickness(1:2) EsMeanThickness(12:12:156) EsMeanThickness(13:12:156) EsMeanThickness(14:12:156)];
% vecS = [EsMeanThickness(6:12:156) EsMeanThickness(7:12:156) EsMeanThickness(8:12:156)];
%mean(vecW) = 1.954 km, 98.85 km height.
%mean(vecS) = 1.893 km, 98.12 km height.


%---------If I want to smooth the data--------
%doing a 3 month running mean because altitude and thickness varies with
%seasons.
EsMeanLatitudeMorning = movmean(EsMeanLatitudeMorning,3);
EsMeanLatitudeEvening = movmean(EsMeanLatitudeEvening,3);
EsMeanLongitudeMorning = movmean(EsMeanLongitudeMorning,3);
EsMeanLongitudeEvening = movmean(EsMeanLongitudeEvening,3);
EsMeanLocalTimeMorning = movmean(EsMeanLocalTimeMorning,3);
EsMeanLocalTimeEvening = movmean(EsMeanLocalTimeEvening,3);
%-----------------------------

%for morning times 8-12 and evening times 19-23
figure(26);
subplot(3,1,1); plot(unique(TotalRO_MonthNumber),EsMeanLatitudeMorning,'b.-')
hold on;
plot(unique(TotalRO_MonthNumber),EsMeanLatitudeEvening,'k.-')
xlabel('Year')
ylabel('Es Latitude [°]')
xticks = 0:12:156;
xlabels = 2007:2020;
set(gca, 'XTick', xticks, 'XTickLabel', xlabels);
xtickangle(45)
ylim([-40 40])
hold on;
SummerStart = 5:12:156;
SummerStop = 8:12:156;
WinterStart = 11:12:154; %excluding first winter in 2007
WinterStop = 14:12:156;
for i = 1:length(SummerStart)
    patch([SummerStart(i) SummerStart(i) SummerStop(i) SummerStop(i)],[-40 40 40 -40],'y','EdgeColor','y','FaceAlpha',0.25,'EdgeAlpha',0.25,'Handlevisibility','Off')
end
for j = 1:length(WinterStart)
    patch([WinterStart(j) WinterStart(j) WinterStop(j) WinterStop(j)],[-40 40 40 -40],'c','EdgeColor','c','FaceAlpha',0.25,'EdgeAlpha',0.25,'Handlevisibility','Off')
end
patch([0 0 2 2],[-40 40 40 -40],'c','EdgeColor','c','FaceAlpha',0.25,'EdgeAlpha',0.25,'Handlevisibility','Off')
legend('Morning 8-12', 'Evening 19-23')

subplot(3,1,2); plot(unique(TotalRO_MonthNumber),EsMeanLongitudeMorning,'b.-')
hold on;
plot(unique(TotalRO_MonthNumber),EsMeanLongitudeEvening,'k.-')
xlabel('Year')
ylabel('Es longitude [°]')
xticks = 0:12:156;
xlabels = 2007:2020;
set(gca, 'XTick', xticks, 'XTickLabel', xlabels);
xtickangle(45)
hold on;
%adding scaled monthly SNnumber (12 month running mean) which denotes solar activity
% hold on;
% plot(1:length(RM),RM*0.1 + mean(EsMeanHeight)-4,'r-')

%adding patches
ylim([-40 40])


SummerStart = 5:12:156;
SummerStop = 8:12:156;
WinterStart = 11:12:154; %excluding first winter in 2007
WinterStop = 14:12:156;
for i = 1:length(SummerStart)
    patch([SummerStart(i) SummerStart(i) SummerStop(i) SummerStop(i)],[-40 40 40 -40],'y','EdgeColor','y','FaceAlpha',0.25,'EdgeAlpha',0.25)
end
for j = 1:length(WinterStart)
    patch([WinterStart(j) WinterStart(j) WinterStop(j) WinterStop(j)],[-40 40 40 -40],'c','EdgeColor','c','FaceAlpha',0.25,'EdgeAlpha',0.25)
end
patch([0 0 2 2],[-40 40 40 -40],'c','EdgeColor','c','FaceAlpha',0.25,'EdgeAlpha',0.25)


subplot(3,1,3); plot(unique(TotalRO_MonthNumber),EsMeanLocalTimeMorning,'b.-')
hold on;
plot(unique(TotalRO_MonthNumber),EsMeanLocalTimeEvening,'k.-')
xlabel('Year')
ylabel('Es local time [hours]')
xticks = 0:12:156;
xlabels = 2007:2020;
set(gca, 'XTick', xticks, 'XTickLabel', xlabels);
xtickangle(45)
ylim([8 22])
hold on;
SummerStart = 5:12:156;
SummerStop = 8:12:156;
WinterStart = 11:12:154; %excluding first winter in 2007
WinterStop = 14:12:156;
for i = 1:length(SummerStart)
    patch([SummerStart(i) SummerStart(i) SummerStop(i) SummerStop(i)],[8 22 22 8],'y','EdgeColor','y','FaceAlpha',0.25,'EdgeAlpha',0.25)
end
for j = 1:length(WinterStart)
    patch([WinterStart(j) WinterStart(j) WinterStop(j) WinterStop(j)],[8 22 22 8],'c','EdgeColor','c','FaceAlpha',0.25,'EdgeAlpha',0.25)
end
patch([0 0 2 2],[8 22 22 8],'c','EdgeColor','c','FaceAlpha',0.25,'EdgeAlpha',0.25)
suptitle('Monthly means of Es latitude, longitude, and local time, for morninng and evening')


%% Local time monthly time series



minlat = 0;     %adjust these if needed. They determine where the bins start
mintime = 1;
gridspacing_lat = 0.25;
gridspacing_month = 1; %one month
latidx = 1 + floor((TotalEventLocalTime - minlat) ./ gridspacing_lat);
timedx = 1 + floor((TotalEs_MonthNumber- mintime) ./ gridspacing_month);

A = accumarray( [latidx(:), timedx(:)], 1); 
LocalTimeOcc = zeros(size(A));
for i = 1:length(A(1,:))-12
   LocalTimeOcc(:,i) = A(:,i)./sum(A(:,i)); 
end

%IF I want to smoothen data
K = (1/9)*ones(3); %defining kernel
LocalTimeOcc = conv2(LocalTimeOcc,K,'same');

figure(15); imagesc(LocalTimeOcc) %imagesc(Zsmooth1)
h = colorbar;
ylabel(h, '%')
colormap(jet)
%caxis([0.3 1.3])
xlabel('Year')
ylabel('Local time [hr]')

yticks = 0:length(A(:,1))/6:length(A(:,1)); %adjust as appropriate, positive integers only
yticks(1) = 1;
ylabels = {'0' '4' '8' '12' '16' '20' '24'};  %time labels
set(gca, 'YTick', yticks, 'YTickLabel', ylabels);
set(gca,'YDir','normal')
%xticks = [1 4:4:84]; %adjust as appropriate, positive integers only
xticks = 1:12:156;
xlabels = 2007:2019; 
set(gca, 'XTick', xticks, 'XTickLabel', xlabels);
xtickangle(45)
title(['Local time relative occurrence 2007-2019'])

%adding scaled monthly SNnumber which denotes solar activity
hold on
plot(1:length(RM),RM*0.8,'w-')


%% plotting Es height and thickness distribution split into local times

LocalTimeMorningIndex = find(TotalEventLocalTime < 12 & TotalEventLocalTime > 8);
LocalTimeEveningIndex = find(TotalEventLocalTime < 23 & TotalEventLocalTime > 19);

TotalDelta_H_Morning = TotalDelta_H(LocalTimeMorningIndex);
TotalDelta_H_Evening = TotalDelta_H(LocalTimeEveningIndex);
TotalEventHeight_Morning = TotalEventHeight(LocalTimeMorningIndex);
TotalEventHeight_Evening = TotalEventHeight(LocalTimeEveningIndex);

TotalEventHeight_MorningBelow80km = TotalEventHeight_Morning(TotalEventHeight_Morning < 80);
TotalDelta_H_MorningBelow80km = TotalDelta_H_Morning(TotalEventHeight_Morning < 80);
TotalDelta_H_MorningAbove80km = TotalDelta_H_Morning(TotalEventHeight_Morning >= 80 & TotalEventHeight_Morning <= 125);

TotalEventHeight_EveningBelow80km = TotalEventHeight_Evening(TotalEventHeight_Evening < 80);
TotalDelta_H_EveningBelow80km = TotalDelta_H_Evening(TotalEventHeight_Evening < 80);
TotalDelta_H_EveningAbove80km = TotalDelta_H_Evening(TotalEventHeight_Evening >= 80 & TotalEventHeight_Evening <= 125);


figure(100)
subplot(1,2,1)
yyaxis left
histogram(TotalDelta_H_MorningAbove80km, 'DisplayStyle', 'stairs','Binwidth',0.1) %'BinWidth',1  ,for 1 km bins
ylabel('Morning 08-12 UTC Above 80 km')
hold on;
yyaxis right
histogram(TotalDelta_H_MorningBelow80km, 'DisplayStyle', 'stairs','Binwidth',0.1)
ylabel('Morning 08-12 UTC Below 80 km')
xlabel('Height [km]')
ylim([0 1400])
title(['Mean thickness above 80 km: ' num2str(mean(TotalDelta_H_MorningAbove80km),3) ' km. Below 80 km: ' num2str(mean(TotalDelta_H_MorningBelow80km),3) ' km.'])

subplot(1,2,2)
yyaxis left
histogram(TotalDelta_H_EveningAbove80km,'Displaystyle','stairs','BinWidth',0.1)
ylabel('Evening 19-23 UTC Above 80 km')
ylim([0 3000])
hold on;
yyaxis right
histogram(TotalDelta_H_EveningBelow80km,'Displaystyle','stairs','BinWidth',0.1)
ylabel('Evening 19-23 UTC Below 80 km')
xlabel('Es cloud thickness [km]')
ylim([0 600])
title(['Mean thickness above 80 km: ' num2str(mean(TotalDelta_H_EveningAbove80km),3) ' km. Below 80 km: ' num2str(mean(TotalDelta_H_EveningBelow80km),3) ' km.'])





figure(101)
subplot(1,2,1)
histogram(TotalEventHeight_Morning, 'DisplayStyle', 'stairs','Binwidth',1) %'BinWidth',1  ,for 1 km bins
hold on;
histogram(TotalEventHeight_Evening, 'DisplayStyle', 'stairs','Binwidth',1)
xlabel('Height [km]')
legend('Morning 08-12 UTC', 'Evening 19-23 UTC')
title({['Total detections at: '], ['Morning ' num2str(length(TotalEventHeight_Morning)) '. Evening ' num2str(length(TotalEventHeight_Evening)) '.']})

subplot(1,2,2)
histogram(TotalDelta_H_Morning,'Displaystyle','stairs','BinWidth',0.1)
hold on;
histogram(TotalDelta_H_Evening,'Displaystyle','stairs','BinWidth',0.1)
ylabel('Thickness [km]')
xlabel('Es cloud thickness [km]')
legend('Morning 08-12 UTC', 'Evening 19-23 UTC')
title({['Mean thickness at:'], ['Morning ' num2str(mean(TotalDelta_H_Morning)) ' km. Evening ' num2str(mean(TotalDelta_H_Evening)) ' km.']})


%adding labels to figure in latex
figure(111)
subplot(1,2,1)
yyaxis left
histogram(TotalDelta_H_EveningAbove80km2, 'DisplayStyle', 'stairs','Binwidth',0.1) %'BinWidth',1  ,for 1 km bins
ylabel('Evening 19-23 UTC Above 80 km')
hold on;
yyaxis right
histogram(TotalDelta_H_EveningBelow80km2, 'DisplayStyle', 'stairs','Binwidth',0.1)
ylabel('Evening 19-23 UTC Below 80 km')
xlabel('Height [km]')
ylim([0 600])
title({['Solar minimum.'],['Mean thickness above 80 km: ' num2str(mean(TotalDelta_H_EveningAbove80km2),3) ' km. Below 80 km: ' num2str(mean(TotalDelta_H_EveningBelow80km2),3) ' km.']})

subplot(1,2,2)
yyaxis left
histogram(TotalDelta_H_EveningAbove80km,'Displaystyle','stairs','BinWidth',0.1)
ylabel('Evening 19-23 UTC Above 80 km')
ylim([0 3000])
hold on;
yyaxis right
histogram(TotalDelta_H_EveningBelow80km,'Displaystyle','stairs','BinWidth',0.1)
ylabel('Evening 19-23 UTC Below 80 km')
xlabel('Es cloud thickness [km]')
ylim([0 600])
title({['Solar maximum.'],['Mean thickness above 80 km: ' num2str(mean(TotalDelta_H_EveningAbove80km),3) ' km. Below 80 km: ' num2str(mean(TotalDelta_H_EveningBelow80km),3) ' km.']})






%% Doing same plots thickness and height, but not for morning and evening, but split in height.


%EsMeanThickness = movmean(EsMeanThickness,3); %3 month running mean
%EsMeanHeight = movmean(EsMeanHeight,3); %Good to look at the 3 month running mean (as it changes with seasons, thus 3 months), shows more clearly the effects.


figure(20);
subplot(2,1,1); plot(unique(TotalRO_MonthNumber),EsMeanThickness,'b.-')
xlabel('Year')
ylabel('Es thickness [km]')
xticks = 0:12:156
xlabels = 2007:2020;
set(gca, 'XTick', xticks, 'XTickLabel', xlabels);
xtickangle(45)
ylim([1.8 2.2])
hold on;
SummerStart = 5:12:156;
SummerStop = 8:12:156;
WinterStart = 11:12:154; %excluding first winter in 2007
WinterStop = 14:12:156;
for i = 1:length(SummerStart)
    patch([SummerStart(i) SummerStart(i) SummerStop(i) SummerStop(i)],[1.4 2.2 2.2 1.4],'y','EdgeColor','y','FaceAlpha',0.25,'EdgeAlpha',0.25)
end
for j = 1:length(WinterStart)
    patch([WinterStart(j) WinterStart(j) WinterStop(j) WinterStop(j)],[1.4 2.2 2.2 1.4],'c','EdgeColor','c','FaceAlpha',0.25,'EdgeAlpha',0.25)
end
patch([0 0 2 2],[1.41 2.19 2.19 1.41],'c','EdgeColor','c','FaceAlpha',0.25,'EdgeAlpha',0.25)


subplot(2,1,2); plot(unique(TotalRO_MonthNumber),EsMeanHeight,'b.-')
xlabel('Year')
ylabel('Es altitude [km]')
xticks = 0:12:156;
xlabels = 2007:2020;
set(gca, 'XTick', xticks, 'XTickLabel', xlabels);
xtickangle(45)
hold on;
%adding scaled monthly SNnumber (12 month running mean) which denotes solar activity
hold on;
plot(1:length(RM),RM*0.1 + mean(EsMeanHeight)-4,'r-')

%adding patches
ylim([94.62,107])
suptitle('Monthly means of Es thickness and altitude')

SummerStart = 5:12:156;
SummerStop = 8:12:156;
WinterStart = 11:12:154; %excluding first winter in 2007
WinterStop = 14:12:156;
for i = 1:length(SummerStart)
    patch([SummerStart(i) SummerStart(i) SummerStop(i) SummerStop(i)],[94.63 106.99 106.99 94.63],'y','EdgeColor','y','FaceAlpha',0.25,'EdgeAlpha',0.25)
end
for j = 1:length(WinterStart)
    patch([WinterStart(j) WinterStart(j) WinterStop(j) WinterStop(j)],[94.63 106.99 106.99 94.63],'c','EdgeColor','c','FaceAlpha',0.25,'EdgeAlpha',0.25)
end
patch([0 0 2 2],[94.63 106.99 106.99 94.63],'c','EdgeColor','c','FaceAlpha',0.25,'EdgeAlpha',0.25)


%% Thickness and height plots monthly means, at morning and evening local times.

% %---------If I want to smooth the data--------
% %doing a 3 month running mean because altitude and thickness varies with
%seasons.
% EsMeanHeightMorning = movmean(EsMeanHeightMorning,3);
% EsMeanHeightEvening = movmean(EsMeanHeightEvening,3);
% EsMeanThicknessMorning = movmean(EsMeanThicknessMorning,3);
% EsMeanThicknessEvening = movmean(EsMeanThicknessEvening,3);
% %-----------------------------



figure(30);
%subplot(2,1,1);
%yyaxis left
plot(unique(TotalRO_MonthNumber),EsMeanHeightMorning,'b.-')
xlabel('Year')
ylabel('Es < 80 km altitude [km]')
xticks = 0:12:156;
xlabels = 2007:2020;
set(gca, 'XTick', xticks, 'XTickLabel', xlabels);
xtickangle(45)
hold on;
%adding scaled monthly SNnumber (12 month running mean) which denotes solar activity
plot(1:length(RM),RM*0.10 + mean(EsMeanHeightMorning(1:50))-8,'r-','Handlevisibility','Off')

%adding patches
%ylim([98 110])
%ylim([88 108])
ylim([50 70])

SummerStart = 5:12:156;
SummerStop = 8:12:156;
WinterStart = 11:12:154; %excluding first winter in 2007
WinterStop = 14:12:156;
for i = 1:length(SummerStart)
    patch([SummerStart(i) SummerStart(i) SummerStop(i) SummerStop(i)],[50 70 70 50],'y','EdgeColor','y','FaceAlpha',0.25,'EdgeAlpha',0.25,'Handlevisibility','Off')
end
for j = 1:length(WinterStart)
    patch([WinterStart(j) WinterStart(j) WinterStop(j) WinterStop(j)],[50 70 70 50],'c','EdgeColor','c','FaceAlpha',0.25,'EdgeAlpha',0.25,'Handlevisibility','Off')
end
patch([0 0 2 2],[50 70 70 50],'c','EdgeColor','c','FaceAlpha',0.25,'EdgeAlpha',0.25,'Handlevisibility','Off')
% yyaxis right
% ylabel('Es < 80 km altitude [km]')
% ax = gca;
% ax.YAxis(2).Color = 'k';
% %ylim([54 66])
% ylim([50 70])
plot(unique(TotalRO_MonthNumber),EsMeanHeightEvening,'k.-')
legend('Morning','Evening')
%legend('Es > 80 km','Es < 80 km')
 
% 
% %-----------------------thickness----------------
figure(31);
%subplot(2,1,1); 
plot(unique(TotalRO_MonthNumber),EsMeanThicknessMorning,'b.-')
xlabel('Year')
ylabel('Es < 80 km thickness [km]')
xticks = 0:12:156;
xlabels = 2007:2020;
set(gca, 'XTick', xticks, 'XTickLabel', xlabels);
xtickangle(45)
hold on;
%adding scaled monthly SNnumber (12 month running mean) which denotes solar activity
plot(1:length(RM),RM*0.005 + mean(EsMeanThicknessMorning(1:50)) - 0.2,'r-','Handlevisibility','Off')
plot(unique(TotalRO_MonthNumber),EsMeanThicknessEvening,'k.-')
legend('Morning','Evening')
%legend('Es > 80 km','Es < 80 km')

%adding patches
ylim([1.4 2.4])
%ylim([1 3])

SummerStart = 5:12:156;
SummerStop = 8:12:156;
WinterStart = 11:12:154; %excluding first winter in 2007
WinterStop = 14:12:156;
for i = 1:length(SummerStart)
    patch([SummerStart(i) SummerStart(i) SummerStop(i) SummerStop(i)],[1 3 3 1],'y','EdgeColor','y','FaceAlpha',0.25,'EdgeAlpha',0.25,'Handlevisibility','Off')
end
for j = 1:length(WinterStart)
    patch([WinterStart(j) WinterStart(j) WinterStop(j) WinterStop(j)],[1 3 3 1],'c','EdgeColor','c','FaceAlpha',0.25,'EdgeAlpha',0.25,'Handlevisibility','Off')
end
patch([0 0 2 2],[1 3 3 1],'c','EdgeColor','c','FaceAlpha',0.25,'EdgeAlpha',0.25,'Handlevisibility','Off')



