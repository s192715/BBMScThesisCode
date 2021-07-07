%Main script that analyses a COSMIC-1 profile and outputs if S4
%scintillations are detected and if Es is detected, and extracts height and
%thickness of each Es detection as well as other variables describing the
%measurement such as time and geohgraphical coordinates.
%It calls the secondary script MScThesisScriptCombination.m, and outputs variables in a .mat file.

%-Discards "bad" RO measurements with 2ndTryCriteria, so profiles with mean A20km (above 30km height) below 300 V/V are cut.
%-Profile height cutoff at 30 km to avoid tropospheric effects.
%-Profiles only having heights lower than 30 km are removed.
%-Es multidetections are changed to single detections.
%-Occultation point of the profiles is at tangent height 30 km

------------------------------------------
MonthName = "2019_6";

%Path to directory that contains the COSMIC-1 data:
MyDataDir = dir(["/work2/s192715/Data/Months/" + MonthName]);
MyDataDirString = ["/work2/s192715/Data/Months/" + MonthName];

%Variables that are extracted:
TotalDelta_H = [];
TotalEventHeight = [];
TotalEventLat = [];
TotalEventLon = [];
TotalEventTime = [];
TotalEventAbsTime = [];
n = []; % a number that tracks from which file an event occured
TotalEventAzi = [];

%allocating for RO measurement statistics
TotalEs_present = [];
TotalNumberOfEvents = []; %shows number of Es for each RO measurement
TotalS4_present = [];
TotalBoth_present = [];
TotalNeither_present = [];
TotalOnlyEs_present = [];
TotalOnlyS4_present = [];

TotalRO_lat = [];
TotalRO_lon = [];
TotalRO_Time = [];
TotalRO_MeanS4 = [];
TotalRO_MaxS4 = [];
TotalRO_HeightOfS4Max = [];

TotalMeanSNR60_80 = [];
TotalMeanSNR = [];
TotalMeanSNR_LargeScale = [];
TotalMaxHeight = [];

FilestampIndex = 1;
FigureNumber = 1;
for k = 3:length(MyDataDir) %for each directory in MyDataDir (skips first 2 elements because they are just . and ..)
    tic
    MyDir = dir(MyDataDirString + "/" + MyDataDir(k).('name'));
    MyDirString = MyDataDirString + "/" + MyDataDir(k).('name');
    
    disp(k-2) %displays the number of day that is being calculated
    
    number = 3:length(MyDir);

    for i = 1:length(number)
        %tic
        fileName = MyDir(number(i)).('name');
        fileLocation = MyDirString + "/" + fileName;
        ncid = netcdf.open(fileLocation);
       
        [Delta_H, EventHeight, EventTime, EventAbsTime, EventLat, EventLon, EventAzi, Es_present, NumberOfEvents, S4_present, Both_present, Neither_present, OnlyEs_present, OnlyS4_present, RO_lat, RO_lon, RO_Time, RO_MaxS4, RO_MeanS4, RO_HeightOfS4Max, Filestamp, MeanSNR60_80, MeanSNR, MeanSNR_LargeScale, MaxHeight] = MScThesisScriptCombination4(ncid,0);
        FigureNumber = FigureNumber + 1;
        
        %append to vectors here
        TotalDelta_H = [TotalDelta_H,Delta_H]; %one column for each event
        TotalEventHeight = [TotalEventHeight,EventHeight];
        TotalEventLat = [TotalEventLat,EventLat];
        TotalEventLon = [TotalEventLon,EventLon];
        TotalEventTime = [TotalEventTime,EventTime];
        TotalEventAbsTime = [TotalEventAbsTime,EventAbsTime]; %these are in julian dates (not modified).
        n = [n ones(1,length(Delta_H))*number(i)]; %a number that tracks from which file an event occurred
        TotalEventAzi = [TotalEventAzi, EventAzi];
        
        %statistics for each RO measurement:
        TotalEs_present = [TotalEs_present, Es_present]; %1 for Es present, 0 for not present
        TotalNumberOfEvents = [TotalNumberOfEvents, NumberOfEvents]; %shows how many Es events there were in each RO measurement
        TotalS4_present = [TotalS4_present, S4_present];
        TotalBoth_present = [TotalBoth_present, Both_present];
        TotalNeither_present = [TotalNeither_present, Neither_present];
        TotalOnlyEs_present = [TotalOnlyEs_present, OnlyEs_present];
        TotalOnlyS4_present = [TotalOnlyS4_present, OnlyS4_present];

        TotalRO_lat = [TotalRO_lat, RO_lat];
        TotalRO_lon = [TotalRO_lon, RO_lon];
        TotalRO_Time = [TotalRO_Time, RO_Time];
        TotalRO_MeanS4 = [TotalRO_MeanS4, RO_MeanS4];
        TotalRO_MaxS4 = [TotalRO_MaxS4, RO_MaxS4];    
        TotalRO_HeightOfS4Max = [TotalRO_HeightOfS4Max, RO_HeightOfS4Max];
        
        TotalMeanSNR60_80 = [TotalMeanSNR60_80 MeanSNR60_80];
        TotalMeanSNR = [TotalMeanSNR MeanSNR];
        TotalMeanSNR_LargeScale = [TotalMeanSNR_LargeScale MeanSNR_LargeScale];
        TotalMaxHeight = [TotalMaxHeight MaxHeight];
        

        netcdf.close(ncid)
        
        %collects profiles that were cut out with 2ndTryCriteria
        if isempty(Filestamp) == 0 %if Filestamp is not empty
            FilestampCell{FilestampIndex} = Filestamp;
            FilestampIndex = FilestampIndex + 1;
        end
        
    end
    toc
end

%Percentages in RO measurements:
TotalNumberROmeasurements = length(TotalEs_present);

Es_perc = 100*sum(TotalEs_present)/TotalNumberROmeasurements;
S4_perc = 100*sum(TotalS4_present)/TotalNumberROmeasurements;
Both_perc = 100*sum(TotalBoth_present)/TotalNumberROmeasurements;
Neither_perc = 100*sum(TotalNeither_present)/TotalNumberROmeasurements;
OnlyEs_perc = 100*sum(TotalOnlyEs_present)/TotalNumberROmeasurements;
OnlyS4_perc = 100*sum(TotalOnlyS4_present)/TotalNumberROmeasurements;

% %k = 33; %for July2009
% k = 30; %for February2009
% X = categorical({'Es','Only Es','S4','Only S4','Both','Neither'});
% X = reordercats(X,{'Es','Only Es','S4','Only S4', 'Both','Neither'});
% figure(1); bar(X,[Es_perc OnlyEs_perc S4_perc OnlyS4_perc Both_perc Neither_perc])
% title({['2009 February statistics.'],[num2str(TotalNumberROmeasurements) ' RO measurements over ' num2str(k-2) ' days.']})
% ylabel('%')
% ylim([0 100])


%--------------------------------------------
%Extracting local time in hours from TotalEventAbsTime

%Extracting UTC time of day from juliandate in TotalEventAbsTime:
TotalEventTimeUTC  = datevec(datetime(TotalEventAbsTime,'ConvertFrom','juliandate'));
TotalEventTimeUTC_HourOfDay = 24*datenum(hours(TotalEventTimeUTC(:,4)) + minutes(TotalEventTimeUTC(:,5)) + seconds(TotalEventTimeUTC(:,6)));

%converting utc time to local time, depending on time zones:
[zd] = timezone(TotalEventLon);
TotalEventLocalTime = (TotalEventTimeUTC_HourOfDay - zd)';


%making all values positive between 0 <= t < 24 i.e. displays the local time of the
%previous date instead for negative value:
for i = 1:length(TotalEventLocalTime)
    if TotalEventLocalTime(i) < 0
        TotalEventLocalTime(i) = 24 + TotalEventLocalTime(i);
    elseif TotalEventLocalTime(i) >= 24
        TotalEventLocalTime(i) = TotalEventLocalTime(i) - 24;
    end
end

%--------------------------------------------------------------------------
% %For plotting distributions for Es events:
% figure(2); subplot(2,3,1)
% plot(TotalEventHeight,TotalEventLat,'k.','Markersize',2.5)
% xlabel('Height [km]')
% ylabel('Latitude [°]')
% ylim([-90 90])
% yticks([-90 -45 0 45 90])
% grid on;
% 
% subplot(2,3,2)
% plot(TotalEventLocalTime,TotalEventLat,'k.','Markersize',2.5)
% xlabel('Local time [Hours]')
% ylabel('Latitude [°]')
% xlim([0 24])
% xticks([0 4 8 12 16 20 24])
% ylim([-90 90])
% yticks([-90 -45 0 45 90])
% grid on;
% 
% subplot(2,3,3)
% plot(TotalEventHeight,TotalDelta_H,'k.','Markersize',2.5)
% xlabel('Height [km]')
% ylabel('Es cloud thickness [km]')
% grid on;
% 
% subplot(2,3,4)
% histogram(TotalEventHeight, 'DisplayStyle', 'stairs','Binwidth',1) %'BinWidth',1  ,for 1 km bins
% xlabel('Height [km]')
% ylabel('Number')
% %ylim([0 800]) %For July2009 comparison
% ylim([0 400]) %For February2009 comparison
% 
% % % 
% % % [values, edges] = histcounts(TotalEventHeight, 'Binwidth',1);
% % % centers = (edges(1:end-1)+edges(2:end))/2;
% % % figure
% % % plot(centers, values, 'k-')
% 
% subplot(2,3,5)
% TotalDelta_H_80_125 = TotalDelta_H(80 <= TotalEventHeight & TotalEventHeight <= 125);
% TotalDelta_H_under80 = TotalDelta_H(TotalEventHeight < 80);
% yyaxis left
% histogram(TotalDelta_H_80_125,'Displaystyle','stairs','BinWidth',0.1)
% ylabel('80 km ≤ h ≤ 125 km')
% %ylim([0 1000]) %For July2009 comparison
% ylim([0 400]) %For February2009 comparison
% %xlim([0 7])
% hold on;
% yyaxis right
% histogram(TotalDelta_H_under80,'Displaystyle','stairs','BinWidth',0.1)
% ylabel('h < 80 km')
% %ylim([0 200]) %For July2009 comparison
% ylim([0 50]) %For February2009 comparison
% %histogram(TotalDelta_H,'Displaystyle','stairs','BinWidth',0.1) % could do 0.1 km bins
% %ylabel('Number')
% xlabel('Es cloud thickness [km]')
% 
% 
% subplot(2,3,6)
% geoplot(TotalEventLat,TotalEventLon,'k.','Markersize',1)
% geolimits([-90 90],[-180 180])
% geobasemap landcover
% suptitle(['2009 February Total RO measurements: ' num2str(TotalNumberROmeasurements) '. Total Es Detections: ' num2str(length(TotalDelta_H)) '. Es occurance in RO measurement: ' num2str(Es_perc,'%4.2f') '%.'])
% 

%-------------------------------------------------------------------------
%path for output variable:
save(["/zhome/e8/9/144512/Desktop/MScThesis/MatlabScripts/SavedMatlabVariables/MonthsMainRunAll2Test/" + MonthName + ".mat"],'TotalDelta_H','TotalEventHeight','TotalEventLat','TotalEventLon','TotalEventTime','TotalEventAbsTime','n','TotalEventLocalTime','TotalEventAzi', 'TotalEs_present', 'TotalNumberOfEvents', 'TotalS4_present', 'TotalBoth_present', 'TotalNeither_present', 'TotalOnlyEs_present', 'TotalOnlyS4_present', 'TotalRO_lat', 'TotalRO_lon', 'TotalRO_Time', 'TotalRO_MaxS4', 'TotalRO_MeanS4', 'TotalRO_HeightOfS4Max', 'FilestampCell', 'TotalMeanSNR60_80', 'TotalMeanSNR', 'TotalMeanSNR_LargeScale', 'TotalMaxHeight')

%load(["/zhome/e8/9/144512/Desktop/MScThesis/MatlabScripts/SavedMatlabVariables/MonthsMainRunAll/CombineMonths/" + MonthName + ".mat"])
%--------------------------------------------------------------------------







