%When using data from scratch folder:
MyDataDir = dir("/work1/s192715/Data/test");
MyDataDirString = "/work1/s192715/Data/test";

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

FigureNumber = 1;
for k = 3:length(MyDataDir) %for each directory in MyDataDir (skips first 2 elements because they are just . and ..)
    tic
    k = k;
    MyDir = dir(MyDataDirString + "/" + MyDataDir(k).('name'));
    MyDirString = MyDataDirString + "/" + MyDataDir(k).('name');
    
    disp(k-2) %displayes the number of day that is being calculated
    
    number = 3:length(MyDir);
    
    n = []; %starts over after every day
    TotalDelta_H = [];

    for i = 1:length(number)
        %tic
        fileName = MyDir(number(i)).('name');
        fileLocation = MyDirString + "/" + fileName;
        ncid = netcdf.open(fileLocation);
        
        %disp(FigureNumber)
        [Delta_H, EventHeight, EventTime, EventAbsTime, EventLat, EventLon, EventAzi, Es_present, NumberOfEvents, S4_present, Both_present, Neither_present, OnlyEs_present, OnlyS4_present, RO_lat, RO_lon, RO_Time] = MScThesisScriptCombination2(ncid,0);
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
        
        nn{k} = n; %number that tracks which day the file was in
        TotalDelta_H_Cell{k} = TotalDelta_H;

        
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
                
        netcdf.close(ncid)
        
        %toc
    end
    toc
end

 
% [counts,value] = groupcounts(TotalDelta_H');
% SameThickness = sum(counts)-length(counts); %Number of events that have the same thickness
% MaxCounts = max(counts);
% ThicknessAtMaxCounts = value(counts==MaxCounts);
% SortedCounts = sort(counts); %lowest to highest

%----
% d = [true, diff(nn{3}) ~= 0, true];  % TRUE if values change
% n = diff(find(d));               % Number of repetitions
% Y = repelem(n, n); %outputs vector of same length as input vector
% 
% nn{3}(find(Y==8)) %9 was manually selected, it is max(Y) i.e. in one file 9 detections were made with same thickness
%----
CellIndex = 8;
d = [true, diff(TotalDelta_H_Cell{CellIndex}) ~= 0, true];  % TRUE if values change
n = diff(find(d));               % Number of repetitions
Y = repelem(n, n) %outputs vector of same length as input vector

TotalDelta_H_Cell{CellIndex}(find(Y==max(Y)))
nn{CellIndex}(find(Y==max(Y))) %profile number in the given day
Ymax = max(Y)
%----


%% 

%-----to plot manually profile of a file in a given day
MyDir = dir("/work1/s192715/Data/Months/2009_2/atmPhs_repro2013_2009_035");
MyDirString = "/work1/s192715/Data/Months/2009_2/atmPhs_repro2013_2009_035";
%number = 123;
number = 1267;

FigureNumber = 1;
for i = number
tic
fileName = MyDir(number).('name');
fileLocation = MyDirString + "/" + fileName;
ncid = netcdf.open(fileLocation);

disp(FigureNumber)
PlotIndex = 1; %1 for plot, 0 for no plot.
fail = PlotFailedProfiles(ncid,PlotIndex);
FigureNumber = FigureNumber + 1;

end