%Finding magnetic equator coordinates so it can be plotted (used in
%function GlobalMaps.m)
lat = -15:0.01:15;
lon = -179:0.1:180;
TotalDIP = zeros(length(lat),length(lon));

IterationNumber = 1;
LatPlacement = 1;
LonPlacement = 1;
for i = lat
    disp(IterationNumber)
    tic
    for j = lon
        [~, ~, ~, DIP] = igrfmagm(60000, i, j, decyear(2009,1,1),12);
        TotalDIP(LatPlacement,LonPlacement) = DIP;
        LonPlacement = LonPlacement + 1;
    end
    toc
    IterationNumber = IterationNumber + 1;
    LatPlacement = LatPlacement + 1;
    LonPlacement = 1;
end

[MinDipAngle, LatIndex] = min(abs(TotalDIP));
closestValue1 = TotalDIP(LatIndex);

plot(lon,lat(LatIndex))

save('/zhome/e8/9/144512/Desktop/MScThesis/MatlabScripts/SavedMatlabVariables/MagneticEquator.mat','lon','lat','LatIndex')
%load('/zhome/e8/9/144512/Desktop/MScThesis/MatlabScripts/SavedMatlabVariables/MagneticEquator.mat')

%% Plotting geomagnetic field components

model_epoch = '2015';
decimal_year = 2015;


height = 1000000;

% Geodetic Longitude value in degrees to use for latitude sweep. 
geod_lon = -180:1:180;      %degrees

% Geodetic Latitude values to sweep. 
geod_lat = -89.5:.5:89.5; %degrees

% Loop through longitude values for each array of latitudes -89.5:89.5.
for lonIdx = size(geod_lon,2):-1:1
    tic
    for latIdx = size(geod_lat,2):-1:1

    % Use WRLDMAGM function to obtain magnetic parameters for each lat/lon value. 
    [xyz, h, dec, dip, f] = wrldmagm(height, geod_lat(latIdx),geod_lon(lonIdx), decimal_year, model_epoch);

    %Store results 
    WMMResults(latIdx,1:7,lonIdx) = [xyz' h dec dip f];

    end
    toc
end


% WMMFileName = 'WMMResults_Epoch_2015v2_decyear_2015.mat';
% load(WMMFileName);


landAreas = shaperead('landareas.shp','UseGeoCoords',true);
plotWMM = load('astPlotWMM.mat');

hH = figure;
set(hH,'Position',[0 0 827 620],'Color','white')
astPlotWMMContours( WMMResults, plotWMM, 4, landAreas, geod_lat, geod_lon, decimal_year, model_epoch)
fullFileName = '/zhome/e8/9/144512/Desktop/MScThesis/Figures/MagneticFieldPlots/hH_1000km.png';
saveas(hH,fullFileName);

hDIP = figure;
set(hDIP,'Position',[0 0 827 620],'Color','white')
astPlotWMMContours( WMMResults, plotWMM, 6, landAreas, geod_lat, geod_lon, decimal_year, model_epoch)
fullFileName = '/zhome/e8/9/144512/Desktop/MScThesis/Figures/MagneticFieldPlots/hDIP_1000km.png';
saveas(hDIP,fullFileName);

hF = figure;
set(hF,'Position',[0 0 827 620],'Color','white')
astPlotWMMContours( WMMResults, plotWMM, 7, landAreas, geod_lat, geod_lon, decimal_year, model_epoch)
fullFileName = '/zhome/e8/9/144512/Desktop/MScThesis/Figures/MagneticFieldPlots/hF_1000km.png';
saveas(hF,fullFileName);


 H = fspecial('gaussian',[3 3],0.5);
 A = [1 2 3; 4 5 6; 7 8 9];
 out = imfilter(A,H,'replicate');


out = conv2(H,A,'same')
