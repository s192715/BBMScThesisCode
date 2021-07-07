function p = PlotFailedProfiles(ncid,PlotIndex)

%same as MScThesisScriptCombination.m but for tangent point calculation the
%Fortran routine has been implemented, making this script much faster.

%FigureNumber input can be set to zero for no plot

%for selected RO measurements in netCDF files, can plot:
%SNR vs time
%SNR vs h_t
%S4 vs time
%marks lat and lon on a map

%Also: Finds all sporadicE events for the selected RO measurement
%returns Es cloud thickness, height, time from measurement start, absolute time, lat,
%and lon of each event.


%-Can handle zero events, i.e. if event has FirstCriteriaIndex zero or 1.
%-Excludes first and last 1.1 km of the height interval, otherwise criteria
%in condition 2 can not be evaluated.
%-When calculating the height, if two or more height elements are the minimum, the
%first one gets selected.

%SNR vs time
caL1Snr = netcdf.getVar(ncid,1);
time = netcdf.getVar(ncid,0);
fileStamp = netcdf.getAtt(ncid,netcdf.getConstant('NC_GLOBAL'),'fileStamp');

%--------------------------------------------------------------
%SNR vs h_t

%first, get datetime from netCDF and then convert to Julian (J2000)
year = double(netcdf.getAtt(ncid,netcdf.getConstant('NC_GLOBAL'),'year'));
month = double(netcdf.getAtt(ncid,netcdf.getConstant('NC_GLOBAL'),'month'));
day = double(netcdf.getAtt(ncid,netcdf.getConstant('NC_GLOBAL'),'day'));
hour = double(netcdf.getAtt(ncid,netcdf.getConstant('NC_GLOBAL'),'hour'));
minute = double(netcdf.getAtt(ncid,netcdf.getConstant('NC_GLOBAL'),'minute'));
second = double(netcdf.getAtt(ncid,netcdf.getConstant('NC_GLOBAL'),'second'));

%getting the ECI coordinates of LEO and GPS
xLeo_ECI = netcdf.getVar(ncid,4); %units km
yLeo_ECI = netcdf.getVar(ncid,5);
zLeo_ECI = netcdf.getVar(ncid,6);
%xdLeo_ECI = netcdf.getVar(ncid,7); %units km/s
%ydLeo_ECI = netcdf.getVar(ncid,8);
%zdLeo_ECI = netcdf.getVar(ncid,9);

xGPS_ECI = netcdf.getVar(ncid,10);
yGPS_ECI = netcdf.getVar(ncid,11);
zGPS_ECI = netcdf.getVar(ncid,12);
%xdGPS_ECI = netcdf.getVar(ncid,13);
%ydGPS_ECI = netcdf.getVar(ncid,14);
%zdGPS_ECI = netcdf.getVar(ncid,15);

%--------------------------------------------------------------------------
%FORTRAN Transforming ECI into ECEF
r_LEO_ECI = [xLeo_ECI'; yLeo_ECI'; zLeo_ECI']; %format: 3xN
r_GPS_ECI = [xGPS_ECI'; yGPS_ECI'; zGPS_ECI']; %format: 3xN

jdf = juliandate([year month day 0 0 0]); %fractional Julian day at 00:00:00 UTC
utc_seconds = hour*3600 + minute*60 + second*1 + time; %seconds since 00:00:00 UTC

jdf2000 = 2451545;
tu = (jdf - jdf2000)./36525; %julian centuries since J2000

%This is only to have utc time format, used later in retrieving absolute time:
utc_start_time = [year, month, day, hour, minute, second];
utc = datevec(datetime(utc_start_time) + seconds(time));

%-------------------
%From eci2eci_nosofa
%trying first eci2eci, transformation from eci to intermediate eci

%Formula for precession from Hoffman-Wellenhof et al., 2008.
zet = (2306.2181*tu + 0.30188*tu^2 + 0.017998*tu^3)*pi/648000;
z   = (2306.2181*tu + 1.09468*tu^2 + 0.018203*tu^3)*pi/648000;
ups = (2004.3109*tu - 0.42665*tu^2 - 0.041833*tu^3)*pi/648000;
szet = sin(zet);
czet = cos(zet);
sz   = sin(z);
cz   = cos(z);
sups = sin(ups);
cups = cos(ups);                 
RP = [cz*cups*czet-sz*szet, -cz*cups*szet-sz*czet, -cz*sups; sz*cups*czet+cz*szet , -sz*cups*szet+cz*czet , -sz*sups; sups*czet, -sups*szet, cups];
Qt = RP; %this is the rotation matrix

%intermediate ECI coordinates
r_LEO_ECI_Interm = Qt*r_LEO_ECI;
r_GPS_ECI_Interm = Qt*r_GPS_ECI;

%--------------------
%From eci2ecf_nosofa
%formulas from Astronomical Almanac, 1993.
gmst0 = 24110.54841 + 8640184.812866*tu + 0.093104*tu.^2 - (6.2*10^(-6))*tu.^3;
gmst  = mod(gmst0 + utc_seconds*1.0027379093, 86400); %mod returns remainder after division
 
phi = gmst*2*pi/86400; %in radians

pa = repmat([0 0 1], [length(time),1]); %polar axis (in format to match each time step coordinate)
r_LEO_ECEF_test = vector_rotate(r_LEO_ECI_Interm', pa, -phi); %ECI -> ECF, negative phi.
r_GPS_ECEF_test = vector_rotate(r_GPS_ECI_Interm', pa, -phi);

r_LEO_ECEF = r_LEO_ECEF_test'; %format 3xN
r_GPS_ECEF = r_GPS_ECEF_test';

%--------------------------------------------------------------------------
%FORTRAN, calculating the tangent height h_t, latitude, and longitude as a function of time step:

%1. find the perigee point i.e. tangent point in ECEF(for all timesteps).
%2. find the occultation point (perigee point that is closest to the wgs84 ellipsoid).
%3. find radius of curvature r_coc in ECEF coordinates by travelling from the occultation point
    %one Re in direction opposite to the normal of the surface.
%4. find the perigee points in ECEF and corresponding lat lon h w.r.t. r_coc.

wgs84 = wgs84Ellipsoid('kilometer');

r_LEO_ECEF2 = r_LEO_ECEF'; %format Nx3
r_GPS_ECEF2 = r_GPS_ECEF';

%height of the tangent point
slta = impact_paramter(r_LEO_ECEF2, r_GPS_ECEF2);

ro = sqrt(dot(r_LEO_ECEF2, r_LEO_ECEF2,2)); %length of LEO vector
alpha = acos(slta./ro);

%coordinates of the tangent point in ECEF
perigee = vector_rotate(r_LEO_ECEF2, cross(r_LEO_ECEF2, r_GPS_ECEF2), alpha).*(slta./ro);

%----------------------------------------
%Finding r_coc, the ECEF coordinates of the center of curvature.
%transforming the perigee in ECEF to geodetic
[la lo h] = ecef2geodetic(wgs84,perigee(:,1),perigee(:,2),perigee(:,3));

%iocc is index of the occultation point, finding the index of lowest perigee:
iocc = find(abs(h) == abs(min(abs(h))));

%starting values, lat and lon of perigee point
latF = la(iocc);
lonF = lo(iocc);

%surface normal:
n = [cosd(latF).*cosd(lonF); cosd(latF).*sind(lonF); sind(latF)];

%equatorial radius in km
Re = 6378.137;

%ECEF coordinates of center of curvature:
r_coc = perigee(iocc,:) - Re*n'; 

%----------------------------------------
%second iteration, now with r_coc fix:

r_LEO_ECEF_new = r_LEO_ECEF2 - r_coc; %format Nx3
r_GPS_ECEF_new = r_GPS_ECEF2 - r_coc;

%height of the tangent point
slta = impact_paramter(r_LEO_ECEF_new, r_GPS_ECEF_new);

ro = sqrt(dot(r_LEO_ECEF_new, r_LEO_ECEF_new,2)); %length of LEO vector
alpha = acos(slta./ro);

%coordinates of the tangent point in ECEF
perigee = vector_rotate(r_LEO_ECEF_new, cross(r_LEO_ECEF_new, r_GPS_ECEF_new), alpha).*(slta./ro);

lat_tp_new = 90 - rad2deg(acos(perigee(:,3)./sqrt(sum(perigee.^2,2))));
lon_tp_new = rad2deg(atan2(perigee(:,2), perigee(:,1)));
h_t_new = slta-Re;

h_t = h_t_new'; %format 1xN, same format as from old MyRoutine.
latitude = lat_tp_new';
longitude = lon_tp_new';
%--------------------------------------------------------------------------

%Know if it's setting or rising occultation
setting = netcdf.getAtt(ncid,netcdf.getConstant('NC_GLOBAL'),'setting');


%Calculating azimuth:
%1. Unit vector from LEO to GPS
LEOtoGPS = (r_GPS_ECEF-r_LEO_ECEF)'./vecnorm((r_GPS_ECEF-r_LEO_ECEF)',2,2);
%GPStoLEO = (r_LEO_ECEF-r_GPS_ECEF)'./vecnorm((r_LEO_ECEF-r_GPS_ECEF)',2,2)
%2. ECEF of tangent point
[Xtangent, Ytangent, Ztangent] = geodetic2ecef(wgs84,latitude,longitude,h_t);
%3. ECEF of the point that has traveled the unit vector of GPStoLEO from tangent point:
Xpoint = Xtangent' + LEOtoGPS(:,1);
Ypoint = Ytangent' + LEOtoGPS(:,2);
Zpoint = Ztangent' + LEOtoGPS(:,3);
%4. Transforming the ECEF of the point to geodetic
[LatPoint, LonPoint, ~] = ecef2geodetic(wgs84, Xpoint, Ypoint, Zpoint);
%5. Finding azimuth between two points in geodetic coordintates
Azi = azimuth(latitude',longitude',LatPoint,LonPoint,wgs84);

%-------------------------------------------------------------------------------------------------
%Finding SporadicE events

%sort so h_t contains ascending values:
%flipping all vectors so each element corresponds to the correct height
reversed_time = 0;
if h_t(end) < h_t(1)
    h_t = flip(h_t);
    caL1Snr = flip(caL1Snr);
    latitude = flip(latitude);
    longitude = flip(longitude);
    time = flip(time);
    utc = flip(utc); %first line now corresponds to first element in time
    Azi = flip(Azi);
    reversed_time = 1; %to indicate the in time_pos, the last value actually corresponds to the first measurement
end

%Extracting only the values corresponding to positive h_t values
NegativeValuesIndex = find(h_t < 0);

%taking into account if measurement actually doesn't contain any negative
%height values (happened at least for one RO measurement, n=1172 in 2014_001).
if isempty(NegativeValuesIndex)
    NegativeValuesIndex = 0;
end

NumberOfPositiveValues = length(h_t)-NegativeValuesIndex(end);
h_t_pos = h_t(end-NumberOfPositiveValues+1:end);
caL1Snr_pos = caL1Snr(end-NumberOfPositiveValues+1:end);
latitude_pos = latitude(end-NumberOfPositiveValues+1:end);
longitude_pos = longitude(end-NumberOfPositiveValues+1:end);
time_pos = time(end-NumberOfPositiveValues+1:end);
utc_pos = utc(end-NumberOfPositiveValues+1:end,:);
Azi_pos = Azi(end-NumberOfPositiveValues+1:end,:);


%calculating the 1km and 20km sliding averages
A1 = movmean(caL1Snr_pos,1,'SamplePoints',h_t_pos);
%plot(h_t_pos,A1/10,'r-')

A20 = movmean(caL1Snr_pos,20,'SamplePoints',h_t_pos);
%plot(h_t_pos,A20/10,'y-')

%Index for h_t where Es has been identified with criteria 1 (in paper 2).
FirstCriteriaIndex = find(A1 < 0.6*A20);

% if PlotIndex == 1
%     % %plot to analyze each profile better:
%     plot(h_t_pos,caL1Snr_pos/10,'b-','HandleVisibility','off')
%     hold on;
%     plot(h_t_pos,A1/10,'r-')
%     plot(h_t_pos,A20/10,'y-')
%     plot(h_t_pos(FirstCriteriaIndex(1:3:end)),caL1Snr_pos(FirstCriteriaIndex(1:3:end)),'k.','HandleVisibility','off')
%     %plots every third element of FirstCriteriaIndex just for less density, but
%     %shows where the algorithm detects an event.
%     xlabel('Straigt-line-tangent height [km]')
%     ylabel('SNR [V/V]')
%     legend('1km mean','20km mean')
% end


%This piece of code identifies the number of sporadicE events in a RO measurement (given its indexes are continuous). 
%Also, creates variable EventSpan with the indexes of each event (each line is one event).
NumberOfEvents = 0;
if isempty(FirstCriteriaIndex) == 0 %if the FirstCriteriaIndex is not empty
    NumberOfEvents = 1; %number of Es events
    EventSpan = zeros(length(FirstCriteriaIndex),length(FirstCriteriaIndex)); %EventSpan are the FirstCriteriaIndex of each event
    
    if length(FirstCriteriaIndex) == 1 %accounts for if FirstCriteriaIndex is only of length one (then lower loop won't run).
        EventSpan = FirstCriteriaIndex;
    end
    
    for i = 1:length(FirstCriteriaIndex)-1
        if FirstCriteriaIndex(i+1) - FirstCriteriaIndex(i) == 1
            EventSpan(NumberOfEvents,i) = FirstCriteriaIndex(i);
        
        elseif FirstCriteriaIndex(i+1)-FirstCriteriaIndex(i) ~= 1
            EventSpan(NumberOfEvents,i) = FirstCriteriaIndex(i);
            NumberOfEvents = NumberOfEvents + 1;
        end
        
        %to include the last element of FirstCriteriaIndex:
        if FirstCriteriaIndex(i) == FirstCriteriaIndex(end-1)
            EventSpan(NumberOfEvents,i+1) = FirstCriteriaIndex(end);
        end
    end
    
    %removing all zero rows:
    EventSpan(sum(EventSpan,2) == 0,:) = []; %Removing lines containing only zero.
    %EventSpan is now matrix where each row is an event, containing the indexes of h_t_pos (and thus also caL1Snr) where the event occurs.
    
    %This isn't necessary after thickness condition fixed, FirstCriteriaIndex (i.e. Hspan) can be of length 1
%     RemoveSingleElementIndex = []; %if only one height element, then remove that as an event:
%     for i = 1:NumberOfEvents
%         if nnz(EventSpan(i,:)) == 1 %nnz returns number of non-zero elements
%             RemoveSingleElementIndex = [RemoveSingleElementIndex i];
%         end
%     end
%     EventSpan(RemoveSingleElementIndex,:) = [];
%     NumberOfEvents = NumberOfEvents - length(RemoveSingleElementIndex);


    %if it were to happen that all NumberOfEvents would be removed in the
    %above, then these would be empty, so all fine.
    Delta_H = zeros(1,NumberOfEvents); %allocating vector for thickness values (first value corresponds to event 1 etc.)
    EventHeight = zeros(1,NumberOfEvents);
    EventTime = zeros(1,NumberOfEvents); %time since start of occultation
    EventAbsTime = zeros(1,NumberOfEvents); %these will be modified julian dates
    EventLat = zeros(1,NumberOfEvents);
    EventLon = zeros(1,NumberOfEvents);
    EventAzi = zeros(1,NumberOfEvents);
   
        
    %Obtaining Es layer thickness of each Es event (new vector Hspan always created, contains h_t_pos index for each event) 
    Hspan = zeros(1,length(EventSpan(1,:)));
    for i = 1:NumberOfEvents %for each event i.e. each line in EventSpan
        Hspan(1:length(EventSpan(i,:))) = EventSpan(i,:); %Hspan is made so line in EventSpan matrix can be worked with as vector.
        Hspan(Hspan == 0) = []; %removing zeros in the EventSpan line

        LowerSideExtension = []; %lower side extension of the FirstCriteriaIndex, spans the lower end of (h1,h2)
        UpperSideExtension = [];
        p = 1;
        k = 1;
        
        if Hspan(1) == 1 %makes sure if Hspan(1) is the first element, then LowerSideExtension will be empty.
            LowerSideExtension = [];
        else
            while A1(Hspan(1)-p) < 0.95*A20(Hspan(1)-p)       
                LowerSideExtension = [Hspan(1)-p, LowerSideExtension];
                p = p+1;
                if Hspan(1)-p <= 0 %lower bound assigned to first element, this will be disregarded as an event later
                    break          %when evaluating the h1-1km condition because it can't evaluate that.
                end     
            end
        end
                
        if Hspan(end) == length(h_t_pos) %makes sure if Hspan(end) is the last element, then UpperSideExtension will be empty.
            UpperSideExtension = [];
        else
            while A1(Hspan(end)+k) < 0.95*A20(Hspan(end)+k)         
                UpperSideExtension = [UpperSideExtension, Hspan(end)+k];
                k = k+1;
                if Hspan(end)+k >= length(h_t_pos)
                    break
                end  
            end
        end

        %The FirstCriteriaIndexCondition spans the (h1,h2) interval in paper 2, and
        %corresponds to the thickness of the cloud.
        FirstCriteriaIndexCondition = [LowerSideExtension, Hspan, UpperSideExtension];  
        
        
        %---------------------
        %Might need to add some condition here if
        %FirstCriteriaIndexCondition can be 1:
        if length(FirstCriteriaIndexCondition) == 1
            Delta_H(i) = [];
            EventHeight(i)= [];
            EventTime(i) = [];
            EventLat(i) = [];
            EventLon(i) = [];
            EventAbsTime(i) = [];
            EventAzi(i) = [];
            continue
            %break
        end
        %---------------------
               
        %FirstCriteriaIndexCondition corresponds to the interval (h_1,h_2) in paper 2
        if isempty(FirstCriteriaIndexCondition) == 0
            Delta_H(i) = h_t_pos(FirstCriteriaIndexCondition(end)) - h_t_pos(FirstCriteriaIndexCondition(1));
        end

        %middle element of FirstCriteriaIndexCondition found, which is approximately the middle height (if height
        %approximated as continuously increasing at each time step over the interval).
        FirstCriteriaIndexCondition_Middle = ceil((FirstCriteriaIndexCondition(1) + FirstCriteriaIndexCondition(end))/2);

        %extracting other parameters associated with each Es event:
        %These are the value in the middle of the cloud of each event.
        %placed here to use FirstCriteriaIndexCondition_Middle.
        %pos values are only the values at indexes that had positive height.
        EventHeight(i)= h_t_pos(FirstCriteriaIndexCondition_Middle);
        EventTime(i) = time_pos(FirstCriteriaIndexCondition_Middle);
        EventLat(i) = latitude_pos(FirstCriteriaIndexCondition_Middle);
        EventLon(i) = longitude_pos(FirstCriteriaIndexCondition_Middle);
        EventAbsTime(i) = juliandate(utc_pos(FirstCriteriaIndexCondition_Middle,:)); 
        EventAzi(i) = Azi_pos(FirstCriteriaIndexCondition_Middle);

    %Criteria2:
        %finding the std of the middle interval
        if Delta_H(i) < 1
            A_std_1km_middle(i) = std(caL1Snr_pos(FirstCriteriaIndexCondition));
        else

            %finding the first and last index of the thickness interval
            FirstCriteriaIndexCondition_First = FirstCriteriaIndexCondition_Middle;
            FirstCriteriaIndexCondition_End = FirstCriteriaIndexCondition_Middle;
            while h_t_pos(FirstCriteriaIndexCondition_Middle) - h_t_pos(FirstCriteriaIndexCondition_First) <= 0.5 %finding lower bound of 1 km interval from center
                FirstCriteriaIndexCondition_First = FirstCriteriaIndexCondition_First-1;
                if FirstCriteriaIndexCondition_First == 1 %our height resolution is about 50m, this can only happen if Delta_H is just above 1
                    break                                 %and situated close to zero height (would be cut later), but we'd then be estimating the thickness wrong by max ~50 meters. 
                end
            end
            while h_t_pos(FirstCriteriaIndexCondition_End) - h_t_pos(FirstCriteriaIndexCondition_Middle) <= 0.5 %finding upper bound of 1 km interval from center
                FirstCriteriaIndexCondition_End = FirstCriteriaIndexCondition_End + 1;
                if FirstCriteriaIndexCondition_End == length(h_t_pos) %if the U-shape is at the end and index reaches end of heights
                    break
                end
            end

            FirstCriteriaIndexCondition_1kmInterval = FirstCriteriaIndexCondition_First:FirstCriteriaIndexCondition_End;
            A_std_1km_middle(i) = std(caL1Snr_pos(FirstCriteriaIndexCondition_1kmInterval));

        end

        
        %Excludes event at first and last 1.1 km of the height interval:
        %(because then it can't evaluate the std at the side intervals, for criteria 2)
        if h_t_pos(FirstCriteriaIndexCondition(1)) - h_t_pos(1) < 1.1 || h_t_pos(end) - h_t_pos(FirstCriteriaIndexCondition(end)) < 1.1
            Delta_H(i) = [];
            EventHeight(i)= [];
            EventTime(i) = [];
            EventLat(i) = [];
            EventLon(i) = [];
            EventAbsTime(i) = [];
            EventAzi(i) = [];
            continue            
        else
            %finding indexes of lower side interval (h_1,h_1 - 1km)
            k = 1;
            while h_t_pos(FirstCriteriaIndexCondition(1)) - h_t_pos(FirstCriteriaIndexCondition(1) - k) < 1
                k = k+1;
            end
            LowerSideIntervalIndex = (FirstCriteriaIndexCondition(1)-k):FirstCriteriaIndexCondition(1);
            
            %finding indexes of upper side interval (h_2,h_2 + 1 km)
            p = 1;
            while h_t_pos(FirstCriteriaIndexCondition(end) + p) - h_t_pos(FirstCriteriaIndexCondition(end)) < 1
                p = p+1;
            end
            UpperSideIntervalIndex = FirstCriteriaIndexCondition(end):(FirstCriteriaIndexCondition(end)+k);
            
            A_std_lower(i) = std(caL1Snr_pos(LowerSideIntervalIndex));
            A_std_upper(i) = std(caL1Snr_pos(UpperSideIntervalIndex));            
        end
            
        %Removing events that do not fulfill criteria 2:
        if (A_std_1km_middle(i) < A_std_lower(i)) && (A_std_1km_middle(i) < A_std_upper(i))
            if (A_std_1km_middle(i) < (1/3)*A_std_lower(i)) || (A_std_1km_middle(i) < (1/3)*A_std_upper(i))
                  %nothing happens if these criteria are fulfilled
            else
                Delta_H(i) = 0;%Delta_H(i) had already been calculated from criteria 1, now if criteria 2 is not fulfilled then it will be given zero value.
            end
        else
            Delta_H(i) = 0;
        end

    end
    
else
    %so, if FirstCriteriaIndex is empty, every output is empty i.e. no event.
    Delta_H = [];
    EventHeight = [];
    EventTime = [];
    EventAbsTime = [];
    EventLat = [];
    EventLon = [];
    EventAzi = [];
    
end

%clean up vectors from zero
ZeroInd = find(~Delta_H);

Delta_H(ZeroInd) = [];
EventHeight(ZeroInd) = [];
EventLat(ZeroInd) = [];
EventLon(ZeroInd) = [];
EventTime(ZeroInd) = [];
EventAbsTime(ZeroInd) = [];
EventAzi(ZeroInd) = [];


%----------------------------
%clean up vectors with h < 20 km
HeightInd = find(EventHeight < 20);

Delta_H(HeightInd) = [];
EventHeight(HeightInd) = [];
EventLat(HeightInd) = [];
EventLon(HeightInd) = [];
EventTime(HeightInd) = [];
EventAbsTime(HeightInd) = [];
EventAzi(HeightInd) = [];

if isempty(Delta_H) == 1 %if Delta_H is empty we have zero Es events
    NumberOfEvents = 0;
end


if PlotIndex == 1
    % %plot to analyze each profile better:
    plot(h_t_pos,caL1Snr_pos/10,'b-','HandleVisibility','off')
    hold on;
    plot(h_t_pos,A1/10,'r-')
    plot(h_t_pos,A20/10,'y-')
    plot(h_t_pos(FirstCriteriaIndex(1:3:end)),caL1Snr_pos(FirstCriteriaIndex(1:3:end)),'k.','HandleVisibility','off')
    %plots every third element of FirstCriteriaIndex just for less density, but
    %shows where the algorithm detects an event.
    xlabel('Straigt-line-tangent height [km]')
    ylabel('SNR [V/V]')
    legend('1km mean','20km mean')
    title({['Thickness [km]: ' num2str(Delta_H,3)],['Height [km]: ' num2str(EventHeight,3)]})
end




end