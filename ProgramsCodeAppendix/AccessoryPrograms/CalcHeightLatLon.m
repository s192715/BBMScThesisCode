function [h_t_pos latitude_pos longitude_pos caL1Snr_pos A1 A20] = CalcHeightLatLon(ncid);

%The function is used in PlotFilestampCell.m
%function that returns h_t_pos, latitude_pos, longitude_pos, and caL1Snr_pos of a RO profile
%for h_t above 20 km. Also returns A1 and A20 to include in PlotFilestampCell.m

time = netcdf.getVar(ncid,0);
caL1Snr = netcdf.getVar(ncid,1);

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
% zet = (2306.2181*tu + 0.30188*tu^2 + 0.017998*tu^3)*pi/648000;
% z   = (2306.2181*tu + 1.09468*tu^2 + 0.018203*tu^3)*pi/648000;
% ups = (2004.3109*tu - 0.42665*tu^2 - 0.041833*tu^3)*pi/648000;
% szet = sin(zet);
% czet = cos(zet);
% sz   = sin(z);
% cz   = cos(z);
% sups = sin(ups);
% cups = cos(ups);                 
% RP = [cz*cups*czet-sz*szet, -cz*cups*szet-sz*czet, -cz*sups; sz*cups*czet+cz*szet , -sz*cups*szet+cz*czet , -sz*sups; sups*czet, -sups*szet, cups];
% Qt = RP; %this is the rotation matrix
% 
% %intermediate ECI coordinates
% r_LEO_ECI_Interm = Qt*r_LEO_ECI;
% r_GPS_ECI_Interm = Qt*r_GPS_ECI;

%--------------------
%From eci2ecf_nosofa
%formulas from Astronomical Almanac, 1993.
gmst0 = 24110.54841 + 8640184.812866*tu + 0.093104*tu.^2 - (6.2*10^(-6))*tu.^3;
gmst  = mod(gmst0 + utc_seconds*1.0027379093, 86400); %mod returns remainder after division
 
phi = gmst*2*pi/86400; %in radians

pa = repmat([0 0 1], [length(time),1]); %polar axis (in format to match each time step coordinate)
r_LEO_ECEF_test = vector_rotate(r_LEO_ECI', pa, -phi); %ECI -> ECF, negative phi.
r_GPS_ECEF_test = vector_rotate(r_GPS_ECI', pa, -phi);

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

if isempty(find(h_t >= 20)) == 1
    caL1Snr_pos = caL1Snr;
    h_t_pos = h_t;
    A1 = movmean(caL1Snr_pos,1,'SamplePoints',h_t_pos);
    A20 = movmean(caL1Snr_pos,20,'SamplePoints',h_t_pos);
    latitude_pos = latitude;
    longitude_pos = longitude;
    
    return
end

%Extracting only the values corresponding to above h_t 20km
NegativeValuesIndex = find(h_t < 0);

%taking into account if measurement actually doesn't contain any
%height values above zero (happened at least for one RO measurement, n=1172 in 2014_001).
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


end

