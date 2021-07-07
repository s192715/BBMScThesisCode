%This script is made to analyze the FilestampCell variable that is
%outputted in version 3 of MScThesisMainScript3.m and
%MScThesisScriptCombination3.m, where a cutoff for the low SNR has been
%applied.

%loading in data which includes the FilestampCell variable.
%load('/zhome/e8/9/144512/Desktop/MScThesis/MatlabScripts/SavedMatlabVariables/ROcriteriaTest/2009July_FirstTry.mat')
year = '2009';

%Directory of the data with the filestamps
MyDataDir = dir("/work1/s192715/Data/Months/2009_2");
MyDataDirString = "/work1/s192715/Data/Months/2009_2";
 
%MeanSNR = [];
NewMeanSNR = [];
for i = 1:1:length(FilestampCell) %plotting every 100th file
    FilestampString = FilestampCell{i};
    FilestampDay = FilestampString(11:13);
    
    FileDir = ['atmPhs_repro2013_' year '_' FilestampDay];
    FileName = ['atmPhs_' FilestampString '_2013.3520_nc'];
    FileNameFullPath = [MyDataDirString + '/' + FileDir + '/' + FileName];
    
    
    ncid = netcdf.open(FileNameFullPath);
    [h_t_pos latitude_pos longitude_pos caL1Snr_pos A1 A20] = CalcHeightLatLon(ncid);
       
    if isempty(find(h_t_pos >= 20)) == 1
        plot(h_t_pos,caL1Snr_pos/10,'b-')
        hold on;
        xlabel('Straight-line-tangent height [km]')
        ylabel('SNR [V/V]')
        title(FilestampString)
        plot(h_t_pos,A1/10,'r-')
        plot(h_t_pos,A20/10,'y-')
        legend('caL1SNR','1km mean','20 km mean')
        
        pause
        close
        
    else
        
        Below20kmIndex = find(h_t_pos < 20);

        A1_Above20km = A1;
        A1_Above20km(Below20kmIndex) = [];
        A20_Above20km = A20;
        A20_Above20km(Below20kmIndex) = [];
        h_t_pos_Above20km = h_t_pos;
        h_t_pos_Above20km(Below20kmIndex) = [];
        caL1Snr_pos_Above20km = caL1Snr_pos;
        caL1Snr_pos_Above20km(Below20kmIndex) = [];

        plot(h_t_pos_Above20km,caL1Snr_pos_Above20km/10,'b-')
        hold on;
        xlabel('Straight-line-tangent height [km]')
        ylabel('SNR [V/V]')
        title(FilestampString)
        plot(h_t_pos_Above20km,A1_Above20km/10,'r-')
        plot(h_t_pos_Above20km,A20_Above20km/10,'y-')
        legend('caL1SNR','1km mean','20 km mean')

        pause
        close
    end
    
    
%     Below20kmIndex = find(h_t_pos < 20);
% 
%     A1_Above20km = A1;
%     A1_Above20km(Below20kmIndex) = [];
%     A20_Above20km = A20;
%     A20_Above20km(Below20kmIndex) = [];
%     h_t_pos_Above20km = h_t_pos;
%     h_t_pos_Above20km(Below20kmIndex) = [];
%     caL1Snr_pos_Above20km = caL1Snr_pos;
%     caL1Snr_pos_Above20km(Below20kmIndex) = [];
%     
%     plot(h_t_pos_Above20km,caL1Snr_pos_Above20km/10,'b-')
%     hold on;
%     xlabel('Straight-line-tangent height [km]')
%     ylabel('SNR [V/V]')
%     title(FilestampString)
%     plot(h_t_pos_Above20km,A1_Above20km/10,'r-')
%     plot(h_t_pos_Above20km,A20_Above20km/10,'y-')
%     legend('caL1SNR','1km mean','20 km mean')
%     
%     %minimumA20 = min(A20_Above20km./10)
%     %MeanSNR = [MeanSNR mean(caL1Snr_pos(60 < h_t_pos < 80))];
%     %NewMeanSNR = [NewMeanSNR mean(caL1Snr_pos_Above20km)];
%     
%     pause
%     close
end
 

% histogram(NewMeanSNR/10,[0:4:200], 'Normalization','probability')
% xlabel('SNR [V/V]')
% ylabel('%')
% title('July 2009. Mean SNR of excluded profiles. Bin size: 4 V/V.')

