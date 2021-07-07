%plots output of multiple detections in february 2009

load('/zhome/e8/9/144512/Desktop/MScThesis/MatlabScripts/SavedMatlabVariables/Misc/2009_2MultipleDetectionTest2.mat')     
        

MultipleDetectionsFilenames = string(MultipleDetectionsFilestamps);
%[UniqueFilestamps ia ic] = unique(MultipleDetectionsFilestamps);

[UniqueHeight ia ic] = unique(cell2mat(MultipleDetectionsHeight));
for i = 1:length(ia);
    id(i) = nnz(ic == i); %id is the number of how often the same height appears, denotes multiple detections.
end
  


%Find the filestamp and thickness of each corresponding element in unique height
MultipleDetectionsThickness = cell2mat(MultipleDetectionsThickness);
UniqueFilenames = [];
UniqueThickness = [];
for p = 1:length(UniqueHeight)
    index = find(UniqueHeight(p) == cell2mat(MultipleDetectionsHeight)); %index elements contain same height values
    UniqueFilenames = [UniqueFilenames MultipleDetectionsFilenames(index(1))];
    UniqueThickness = [UniqueThickness MultipleDetectionsThickness(index(1))];
end



%%

%edited to analyse all profiles in a month to extract valid single detections to show plots in thesis of how the U-shape algorithm works
%load('/zhome/e8/9/144512/Desktop/MScThesis/MatlabScripts/SavedMatlabVariables/MonthsMainRunAll/CombineMonths/2009_2.mat')     
      
year = '2009';
% %UniqueFilenames = FilestampCellEdetections;
% %UniqueFilenames = ["C001.2009.032.18.46.G17", "C001.2009.032.00.36.G25","C001.2009.032.00.36.G32"];
% UniqueFilenames{1} = 'C002.2009.032.10.56.G06';
% UniqueFilenames{2} = 'C002.2009.032.19.53.G06';
% UniqueFilenames{3} = 'C001.2009.032.00.36.G32';

%Directory of the data with the filestamps
MyDataDir = dir("/work1/s192715/Data/Months/2009_2");
MyDataDirString = "/work1/s192715/Data/Months/2009_2";


%To only look at profiles with more than 2-fold detections, run:
% Index2fold = find(id==2);
% UniqueHeight(Index2fold) = [];
% UniqueThickness(Index2fold) = [];
% UniqueFilenames(Index2fold) = [];
% id(Index2fold) = [];

for i = 1:length(UniqueFilenames)

    FilestampString = UniqueFilenames{i};
    FilestampDay = FilestampString(11:13);
    
    FileDir = ['atmPhs_repro2013_' year '_' FilestampDay];
    FileName = ['atmPhs_' FilestampString '_2013.3520_nc'];
    FileNameFullPath = [MyDataDirString + '/' + FileDir + '/' + FileName];
    
    
    
    ncid = netcdf.open(FileNameFullPath);
    [h_t_pos latitude_pos longitude_pos caL1Snr_pos A1 A20] = CalcHeightLatLon(ncid);
    [Delta_H, EventHeight, EventTime, EventAbsTime, EventLat, EventLon, EventAzi, Es_present, NumberOfEvents, S4_present, Both_present, Neither_present, OnlyEs_present, OnlyS4_present, RO_lat, RO_lon, RO_Time, RO_MaxS4, RO_MeanS4, RO_HeightOfS4Max, Filestamp] = MScThesisScriptCombination4(ncid,0);
    fileStamp = netcdf.getAtt(ncid,netcdf.getConstant('NC_GLOBAL'),'fileStamp');
    
    
    
    if isempty(find(h_t_pos >= 30)) == 1 %if don't have any heights above 30km
        plot(h_t_pos,caL1Snr_pos/10,'b-')
        hold on;
        xlabel('Straight-line-tangent height [km]')
        ylabel('SNR [V/V]')
        title({[FilestampString '. NoOfDetections: ' num2str(id(i)) '.'],['Height: ' num2str(UniqueHeight(i)) '. Thickness: ' num2str(UniqueThickness(i)) '.']})
        %title({[FilestampString '.'],['Height: ' num2str(EventHeight) '. Thickness: ' num2str(Delta_H) '.']})
        plot(h_t_pos,A1/10,'r-','Linewidth',2)
        plot(h_t_pos,A20/10,'y-','Linewidth',2)
        plot(h_t_pos,0.6.*A20/10,'k-','Linewidth',2)
        legend('caL1SNR','1km mean','20 km mean','0.6*20km mean')
        
        pause
        close
        
    else
        
        Below20kmIndex = find(h_t_pos < 30);

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
        title({[FilestampString '. NoOfDetections: ' num2str(id(i)) '.'],['Height: ' num2str(UniqueHeight(i)) '. Thickness: ' num2str(UniqueThickness(i)) '.']})
        %title({[FilestampString '.'],['Height: ' num2str(EventHeight) '. Thickness: ' num2str(Delta_H) '.']})
        plot(h_t_pos_Above20km,A1_Above20km/10,'r-','LineWidth',2)
        plot(h_t_pos_Above20km,A20_Above20km/10,'y-','LineWidth',2)
        plot(h_t_pos_Above20km,0.6.*A20_Above20km/10,'k-','LineWidth',2)
        legend('caL1SNR','1km mean','20 km mean','0.6*20km mean')
        
%         if i == 1
%         legend('caL1SNR','1km mean','20 km mean','0.6*20km mean')
%         han=axes(fig,'visible','off'); 
%         han.XLabel.Visible='on';
%         han.YLabel.Visible='on';
%         ylabel(han,'L1 C/A SNR [V/V]');
%         xlabel(han,'Straight-line tangent height [km]');
%         end
        
        pause
        close 
    end
    
end
 