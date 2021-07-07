%finding double detections in the Es algorithm, seeing if they are common
%to low h_t height (could then typically be because of atmosphere, if so maybe then increase 20km cutoff)

%CONCLUSION:
%mean(MultipleDetectionsHeight) is 100 km, so double or multiple detections
%do not seem to be because atmosphere at low tangent height.

%This profile had 5-fold detection.
%load('/zhome/e8/9/144512/Desktop/MScThesis/MatlabScripts/SavedMatlabVariables/Misc/SingleSavedProfileTest.mat')

%loading February 2009
load('/zhome/e8/9/144512/Desktop/MScThesis/MatlabScripts/SavedMatlabVariables/MonthsMain2/CombineMonths/2009_2.mat')
Month = 'February'
year = 2009;

%testing to detect values of same thickness.
au        = unique(TotalDelta_H); %finds unique thickness values
[N, Edge] = histcounts(TotalDelta_H, au); %if f.ex. thickness value is repeated 2 times, N will be 2.
ThicknessValuesRepeated  = au(N > 1); %thickness values that are repeared
HowOftenRepeated = N(N > 1);

k = 1;
for i = 1:length(ThicknessValuesRepeated)
    index = find(TotalDelta_H == ThicknessValuesRepeated(i));
    if length(index) >= 2
        MarkIndex{1,k} = index;
        k = k+1;
    end
end
%Each cell in MarkIndex contains the multiple detection Event indices,
%i.e. we can see at which thickness and height we had multiple detections by TotalDelta_H(MarkIndex{1,1})



for j = 1:length(MarkIndex)
    MultipleDetectionsThickness(j) = mean(TotalDelta_H(MarkIndex{1,j})); %this is a vector of same values, so the mean is just to extract one of the same value numbers.
    MultipleDetectionsHeight(j) = mean(TotalEventHeight(MarkIndex{1,j})); %heights of the events that had multiple detections
    MultipleDetectionsLat(j) = mean(TotalEventLat(MarkIndex{1,j}));
    MultipleDetectionsLocalTime(j) = mean(TotalEventLocalTime(MarkIndex{1,j}));
end

figure(10); subplot(2,2,1)
histogram(MultipleDetectionsThickness,'BinWidth',1)
xlabel('Thickness [km]')

subplot(2,2,2)
histogram(MultipleDetectionsHeight,'BinWidth',10)
xlabel('Height [km]')

subplot(2,2,3)
histogram(MultipleDetectionsLat,'BinWidth',10)
xlabel('Latitude [°]')

subplot(2,2,4)
histogram(MultipleDetectionsLocalTime,'BinWidth',1)
xlabel('Local time [Hours]')
suptitle({[Month ' ' num2str(year) '.' ' Total detections: ' num2str(length(TotalDelta_H)) '.'],['Multi-detections: ' num2str(length(MarkIndex)) '. Double-detections: ' num2str(sum(HowOftenRepeated == 2)) '.']})



%  a = ['aaa'; 'aaa'; 'bbb'; 'bbb'; 'bbb'; 'bbb'; 'bbb'; 'bbb'; 'ccc'; 'ccc'; 'ccc'; 'ccc'; 'ddd'; 'ddd'; 'ddd']
%  [C ia ic] = unique(a,'rows') %gefur rétt
% 
%  for i = 1:length(ia);
%      id(i) = nnz(ic == i); %id is the number of how often the same string appears, denotes multiple detections.
%  end
%  %C contains the name of filenames having the profile of corresponding
%  %number of multiple detections in id.
 