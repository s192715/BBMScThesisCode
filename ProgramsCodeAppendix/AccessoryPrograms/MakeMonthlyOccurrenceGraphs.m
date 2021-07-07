%this was just to make simple line graphs of monthly occurrence, not
%significant script.

MyDataDir = dir("/zhome/e8/9/144512/Desktop/MScThesis/MatlabScripts/SavedMatlabVariables/MonthsMain2/MakeMonthlyOccurrenceGraphs");
MyDataDirString = "/zhome/e8/9/144512/Desktop/MScThesis/MatlabScripts/SavedMatlabVariables/MonthsMain2/MakeMonthlyOccurrenceGraphs";

S4_Occ = zeros(1,length(MyDataDir)-2);
Es_Occ = zeros(1,length(MyDataDir)-2);
S4_Amp = zeros(1,length(MyDataDir)-2);
for k = 3:length(MyDataDir) %for each directory in MyDataDir (skips first 2 elements because they are just . and ..)

    Month = dir(MyDataDirString + "/" + MyDataDir(k).('name'));
    MonthString = MyDataDirString + "/" + MyDataDir(k).('name');
    load(MonthString)
    
    S4_Occ(k-2) = (100*sum(TotalS4_present)/length(TotalS4_present));
    Es_Occ(k-2) = (100*sum(TotalEs_present)/length(TotalEs_present));
    
    Index_Above03 = TotalRO_MaxS4 > 0.3;
    S4_Amp(k-2) = (sum(TotalRO_MaxS4(Index_Above03))/length(TotalRO_MaxS4(Index_Above03)));
       
    
end

%to make months appear in correct order
Index = [1 5 6 7 8 9 10 11 12 2 3 4 13 17 18 19 20 21 22 23 24 14 15 16 25 29 30 31 32 33 34 35 36 26 27 28 37 41 42 43 44 45 46 47 48 38 39 40 49 53 54 55 56 57 58 59 60 50 51 52 61 65 66 67 68 69 70 71 72 62 63 64 73 77 78 79 80 81 82 83 84 74 75 76];                                                  
S4_Occ = S4_Occ(Index);
Es_Occ = Es_Occ(Index);
S4_Amp = S4_Amp(Index);


figure(1); plot(1:length(Index),S4_Occ, 'b-*')
xlabel('Months')
xticks = 1:12:84;
xlabels = {"Jan '07" "Jan '08" "Jan '09" "Jan '10" "Jan '11" "Jan '12" "Jan '13"}; 
set(gca, 'XTick', xticks, 'XTickLabel', xlabels);
xtickangle(45)
ylabel('%')
title('Monthly S4 Global Occurrence')

figure(2); plot(1:length(Index),Es_Occ, 'b-*')
xlabel('Months')
xticks = 1:12:84;
xlabels = {"Jan '07" "Jan '08" "Jan '09" "Jan '10" "Jan '11" "Jan '12" "Jan '13"}; 
set(gca, 'XTick', xticks, 'XTickLabel', xlabels);
xtickangle(45)
ylabel('%')
title('Monthly Es Global Occurrence')

figure(3); plot(1:length(Index),S4_Amp, 'b-*')
xlabel('Months')
xticks = 1:12:84;
xlabels = {"Jan '07" "Jan '08" "Jan '09" "Jan '10" "Jan '11" "Jan '12" "Jan '13"}; 
set(gca, 'XTick', xticks, 'XTickLabel', xlabels);
xtickangle(45)
ylabel('S4')
title('Monthly S4 Global Mean Amplitude')

