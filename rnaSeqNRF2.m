%% clear
clear; 
clc;
%% Multiple myeloma Cells (MM) RPMI-8226
%Getting the data files 
    excelDataMM = '/Users/harshitakumar/Documents/Research/Letterio Research/Letterio Data/RNA Seq/RNA-Seq MM 2.xlsx';
    [~, text, dataMM] = xlsread(excelDataMM);
%grabbing the specific data for each concentration group 
    %[controlDataMM, controlData2MM] = deal(cell2mat(dataMM(2:end, 8:10)));
    [controlDataMM, controlData2MM] = deal(cell2mat(dataMM(2:end, 9:10)));
    controlDataPCMM = deal(cell2mat(dataMM(2:end, 8:10)));
    [lowConcDataMM, lowConcData2MM, lowConcData3MM] = deal(cell2mat(dataMM(2:end, 11:13)));
    [highConcDataMM, highConcData2MM]= deal(cell2mat(dataMM(2:end, 14:16)));
%calculating the mean
    controlMeanMM = mean(controlDataMM, 2);
    lowConcMeanMM = mean(lowConcDataMM, 2);
    highConcMeanMM = mean(highConcDataMM, 2);
%deleting rows in which the mean values for each group is less than 1
    %group1 = control + low; group 2 = control + high
    toDelete1 = (controlMeanMM < 1) & (lowConcMeanMM < 1);
    controlDataMM(toDelete1,:) = [];
    lowConcDataMM(toDelete1,:) = [];
    toDelete2 = (controlMeanMM < 1) & (highConcMeanMM < 1);
    controlData2MM(toDelete2,:) = [];
    highConcDataMM(toDelete2,:) = [];
    controlDataPCMM(toDelete2,:) = [];
    lowConcData3MM(toDelete2,:) = [];
    toDelete3 = (lowConcMeanMM < 1) & (highConcMeanMM < 1);
    lowConcData2MM(toDelete3,:) = [];
    highConcData2MM(toDelete3,:) = [];
%creating labels without the irrelevant data points
    [labelLowMM, labelHighMM, labelBothMM] = deal(text(2:end,6));
    labelLowMM(toDelete1, :) = [];
    labelHighMM(toDelete2, :) = [];
    labelBothMM(toDelete3, :) = [];
%mean 
    controlMeanMM = mean(controlDataMM, 2);
    controlMean2MM = mean(controlData2MM, 2); 
    controlMean2PCMM = mean(controlDataPCMM, 2);
    lowConcMeanMM = mean(lowConcDataMM, 2);
    lowConcMean2MM = mean(lowConcData2MM, 2);
    lowConcMean3MM = mean(lowConcData3MM, 2);
    highConcMeanMM = mean(highConcDataMM, 2);
    highConcMean2MM = mean(highConcData2MM, 2);
%standard deviations
    controlSDMM = std(controlDataMM, 0, 2);
    controlSD2MM = std(controlData2MM, 0, 2);
    controlSD2PCMM = std(controlDataPCMM, 0, 2);
    lowConcSDMM = std(lowConcDataMM, 0, 2);
    lowConcSD2MM = std(lowConcData2MM, 0, 2);
    lowConcSD3MM = std(lowConcData3MM, 0, 2);
    highConcSDMM = std(highConcDataMM, 0, 2);    
    highConcSD2MM = std(highConcData2MM, 0, 2);
%calculating the P-values
    p_Value1MM = mattest(controlDataMM, lowConcDataMM, 'VarType', 'equal',  'Labels', labelLowMM);
    controlHighp_ValueMM  = mattest(controlData2MM, highConcDataMM, 'VarType', 'equal', 'Labels', labelHighMM);
    controlHighp_ValuePCMM  = mattest(controlDataPCMM, highConcDataMM, 'VarType', 'equal', 'Labels', labelHighMM);
    controlLowp_ValueMM =  mattest(controlData2MM, lowConcData3MM, 'VarType', 'equal', 'Labels', labelHighMM);
    controlLowp_ValuePCMM =  mattest(controlDataPCMM, lowConcData3MM, 'VarType', 'equal', 'Labels', labelHighMM);
    lowHighp_ValueMM = mattest(lowConcMean3MM, highConcDataMM, 'VarType', 'equal', 'Labels', labelHighMM);
    p_Value2MM  = mattest(lowConcData2MM, highConcData2MM, 'VarType', 'equal', 'Labels', labelHighMM);
%tables 
    %controlLowTableMM = table(labelLowMM, controlMeanMM, controlSDMM, lowConcMeanMM, lowConcSDMM, p_Value1MM, lowConcMeanMM./controlMeanMM, controlDataMM, lowConcDataMM);
    controlHighTableMM = table(labelHighMM, controlMean2MM, lowConcMean3MM, highConcMeanMM, controlSD2MM, lowConcSD3MM, highConcSDMM, controlLowp_ValueMM, lowHighp_ValueMM, controlHighp_ValueMM, lowConcMean3MM./controlMean2MM, highConcMeanMM./lowConcMean3MM,  highConcMeanMM./controlMean2MM,controlData2MM, lowConcData3MM, highConcDataMM);
    controlHighTable2MM = table(labelHighMM, controlMean2PCMM, lowConcMean3MM, highConcMeanMM, controlSD2PCMM, lowConcSD3MM, highConcSDMM, controlLowp_ValuePCMM, lowHighp_ValueMM, controlHighp_ValuePCMM, lowConcMean3MM./controlMean2PCMM, highConcMeanMM./lowConcMean3MM,  highConcMeanMM./controlMean2PCMM,controlDataPCMM, lowConcData3MM, highConcDataMM);
    %controlBothTableMM = table(labelBothMM, lowConcMean2MM, lowConcSD2MM, highConcMean2MM, highConcSD2MM, p_Value2MM, highConcMean2MM./lowConcSD2MM, lowConcData2MM, highConcData2MM);
    [controlLowTableMM.Properties.VariableNames{9}, controlBothTableMM.Properties.VariableNames{9}] = deal('LowHighFoldChangeMM');
    controlHighTableMM.Properties.VariableNames{11} = 'ControlLowFoldchangeMM';
    controlHighTableMM.Properties.VariableNames{12} = 'LowHighFoldChangeMM';
    controlHighTableMM.Properties.VariableNames{13} = 'ControlHighFoldchangeMM';
% %volcano calculations and plots
%     header = {'MM Gene Names'; 'MM P-Value'; 'MM Absolute Fold Change'};
%     %control vs low
%         volcanoLowMM = mavolcanoplot(lowConcMeanMM , controlMeanMM, p_Value1MM, 'Labels', labelLowMM, 'LogTrans', true, 'Plotonly', true);
%         title('Multiple Myeloma cells comparison of Control and 100nM (low Concentration) 2P')
%         lowTableMM = table(volcanoLowMM.GeneLabels, volcanoLowMM.PValues, negToAbs(volcanoLowMM.FoldChanges), 'VariableNames', header); 
%     %control vs high
%         volcanoHighMM = mavolcanoplot(highConcMeanMM, controlMean2MM, controlHighp_ValueMM, 'Labels', labelHighMM, 'LogTrans', true);
%         title('Multiple Myeloma cells comparison of Control and 400nM (high Concentration) 2P')
%         highTableMM = table(volcanoHighMM.GeneLabels, volcanoHighMM.PValues, negToAbs(volcanoHighMM.FoldChanges), 'VariableNames', header);
%         [~, posm, ~] = intersect(labelHighMM, volcanoHighMM.GeneLabels);
%         controlHighDiffMM = controlHighTableMM(posm, :);
%     %low vs high  
%         volcanoBothMM = mavolcanoplot(highConcMean2MM, lowConcMean2MM, p_Value2MM, 'Labels', labelBothMM, 'LogTrans', true);
%         title('Multiple Myeloma cells comparison of 100nM (low Concentration) and 400nM (high Concentration) 2P')
%         bothTableMM = table(volcanoBothMM.GeneLabels, volcanoBothMM.PValues, negToAbs(volcanoBothMM.FoldChanges), 'VariableNames', header);
clear excelDataMM dataMM 
clear controlDataMM controlData2MM lowConcDataMM lowConcData2MM lowConcData3MM highConcDataMM highConcData2MM 
clear toDelete1 toDelete2 toDelete3 text labelLowMM labelHighMM labelBothMM
clear p_Value1MM p_Value2MM controlHighp_ValueMM controlLowp_ValueMM lowHighp_ValueMM
clear controlMeanMM controlMean2MM lowConcMeanMM lowConcMean2MM lowConcMean3MM highConcMeanMM highConcMean2MM
clear controlSDMM controlSD2MM lowConcSDMM lowConcSD2MM lowConcSD3MM highConcSDMM highConcSD2MM
%% DIPG Cells SF8628 
%Getting the data files 
    excelDataDIPG = '/Users/harshitakumar/Documents/Research/Letterio Research/Letterio Data/RNA Seq/4_genes_fpkm_expression_DIPG.xlsx';
    [~, ~, dataDIPG] = xlsread(excelDataDIPG);
%grabbing the specific data for each concentration group
    labelLowDIPG = dataDIPG(2:end,6);
    clear controlDataDIPG controlData2DIPG
    [controlDataDIPG, controlData2DIPG] = deal(cell2mat(dataDIPG(2:end, 14:16)));
    [lowConcDataDIPG, lowConcData2DIPG, lowConcData3DIPG] = deal(cell2mat(dataDIPG(2:end, 17:19)));
    [highConcDataDIPG,highConcData2DIPG] = deal(cell2mat(dataDIPG(2:end, 20:22)));
%calculating the mean
    controlMeanDIPG = mean(controlDataDIPG, 2);
    lowConcMeanDIPG = mean(lowConcDataDIPG, 2);
    highConcMeanDIPG = mean(highConcDataDIPG, 2);
%deleting rows in which the mean values for each group is less than 1
    %group1 = control + low; group 2 = control + high
    toDelete4 = (controlMeanDIPG < 1) & (lowConcMeanDIPG < 1);
    controlDataDIPG(toDelete4,:) = [];
    lowConcDataDIPG(toDelete4,:) = [];
    toDelete5 = (controlMeanDIPG < 1) & (highConcMeanDIPG < 1);
    controlData2DIPG(toDelete5,:) = [];
    highConcDataDIPG(toDelete5,:) = [];
    lowConcData3DIPG(toDelete5,:) = [];
    toDelete6 = (lowConcMeanDIPG < 1) & (highConcMeanDIPG < 1);
    lowConcData2DIPG(toDelete6,:) = [];
    highConcData2DIPG(toDelete6,:) = [];
%recalculating the mean
    controlMeanDIPG = mean(controlDataDIPG, 2);
    controlMean2DIPG = mean(controlData2DIPG, 2);
    lowConcMeanDIPG = mean(lowConcDataDIPG, 2);
    lowConcMean2DIPG = mean(lowConcData2DIPG, 2);
    lowConcMean3DIPG = mean(lowConcData3DIPG, 2);
    highConcMeanDIPG = mean(highConcDataDIPG, 2);
    highConcMean2DIPG = mean(highConcData2DIPG, 2);
%calculating the standard deviation
    controlSDDIPG = std(controlDataDIPG, 0, 2);
    controlSD2DIPG = std(controlData2DIPG, 0, 2);
    lowConcSDDIPG = std(lowConcDataDIPG, 0, 2);
    lowConcSD2DIPG = std(lowConcData2DIPG, 0, 2);
    lowConcSD3DIPG = std(lowConcData3DIPG, 0, 2);
    highConcSDDIPG = std(highConcDataDIPG, 0, 2);
    highConcSD2DIPG = std(highConcData2DIPG, 0, 2);
%creating labels without the irrelevant data points
    labelLowhDIPG = dataDIPG(2:end, 6);
    labelLowDIPG(toDelete4, :) = [];
    labelHighDIPG = dataDIPG(2:end, 6);
    labelHighDIPG(toDelete5, :) = [];
    labelBothDIPG = dataDIPG(2:end, 6);
    labelBothDIPG(toDelete6, :) = [];
    clear toDelete4 toDelete5 toDelete6
%calculating p-value
    p_ValueDIPG2 = mattest(controlDataDIPG, lowConcDataDIPG, 'VarType', 'equal',  'Labels', labelLowDIPG);
    controlHighp_ValueDIPG  = mattest(controlData2DIPG, highConcDataDIPG, 'VarType', 'equal',  'Labels', labelHighDIPG);
    controlLowp_ValueDIPG  = mattest(controlData2DIPG, lowConcData3DIPG, 'VarType', 'equal',  'Labels', labelHighDIPG);
    lowHighp_ValueDIPG  = mattest(lowConcData3DIPG, highConcDataDIPG, 'VarType', 'equal',  'Labels', labelHighDIPG);
    p_ValueDIPG3  = mattest(lowConcData2DIPG, highConcData2DIPG, 'VarType', 'equal', 'Labels', labelHighDIPG);
%tables 
    controlLowTableDIPG = table(labelLowDIPG, controlMeanDIPG, controlSDDIPG, lowConcMeanDIPG, lowConcSDDIPG, p_ValueDIPG2, lowConcMeanDIPG./controlMeanDIPG, controlDataDIPG, lowConcDataDIPG);
    controlHighTableDIPG = table(labelHighDIPG, controlMean2DIPG, lowConcMean3DIPG, highConcMeanDIPG, controlSD2DIPG,lowConcSD3DIPG, highConcSDDIPG, controlLowp_ValueDIPG,lowHighp_ValueDIPG, controlHighp_ValueDIPG,lowConcMean3DIPG./controlMean2DIPG, highConcMeanDIPG./lowConcMean3DIPG, highConcMeanDIPG./controlMean2DIPG,controlData2DIPG, lowConcData3DIPG, highConcDataDIPG);
    controlBothTableDIPG = table(labelBothDIPG, lowConcMean2DIPG, lowConcSD2DIPG, highConcMean2DIPG, highConcSD2DIPG, p_ValueDIPG3, highConcMean2DIPG./lowConcSD2DIPG, lowConcData2DIPG, highConcData2DIPG);
    [controlLowTableDIPG.Properties.VariableNames{9}, controlBothTableDIPG.Properties.VariableNames{9}] = deal('LowHighFoldChangeDIPG');
    controlHighTableDIPG.Properties.VariableNames{11} = 'ControlLowFoldchangeDIPG';
    controlHighTableDIPG.Properties.VariableNames{12} = 'LowHighFoldChangeDIPG';
    controlHighTableDIPG.Properties.VariableNames{13} = 'ControlHighFoldchangeDIPG';
% %volcano calculations and plots
%     header2 = {'DIPG Gene Names'; 'DIPG P-Value'; 'DIPG Absolute Fold Change'};
%     %control vs low
%         volcanoLowDIPG = mavolcanoplot(lowConcMeanDIPG, controlMeanDIPG, p_ValueDIPG, 'Labels', labelLowDIPG, 'LogTrans', true);
%         title('DIPG Cells comparison of Control and Low Concentration 2P')
%         lowTableDIPG = table(volcanoLowDIPG.GeneLabels, volcanoLowDIPG.PValues, negToAbs(volcanoLowDIPG.FoldChanges), 'VariableNames', header2); 
%     %control vs high
%         volcanoHighDIPG = mavolcanoplot(highConcMeanDIPG, controlMean2DIPG, controlHighp_ValueDIPG, 'Labels', labelHighDIPG, 'LogTrans', true);
%         title('DIPG Cells comparison of Control and High Concentration 2P')
%         highTableDIPG = table(volcanoHighDIPG.GeneLabels, volcanoHighDIPG.PValues, negToAbs(volcanoHighDIPG.FoldChanges), 'VariableNames', header2); 
%     %low vs high   
%         volcanoBothDIPG = mavolcanoplot(highConcMean2DIPG, lowConcMean2DIPG, p_ValueDIPG3, 'Labels', labelBothDIPG, 'LogTrans', true);
%         title('DIPG Cells comparison of low and High Concentration 2P')
%         bothTableDIPG = table(volcanoBothDIPG.GeneLabels, volcanoBothDIPG.PValues, negToAbs(volcanoBothDIPG.FoldChanges), 'VariableNames', header2); 
clear excelDataDIPG dataDIPG p_ValueDIPG3 p_ValueDIPG2 labelLowhDIPG
clear controlDataDIPG controlData2DIPG lowConcDataDIPG lowConcData2DIPG lowConcData3DIPG highConcDataDIPG highConcData2DIPG 
clear toDelete1 toDelete2 toDelete3 text labelLowDIPG labelHighDIPG labelBothDIPG
clear p_Value1DIPG p_Value2DIPG controlHighp_ValueDIPG controlLowp_ValueDIPG lowHighp_ValueDIPG
clear controlMeanDIPG controlMean2DIPG lowConcMeanDIPG lowConcMean2DIPG lowConcMean3DIPG highConcMeanDIPG highConcMean2DIPG
clear controlSDDIPG controlSD2DIPG lowConcSDDIPG lowConcSD2DIPG lowConcSD3DIPG highConcSDDIPG highConcSD2DIPG
%% Ovarian OVCAR2
% Cancer cells
%Getting the data files 
    excelDataOv = '/Users/harshitakumar/Documents/Research/Letterio Research/Letterio Data/RNA Seq/RNA-Seq Ovarian 2.xlsx';
    [~, ~, dataOv] = xlsread(excelDataOv);
    lengthOv = 60612;
    clear excelDataOv
%getting the specific data for each group
    [controlDataOv, controlData2Ov] = deal(cell2mat(dataOv(2:end, 3:5)));
    [lowConcDataOv, lowConcData2Ov, lowConcData3Ov] = deal(cell2mat(dataOv(2:end, 6:8)));
    [highConcDataOv,highConcData2Ov] = deal(cell2mat(dataOv(2:end, 9:11)));
%calculating the mean
    controlMeanOv = mean(controlDataOv, 2);
    lowConcMeanOv = mean(lowConcDataOv, 2);
    highConcMeanOv = mean(highConcDataOv, 2);
%deleting rows in which the mean values for each group is less than 1
    %group1 = control + low; group 2 = control + high
    toDelete7 = (controlMeanOv < 1) & (lowConcMeanOv < 1);
    controlDataOv(toDelete7,:) = [];
    lowConcDataOv(toDelete7,:) = [];
    toDelete8 = (controlMeanOv < 1) & (highConcMeanOv < 1);
    controlData2Ov(toDelete8,:) = [];
    highConcDataOv(toDelete8,:) = [];
    lowConcData3Ov(toDelete8,:) = [];
    toDelete9 = (lowConcMeanOv < 1) & (highConcMeanOv < 1);
    lowConcData2Ov(toDelete9,:) = [];
    highConcData2Ov(toDelete9,:) = [];
%recalculating the mean
    controlMeanOv = mean(controlDataOv, 2);
    controlMean2Ov = mean(controlData2Ov, 2);
    lowConcMeanOv = mean(lowConcDataOv, 2);
    lowConcMean2Ov = mean(lowConcData2Ov, 2);
    lowConcMean3Ov = mean(lowConcData3Ov, 2);
    highConcMeanOv = mean(highConcDataOv, 2);
    highConcMean2Ov = mean(highConcData2Ov, 2);
%calculating the standard deviation
    controlSDOv = std(controlDataOv, 0, 2);
    controlSD2Ov = std(controlData2Ov, 0, 2);
    lowConcSDOv = std(lowConcDataOv, 0, 2);
    lowConcSD2Ov = std(lowConcData2Ov, 0, 2);
    lowConcSD3Ov = std(lowConcData3Ov, 0, 2);
    highConcSDOv = std(highConcDataOv, 0, 2);
    highConcSD2Ov = std(highConcData2Ov, 0, 2);
%creating labels without the irrelevant data points
    [labelLowOv, labelHighOv, labelBothOv]= deal(dataOv(2:end,1));
    labelLowOv(toDelete7, :) = [];
    labelHighOv(toDelete8, :) = [];
    labelBothOv(toDelete9, :) = [];
    clear toDelete7 toDelete8 toDelete9
%new lengths
    lengthLowOv = length(lowConcDataOv);
    lengthHighOv = length(highConcDataOv);
    lengthBothOv = length(highConcData2Ov);
%calculating the p-value
    p_ValueOv1 = mattest(controlDataOv, lowConcDataOv,'VarType', 'equal',  'Labels', labelLowOv);
    controlHighp_ValueOv  = mattest(controlData2Ov, highConcDataOv, 'VarType', 'equal',  'Labels', labelHighOv);
    controlLowp_ValueOv  = mattest(controlData2Ov, lowConcData3Ov, 'VarType', 'equal',  'Labels', labelHighOv);
    lowHighp_ValueOv  = mattest(lowConcData3Ov, highConcDataOv, 'VarType', 'equal',  'Labels', labelHighOv);
    p_ValueOv2  = mattest(lowConcData2Ov, highConcData2Ov, 'VarType', 'equal',  'Labels', labelHighOv);
%tables
    controlLowTableOv = table(labelLowOv, controlMeanOv, controlSDOv, lowConcMeanOv, lowConcSDOv, p_ValueOv1, lowConcMeanOv./controlMeanOv,controlDataOv, lowConcDataOv);
    controlHighTableOv = table(labelHighOv, controlMean2Ov, lowConcMean3Ov, highConcMeanOv, controlSD2Ov, lowConcSD3Ov, highConcSDOv, controlLowp_ValueOv, lowHighp_ValueOv, controlHighp_ValueOv, lowConcMean3Ov./controlMean2Ov, highConcMeanOv./lowConcMean3Ov, highConcMeanOv./controlMean2Ov,controlData2Ov, lowConcData3Ov, highConcDataOv);
    controlBothTableOv = table(labelBothOv, lowConcMean2Ov, lowConcSD2Ov, highConcMean2Ov, highConcSD2Ov, p_ValueOv2, highConcMean2Ov./lowConcSD2Ov, lowConcData2Ov, highConcData2Ov);
    [controlLowTableOv.Properties.VariableNames{9}, controlHighTableOv.Properties.VariableNames{12}, controlBotbTableOv.Properties.VariableNames{9}] = deal('LowHighFoldChangeOv');
    controlHighTableOv.Properties.VariableNames{11} = 'ControlLowFoldchangeOv';
    controlHighTableOv.Properties.VariableNames{13} = 'ControlHighFoldchangeOv';

% %volcano calculations and plots
%     header3 = {'Ov Gene Names'; 'Ov P-Value'; 'Ov Absolute Fold Change'};
%     volcanoLowOv = mavolcanoplot(lowConcDataOv, controlMeanOv, p_ValueOv1, 'Labels', labelLowOv, 'LogTrans', true);
%     title('Ovarian cancer Cells comparison of Control and Low Concentration 2P')
%     lowTableOv = table(volcanoLowOv.GeneLabels, volcanoLowOv.PValues, negToAbs(volcanoLowOv.FoldChanges), 'VariableNames', header3); 
%     volcanoHighOv = mavolcanoplot(highConcMeanOv, controlMean2Ov, controlHighp_ValueOv, 'Labels', labelHighOv, 'LogTrans', true);
%     title('Ovarian cancer Cells comparison of Control and High Concentration 2P')
%     highTableOv = table(volcanoHighOv.GeneLabels, volcanoHighOv.PValues, negToAbs(volcanoHighOv.FoldChanges), 'VariableNames', header3); 
%     volcanoBothOv = mavolcanoplot(highConcMean2Ov, lowConcMean2Ov, lowHighp_ValueOv2, 'Labels', labelBothOv, 'LogTrans', true);
%     title('Ovarian cancer Cells comparison of low and High Concentration 2P')
%     bothTableOv = table(volcanoBothOv.GeneLabels, volcanoBothOv.PValues, negToAbs(volcanoBothOv.FoldChanges), 'VariableNames', header3); 
clear excelDataOv dataOv p_ValueOv3 p_ValueOv2 labelLowhOv
clear controlDataOv controlData2Ov lowConcDataOv lowConcData2Ov lowConcData3Ov highConcDataOv highConcData2Ov 
clear toDelete1 toDelete2 toDelete3 text labelLowOv labelHighOv labelBothOv
clear p_ValueOv1 p_Value2Ov controlHighp_ValueOv controlLowp_ValueOv lowHighp_ValueOv
clear controlMeanOv controlMean2Ov lowConcMeanOv lowConcMean2Ov lowConcMean3Ov highConcMeanOv highConcMean2Ov
clear controlSDOv controlSD2Ov lowConcSDOv lowConcSD2Ov lowConcSD3Ov highConcSDOv highConcSD2Ov
clear lengthBothOv lengthHighOv lengthLowOv lengthOv
%% Fibroblast PRIMARY DERMAL FIBROBLAST
%Getting the data files 
    excelDataFib = '/Users/harshitakumar/Documents/Research/Letterio Research/Letterio Data/RNA Seq/4_genes_fpkm_expression F+KO.xlsx';
    [~, text, dataFib] = xlsread(excelDataFib);
%grabbing the specific data for each concentration group 
    [controlDataFib, controlData2Fib] = deal(cell2mat(dataFib(2:end, 14:16)));
    [lowConcDataFib, lowConcData2Fib, lowConcData3Fib] = deal(cell2mat(dataFib(2:end, 20:22)));
    [highConcDataFib, highConcData2Fib]= deal(cell2mat(dataFib(2:end, 17:19)));
%calculating the mean
    controlMeanFib = mean(controlDataFib, 2);
    lowConcMeanFib = mean(lowConcDataFib, 2);
    highConcMeanFib = mean(highConcDataFib, 2);
%deleting rows in which the mean values for each group is less than 1
    %group1 = control + low; group 2 = control + high
    toDelete1 = (controlMeanFib < 1) & (lowConcMeanFib < 1);
    controlDataFib(toDelete1,:) = [];
    lowConcDataFib(toDelete1,:) = [];
    toDelete2 = (controlMeanFib < 1) & (highConcMeanFib < 1);
    controlData2Fib(toDelete2,:) = [];
    highConcDataFib(toDelete2,:) = [];
    lowConcData3Fib(toDelete2,:) = [];
    toDelete3 = (lowConcMeanFib < 1) & (highConcMeanFib < 1);
    lowConcData2Fib(toDelete3,:) = [];
    highConcData2Fib(toDelete3,:) = [];
%creating labels without the irrelevant data points
    [labelLowFib, labelHighFib, labelBothFib] = deal(text(2:end,6));
    labelLowFib(toDelete1, :) = [];
    labelHighFib(toDelete2, :) = [];
    labelBothFib(toDelete3, :) = [];
%mean 
    controlMeanFib = mean(controlDataFib, 2);
    controlMean2Fib = mean(controlData2Fib, 2);
    lowConcMeanFib = mean(lowConcDataFib, 2);
    lowConcMean2Fib = mean(lowConcData2Fib, 2);
    lowConcMean3Fib = mean(lowConcData3Fib, 2);
    highConcMeanFib = mean(highConcDataFib, 2);
    highConcMean2Fib = mean(highConcData2Fib, 2);
%standard deviations
    controlSDFib = std(controlDataFib, 0, 2);
    controlSD2Fib = std(controlData2Fib, 0, 2);
    lowConcSDFib = std(lowConcDataFib, 0, 2);
    lowConcSD2Fib = std(lowConcData2Fib, 0, 2);
    lowConcSD3Fib = std(lowConcData3Fib, 0, 2);
    highConcSDFib = std(highConcDataFib, 0, 2);    
    highConcSD2Fib = std(highConcData2Fib, 0, 2);
%calculating the P-values
    p_Value1Fib = mattest(controlDataFib, lowConcDataFib, 'VarType', 'equal',  'Labels', labelLowFib);
    controlHighp_ValueFib  = mattest(controlData2Fib, highConcDataFib, 'VarType', 'equal', 'Labels', labelHighFib);
    controlLowp_ValueFib =  mattest(controlData2Fib, lowConcData3Fib, 'VarType', 'equal', 'Labels', labelHighFib);
    lowHighp_ValueFib = mattest(lowConcMean3Fib, highConcDataFib, 'VarType', 'equal', 'Labels', labelHighFib);
    p_Value2Fib  = mattest(lowConcData2Fib, highConcData2Fib, 'VarType', 'equal', 'Labels', labelHighFib);
%tables 
    %controlLowTableMM = table(labelLowMM, controlMeanMM, controlSDMM, lowConcMeanMM, lowConcSDMM, p_Value1MM, lowConcMeanMM./controlMeanMM, controlDataMM, lowConcDataMM);
    controlHighTableFib = table(labelHighFib, controlMean2Fib, lowConcMean3Fib, highConcMeanFib, controlSD2Fib, lowConcSD3Fib, highConcSDFib, controlLowp_ValueFib, lowHighp_ValueFib, controlHighp_ValueFib, lowConcMean3Fib./controlMean2Fib, highConcMeanFib./lowConcMean3Fib,  highConcMeanFib./controlMean2Fib,controlData2Fib, lowConcData3Fib, highConcDataFib);
    %controlBothTableMM = table(labelBothMM, lowConcMean2MM, lowConcSD2MM, highConcMean2MM, highConcSD2MM, p_Value2MM, highConcMean2MM./lowConcSD2MM, lowConcData2MM, highConcData2MM);
    [controlLowTableFib.Properties.VariableNames{9}, controlBothTableFib.Properties.VariableNames{9}] = deal('LowHighFoldChangeFib');
%     fibGSR = cell2table({'GSR', 0, 0, 0, 0, 0, 0, 0, 0,0, 0, 0, 0,0, 0, 0});
%     controlHighTableFib = [controlHighTableFib; fibGSR];
    controlHighTableFib.Properties.VariableNames{11} = 'ControlLowFoldchangeFib';
    controlHighTableFib.Properties.VariableNames{12} = 'LowHighFoldChangeFib';
    controlHighTableFib.Properties.VariableNames{13} = 'ControlHighFoldchangeFib';
    [~, p, ~] = intersect(controlHighTableOv.labelHighOv, 'UNKL');    
    controlHighTableFib(end+1, :) = controlHighTableOv(p, :);
    [~, p2, ~] = intersect(controlHighTableFib.labelHighFib, 'UNKL');    
    controlHighTableFib(p2, 2:end) = table(0);
% % %volcano calculations and plots
%     header = {'Fib Gene Names'; 'Fib P-Value'; 'Fib Absolute Fold Change'};
%     %control vs low
%         volcanoLowFib = mavolcanoplot(lowConcMeanFib , controlMeanFib, p_Value1Fib, 'Labels', labelLowFib, 'LogTrans', true, 'Plotonly', true);
%         title('Multiple Myeloma cells comparison of Control and 100nM (low Concentration) 2P')
%         lowTableFib = table(volcanoLowFib.GeneLabels, volcanoLowFib.PValues, negToAbs(volcanoLowFib.FoldChanges), 'VariableNames', header); 
%     % control vs high
%         volcanoHighFib = mavolcanoplot(highConcMeanFib, controlMean2Fib, controlHighp_ValueFib, 'Labels', labelHighFib, 'LogTrans', true);
%         title('Multiple Myeloma cells comparison of Control and 400nM (high Concentration) 2P')
%         highTableFib = table(volcanoHighFib.GeneLabels, volcanoHighFib.PValues, negToAbs(volcanoHighFib.FoldChanges), 'VariableNames', header);
%         [~, posm, ~] = intersect(labelHighFib, volcanoHighFib.GeneLabels);
%         controlHighDiffFib = controlHighTableFib(posm, :);
%     %low vs high  
%         volcanoBothMM = mavolcanoplot(highConcMean2MM, lowConcMean2MM, p_Value2MM, 'Labels', labelBothMM, 'LogTrans', true);
%         title('Multiple Myeloma cells comparison of 100nM (low Concentration) and 400nM (high Concentration) 2P')
%         bothTableMM = table(volcanoBothMM.GeneLabels, volcanoBothMM.PValues, negToAbs(volcanoBothMM.FoldChanges), 'VariableNames', header);
clear excelDataFib dataFib 
clear controlDataFib controlData2Fib lowConcDataFib lowConcData2Fib lowConcData3Fib highConcDataFib highConcData2Fib 
clear toDelete1 toDelete2 toDelete3 text labelLowFib labelHighFib labelBothFib
clear p_Value1Fib p_Value2Fib controlHighp_ValueFib controlLowp_ValueFib lowHighp_ValueFib
clear controlMeanFib controlMean2Fib lowConcMeanFib lowConcMean2Fib lowConcMean3Fib highConcMeanFib highConcMean2Fib
clear controlSDFib controlSD2Fib lowConcSDFib lowConcSD2Fib lowConcSD3Fib highConcSDFib highConcSD2Fib
%% KEAP Knockout Cells KEAP1 KO ARH-77
%Getting the data files 
    excelDataK = '/Users/harshitakumar/Documents/Research/Letterio Research/Letterio Data/RNA Seq/4_genes_fpkm_expression F+KO.xlsx';
    [~, text, dataK] = xlsread(excelDataK);
%grabbing the specific data for each concentration group 
    knockoutDataK = cell2mat(dataK(2:end, 23:25));
    wildTypeDataPCK = cell2mat(dataK(2:end, 32:34));
    wildTypeDataK = cell2mat(dataK(2:end, 32:33));
%calculating the mean
    wildTypeMeanK = mean(wildTypeDataK, 2);
    wildTypeMeanPCK = mean(wildTypeDataPCK, 2);
    knockoutMeanK = mean(knockoutDataK, 2);
%deleting rows in which the mean values for each group is less than 1
    %group1 = control + low; group 2 = control + high
    toDelete1 = (wildTypeMeanK < 1) & (knockoutMeanK < 1);
    wildTypeDataK(toDelete1,:) = [];
    knockoutDataK(toDelete1,:) = [];
    wildTypeDataPCK(toDelete1,:) = [];
%creating labels without the irrelevant data points
    labelK = text(2:end,6);
    labelK(toDelete1, :) = [];
%mean 
    wildTypeMeanK = mean(wildTypeDataK, 2);
    wildTypeMeanPCK = mean(wildTypeDataPCK, 2);
    knockoutMeanK = mean(knockoutDataK, 2);
%standard deviations
    wildTypeSDK = std(wildTypeDataK, 0, 2);
    wildTypeSDPCK = std(wildTypeDataPCK, 0, 2);
    knockoutSDK = std(knockoutDataK, 0, 2);
%calculating the P-values
    p_Value1K = mattest(wildTypeDataK, knockoutDataK, 'VarType', 'equal',  'Labels', labelK);
    p_Value1PCK = mattest(wildTypeDataPCK, knockoutDataK, 'VarType', 'equal',  'Labels', labelK);
%tables 
    controlTableK = table(labelK, wildTypeMeanK, knockoutMeanK, wildTypeSDK, knockoutSDK, p_Value1K, knockoutMeanK./wildTypeMeanK, wildTypeDataK, knockoutDataK);
    controlTablePCK = table(labelK, wildTypeMeanPCK, knockoutMeanK, wildTypeSDPCK, knockoutSDK, p_Value1PCK, knockoutMeanK./wildTypeMeanPCK, wildTypeDataPCK, knockoutDataK);    
    controlTableK.Properties.VariableNames{7} = 'knockoutWildtypeFoldChange';
% % %volcano calculations and plots
%     header = {'K Gene Names'; 'K P-Value'; 'K Absolute Fold Change'};
%     %control vs low
%         volcanoLowK = mavolcanoplot(lowConcMeanK , controlMeanK, p_Value1K, 'Labels', labelLowK, 'LogTrans', true, 'Plotonly', true);
%         title('Multiple Myeloma cells comparison of Control and 100nM (low Concentration) 2P')
%         lowTableFib = table(volcanoLowK.GeneLabels, volcanoLowK.PValues, negToAbs(volcanoLowK.FoldChanges), 'VariableNames', header); 
%     %control vs high
%         volcanoHighK = mavolcanoplot(highConcMeanK, controlMean2K, controlHighp_ValueK, 'Labels', labelHighK, 'LogTrans', true);
%         title('Multiple Myeloma cells comparison of Control and 400nM (high Concentration) 2P')
%         highTableK = table(volcanoHighK.GeneLabels, volcanoHighK.PValues, negToAbs(volcanoHighK.FoldChanges), 'VariableNames', header);
%         [~, posm, ~] = intersect(labelHighK, volcanoHighK.GeneLabels);
%         controlHighDiffK = controlHighTableK(posm, :);
%     %low vs high  
%         volcanoBothMM = mavolcanoplot(highConcMean2MM, lowConcMean2MM, p_Value2MM, 'Labels', labelBothMM, 'LogTrans', true);
%         title('Multiple Myeloma cells comparison of 100nM (low Concentration) and 400nM (high Concentration) 2P')
%         bothTableMM = table(volcanoBothMM.GeneLabels, volcanoBothMM.PValues, negToAbs(volcanoBothMM.FoldChanges), 'VariableNames', header);
clear excelDataK dataK toDelete1 text p p2
clear labelK wildTypeMeanK knockoutMeanK wildTypeSDK knockoutSDK p_Value1K wildTypeDataK knockoutDataK
%% TCGA + C CLE Data LUAD TCGA, NSCLC
%loading data files
excelCCLEData = '/Users/harshitakumar/Documents/Research/Letterio Research/Letterio Data/RNA Seq/KEAP1_filered_lung_logRNA2.xlsx';
[~, ~, dataFile1] = xlsread(excelCCLEData);
dataCCLE = cell2mat(dataFile1(2:end, [1 3]));
labelCCLE = dataFile1(2:end, 4);
excelTCGAData = '/Users/harshitakumar/Documents/Research/Letterio Research/Letterio Data/RNA Seq/TCGA_KEAP1_NFE2L2.xlsx';
[~, ~, dataFile] = xlsread(excelTCGAData);
dataTCGA = cell2mat(dataFile(2:end, [7 8]));
labelTCGA = dataFile(2:end, 1);
hmTCGAData = [];
i = 1;
%grapping the TCGA data that is upregulated
while (length(hmTCGAData) < 500)
    if (dataTCGA(i, 1) > 0)
        hmTCGAData = [hmTCGAData; labelTCGA(i, :), dataTCGA(i, 1), dataTCGA(i, 2) ];
    end
    i = i + 1;
end
hmCLData = [];
i = 1;
%grabbing the Cell line data that is upregulated
while (length(hmCLData) < 500)
    if (dataCCLE(i, 1) > 0)
        hmCLData = [hmCLData; labelCCLE(i, :), dataCCLE(i, 1), dataCCLE(i, 2) ];
    end
    i = i + 1;
end
%grabbing TCGA data that is downregulated
i = 1;
downTCGA = [];
while (length(downTCGA) < 200)
    if (dataTCGA(i, 1) < 0)
        downTCGA = [downTCGA; labelTCGA(i, :), dataTCGA(i, 1), dataTCGA(i, 2) ];
    end
    i = i + 1;
end
%grabbing CCLE data that is downregulated
i = 1;
downCCLE = [];
while (length(downCCLE) < 200)
    if (dataCCLE(i, 1) < 0)
        downCCLE = [downCCLE; labelCCLE(i, :), dataCCLE(i, 1), dataCCLE(i, 2) ];
    end
    i = i + 1;
end
%creating the tables
hmTCGAData = table(hmTCGAData(:,1), 2.^cell2mat(hmTCGAData(:,2)), cell2mat(hmTCGAData(:,3)));
hmCLData = table(hmCLData(:,1), 2.^cell2mat(hmCLData(:,2)), 10.^cell2mat(hmCLData(:,3)));
downTCGA = table(downTCGA(:,1), 2.^cell2mat(downTCGA(:,2)), cell2mat(downTCGA(:,3)));
downCCLE = table(downCCLE(:,1), 2.^cell2mat(downCCLE(:,2)), 10.^cell2mat(downCCLE(:,3)));
dataTCGA = table(labelTCGA, 2.^dataTCGA(:,1), dataTCGA(:,2));
dataCCLE = table(labelCCLE, 2.^dataCCLE(:,1), 10.^dataCCLE(:,2));
%Volcano plots
    volcanoTCGA = mavolcanoplot(hmTCGAData.(2), ones(height(hmTCGAData), 1), hmTCGAData.(3), 'Labels', hmTCGAData.(1), 'LogTrans', true, 'Plotonly', true);
    title('LUAD TCGA Differentially Expressed Genes', 'FontSize', 14)
    volcanoCCLE = mavolcanoplot(hmCLData.(2), ones(height(hmCLData), 1), hmCLData.(3), 'Labels', hmCLData.(1), 'LogTrans', true, 'Plotonly', true);
    title('NSCLC (CCLE) Differentially Expressed Genes', 'FontSize', 14)
clear i excelCCLEData excelTCGAData dataFile dataFile1 ans labelTCGA labelCCLE
%% Filtering for nRF2 Pathway Genes
    pvalue = 0.1; %selecting the p-value for table
    %selecting the genes for each dataset
    nRF2GenesMM = [];
    for i = 1:height(controlHighTableMM)
        %filtering all p-values 
        if (controlHighTableMM.controlLowp_ValueMM(i) < pvalue)
            %genes need to increase from control to low then decrease, concentrations at high to be greater than control
            if ((controlHighTableMM.ControlLowFoldchangeMM(i)> 1.25) )
                %& (controlHighTableMM.LowHighFoldChangeMM(i) < 5) && (controlHighTableMM.ControlHighFoldchangeMM(i) > 0.75)
                nRF2GenesMM = [nRF2GenesMM; controlHighTableMM(i, :)];
            end 
        end
    end
    nRF2GenesDIPG = [];
    for i = 1:height(controlHighTableDIPG)
        %filtering all p-values 
        if (controlHighTableDIPG.controlLowp_ValueDIPG(i) < pvalue)
            %genes need to double from control to low then decrease
            if ((controlHighTableDIPG.ControlLowFoldchangeDIPG(i)> 1.25))
                % && (controlHighTableDIPG.LowHighFoldChangeDIPG(i) < 5) && (controlHighTableDIPG.ControlHighFoldchangeDIPG(i) > 0.75)
                nRF2GenesDIPG = [nRF2GenesDIPG; controlHighTableDIPG(i, :)];
            end 
        end
    end
    nRF2GenesOv = [];
    for i = 1:height(controlHighTableOv)
        %filtering all p-values
        if (controlHighTableOv.controlLowp_ValueOv(i) < pvalue) 
            %genes need to double from control to low then decrease
            if ((controlHighTableOv.ControlLowFoldchangeOv(i)> 1.25) )
                %&& (controlHighTableOv.LowHighFoldChangeOv(i) < 5)&& (controlHighTableOv.ControlHighFoldchangeOv(i) > 0.75)
                nRF2GenesOv = [nRF2GenesOv; controlHighTableOv(i, :)];
            end 
        end
    end
    nRF2GenesFib = [];
    for i = 1:height(controlHighTableFib)
        %filtering all p-values 
        if (controlHighTableFib.controlLowp_ValueFib(i) < pvalue) 
            %genes need to double from control to low then decrease
            if ((controlHighTableFib.ControlLowFoldchangeFib(i)> 1.25) )
                %&& (controlHighTableFib.LowHighFoldChangeFib(i) < 5) && (controlHighTableFib.ControlHighFoldchangeFib(i) > 0.75)
                nRF2GenesFib = [nRF2GenesFib; controlHighTableFib(i, :)];
            end 
        end
    end
    nRF2GenesK = [];
    for i = 1:height(controlTableK)
        %filtering all p-values and filtering and fold change more than 1.25
        if ((controlTableK.p_Value1K(i) < pvalue) && (controlTableK.knockoutWildtypeFoldChange(i)> 1.25))
            nRF2GenesK = [nRF2GenesK; controlTableK(i, :)];
        end
    end
    %nrf2 genes that show up in 4 cell lines out of 5
        nRF2MDOF = intersect(intersect(intersect(nRF2GenesOv.labelHighOv, nRF2GenesMM.labelHighMM), nRF2GenesDIPG.labelHighDIPG),nRF2GenesFib.labelHighFib);
        nRF2MDOK = intersect(intersect(intersect(nRF2GenesOv.labelHighOv, nRF2GenesMM.labelHighMM), nRF2GenesDIPG.labelHighDIPG),nRF2GenesK.labelK);
        nRF2MDKF = intersect(intersect(intersect(nRF2GenesFib.labelHighFib, nRF2GenesMM.labelHighMM), nRF2GenesDIPG.labelHighDIPG),nRF2GenesK.labelK);
        nRF2MKOF = intersect(intersect(intersect(nRF2GenesOv.labelHighOv, nRF2GenesMM.labelHighMM), nRF2GenesFib.labelHighFib),nRF2GenesK.labelK);
        nRF2KDOF = intersect(intersect(intersect(nRF2GenesOv.labelHighOv, nRF2GenesFib.labelHighFib), nRF2GenesDIPG.labelHighDIPG),nRF2GenesK.labelK);
        nRF24of5 = union(union(union(union(nRF2MDOF, nRF2MDOK), nRF2MDKF), nRF2MKOF), nRF2KDOF);
clear i nRF2MDOF nRF2MDOK nRF2MDKF nRF2MKOF nRF2KDOF pvalue
%% Filtering and expression table for downregulated genes
    downUnion = union(downCCLE.(1), downTCGA.(1));
    downIntersect = intersect(downCCLE.(1), downTCGA.(1));
    chosenList = downUnion;
    data = [];
    dataDownTable = [];
    data2 = [];
    dataDownTable2 = [];
    downList = [];
    for i = 1:length(chosenList)
        gene = chosenList(i);
        data = [];
        data2 = [];
        for j = 1:7
            switch j
                case 1
                    cellLine = controlHighTableMM;
                    f = 11;
                    p = 8;
                case 2
                    cellLine = controlHighTableDIPG; 
                    f = 11;
                    p = 8;
                case 3
                    cellLine = controlHighTableOv;
                    f = 11;
                    p = 8;
                case 4
                    cellLine = controlHighTableFib;
                    f = 11;
                    p = 8;
                case 5
                    cellLine = controlTableK;
                    f = 7;
                    p = 6;
                case 6
                    cellLine = downTCGA;
                    f = 2;
                    p = 3;
                case 7
                    cellLine = downCCLE;
                    f = 2;
                    p = 3;
            end
            location = find(table2array(cellLine(:,1)) == string(gene));
            if (j <= 5)
                if (isempty(location) == 0)
                    if ((cellLine{location(1), p} <0.1) && (cellLine{location(1), f} <1))
                        data = [data, [cellLine{location(1), f}, cellLine{location(1), p}]];
                        data2 = [data2, [cellLine{location(1), f}, cellLine{location(1), p}]];
                        %data(i, j+1) = cellLine{i, f};
                        %data(i, j+6) =  cellLine{i, p};
                        downList = [downList; gene];
                    else
                        data = [data, [0, 0]];
                        data2 = [data2, [cellLine{location(1), f}, cellLine{location(1), p}]];
                    end
                else
                    data = [data, [NaN, NaN]];
                    data2 = [data2, [NaN, NaN]];
                end
            elseif (j == 6)
                location2 = find(table2array(dataTCGA(:,1)) == string(gene));
                if (isempty(location) == 0)
                    data = [data, [cellLine{location(1), f}, cellLine{location(1), p}]];
                elseif (isempty(location2) == 0)
                    data = [data, [0, 0]];
                else
                    data = [data, [NaN, NaN]];
                end
            elseif (j == 7)
                location2 = find(table2array(dataCCLE(:,1)) == string(gene));
                if (isempty(location) == 0)
                    data = [data, [cellLine{location(1), f}, cellLine{location(1), p}]];
                elseif (isempty(location2) == 0)
                    data = [data, [0, 0]];
                else
                    data = [data, [NaN, NaN]];
                end
            end
        end
        dataDownTable = [dataDownTable; [gene, data(1), data(2), data(3), data(4), data(5), data(6), data(7), data(8), data(9), data(10), data(11), data(12), data(13), data(14)]];
        dataDownTable2 = [dataDownTable2; [gene, data2(1), data2(2), data2(3), data2(4), data2(5), data2(6), data2(7), data2(8), data2(9), data2(10)]];
    end
    downList = table(downList);
    downList.Properties.VariableNames{1} = 'Labels';
    downCount = groupcounts(downList, "Labels");
    expressedIn = zeros(length(dataDownTable), 1);
    downRegIn = zeros(length(dataDownTable), 1);
    avFoldChange = zeros(length(dataDownTable), 1);
    geoMeanPValue = ones(length(dataDownTable), 1);
    for i = 1:length(dataDownTable)
        for j = 2:2:14
            if (cell2mat(dataDownTable(i, j)) == 0)
                expressedIn(i) = expressedIn(i) + 1;
            elseif ((cell2mat(dataDownTable(i, j)) > 0) || (cell2mat(dataDownTable(i, j)) < 0))
                expressedIn(i) = expressedIn(i) + 1;
                downRegIn(i) = downRegIn(i) + 1;
                avFoldChange(i) = avFoldChange(i) + cell2mat(dataDownTable(i, j));
                geoMeanPValue(i) = geoMeanPValue(i) * cell2mat(dataDownTable(i, j+1));
            end
        end
    end
    percentExpressed = downRegIn./expressedIn * 100;
    avFoldChange = avFoldChange./downRegIn;
    geoMeanPValue = geoMeanPValue.^(1./downRegIn);
    percentDownTable = [table(dataDownTable(:,1)), table(expressedIn), table(downRegIn), table(percentExpressed), table(avFoldChange), table(geoMeanPValue)];
    percentDownTable.Properties.VariableNames{1} = 'Labels';
    clear expressedIn downRegIn avFoldChange geoMeanPValue i data percentExpressed location cellLine f p j gene chosenList downList location2 data2 downUnion downIntersect
%% Expression table for upregulated
% putting all the genes with all the cancers in one table
    upUnion = union(union(union(union(union(union(nRF2GenesMM.(1), nRF2GenesDIPG.(1)), nRF2GenesK.(1)), nRF2GenesFib.(1)), nRF2GenesOv.(1)), hmCLData.(1)), hmTCGAData.(1));
    data = [];
    dataUpTable = [];
    data2 = [];
    dataUpTable2 = [];
    heat = [];
    for i = 1:length(upUnion)
        gene = upUnion(i);
        data = [];
        for j = 1:7
            switch j
                case 1
                    cellLine = controlHighTableMM;
                    nLine = nRF2GenesMM;
                    f = 11;
                    p = 8;
                case 2
                    cellLine = controlHighTableDIPG; 
                    nLine = nRF2GenesDIPG;
                    f = 11;
                    p = 8;
                case 3
                    cellLine = controlHighTableOv;
                    nLine = nRF2GenesOv;
                    f = 11;
                    p = 8;
                case 4
                    cellLine = controlHighTableFib;
                    nLine = nRF2GenesFib;
                    f = 11;
                    p = 8;
                case 5
                    cellLine = controlTableK;
                    nLine = nRF2GenesK;
                    f = 7;
                    p = 6;
                case 6
                    nLine = hmTCGAData;
                    cellLine = dataTCGA;
                    f = 2;
                    p = 3;
                case 7
                    nLine = hmCLData;
                    cellLine = dataCCLE;
                    f = 2;
                    p = 3;
            end
            locationCL = find(table2array(nLine(:,1)) == string(gene));
            locationNL = find(table2array(cellLine(:,1)) == string(gene));
            if (isempty(locationCL) == 0) %if in NRFlist
                data = [data, [nLine{locationCL(1), f}, nLine{locationCL(1), p}]];
                data2 = [data2, [nLine{locationCL(1), f}, nLine{locationCL(1), p}]];
            elseif (isempty(locationNL) == 0) %if expressed
                data = [data, [0, 0]];

            else %if not expressed
                data = [data, [NaN, NaN]];
            end
        end
        dataUpTable = [dataUpTable; [gene, data(1), data(2), data(3), data(4), data(5), data(6), data(7), data(8), data(9), data(10), data(11), data(12), data(13), data(14)]];
    end
    expressedIn = zeros(length(dataUpTable), 1);
    upRegIn = zeros(length(dataUpTable), 1);
    avFoldChange = zeros(length(dataUpTable), 1);
    geoMeanPValue = ones(length(dataUpTable), 1);
    log2FoldChange = zeros(length(dataUpTable), 1);
    rangeFold = zeros(length(dataUpTable), 1);
    minFold = zeros(length(dataUpTable), 1);
    maxFold = zeros(length(dataUpTable), 1);
    for i = 1:length(dataUpTable)
        r = [];
        for j = 2:2:14
            if (cell2mat(dataUpTable(i, j)) == 0)
                expressedIn(i) = expressedIn(i) + 1;
            elseif (cell2mat(dataUpTable(i, j)) > 0)
                expressedIn(i) = expressedIn(i) + 1;
                upRegIn(i) = upRegIn(i) + 1;
                avFoldChange(i) = avFoldChange(i) + cell2mat(dataUpTable(i, j));
                log2FoldChange(i) = log2FoldChange(i) + log2(cell2mat(dataUpTable(i, j)));
                geoMeanPValue(i) = geoMeanPValue(i) * cell2mat(dataUpTable(i, j+1));
                r = [r, cell2mat(dataUpTable(i, j))];
            end
        end
        rangeFold(i) = range(r);
        minFold(i) = min(r);
        maxFold(i) = max(r);
    end
    percentExpressed = upRegIn./expressedIn * 100;
    avFoldChange = avFoldChange./upRegIn;
    log2FoldChange = log2FoldChange./upRegIn;
    geoMeanPValue = geoMeanPValue.^(1./upRegIn);
    percentUpTable = [table(dataUpTable(:,1)), table(expressedIn), table(upRegIn), table(percentExpressed), table(avFoldChange), table(geoMeanPValue), table(log2FoldChange), table(rangeFold), table(minFold), table(maxFold)];
    percentUpTable.Properties.VariableNames{1} = 'Labels';
    [percentUpTable, idx] = sortrows(percentUpTable,"upRegIn","descend");
    dataUpTable = dataUpTable(idx, :);
    clear expressedIn upRegIn avFoldChange geoMeanPValue i data percentExpressed location cellLine f p j gene chosenList downList nLine
    clear rangeFold minFold r maxFold log2FoldChange heat upUnion locationCL locationNL idx
%% nRF2 processing using online data as well
    allNRF2list = table([nRF2GenesOv.labelHighOv; nRF2GenesMM.labelHighMM; nRF2GenesDIPG.labelHighDIPG; nRF2GenesFib.labelHighFib; nRF2GenesK.labelK; hmCLData.(1); hmTCGAData.(1)]);
    allNRF2list.Properties.VariableNames{1} = 'Labels';
    geneCounts = groupcounts(allNRF2list, "Labels");
    nRF27of7 = [];
    nRF26of7 = [];
    nRF25of7 = [];
    nRF24of7 = [];
    for i = 1:size(geneCounts, 1)
        if (table2array(geneCounts(i, 2)) >= 7)
            nRF27of7 = [nRF27of7; geneCounts(i, 1)];
        end
        if (table2array(geneCounts(i, 2)) == 6)
            nRF26of7 = [nRF26of7; geneCounts(i, 1)];
        end
        if (table2array(geneCounts(i, 2)) == 5)
            nRF25of7 = [nRF25of7; geneCounts(i, 1)];
        end
        if (table2array(geneCounts(i, 2)) >= 4)
            nRF24of7 = [nRF24of7; geneCounts(i, 1)];
        end
    end
    nRF27of7 = nRF27of7([1:end], :);
    %fold Change + pValue table
        fold = [];
        foldChange = [];
        pValue = [];
        pValueChange = [];
        for i = 1:length(nRF27of7.(1))
            gene = string(nRF27of7{i,1});
            fold = [];
            pValue = [];
            for j = 1:7
                switch j
                    case 1
                        cellLine = nRF2GenesMM;
                        f = 11;
                        p = 8;
                    case 2
                        cellLine = nRF2GenesDIPG; 
                        f = 11;
                        p = 8;
                    case 3
                        cellLine = nRF2GenesOv;
                        f = 11;
                        p = 8;
                    case 4
                        cellLine = nRF2GenesFib;
                        f = 11;
                        p = 8;
                    case 5
                        cellLine = nRF2GenesK;
                        f = 7;
                        p = 6;
                    case 6
                        cellLine = hmTCGAData;
                        f = 2;
                        p = 3;
                    case 7
                        cellLine = hmCLData;
                        f = 2;
                        p = 3;
                end
                location = find(cellLine.(1) == gene);
                if isempty(location)
                    fold = [fold, NaN];
                    pValue = [pValue, NaN];
                else                
                    location = location(1);
                    fold = [fold, cellLine{location,f}];
                    pValue = [pValue, cellLine{location,p}];
                end
            end
            foldChange = [foldChange; fold];
            pValueChange = [pValueChange; pValue];
        end
        avFoldChange = mean(foldChange(:,1:5), 2, 'omitnan');
        avPValue = geomean(pValueChange(:,1:5), 2, 'omitnan');
        nRF27of7 = [nRF27of7, table(foldChange), table(pValueChange), table(avFoldChange), table(avPValue)];
    nRF25of7 = nRF25of7([1:end], :);
    %fold Change + pValue table
        fold = [];
        foldChange = [];
        pValue = [];
        pValueChange = [];
        for i = 1:length(nRF25of7.(1))
            gene = string(nRF25of7{i,1});
            fold = [];
            pValue = [];
            for j = 1:7
                switch j
                    case 1
                        cellLine = nRF2GenesMM;
                        f = 11;
                        p = 8;
                    case 2
                        cellLine = nRF2GenesDIPG; 
                        f = 11;
                        p = 8;
                    case 3
                        cellLine = nRF2GenesOv;
                        f = 11;
                        p = 8;
                    case 4
                        cellLine = nRF2GenesFib;
                        f = 11;
                        p = 8;
                    case 5
                        cellLine = nRF2GenesK;
                        f = 7;
                        p = 6;
                    case 6
                        cellLine = hmTCGAData;
                        f = 2;
                        p = 3;
                    case 7
                        cellLine = hmCLData;
                        f = 2;
                        p = 3;
                end
                location = find(cellLine.(1) == gene);
                if isempty(location)
                    fold = [fold, NaN];
                    pValue = [pValue, NaN];
                else                
                    location = location(1);
                    fold = [fold, cellLine{location,f}];
                    pValue = [pValue, cellLine{location,p}];
                end
            end
            foldChange = [foldChange; fold];
            pValueChange = [pValueChange; pValue];
        end
        avFoldChange = mean(foldChange(:,1:5), 2, 'omitnan');
        avPValue = geomean(pValueChange(:,1:5), 2, 'omitnan');
        nRF25of7 = [nRF25of7, table(foldChange), table(pValueChange), table(avFoldChange), table(avPValue)];
        clear fold pValue p i f j location avFoldChange avPValue gene allNRF2list foldChange cellLine  pValueChange geneCounts 
%% Plotting Core genes Relative Fold Change
pValueCore = [];
for j = 2:2
    gene = nRF27of7(j, :);
    barData = [];
    errData = [];
    actualDataPoints = [];
    textt = [];
    labels = {};
    pcor = [];
    for i = 1:4
        switch i
            case 4
                cellLine = controlHighTableMM;
                titles = ('MM Cells');
            case 1
                cellLine = controlHighTableDIPG;  
                titles = ('DIPG Cells');
            case 2
                cellLine = controlHighTableOv;
                titles = ('Ovarian Cells');
            case 3
                cellLine = controlHighTableFib;
                titles = ('Fibroblast Cells');
        end
        [~, p, ~] = intersect(cellLine.(1), gene.(1));
        data = table2array(cellLine(p, 2:end));
        if (data(1,7) < 0.01)
                if (data(1, 9) < 0.01)
                    g = {'**', '**'};
                elseif (data(1, 9) < 0.05)
                    g = {'**', '*'};
                else
                    g = {'**', '^{NS}'};
                end
            elseif (data(1,7) < 0.05)
                if (data(1, 9) < 0.01)
                    g = {'*', '**'};    
                elseif (data(1, 9) < 0.05)
                    g = {'*', '*'};
                else
                    g = {'*', '^{NS}'};
                end
        else 
                if (data(1, 9) < 0.01)
                    g = {'^{0.056}', '**'};
                elseif (data(1, 9) < 0.05)
                    g = {'^{0.072}', '*'};
                else
                    g = {'^NS', '^{NS}'};
                end
            end
        barData = [barData, [data(1, 2)/data(1,1), data(1, 3)/data(1,1)]];
        actualDataPoints = [actualDataPoints; table2array(cellLine(p,15))/data(1,1);table2array(cellLine(p,16))/data(1,1)];
        errData = [errData, [data(1,5)/data(1,1), data(1,6)/data(1,1)]];
        textt = [textt, [data(1, 2)/data(1,1), data(1, 3)/data(1,1)] + [data(1,5)/data(1,1), data(1,6)/data(1,1)]];
        labels = [labels, g];
        pcor = [pcor, data(1,7), data(1,9)];
    end
    cellLine = controlTableK;
    [~, p, ~] = intersect(cellLine.(1), gene.(1));
    data = table2array(cellLine(p, 2:end));
    if (data(1, 5) < 0.01)
        g = {'**'};
    elseif (data(1, 5) < 0.05)
        g = {'*'};
    else
        g = {'^{NS}'};
    end
    barData = [barData, [data(1, 2)/data(1,1)]];
    errData = [errData, data(1, 4)/data(1,1)];
    textt = [textt, [data(1, 2)/data(1,1)] + [data(1,4)/data(1,1)]];
    labels = [labels, g];
    pValueCore = [pValueCore; [pcor,data(1,5)]];
    actualDataPoints = [actualDataPoints; table2array(cellLine(p,9))/data(1,1)];
     %FIGURE
    figure()
    a = 1:9;
    bar(1:2:a(end-2), barData(1:2:a(end-2)), 0.2,  'FaceColor', [0.6 0.6 0.6]);
    hold on
    bar(0, 1, 0.4, 'w')
    bar(a(end), barData(a(end)), 0.4,  'FaceColor', [0.4 0.4 0.4]);
    bar(2:2:a(end)-1, barData(2:2:a(end-1)), 0.2, 'FaceColor', [0.8 0.8 0.8]);
    bar(a(end)+0.2, 0)
    bar(-0.2, 0)
    swarmchart(1:a(end), actualDataPoints', 30 ,"black", 'LineWidth', 1)
    %er = errorbar(a, barData, zeros(1, length(errData)), errData);
    %er.Color = [0 0 0]; er.LineStyle = 'none';
    text(1:a(end)-1, 1.1*textt(1:a(end-1)), labels(1:a(end-1)), 'FontSize', 14, 'FontWeight', 'bold', 'HorizontalAlignment', 'center')
    text(a(end), 1.1*textt(a(end)), labels(a(end)), 'FontSize', 14, 'FontWeight', 'bold', 'HorizontalAlignment', 'center')
    xticks(0:a(end))
    xticklabels({'C', 'S_{L}', 'S_{H}', 'O_{L}', 'O_{H}', 'F_{L}', 'F_{H}','R_{L}','R_{H}', 'A_{K}'})    
    set(gcf,'color','w');
    set(gcf, 'Units', 'normalized', 'Position', [0 1 0.24 0.16]);
    title(string(gene.(1)), 'FontWeight', 'bold', 'FontSize', 16, 'FontAngle', 'italic')
    %title('~CBR3', 'FontWeight', 'bold', 'FontSize', 16, 'FontAngle', 'italic')
    ylabel('Fold Change', 'FontWeight', 'bold', 'FontSize', 12)
    xlabel('Treatment', 'FontWeight', 'bold', 'FontSize', 13)
    axis tight
    ylim([0 max(textt)*1.25])
    hold off
%     if (j == 2)
%         ylim([0 240]);
%         breakyaxis([30 220]);
%     end
end
    clear cellLine k bar7of7Genes er p titles tit  gene i  g barData errData textt p titles gene
    clear pcor a
%% Plotting Conditional genes Relative Fold Change
geneList = {'HMOX1', 'AKIRIN2', 'G6PD', 'TKT', 'TALDO1', 'PGD', 'B4GALNT1', 'POPDC3', 'PLAT', 'CDC42EP1', 'EVA1A', 'AKR1C1', 'AKR1C2', 'CDC42EP3'};
EMTgenes = {'NFE2L2',  'ZEB1', 'SNAI1', 'CDH1', 'NOTCH1', 'DLL4'}
pValueOther = [];
for j = 1:length(EMTgenes)
   gene = EMTgenes(j);
    barData = [];
    actualDataPoints = [];
    errData = [];
    textt = [];
    labels = {};
    num1 = [];
    num2 = [];
    pOr = [];
    for i = [1 3 5 7]
        switch i
            case 7
                cellLine = controlHighTableMM;
                titles = ('MM Cells');
            case 1
                cellLine = controlHighTableDIPG;  
                titles = ('DIPG Cells');
            case 3
                cellLine = controlHighTableOv;
                titles = ('Ovarian Cells');
            case 5
                cellLine = controlHighTableFib;
                titles = ('Fibroblast Cells');
        end
        [~, p, ~] = intersect(cellLine.(1), gene);
        data = table2array(cellLine(p, 2:end));
        if (isempty(data))
            barData = [barData, [0, 0]];
            actualDataPoints = [actualDataPoints;  [0, 0, 0];[0, 0, 0]];
            errData = [errData, [0, 0]];
            textt = [textt, [0, 0]];
            labels = [labels, {'', ''}];
            pOr = [pOr, 0, 0];
        else
            if (data(1,7) < 0.01)
                if (data(1, 9) < 0.01)
                    g = {'**', '**'};
                elseif (data(1, 9) < 0.05)
                    g = {'**', '*'};
                else
                    g = {'**', '^{NS}'};
                end
            elseif (data(1,7) < 0.05)
                if (data(1, 9) < 0.01)
                    g = {'*', '**'};
                elseif (data(1, 9) < 0.05)
                    g = {'*', '*'};
                else
                    g = {'*', '^{NS}'};
                end
            else 
                if (data(1, 9) < 0.01)
                    g = {'^{NS}', '**'};
                elseif (data(1, 9) < 0.05)
                    g = {'^{NS}', '*'};
                else
                    g = {'^{NS}', '^{NS}'};
                end
            end
            errData = [errData, [data(1,5)/data(1,1), data(1,6)/data(1,1)]];
            actualDataPoints = [actualDataPoints; table2array(cellLine(p,15))/data(1,1);table2array(cellLine(p,16))/data(1,1)];
            barData = [barData, [data(1, 2)/data(1,1), data(1, 3)/data(1,1)]];
            textt = [textt, [data(1, 2)/data(1,1), data(1, 3)/data(1,1)] + [data(1,5)/data(1,1), data(1,6)/data(1,1)]];
            labels = [labels, g];
            pOr = [pOr, data(1,7), data(1,9)];
        end
            num1 = [num1, i];
            num2 = [num2, i+1];
    end
    cellLine = controlTableK;
    [~, p, ~] = intersect(cellLine.(1), gene);
    data = table2array(cellLine(p, 2:end));
    if (isempty(data))
        barData = [barData, [0]];
        actualDataPoints = [actualDataPoints;  [0, 0, 0]];
        errData = [errData, 0];
        textt = [textt, [0]];
        labels = [labels, {''}];
      pValueOther = [pValueOther; [pOr,0]];
    else
        if (data(1, 5) < 0.01)
            g = {'**'};
        elseif (data(1, 5) < 0.05)
            g = {'*'};
        else
            g = {'^{NS}'};
        end
        barData = [barData, [data(1, 2)/data(1,1)]];
        actualDataPoints = [actualDataPoints; table2array(cellLine(p,9))/data(1,1)];
        errData = [errData, data(1, 4)/data(1,1)];
        textt = [textt, [data(1, 2)/data(1,1)] + [data(1,4)/data(1,1)]];
        labels = [labels, g];
        pValueOther = [pValueOther; [pOr,data(1,5)]];
    end
    num1 = [num1, [num1(end) + 2]];
    if (j == 12)
        barData(:, 7:8) = 0;
        errData(:, 7:8) = 0;
        textt(:, 7:8) = 0;
        labels(:, 7:8) = {''};
    end
    %FIGURE
    figure()
    a = 1:9;
    bar(1:2:a(end-2), barData(1:2:a(end-2)), 0.2,  'FaceColor', [0.6 0.6 0.6]);
    hold on
    bar(0, 1, 0.4, 'w')
    bar(a(end), barData(a(end)), 0.4,  'FaceColor', [0.4 0.4 0.4]);
    bar(2:2:a(end)-1, barData(2:2:a(end-1)), 0.2, 'FaceColor', [0.8 0.8 0.8]);
    bar(a(end)+0.2, 0)
    bar(-0.2, 0)
    swarmchart(1:a(end), actualDataPoints', 30 ,"black", 'LineWidth', 1)
    barData
    %er = errorbar(a, barData, zeros(1, length(errData)), errData);
    %er.Color = [0 0 0]; er.LineStyle = 'none';
    text(1:a(end)-1, 1.1*textt(1:a(end-1)), labels(1:a(end-1)), 'FontSize', 14, 'FontWeight', 'bold', 'HorizontalAlignment', 'center')
    text(a(end), 1.1*textt(a(end)), labels(a(end)), 'FontSize', 14, 'FontWeight', 'bold', 'HorizontalAlignment', 'center')
    xticks(0:a(end))
    xticklabels({'C', 'S_{L}', 'S_{H}', 'O_{L}', 'O_{H}', 'F_{L}', 'F_{H}','R_{L}','R_{H}', 'A_{K}'})    
    set(gcf,'color','w');
    set(gcf, 'Units', 'normalized', 'Position', [0 1 0.24 0.16]);
    title(string(EMTgenes(j)), 'FontWeight', 'bold', 'FontSize', 16, 'FontAngle', 'italic')
    ylabel('Fold Change', 'FontWeight', 'bold', 'FontSize', 12)
    xlabel('Treatment', 'FontWeight', 'bold', 'FontSize', 13)
    axis tight
    if (max(textt)*1.2 > 1)
        lim = max(textt)*1.2;
    else
        lim = 1.1;
    end
    ylim([0 lim])
    if (j == 13)
        set(gca, 'YScale', 'log')
        set(gca,'YTick',[0.1, 1, 10, 100, 1000])
        set(gca,'YTickLabel',[0.1 1 10 100 1000])
        ylim([10^-1 10^3.5])
    end
    if (j == 11)
        ylim([0 1.2])
    end
%     if (j == 1)
%         ylim([0 150])
%         breakyaxis([35 40])
%     end
   hold off
end
    clear cellLine k bar7of7Genes er p titles tit  gene i  g  errData textt p titles gene
    clear er ax1 ax2 g gene geneList data p p2 pOr y yt ax 
%% NQO1 Genes
    gene = ("NQO1");
    tiledlayout(1, 4);
    titles = ["SF8628","OVCAR8", "PDF","RPMI-8226"];
    for i = 1:4
        switch i
            case 4
                cellLine = controlHighTableMM;
            case 1
                cellLine = controlHighTableDIPG;  
            case 2
                cellLine = controlHighTableOv;
            case 3
                cellLine = controlHighTableFib;
        end
        [~, p, ~] = intersect(cellLine.(1), gene);
        data = table2array(cellLine(p, 2:end));
        if (data(1,7) < 0.01)
            if (data(1, 9) < 0.01)
                g = {'', '**', '**'};
            elseif (data(1, 9) < 0.05)
                g = {'','**', '*'};
            else
                g = {'','**', '^{NS}'};
            end
        else
            if (data(1, 9) < 0.01)
                g = {'','*', '**'};
            elseif (data(1, 9) < 0.05)
                g = {'','*', '*'};
            else
                g = {'','*', '^{NS}'};
            end
        end
        nexttile
        bar(1:3, data(1:3)/data(1), 0.6,  'FaceColor', [0.7 0.7 0.7]);
        hold on
        if (i <= 3)
            swarmchart(1:3, [table2array(cellLine(p,14))/data(1); table2array(cellLine(p,15))/data(1);table2array(cellLine(p,16))/data(1)], 25, "k" ,'LineWidth', 1.5)
        else
            swarmchart(1:3, [[-5,table2array(cellLine(p,14))/data(1)]; table2array(cellLine(p,15))/data(1);table2array(cellLine(p,16))/data(1)], 25, "k" ,'LineWidth', 1.5)
        end
        %er = errorbar(1:3, data(1:3)/data(1), zeros(1, 3), data(4:6)/data(1));
        %er.Color = [0 0 0]; er.LineStyle = 'none';
        text(1:3, 1.03*(data(1:3)/data(1) + data(4:6)/data(1)), g, 'FontSize', 16, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
        xticklabels({'Control', 'Low', 'High'});
        set(gcf,'color','w');
        title("{\it NQO1}" + {" in "} + titles(i), 'FontWeight', 'bold', 'FontSize', 14);
        ylabel('Relative expression', 'FontSize', 12)
        ylim([0 1.15*max(data(1:3)/data(1) + data(4:6)/data(1))])
    end 
    clear er g i  titles   nRF2 titles
%% Plotting NRF2 genes
    %bar plot for nRF2 likely genes mean over all genes in the list
    %for each cell line separately, and then combined as well
            barNRF2Likely = [mean(table2array(nRF2GenesMM(:, 2:4), 1));mean(table2array(nRF2GenesDIPG(:, 2:4), 1));mean(table2array(nRF2GenesOv(:, 2:4), 1));mean(table2array(nRF2GenesFib(:, 2:4), 1));mean(table2array(nRF2GenesK(:, 2:4), 1))];
            barNRF2Likely = [barNRF2Likely; mean(barNRF2Likely, 1)];
            label = {'Multiple Myeloma', 'DIPG', 'Ovarian', 'Fibroblast', 'KEAP1 Status', 'Combined'};
            clear p1 p2 p3 p4 p5 p6
            figure(1)
            barnRF2 = tiledlayout(2, 3);
            for i = 1:length(barNRF2Likely)
                nexttile
                bar(barNRF2Likely(i,:)')
                set(gca, 'xticklabels', {'Control', 'Low', 'High'})
                title(label(i))
            end
            title(barnRF2, 'NRF2 Genes in Different Cell lines')
            xlabel(barnRF2, '2P Concentrations')
            ylabel(barnRF2, 'FPKM')
clear barNRF2Likely likely i barnRF2 label ans
%% Heatmaps 
    %universalHeatmap 
        universalHeatmap = [];
        ylabelUniversal = [];
        for i = 1:length(dataUpTable)
            heat = [];
            if ((percentUpTable.(3)(i) > 5) && (percentUpTable.(4)(i) > 80))
                ylabelUniversal = [ylabelUniversal; percentUpTable.(1)(i)];
                for j = 2:2:14
                    if (dataUpTable{i,j} > 0)
                        heat = [heat, 1];
                    elseif (dataUpTable{i,j} == 0)
                        heat = [heat, 0];
                    else
                        heat = [heat, -1];
                    end
                end
                universalHeatmap = [universalHeatmap; heat];
            end
        end
        universalHeatmap = universalHeatmap(:, [2 3 4 1 5 6 7]);
        % 1 is upregulated (red), 0 is expressed (white), -1 is not
        % expressed (black)
        figure(2)
        %xlabels2 = char({'DIPG', 'Ovarian', 'Fibroblast', 'MM', 'KEAP1 Knockout', 'TCGA', 'CCLE'});
        xlabels2 = char({'SF8628', 'OVCAR8', 'Primary Dermal Fibroblast', 'RPMI-8226', 'KEAP1 KO ARH-77', 'LUAD TCGA', 'NSCLC'});
        h = heatmap(xlabels2, ylabelUniversal, double(universalHeatmap), 'Title', 'Universal and near universal Upregulated Genes', 'FontSize', 20, 'CellLabelColor', 'none', 'ColorbarVisible', 'on');
        %h.YLabel = (ylabelUniversal, 'FontAngle', 'italic');
        cmap = colormap([0, 0, 0; 1, 1, 1; 1, 0, 0]);
        axs = struct(gca); %ignore warning that this should be avoided
        c = axs.Colorbar;
        c.Ticks = [-1, 0, 1];
        c.TickLabelsMode = 'manual';
        c.TickLabels = {'Not Expressed','Unchanged','Upregulated'};
        set(gcf,'color','w');
        c.FontSize = 12;
    %near universal heat map
        nearUniversalHeatmap = [];
        ylabelNearUniversal = [];
        for i = 1:length(dataUpTable)
            heat = [];
            if ((percentUpTable.(3)(i) < 6) && (percentUpTable.(4)(i) >= 60) && (percentUpTable.(3)(i) > 4))
                if (sum(string(percentUpTable.(1)(i)) == [string({'FTH1P2'}), string({'FTH1P7'}), string({'FTH1P8'})]) >= 1)
                else                
                    ylabelNearUniversal = [ylabelNearUniversal; percentUpTable.(1)(i)];
                    for j = 2:2:14
                        if (dataUpTable{i,j} > 0)
                            heat = [heat, 1];
                        elseif (dataUpTable{i,j} == 0)
                            heat = [heat, 0];
                        else
                            heat = [heat, -1];
                        end
                    end
                    nearUniversalHeatmap = [nearUniversalHeatmap; heat];
                end
            end
        end
        nearUniversalHeatmap = nearUniversalHeatmap(:, [2 3 4 1 5 6 7]);
        % 1 is upregulated (red), 0 is expressed (white), -1 is not
        % expressed (black)
        figure(3)
        h = heatmap(xlabels2, ylabelNearUniversal, double(nearUniversalHeatmap), 'Title','Conditional Upregulated Genes', 'FontSize', 20,  'CellLabelColor', 'none', 'ColorbarVisible', 'on');
        cmap = colormap([0, 0, 0; 1, 1, 1; 1, 0, 0]);
        axs = struct(gca); %ignore warning that this should be avoided
        c = axs.Colorbar;
        c.Ticks = [-1, 0, 1];
        c.TickLabelsMode = 'manual';
        c.TickLabels = {'Not Expressed','Unchanged','Upregulated'};
        set(gcf,'color','w');
        c.FontSize = 12;
    %down regulated genes
        downRegulatedHeatmap = [];
        ylabelDown = [];
        [percentDownTable, idx] = sortrows(percentDownTable,"percentExpressed","descend");
        dataDownTable = dataDownTable(idx, :);
        for i = 1:length(dataDownTable)
            heat = [];
            if ((percentDownTable.(3)(i) >= 4) && (percentDownTable.(4)(i) > 70))
                ylabelDown = [ylabelDown; percentDownTable.(1)(i)];
                for j = 2:2:14
                    if ((dataDownTable{i,j} > 0) || (dataDownTable{i,j} < 0))
                        heat = [heat, 1];
                    elseif (dataDownTable{i,j} == 0)
                        heat = [heat, 0];
                    else
                        heat = [heat, -1];
                    end
                end
                downRegulatedHeatmap = [downRegulatedHeatmap; heat];
            end
        end
        downRegulatedHeatmap = downRegulatedHeatmap(:, [2 3 4 1 5 6 7]);
        figure(4)
        heatmap(xlabels2, ylabelDown, double(downRegulatedHeatmap), 'Title', 'Downregulated Genes','FontSize', 20,  'CellLabelColor', 'none', 'ColorbarVisible', 'on');
        cmap = colormap([0, 0, 0; 1, 1, 1; 0, 0, 1]);
        axs = struct(gca); %ignore warning that this should be avoided
        c = axs.Colorbar;
        c.Ticks = [-1, 0, 1];
        c.TickLabelsMode = 'manual';
        c.TickLabels = {'Not Expressed','Unchanged','Downregulated'};
        set(gcf,'color','w');
        c.FontSize = 12;
% Fold Change Table
    [~, ~, pos] = intersect(union(ylabelNearUniversal, ylabelUniversal), percentUpTable.(1));
    upRegulatedFoldTable = percentUpTable(pos, :);
    [~, ~, pos] = intersect(ylabelDown, percentDownTable.(1));
    downRegulatedFoldTable = percentDownTable(pos, :);
    clear i c xlabels xlabels2 mymap heat downRegulatedHeatmap pos ylabelNearUniversal
    clear axs cmap j dataDownTable2 downCount idx nearUniversalHeatmap universalHeatmap ylabelDown
%% Validation
%Head Neck Cancer + Liver Validation
    excelHeadNeckData = '/Users/harshitakumar/Documents/Research/Letterio Research/Letterio Data/RNA Seq/HN-Nrf2-List.xlsx';
    [~, ~, dataFile1] = xlsread(excelHeadNeckData);
    dataHeadNeck = cell2mat(dataFile1(2:end, [7 8]));
    labelHeadNeck = cell2table(dataFile1(2:end, 1));
    excelLiverData ='/Users/harshitakumar/Documents/Research/Letterio Research/Letterio Data/RNA Seq/liver-TCGA-Nrf2.xlsx';
    [~, ~, dataFile1] = xlsread(excelLiverData);
    dataLiver = cell2mat(dataFile1(2:end, [7 8]));
    labelLiver = cell2table(dataFile1(2:end, 1));
    excelLungData ='/Users/harshitakumar/Documents/Research/Letterio Research/Letterio Data/RNA Seq/LSCC mutations.xlsx';
    [~, ~, dataFile1] = xlsread(excelLungData);
    dataLung = cell2mat(dataFile1(2:end, [7 8]));
    labelLung = cell2table(dataFile1(2:end, 1));
    upRegHeadNeck = [];
    upRegL = [];
    headNeckVal = [];
    liverVal = [];
    upRegLung = [];
    lungVal = [];
    for i = 1:length(ylabelUniversal)
        gene = ylabelUniversal(i);
        location = find(table2array(labelHeadNeck(:,1)) == string(gene));
        if (isempty(location) == 0)
            if ((dataHeadNeck(location, 1) > 0) && (dataHeadNeck(location, 2) < 0.05))
                upRegHeadNeck = [upRegHeadNeck; labelHeadNeck.(1)(location), dataHeadNeck(location, 1), dataHeadNeck(location, 2)];
                headNeckVal = [headNeckVal; [2.^dataHeadNeck(location, 1), dataHeadNeck(location, 2)]];
            else
                headNeckVal = [headNeckVal; [0, 0]];
            end
        else
            headNeckVal = [headNeckVal; [NaN, NaN]];
        end
        location = find(table2array(labelLiver(:,1)) == string(gene));
        if (isempty(location) == 0)
            if ((dataLiver(location, 1) > 0) && (dataLiver(location, 2) < 0.05))
                upRegL = [upRegL; labelLiver.(1)(location), dataLiver(location, 1), dataLiver(location, 2)];
                liverVal = [liverVal; [2.^dataLiver(location, 1), dataLiver(location, 2)]];
            else
                liverVal = [liverVal; [0, 0]];
            end
        else
            liverVal = [liverVal; [NaN, NaN]];
        end
        location = find(table2array(labelLung(:,1)) == string(gene));
        if (isempty(location) == 0)
            if ((dataLung(location, 1) > 0) && (dataLung(location, 2) < 0.05))
                upRegLung = [upRegLung; labelLung.(1)(location), dataLung(location, 1), dataLung(location, 2)];
                lungVal = [lungVal; [2.^dataLung(location, 1), dataLung(location, 2)]];
            else
                lungVal = [lungVal; [0, 0]];
            end
        else
            lungVal = [lungVal; [NaN, NaN]];
        end
    end
    upRegHeadNeck = cell2table(upRegHeadNeck);
    upRegL = cell2table(upRegL);
    upRegLung = cell2table(upRegLung);
   %Validation Heatmaps
%     [~, pos, ~] = intersect(dataUpTable(:,1), ylabelUniversal);
    validationUpTable = [dataUpTable(1:31, :), num2cell(headNeckVal), num2cell(liverVal), num2cell(lungVal)];
    validationHeatMap = [];
    toDel = [];
    sorter = zeros(length(validationUpTable), 1);
    for i = 1:length(validationUpTable)
        heat = [];
        for j = 2:2:10
            if (validationUpTable{i,j} > 2)
                heat = [heat, 2];
                if (j<=14)
                    sorter(i) = sorter(i) + 10;
                end
            elseif ((validationUpTable{i,j} < 2) && (validationUpTable{i,j} > 0))
                heat = [heat, 1];
                if (j<=14)
                    sorter(i) = sorter(i) + 9.5;
                end
            elseif (validationUpTable{i,j} == 0)
                heat = [heat, 0];
                if (j<=14)
                    sorter(i) = sorter(i) + 1;
                end
            else
                heat = [heat, -1];
                if (j<=14)
                    sorter(i) = sorter(i) + 1.5;
                end
            end
        end
        for j = 12:2:20
            switch j
                case 12
                    pVal = 2e-12;
                case 14
                    pVal = 1.91e-4;
                case 16
                    pVal = 1.5e-7;
                case 18
                    pVal = 3.9e-4;
                case 20
                    pVal = 2.5e-19;
            end
            if ((validationUpTable{i,j+1} < pVal) && (validationUpTable{i,j+1} > 0))
                heat = [heat, 2];
                if (j<=14)
                    sorter(i) = sorter(i) + 10;
                end
            elseif (validationUpTable{i,j+1} > pVal)
                heat = [heat, 1];
                if (j<=14)
                    sorter(i) = sorter(i) + 9.5;
                end
            elseif (validationUpTable{i,j+1} == 0)
                heat = [heat, 0];
                if (j<=14)
                    sorter(i) = sorter(i) + 1;
                end
            else
                heat = [heat, -1];
                if (j<=14)
                    sorter(i) = sorter(i) + 1.5;
                end
            end
        end
        validationHeatMap = [validationHeatMap; heat];
    end
    validationHeatMap(toDel, :) = [];
    validationUpTable(toDel, :) = [];
    sorter(toDel, :) = [];
    [sorter, p2] = sort(sorter, 'descend');
    validationUpTable = validationUpTable(p2, :);
    validationHeatMap = validationHeatMap(p2, [2 3 4 1 5 6 7 8 9 10]);
    figure(2)
        xlabels2 = char({'SF8628', 'OVCAR8', 'Primary Dermal Fibroblast', 'RPMI8226', 'KEAP1 KO ARH-77', 'LUAD TCGA', 'NSCLC', 'HNSCC TCGA', 'LIHC TCGA', 'LUSC TCGA'});
        
        heatmap(xlabels2, validationUpTable(:, 1), double(validationHeatMap), 'Title', 'Validation for upregulated Genes', 'FontSize', 14, 'CellLabelColor', 'none', 'ColorbarVisible', 'on');
        cmap = colormap([0, 0, 0; 1, 1, 1; 1, 0.8, 0.8; 1, 0, 0]);
        axs = struct(gca); %ignore warning that this should be avoided
        c = axs.Colorbar;
        c.Ticks = [-1, 0, 1, 2];
        c.TickLabelsMode = 'manual';
        c.TickLabels = {'Not Expressed','Unchanged','Upregulated', 'Highly Upregulated'};
        set(gcf,'color','w');
        c.FontSize = 12;
%% Figure 4 Validation Plots
 core = 1:15;
 nearU = 16:31;
 validationCore = {'\it ABHD4', '\it GLCM', '\it OSGIN1', '\it SLC7A11', '\it SRXN1', '\it AKR1C3', '\it NQO1', '\it ~CBR3', '\it GCLC', '\it PIR', '\it FTH1', '\it GSR', '\it ME1', '\it FTL', '\it EPHX1'};
 validationnearU = {'\it AKR1C1', '\it TRIM16', '\it TRIM16L', '\it ABCB6', '\it SLC48A1', '\it MAFG', '\it PANX2', '\it PGD', '\it SQSTM1', '\it TXRND1', '\it AIFM2', '\it GLA', '\it RIT1', '\it PRDX1', '\it SLC6A6', '\it NQO2'};
figure(3)
        xlabels2 = char({'TCGA HNSC' });    
        heatmap(xlabels2, validationCore, double(validationHeatMap(core, [8])), 'Title', 'TCGA HNSC', 'FontSize', 14, 'CellLabelColor', 'none', 'ColorbarVisible', 'on');
        cmap = colormap([1, 0.8, 0.8; 1, 0, 0]);
        axs = struct(gca); %ignore warning that this should be avoided
        c = axs.Colorbar;
        c.Ticks = [-1, 0, 1, 2];
        c.TickLabelsMode = 'manual';
        c.TickLabels = {'Not Expressed','Unchanged','Upregulated', 'Highly Upregulated'};
        set(gcf,'color','w');
        c.FontSize = 12;
    figure(1)
        xlabels2 = char({'TCGA LUSC'});
        heatmap(xlabels2, validationCore, double(validationHeatMap(core, [10])), 'Title', 'TCGA LUSC', 'FontSize', 14, 'CellLabelColor', 'none', 'ColorbarVisible', 'on');
        cmap = colormap([1, 0.8, 0.8; 1, 0, 0]);
        axs = struct(gca); %ignore warning that this should be avoided
        c = axs.Colorbar;
        c.Ticks = [-1, 0, 1, 2];
        c.TickLabelsMode = 'manual';
        c.TickLabels = {'Not Expressed','Unchanged','Upregulated', 'Highly Upregulated'};
        set(gcf,'color','w');
        c.FontSize = 12;
    figure(4)
        xlabels2 = char({'TCGA LIHC'});
        heatmap(xlabels2, validationCore, double(validationHeatMap(core, [9])), 'Title', 'TCGA LIHC', 'FontSize', 14, 'CellLabelColor', 'none', 'ColorbarVisible', 'on');
        cmap = colormap([1, 1, 1; 1, 0.8, 0.8; 1, 0, 0]);
        axs = struct(gca); %ignore warning that this should be avoided
        c = axs.Colorbar;
        c.Ticks = [-1, 0, 1, 2];
        c.TickLabelsMode = 'manual';
        c.TickLabels = {'Not Expressed','Unchanged','Upregulated', 'Highly Upregulated'};
        set(gcf,'color','w');
        c.FontSize = 12; 
    figure(5)
        xlabels2 = char({'TCGA HNSC'});
        heatmap(xlabels2, validationnearU, double(validationHeatMap(nearU, [8])), 'Title', 'Near Universal Genes in TCGA HNSC', 'FontSize', 14, 'CellLabelColor', 'none', 'ColorbarVisible', 'on');
        cmap = colormap([1, 0.8, 0.8; 1, 0, 0]);
        axs = struct(gca); %ignore warning that this should be avoided
        c = axs.Colorbar;
        c.Ticks = [-1, 0, 1, 2];
        c.TickLabelsMode = 'manual';
        c.TickLabels = {'Not Expressed','Unchanged','Upregulated', 'Highly Upregulated'};
        set(gcf,'color','w');
        c.FontSize = 12;
    figure(2)
        xlabels2 = char({'TCGA LUSC'});
        heatmap(xlabels2, validationnearU, double(validationHeatMap(nearU, [10])), 'Title', 'Near Universal Genes in TCGA LUSC', 'FontSize', 14, 'CellLabelColor', 'none', 'ColorbarVisible', 'on');
        cmap = colormap([1, 0.8, 0.8; 1, 0, 0]);
        axs = struct(gca); %ignore warning that this should be avoided
        c = axs.Colorbar;
        c.Ticks = [-1, 0, 1, 2];
        c.TickLabelsMode = 'manual';
        c.TickLabels = {'Not Expressed','Unchanged','Upregulated', 'Highly Upregulated'};
        set(gcf,'color','w');
        c.FontSize = 12;

    figure(6)
        xlabels2 = char({'TCGA LIHC'});
        heatmap(xlabels2, validationnearU, double(validationHeatMap(nearU, [9])), 'Title', 'Near Universal Genes in TCGA LIHC', 'FontSize', 14, 'CellLabelColor', 'none', 'ColorbarVisible', 'on');
        cmap = colormap([1, 1, 1; 1, 0.8, 0.8; 1, 0, 0]);
        axs = struct(gca); %ignore warning that this should be avoided
        c = axs.Colorbar;
        c.Ticks = [-1, 0, 1, 2];
        c.TickLabelsMode = 'manual';
        c.TickLabels = {'Not Expressed','Unchanged','Upregulated', 'Highly Upregulated'};
        set(gcf,'color','w');
        c.FontSize = 12;
clear h i j c cmap xlabels2 heat dataFile1 dataLiver labelLiver liverVal lungVal headNeckVal gene labelHeadNeck labelLung
clear upRegL upRegLung upRegHeadNeck location axs pVal p2 pos toDel excelHeadNeckData excelLungData excelLiverData core nearU
%     [~, pos, ~] = intersect(labelHeadNeck.(1), downRegulatedFoldTable.(1));
% %     labelHeadNeck3 = labelHeadNeck(pos, :);
% %     dataHeadNeck3 = dataHeadNeck(pos, :);
% %     downRegHN = [];
% %     for i = 1:height(labelHeadNeck3)
% %         if ((dataHeadNeck3(i, 1) < 0) && (dataHeadNeck3(i, 2) < 0.05))
% %             downRegHN = [downRegHN; labelHeadNeck3.(1)(i), dataHeadNeck3(i, 1), dataHeadNeck3(i, 2)];
% %         end
% %     end
% %     downRegHN = cell2table(downRegHN);
% %     [~, pos, ~] = intersect(labelLiver, downRegulatedFoldTable.(1));
% %     labelLiver3 = labelLiver(pos, :);
% %     dataLiver3 = dataLiver(pos, :);
% %     downRegL= [];
% %     for i = 1:height(labelLiver3)
% %         if ((dataLiver3(i, 1) < 0) && (dataLiver3(i, 2) < 0.05))
% %             downRegL = [downRegL; labelLiver3(i), dataLiver3(i, 1), dataLiver3(i, 2)];
% %         end
% %     end
% %     downRegL = cell2table(downRegL);
%     clear excelHeadNeckData dataFile1 dataHeadNeck labelHeadNeck i excelKEAPData labelLiver3 labelLiver2 labelLiver
%     clear dataKEAP labelKEAP posK posK2 pos dataHeadNeck2 dataHeadNeck3 labelHeadNeck3 labelHeadNeck2 dataLiver3 dataLiver2 dataLiver
%     clear location headNeckVal liverVal heat xlabels2 j cmap toDel c upRegHN upRegL  
%     clear gene axs unchangedUpRegL unchangedUpRegHN 
 %Validation Bar Plot
%     upDownUnion = union(ylabelUniversal, ylabelDown);
%     %upDownUnion = ylabelUniversal;
%     [~, pos1, ~] = intersect(dataTCGA.(1), upDownUnion);
%     dataTCGA4 = dataTCGA(pos1, :);
%     dataTCGA4 = sortrows(dataTCGA4,"Var2","descend");
%     [~, pos, ~] = intersect(labelHeadNeck, dataTCGA4.(1));
%     labelHeadNeck4 = categorical(labelHeadNeck(pos, :));
%     dataHeadNeck4 = dataHeadNeck(pos, :);
%     y = [dataHeadNeck4(:, 1), dataTCGA4.(2)];
%     figure(6)
%     bar(labelHeadNeck4, y)
%     ylabel('Log2 Fold Change')
%     xlabel('Genes')
%     legend('Head and Neck TCGA', 'Lung TCGA')
%     excelKEAPData = '/Users/harshitakumar/Documents/Research/Letterio Research/Letterio Data/RNA Seq/TCGA_KEAP1_NFE2L2.xlsx';
%     [~, ~, dataFile1] = xlsread(excelKEAPData);
%     dataKEAP = cell2mat(dataFile1(2:end, [7 8]));
%     labelKEAP = dataFile1(2:end, 1);
%     [~, posK, ~] = intersect(labelKEAP, upRegulatedFoldTable.(1));
%     [~, posK2, ~] = intersect(labelKEAP, downRegulatedFoldTable.(1));
%     pos = union(posK, posK2);
%     labelKEAP = labelKEAP(pos, :);
%     dataKEAP = dataKEAP(pos, :);
%     upRegK = [];
%     downRegK = [];
%     for i = 1:height(labelKEAP)
%         if ((dataKEAP(i, 1) > 0) && (dataKEAP(i, 2) < 0.05))
%             upRegK = [upRegK; labelKEAP(i), dataKEAP(i, 1), dataKEAP(i, 2)];
%         end
%         if ((dataKEAP(i, 1) < 0) && (dataKEAP(i, 2) < 0.05))
%             downRegK = [downRegK; labelKEAP(i), dataKEAP(i, 1), dataKEAP(i, 2)];
%         end
%     end
%     [upRegCommon, posL, posK] = intersect(upRegL(:,1), upRegK(:,1));
%     upRegCommon = table(upRegCommon, upRegL(posL, 2), upRegL(posL, 3), upRegK(posK, 2), upRegK(posK, 3));
%     [downRegCommon, posL2, posK2] = intersect(downRegL(:,1), downRegK(:,1));
%     downRegCommon = table(downRegCommon, downRegL(posL2, 2), downRegL(posL2, 3), downRegK(posK2, 2), downRegK(posK2, 3)); 


%% Universal, almost, conditional, downregulated Heatmaps
%universal =  {'AKR1C3', 'CBR3', 'FTH1', 'GCLC', 'GCLM', 'NQO1', 'OSGIN1', 'PIR', 'SRXN1', 'ABHD4', 'ABCB6', 'AKR1C1', 'B4GALNT1', 'EPHX1', 'FTL', 'ME1', 'NQO2', 'SLC48A1', 'SLC6A6', 'TRIM16', 'TRIM16L', 'SLC7A11', 'AIFM2'};
universalHeatMap = [];
    universal = nRF27of7.(1);
    sorter2 = zeros(length(universal), 1);
    for i = 1:length(universal)
        [~, p, ~] = intersect(validationUpTable(:,1), universal(i));
        universalHeatMap = [universalHeatMap; validationHeatMap(p, 1:7)];
        sorter2(i) = sorter(p);
    end
    [sorter2, p2] = sort(sorter2, 'descend');
    universal = universal(p2);
    universalHeatMap = universalHeatMap(p2, :);
        figure(1)
        xlabels = char({'SF8628', 'OVCAR8', 'Primary Dermal Fibroblast', 'RPMI-8226', 'KEAP1 KO ARH-77', 'TCGA LUAD', 'NSCLC (CCLE)'});
        heatmap(xlabels, validationCore, double(universalHeatMap), 'Title', 'Candidate Core Genes', 'FontSize', 22, 'CellLabelColor', 'none', 'ColorbarVisible', 'on');
        cmap = colormap([1, 0.8, 0.8; 1, 0, 0]);
        axs = struct(gca); %ignore warning that this should be avoided
        c = axs.Colorbar;
        c.Ticks = [-2; -1; 0; 1; 2];
        c.TickLabelsMode = 'manual';
        c.TickLabels = {'Highly Upregulated', 'Upregulated'};
        c.FontSize = 16;
        set(gcf,'color','w');
nearUniversal = nRF26of7.(1);
    nearUniversal(4, :) = [];
    nearUniversalHeatmap = [];
    sorter3 = zeros(length(nearUniversal), 1);
    for i = 1:length(nearUniversal)
        [~, p, ~] = intersect(validationUpTable(:,1), nearUniversal(i));
        nearUniversalHeatmap = [nearUniversalHeatmap; validationHeatMap(p, 1:7)];
        sorter3(i) = sorter(p);
    end
    [sorter3, p2] = sort(sorter3, 'descend');
    nearUniversal = nearUniversal(p2);
    nearUniversalHeatmap = nearUniversalHeatmap(p2, :);
        figure(2)
        heatmap(xlabels, validationnearU, double(nearUniversalHeatmap), 'Title', 'Near Universal Genes', 'FontSize', 22, 'CellLabelColor', 'none', 'ColorbarVisible', 'on');
        cmap = colormap([ 0, 0, 0;  1, 1, 1; 1, 0.8, 0.8; 1, 0, 0]);
        axs = struct(gca); %ignore warning that this should be avoided
        c = axs.Colorbar;
        c.Ticks = [-1, 0, 1, 2];
        c.TickLabelsMode = 'manual';
        c.TickLabels = {'Not Expressed','Unchanged','Upregulated', 'Highly Upregulated'};
        c.FontSize = 16;
        set(gcf,'color','w');
%NOVELDOWNREGULATED
    novelDownregulated = {'TALDO1', 'G6PD', 'HMOX1', 'AKR1C2','PGD' , 'TKT', 'B4GALNT1', 'NR0B1', 'UNKL', 'POPDC3', 'AKIRIN2', 'LRP8'};
    downregulated = {'PLAT', 'EVA1A', 'CDC42EP3'};
    novelDownregulatedHeatmap = [];
    sorter4 = zeros(1, length(novelDownregulated));
    for i = 1:length(novelDownregulated)
        [~, p, ~] = intersect(dataUpTable(:,1), novelDownregulated(i));
        heat = [];
        for j = 2:2:10
            if (dataUpTable{p,j} > 2)
                heat = [heat, 2];
                sorter4(i) = sorter4(i) + 10;
            elseif ((dataUpTable{p,j} < 2) && (dataUpTable{p,j} > 0))
                heat = [heat, 1];
                sorter4(i) = sorter4(i) + 9.5;
            elseif (dataUpTable{p,j} == 0)
                heat = [heat, -1];
                sorter4(i) = sorter4(i) + 1;
            else
                heat = [heat, -2];
                sorter4(i) = sorter4(i) + 1.5;
            end
        end
        for j = 12:2:14
            switch j
                case 12
                    pVal = 2e-12;
                case 14
                    pVal = 1.91e-4;
            end
            if ((dataUpTable{p,j+1} < pVal) && (dataUpTable{p,j+1} > 0))
                heat = [heat, 2];
                sorter4(i) = sorter4(i) + 10;
            elseif (dataUpTable{p,j+1} > pVal)
                heat = [heat, 1];
                sorter4(i) = sorter4(i) + 9.5;
            elseif (dataUpTable{p,j+1} == 0)
                heat = [heat, -1];
                sorter4(i) = sorter4(i) + 1;
            else
                heat = [heat, -2];
                sorter4(i) = sorter4(i) + 1.5;
            end
        end
        novelDownregulatedHeatmap = [novelDownregulatedHeatmap; heat];
    end
    [sorter4, p2] = sort(sorter4, 'descend');
    novelDownregulated = novelDownregulated(p2);
    novelDownregulatedHeatmap = novelDownregulatedHeatmap(p2, :);
    for i = 1:length(downregulated)
        [~, p, ~] = intersect(dataDownTable(:,1), downregulated(i));
        [~, p2, ~] = intersect(dataUpTable(:,1), downregulated(i))
        heat = [];
        for j = 2:2:14
            if (((dataDownTable{p,j} < 1) && (dataDownTable{p,j} > 0)) || (dataDownTable{p,j} < 0))
                heat = [heat, 0];
            elseif (dataDownTable{p,j} == 0)
                heat = [heat, -1];
            else
                heat = [heat, -2];
            end
        end
        novelDownregulatedHeatmap = [novelDownregulatedHeatmap; heat];
    end 
    %novelDownregulated = union(novelDownregulated, downregulated, 'stable');
    novelDownregulated = {'\it PGD'	'\it AKR1C2'	'\it HMOX1'	'\it POPDC3'	'\it B4GALNT1'	'\it UNKL'	'\it NR0B1'	'\it LRP8'	'\it TALDO1'	'\it G6PD'	'\it AKIRIN2'	'\it TKT'	'\it PLAT'	'\it EVA1A'	'\it CDC42EP3'};
    novelDownregulatedHeatmap = novelDownregulatedHeatmap(:, [[2 3 4 1 5 6 7]]);
    novelDownregulatedHeatmap(end, 4) = 1;
    figure(3)
        heatmap(xlabels, novelDownregulated, double(novelDownregulatedHeatmap), 'Title', 'Conditional and Novel Target Genes',  'FontSize', 22, 'CellLabelColor', 'none', 'ColorbarVisible', 'on');
        cmap2 = colormap([ 0, 0, 0;  1, 1, 1; 0, 0, 0.5; 1, 0.8, 0.8; 1, 0, 0]);
        axs2 = struct(gca); %ignore warning that this should be avoided
        c2 = axs2.Colorbar;
        c2.Ticks = [-2, -1, 0, 1, 2];
        c2.TickLabelsMode = 'manual';
        c2.TickLabels = {'Not Expressed','Unchanged','Downregulated','Upregulated', 'Highly Upregulated'};
        c2.FontSize = 16;
        set(gcf,'color','w');
figure(4)
    heatmap(xlabels, validationUpTable(:, 1), double(validationHeatMap(:, 1:7)), 'Title', 'Universal and Near Universal upregulated Genes', 'FontSize', 22, 'CellLabelColor', 'none', 'ColorbarVisible', 'on');
    cmap = colormap([0, 0, 0; 1, 1, 1; 1, 0.8, 0.8; 1, 0, 0]);
    axs = struct(gca); %ignore warning that this should be avoided
    c = axs.Colorbar;
    c.Ticks = [-1, 0, 1, 2];
    c.TickLabelsMode = 'manual';
    c.TickLabels = {'Not Expressed','Unchanged','Upregulated', 'Highly Upregulated'};
    set(gcf,'color','w');
    c.FontSize = 16;
clear c c2  heat axs2 axs cmap cmap2 j b h p2 labels  nearUnversalHeatmap  pVal  
clear downregulated nearUniversal  universal universalHeatMap p i xlabels  sorter2 sorter3 sorter4 xlabels2  
%% Conditional Genes
mutantUnion = table([nRF2GenesK.labelK; hmCLData.(1); hmTCGAData.(1)]);
mutantUnion.Properties.VariableNames{1} = 'Labels';
geneCounts = groupcounts(mutantUnion, "Labels");
nRF23of3 = [];
for i = 1:size(geneCounts, 1)
    if (table2array(geneCounts(i, 2)) == 3)
        nRF23of3 = [nRF23of3; geneCounts(i, 1)];
    end
end
mutantUnion = setdiff(nRF23of3, nRF25of7);
mutantUnion1 = setdiff(nRF23of3, nRF24of7);
clear i geneCounts
%% Fold Change Table
logMinFold = table(log2(percentUpTable.(9)(1:31)), 'VariableNames', "logMinFold");
logMaxFold = table(log2(percentUpTable.(10)(1:31)),  'VariableNames', "logMaxFold");
foldTable = [percentUpTable(1:31, [1 5 9 10 7]), logMinFold, logMaxFold];
universalfoldTable = foldTable(1:15, :);
nearUniversalfoldTable = foldTable(16:end, :);
foldTable = sortrows(foldTable,"avFoldChange","descend");
universalfoldTable = sortrows(universalfoldTable,"avFoldChange","descend");
nearUniversalfoldTable = sortrows(nearUniversalfoldTable,"avFoldChange","descend");
clear logMinFold logMaxFold
%% PCA Plot
   labels = {'Control Data', 'Low Conc Data', 'High Conc Data'};
   cellLabel = {'RPMI 8226 PCA Plot','SF8628 PCA Plot', 'OVCAR8 PCA Plot', 'Primary Dermal Fibroblast PCA Plot', 'KEAP1 KO ARH-77 PCA Plot'};
   for k = 5:5
       switch k
           case 1
               cellLine = controlHighTable2MM;
           case 2
               cellLine = controlHighTableDIPG;
           case 3
               cellLine = controlHighTableOv;
           case 4
               cellLine = controlHighTableFib;
           case 5
               cellLine = controlTablePCK;
       end
       if (k == 5)
            data = log2(table2array(cellLine(:, 8:9)));
       else
           data = log2(table2array(cellLine(:, 14:16)));
       end
       toDelete = [];
       for i = 1:length(data)
           if ((mean(data(i, :)<0)) > 0)
               toDelete = [toDelete; i];
           end
       end
       data(toDelete, :) = [];
       data = data';
       [coeffD, ~, ~, ~, explained] = pca(data);
       samplePC1score = ones(1,height(data));
       samplePC2score = ones(1,height(data));
       for i = 1:height(data)
           for j = 1:width(data)
            samplePC1score(i) = samplePC1score(i)+(data(i,j)*coeffD(j,1));
            %samplePC2score(i) = samplePC2score(i)+(data(i,j)*(perc*coeffD(j,2)+(1-perc)*coeffD(j,3)));
            samplePC2score(i) = samplePC2score(i)+(data(i,j)*coeffD(j,2));
           end
       end
       figure()
       scatter(samplePC1score, samplePC2score)
       hold on
       if (k == 5)
          s1 = scatter(samplePC1score(1:3), samplePC2score(1:3), 60, 'b', "filled");
          s2 = scatter(samplePC1score(4:6), samplePC2score(4:6), 60, 'r', "filled");
           legend([s1, s2], {'Wildtype Data', 'KEAP1 Knockout Data'})
       else
           s1 = scatter(samplePC1score(1:3), samplePC2score(1:3), 60, 'b', "filled");
           s2 = scatter(samplePC1score(4:6), samplePC2score(4:6), 60, 'r', "filled");
           s3 = scatter(samplePC1score(7:9), samplePC2score(7:9), 60, 'k', "filled"); 
            legend([s1, s2, s3], labels)
       end
       title(cellLabel(k), 'FontSize', 14)
       xlabel('PC1 Score')
       ylabel('PC2 Score')
       set(gcf,'color','w');
       hold off
   end

%% Combined PCA Plot
   labels = {'Control Data MM', 'Low Conc Data MM', 'High Conc Data MM', 'Control Data DIPG', 'Low Conc Data DIPG', 'High Conc Data DIPG', 'Control Data Ov', 'Low Conc Data Ov', 'High Conc Data Ov'};
   %Combining data for the three cell lines 
        combinedData = intersect(intersect(controlHighTableOv.labelHighOv, controlHighTableDIPG.labelHighDIPG), controlHighTableMM.labelHighMM);
        %creating tables for high
            %MM
            [~, posMM, ~] = intersect( controlHighTableMM.labelHighMM, combinedData);
            combinedDataMM  = [combinedData, controlHighTableMM(posMM, 14:16)];
            %DIPG
            [~, posDIPG, ~] = intersect( controlHighTableDIPG.labelHighDIPG, combinedData);
            combinedDataDIPG  = [combinedData, controlHighTableDIPG(posDIPG, 14:16)];
            %Ovariance cells
            [~, posOv, ~] = intersect( controlHighTableOv.labelHighOv, combinedData);
            combinedDataOv  = [combinedData, controlHighTableOv(posOv, 14:16)];
            %joining all datasets together
            combinedData = join(join(combinedDataMM, combinedDataDIPG), combinedDataOv);
            combinedData = log2(table2array(combinedData(:, 2:10)));
            toDelete = [];
            for i = 1:length(combinedData)
               if ((mean(combinedData(i, :)<0)) > 0)
                   toDelete = [toDelete; i];
               end
            end
           combinedData(toDelete, :) = [];
            combinedData = combinedData';
   %doing the PCA Analysis
       [coeffCombined, scoreCombined, ~, ~, explained] = pca(combinedData);
       samplePC1score = ones(1,27);
       samplePC2score = ones(1,27);
       for i = 1:width(combinedData)
           for j = 1:height(combinedData)
            samplePC1score(j) = samplePC1score(j)+(combinedData(j,i)*coeffCombined(j,1));
            samplePC2score(j) = samplePC2score(j)+(combinedData(j,i)*coeffCombined(j,2));
           end
       end
       figure(4)
       s1 = scatter(samplePC1score(1:3), samplePC2score(1:3), 'bo');
       hold on
       s2 = scatter(samplePC1score(4:6), samplePC2score(4:6), 'ro');
       s3 = scatter(samplePC1score(7:9), samplePC2score(7:9), 'ko');
       s4 = scatter(samplePC1score(10:12), samplePC2score(10:12), 'b+');
       s5 = scatter(samplePC1score(13:15), samplePC2score(13:15), 'r+');
       s6 = scatter(samplePC1score(16:18), samplePC2score(16:18), 'k+');
       s7 = scatter(samplePC1score(19:21), samplePC2score(19:21), 'b*');
       s8 = scatter(samplePC1score(22:24), samplePC2score(22:24), 'r*');
       s9 = scatter(samplePC1score(25:27), samplePC2score(25:27), 'k*');
       legend([s1, s2, s3, s4, s5, s6, s7, s8, s9], labels)
       title('Combine PCA plot')
       xlabel('PC1 score')
       ylabel('PC2 score')
       hold off
%% Promoter 7of7
% labels = {'TSS', 'ARE-1', 'ARE-3',};
% labels2 = {'ARE-2', 'ARE-4'}; 
labels =  {'TSS', 'ARE-1', 'ARE-2', 'ARE-3', 'ARE-4'};
abhd4 = [0, 1152, 1525, 2357];
gclm = [0, 39];
osgin1 = [0,82, 540];
slc7a11 = [0, 127];
srxn1 = [0,349, 239];
akr1c3 = [0,4805, 1490, 2144];
nqo1 = [0,386];
cbr3 = [0,164, 1389, 2476, 2547];
gclc = [0,3129, 2065, 2755];
pir = [0,2962];
fth1 = [0,4392, 4433, 2057];
gsr  = [0, 860, 4706];
me1 = [0, 2630, 765, 2004, 2779];
ftl = [0, 1368, 3815];
ephx1 = [0, 26, 1718];
titles = {'\it ABHD4', '\it OSGIN1', '\it GCLM', '\it NQO1', '\it SRXN1', '\it SLC7A11', '\it EPHX1', '\it AKR1C3', '\it PIR','\it GCLC','\it CBR3', '\it ME1','\it FTH1', '\it FTL', 'GSR'};
figure(5)
t = tiledlayout(8, 2, 'TileSpacing','compact', 'TileIndexing', 'columnmajor');
title(t, 'Predicted NRF2 Sites in Core Genes', 'FontWeight', 'bold', 'FontSize', 24)
%xlabel(t, 'Base Pairs', 'FontSize', 20)
%m = {'0', '1000', '2000', '3000', '4000', '5000'};
for i = 1:15
    switch i
        case 1
            gene = abhd4;
            l = labels(1:length(gene));
            g = {'', '1152', '1525', '2357'};
        case 2
            gene = osgin1;
            l = {'', 'ARE-1 TSS', 'ARE-2'};
            g = {'', '82', '540'};
        case 3
            gene = gclm;
            l = {'', 'ARE-1 TSS'};
            g = {'', '39'};
        case 4
            gene = nqo1;
            l = {'TSS', 'ARE-1'};
            g = {'', '386'};
        case 5            
            gene = srxn1;
            l = {'TSS', 'ARE-2 ARE-1', ' '};
            g = {'', '349', '239'};
        case 6
            gene = slc7a11;
            l = {'', 'ARE-1 TSS'};
        case 7
            gene = ephx1;
            l = {'', 'ARE-1 TSS', 'ARE-2'};
        case 8           
            gene = akr1c3;
            l = labels(1:length(gene));
            g = {'', '4805', '1490', '2144'};
        case 9
            gene = pir;
            l = labels(1:length(gene));
            g = {'', '2962'};
        case 10            
            gene = gclc;
            l = {'TSS', 'ARE-1',  'ARE-2', 'ARE-3'};
            g = {'', '3129', '2065', '2755'};
        case 11
            gene = cbr3;
            l = {'', 'ARE-1 TSS',  'ARE-2', 'ARE-3 ARE-4', ' '};
            g = {'', '164', '1389', '2476 2547', ' '};
        case 12
            gene = me1;
            l = {'TSS', 'ARE-4 ARE-1',  'ARE-2', 'ARE-3', ''};
        case 13
            gene = fth1;
            l = {'TSS', 'ARE-2 ARE-1', ' ', 'ARE-3'};
            g = {'', '4433 4392', ' ', '2057'};
        case 14
            gene = ftl;
            l = {'TSS', 'ARE-1',  'ARE-2',};
        case 15
            gene = gsr;
            l = labels(1:length(gene));
    end
    y = zeros(1, length(gene));
    nexttile
    plot(gene(2:end), y(2:end)+0.025, 'Marker', "v", 'MarkerSize', 8, 'MarkerFaceColor', "k", 'LineWidth', 1.5, 'LineStyle', "none", 'Color', "k")
    hold on
    plot([0], y(1)+0.025, 'Marker', "v", 'MarkerSize', 8, 'LineWidth', 1.5, 'LineStyle', "none", 'Color', "k")
    text(gene, (y+0.06), l, 'Color', "k", 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center')
    set(gca, 'YTick', [], 'XTick', []);
    set (gca, 'xdir', 'reverse');
    set(gca, 'fontsize', 4)
    set(gca,'XColor','none', 'YColor','none')
    set(gcf,'color','w');
    %axis off;
    if (i <= 7)
        xlim([-100 2600])
        plot([0, 1000, 2000, 2500],  zeros(1, 4), 'Marker', "|", 'MarkerSize', 12, 'LineWidth', 1, 'Color', "k")
        plot((0:2550), zeros(1, 2551) + 1, 'LineWidth', 2, 'Color', "k")
        text([0 1000 2000], (zeros(1, 3)-0.08), {'0', '1000',  '2000'}, 'Color', "k", 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left', 'FontSize', 10)
    elseif (i > 7)
        xlim([-150 5150])
        plot([0, 1000, 2000, 3000, 4000, 5000],  zeros(1, 6), 'Marker', "|", 'MarkerSize', 12, 'LineWidth', 1, 'Color', "k")
        plot((0:5050), zeros(1, 5051) + 1, 'LineWidth', 2, 'Color', "k")
        text([0, 1000, 2000, 3000, 4000, 5000], (zeros(1, 6)-0.08),  {'0', '1000', '2000', '3000', '4000', '5000'}, 'Color', "k", 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left', 'FontSize', 10)
    end
    ylim([-0.15 0.15])
    hold off
    title(string(titles(i)), 'FontSize', 16)
    if (i == 7)
        nexttile
        set(gca, 'YTick', [], 'XTick', []);
        set(gca,'XColor','none', 'YColor','none')
        xlabel('Base Pairs', 'FontSize', 20, 'Color', "k")
    end
    if (i == 15)
        xlabel('Base Pairs', 'FontSize', 20, 'Color', "k")
    end
end
clear j i y title labels maxRound l t g titles gene markings m ftl
clear gclm me1 slc48a1 trim16 trim16L abcb6 akr1c1 nqo1 osgin1 pir srxn1 akr1c3 abhd4 cbr3 fth1 gclc 
%% Promoter 6of7
% labels = {'TSS', 'ARE-1', 'ARE-3',};
% labels2 = {'ARE-2', 'ARE-4'}; 
labels =  {'TSS', 'ARE-1', 'ARE-2', 'ARE-3', 'ARE-4', 'ARE-5', 'ARE-6', 'ARE-7'};
trim16 = [0, 489, 4137, 2616, 985, 4888];
trim16L = [0, 666, 4049, 3901, 56, 313];
sqstm1 = [0, 2303, 4283];
slc48a1 = [0, 1274, 4528];
nqo2 = [0, 2426, 3618];
mafg = [0, 3419];
akr1c1 = [0, 3583, 1717, 4718];
abcb6 = [0, 121, 4511];
txnrd1 = [0, 1505, 2895, 3875];
aifm2 = [0, 1948];
pgd = [0, 3372, 166];
gla = [0, 4406, 2865, 721];
rit1 = [0, 2196, 1062, 1307, 4464, 2032, 2860];
prdx1 = [0, 88, 3290];
titles = {'\it TRIM16', '\it TRIM16L', '\it SQSTM1', '\it SLC48A1', '\it NQO2', '\it MAFG', '\it PGD', '\it AKR1C1','\it ABCB6','\it TXNRD1', '\it GLA','\it AIFM2', '\it RIT1', '\it PRDX1'};
figure(5)
t = tiledlayout(7, 2, 'TileSpacing','compact', 'TileIndexing', 'columnmajor');
title(t, 'Predicted NRF2 Sites in Near Universal Genes', 'FontWeight', 'bold', 'FontSize', 24)
for i = 1:14
    switch i
        case 1
            gene = trim16;
            l = labels(1:length(gene));
        case 2
            gene = trim16L;
            l = {'ARE-5 ARE-4 TSS', 'ARE-1', 'ARE-2 ARE-3', '', '', ''};
        case 3
            gene = sqstm1;
            l = labels(1:length(gene));
        case 4
            gene = slc48a1;
            l = labels(1:length(gene));
        case 5            
            gene = nqo2;
            l = labels(1:length(gene));
        case 6
            gene = mafg;
            l = labels(1:length(gene));
        case 7            
            gene = pgd;
            l = {'ARE-2 TSS', 'ARE-1',''};
        case 8
            gene = akr1c1;
            l = labels(1:length(gene));
        case 9            
            gene = abcb6;
            l = {'ARE-1 TSS', '','ARE-2'};
        case 10
            gene = txnrd1;
            l =  labels(1:length(gene));
        case 11
            gene = gla;
            l = labels(1:length(gene));
        case 12
            gene = aifm2;
            l = labels(1:length(gene));
        case 13
            gene = rit1;
            l = {'TSS', 'ARE-1 ARE-5', 'ARE-3 ARE-2', '', 'ARE-4', '', 'ARE-6'};
        case 14
            gene = prdx1;
            l = {'ARE-1 TSS', '','ARE-2'};
    end
    y = zeros(1, length(gene));
    nexttile
    plot(gene(2:end), y(2:end)+0.025, 'Marker', "v", 'MarkerSize', 8, 'MarkerFaceColor', "k", 'LineWidth', 1.5, 'LineStyle', "none", 'Color', "k")
    hold on
    plot([0], y(1)+0.025, 'Marker', "v", 'MarkerSize', 8, 'LineWidth', 1.5, 'LineStyle', "none", 'Color', "k")
    text(gene, (y+0.06), l, 'Color', "k", 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center')
    set(gca, 'YTick', [], 'XTick', []);
    set (gca, 'xdir', 'reverse');
    set(gca, 'fontsize', 4)
    set(gca,'XColor','none', 'YColor','none')
    set(gcf,'color','w');
    %axis off;
    xlim([-150 5150])
    plot([0, 1000, 2000, 3000, 4000, 5000],  zeros(1, 6), 'Marker', "|", 'MarkerSize', 12, 'LineWidth', 1, 'Color', "k")
    plot((0:5050), zeros(1, 5051) + 1, 'LineWidth', 2, 'Color', "k")
    text([0, 1000, 2000, 3000, 4000, 5000], (zeros(1, 6)-0.08),  {'0', '1000', '2000', '3000', '4000', '5000'}, 'Color', "k", 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left', 'FontSize', 10)
    ylim([-0.15 0.15])
    hold off
    title(string({titles(i)}), 'FontSize', 16)
    if (i == 7) || (i == 14)
        xlabel('Base Pairs', 'FontSize', 20, 'Color', "k")
    end
end
clear j i y title labels maxRound l t g titles gene markings m ftl
clear gclm me1 slc48a1 trim16 trim16L abcb6 akr1c1 nqo1 osgin1 pir srxn1 akr1c3 abhd4 cbr3 fth1 gclc 
clear trim16 trim16L gsr aifm2 nqo2 ephx1 mafg sqstm1 slc6a6 akr1c1 abcb6 txnrd1 
%% Promoter Conditional
% labels = {'TSS', 'ARE-1', 'ARE-3',};
% labels2 = {'ARE-2', 'ARE-4'}; 
labels =  {'TSS', 'ARE-1', 'ARE-2', 'ARE-3', 'ARE-4', 'ARE-5', 'ARE-6', 'ARE-7'};
akr1c2 = [0, 1655, 1689, 2803, 835];
hmox1 = [0, 3912, 3976, 4019, 3134, 4351, 283];
pgd = [0, 3372, 166];
popdc3 = [0, 4853, 529, 452, 303];
b4galnt1 = [0, 2470, 1053, 2276, 678, 2390];
lrp8 = [0, 4268, 3331, 4919, 3291];
unkl = [0, 4038, 3012, 300];
g6pd = [0, 2015];
akirin2 = [0, 4428, 4565, 3162, 1466, 1853];
tkt = [0, 1019, 2302];
taldo1 = [0, 922];
plat = [0, 3326];
eva1a = [0, 3561, 2369, 3001];
titles = {'\it AKR1C2', '\it HMOX1', '\it PGD', '\it POPDC3', '\it B4GALANT1', '\it LRP8', '\it UNKL', '\it G6PD','\it AKIRIN2','\it TKT', '\it TALDO1','\it PLAT', '\it EVA1A'};
figure(5)
t = tiledlayout(7, 2, 'TileSpacing','compact', 'TileIndexing', 'columnmajor');
title(t, 'Predicted NRF2 Sites in Conditional Genes', 'FontWeight', 'bold', 'FontSize', 24)
for i = 1:13
    switch i
        case 1
            gene = akr1c2;
            l = {'TSS', 'ARE-1 ARE-2', '', 'ARE-3', 'ARE-4'};
        case 2
            gene = hmox1;
            l = {'TSS', '', '', 'ARE-5 ARE-3 ARE-2 ARE-1', 'ARE-4', '', 'ARE-6'};
        case 3
            gene = pgd;
            l = {'', 'ARE-1', 'ARE-2 TSS'};
        case 4
            gene = popdc3;
            l =  {'TSS', 'ARE-1', 'ARE-2 ARE-3 ARE-4', '', ''};
        case 5            
            gene = b4galnt1;
            l = {'TSS', 'ARE-1 ARE-3 ARE-5', 'ARE-2', '', 'ARE-4', ''};
        case 6
            gene = lrp8;
            l = {'TSS', 'ARE-1', 'ARE-2 ARE-4', 'ARE-3', ''};
        case 7            
            gene = unkl;
            l = labels(1:length(gene));
        case 8
            gene = g6pd;
            l = labels(1:length(gene));
        case 9            
            gene = akirin2;
            l = {'TSS', '', 'ARE-2 ARE-1', 'ARE-3', 'ARE-4', 'ARE-5'};
        case 10
            gene = tkt;
            l =  labels(1:length(gene));
        case 11
            gene = taldo1;
            l = labels(1:length(gene));
        case 12
            gene = plat;
            l = labels(1:length(gene));
        case 13
            gene = eva1a;
            l = labels(1:length(gene));
    end
    y = zeros(1, length(gene));
    nexttile
    plot(gene(2:end), y(2:end)+0.025, 'Marker', "v", 'MarkerSize', 8, 'MarkerFaceColor', "k", 'LineWidth', 1.5, 'LineStyle', "none", 'Color', "k")
    hold on
    plot([0], y(1)+0.025, 'Marker', "v", 'MarkerSize', 8, 'LineWidth', 1.5, 'LineStyle', "none", 'Color', "k")
    text(gene, (y+0.06), l, 'Color', "k", 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center')
    set(gca, 'YTick', [], 'XTick', []);
    set (gca, 'xdir', 'reverse');
    set(gca, 'fontsize', 4)
    set(gca,'XColor','none', 'YColor','none')
    set(gcf,'color','w');
    %axis off;
    xlim([-150 5150])
    plot([0, 1000, 2000, 3000, 4000, 5000],  zeros(1, 6), 'Marker', "|", 'MarkerSize', 12, 'LineWidth', 1, 'Color', "k")
    plot((0:5050), zeros(1, 5051) + 1, 'LineWidth', 2, 'Color', "k")
    text([0, 1000, 2000, 3000, 4000, 5000], (zeros(1, 6)-0.08),  {'0', '1000', '2000', '3000', '4000', '5000'}, 'Color', "k", 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left', 'FontSize', 10)
    ylim([-0.15 0.15])
    hold off
    title(string(titles(i)), 'FontSize', 16)
    if (i == 6) 
        nexttile
        set(gca, 'YTick', [], 'XTick', []);
        set(gca,'XColor','none', 'YColor','none')
        xlabel('Base Pairs', 'FontSize', 20, 'Color', "k")
    elseif (i == 13)
        xlabel('Base Pairs', 'FontSize', 20, 'Color', "k")
    end
end
clear j i y title labels maxRound l t g titles gene markings m ftl
clear akr1c2 hmox1 pgd popdc3 b4galnt1 lrp8 unkl g6pd akirin2 tkt taldo1 plat eva1a 
%% Rattus 70f7
rattus70f7a = [[-0.08547926	-0.026993275	0.003373623	0.6840873	0.38057947	0.39890814	-0.07499361	0.057988644	-0.040129185	-0.12560415	-0.003373623	0.13019896];...
    [0.001756668	0.13087082	0.028390884	-0.31110573	-0.001756668	-0.35940933	-0.26206112	-0.01307869	0.27232838	0.071320534	-0.123641014	0.16300297];...
    [-0.27887392	-0.17737627	-0.036751747	0.032576084	0.2451973	0.20505047	0.05343628	-0.032576084	0.3279872	0.0922513	-0.22997761	-0.1192565];...
    [-0.065969944	-0.003983021	0.003983021	0.44717073	0.083158016	0.029440403	-0.13136339	0.012561321	-0.06991625	-0.1005702	-0.044097424	0.052853107]; ...
    [0.27639866	0.28521538	0.28268814	1.2726927	2.1871662	2.0407848	-0.5739155	-0.27639866	-0.7423353	-0.7563658	-0.688859	-0.5368242];...
    [0.2985506	0.31991816	0.7535105	1.8996816	2.693192	2.6242642	-0.72481966	-1.088253	-0.3948183	-0.6158204	-0.2985506	-0.8245778];...
    [0.05266261	-0.12314439	0.0996325	3.749725	3.5828378	3.2791593	-1.1544864	-1.3422458	0.045860052	-0.045860052	-0.4043057	-0.20786738]; ...
    [0.59802437	0.3192587	1.3201399	0.63360214	1.3114042	1.4149666	-0.3192587	-1.3511839	-0.9114132	-0.72474384	-1.1477575	-0.45223713]; ...
    [0.36187935	0.26521492	0.57605934	2.5559025	2.8461876	2.9115114	-0.60118294	-0.33578968	-0.8442192	-0.35221386	-0.26521492	-0.31944656]; ...
    [0.6665025	0.8888078	0.32436132	4.7739596	5.0392833	5.38193	-1.2075772	-0.5807185	-0.8259387	-0.5315938	-0.32436132	-0.43174887];...
    [0.10134077	0.07405424	-0.057500362	0.39992	0.03770399	-0.07377577	-0.1454339	-0.00522089	-0.10530138	-0.025844097	0.00522089	0.013363361];...
    [0.8059311	0.6190233	0.30859756	0.47102165	1.0372715	1.6167717	-1.1561909	-0.30859756	-1.3604107	-1.314148	-0.9083805	-0.90011024];...
    [0.6079626	-1.0793166	-0.52570343	0.7108698	4.4091873	4.61061	-1.2035894	1.9811444	-0.17108297	-0.27140093	0.1710825	-3.6877708];...
    [0.19602346	-0.001597881	0.29579306	0.9350333	1.3509383	1.1751475	-0.3841977	-0.39139462	-0.18834066	-0.44583845	-0.20187235	0.001597881];...
    [0.66146374	0.4324093	0.44111633	0.5217581	1.160284	1.3183422	-0.61796665	-0.87683296	-0.4324093	-0.8793154	-0.88236904	-1.0129938]];
rattus70f7b = [[-0.0495553	-0.15832186	-0.06729174	0.44345713	0.30016947	0.2013917	-0.34190512	-0.066009045	0.0495553	0.056107044	-0.061666965	0.19010115];...
    [-0.05632329	0.07536292	0.1763823	0.05632329	-0.1284883	-0.4286325	-0.188267	0.56004405	-0.11716008	0.09253764	-0.34012008	0.21561837];...
    [-0.1365521	-0.013460875	-0.089916945	0.44674754	0.23643231	0.19158053	-0.12608266	-0.35386539	0.013460875	-0.18655658	0.04132819	0.032431364];...
    [0.002337456	-0.085933685	-0.002337456	0.13838196	0.41851044	-0.034259796	-0.31960487	0.049660683	0.039461136	-0.094566345	-0.21005535	0.049420357];...
    [-0.20163584	0.48471117	0.17060423	2.2529025	2.2888713	2.2725673	0.113063335	-0.1352458	-0.40587187	-0.27093935	-0.113063335	-0.19714212];...
    [-0.84990764	0.30753684	-0.50576377	2.1881921	2.8412054	3.238423	-0.09103942	-0.65750337	0.09103942	-0.19520593	0.5005739	-0.18941808];...
    [0.20539379	0.089307785	0.50303555	3.2197104	2.8201094	3.5252218	-1.3895893	-0.86378336	-0.3621807	-0.089307785	-0.29171944	-0.29030037];...
    [0.019182682	0.4005084	0.12930441	2.8165011	1.8799701	2.288179	-0.019182682	-0.2468791	-0.075912	-1.0400968	-0.1648078	-0.7661576];...
    [0.13184595	0.6231618	0.34803343	2.805695	2.6876397	2.9259152	-0.7068305	-0.7855725	-0.5640702	-0.21293211	-0.13184595	-0.30245924];...
    [-0.14057398	0.30631876	0.36156607	4.563579	4.854404	4.8903546	-0.20042181	-0.24891615	-0.06485319	-0.5682645	0.06485319	-0.1321578];
    [-0.007493973	-0.022828102	0.007493973	0.14952946	0.53705883	-0.058659554	-0.14831448	0.11898041	0.041550636	-0.046961784	-0.13361931	0.1397323];...
    [0.23740292	1.3355913	0.9812765	1.6194277	1.0959206	1.773695	-0.9456482	-0.40318298	-0.23740292	-0.6869354	-0.24196815	-0.2775631	];...
    [0.096223354	-0.012323856	0.03261423	0.43440294	0.108043194	0.23175478	-0.32266474	-0.14012575	-0.074615955	-0.103168964	-0.13313341	0.012323856];...
    [-0.18847513	0.4892869	0.04904127	1.5745759	1.5186143	1.8739476	0.011714458	-0.29274702	-0.22419405	-0.2505412	-0.011714458	-0.2065444];...
    [0.25321627	0.983294	0.6660342	1.6951327	1.1540504	2.0333495	-0.3071618	-0.60590315	-0.37934732	-0.32645464	-0.25321627	-0.3258958]];
rattus70f7a = 2.^rattus70f7a;
rattus70f7b = 2.^rattus70f7b;
labels = {'Abhd4', 'Akr1c3', 'Cbr3', 'Fth1', 'Gclc', 'Gclm', 'Nqo1', 'Osgin1', 'Pir', 'Srxn1', 'Ftl', 'Me1', 'Slc7a11', 'Gsr', 'Ephx1'};
rattusTablea = [mean(rattus70f7a(:, 1:3), 2), mean(rattus70f7a(:, 4:6), 2), mean(rattus70f7a(:, 7:9), 2), mean(rattus70f7a(:, 10:12), 2),...
    std(rattus70f7a(:, 1:3), 0, 2), std(rattus70f7a(:, 4:6), 0, 2), std(rattus70f7a(:, 7:9), 0, 2), std(rattus70f7a(:, 10:12),0,  2)];
rattusTableb = [mean(rattus70f7b(:, 1:3), 2), mean(rattus70f7b(:, 4:6), 2), mean(rattus70f7b(:, 7:9), 2), mean(rattus70f7b(:, 10:12), 2),...
    std(rattus70f7b(:, 1:3), 0, 2), std(rattus70f7b(:, 4:6), 0, 2), std(rattus70f7b(:, 7:9), 0, 2), std(rattus70f7b(:, 10:12),0,  2)];
pValueRattusa = [mattest(rattus70f7a(:, 1:3), rattus70f7a(:, 4:6), 'VarType', 'equal',  'Labels', labels),...
    mattest(rattus70f7a(:, 7:9), rattus70f7a(:, 10:12), 'VarType', 'equal',  'Labels', labels),...
    mattest(rattus70f7a(:, 4:6), rattus70f7a(:, 10:12), 'VarType', 'equal',  'Labels', labels)];
pValueRattusb = [mattest(rattus70f7b(:, 1:3), rattus70f7b(:, 4:6), 'VarType', 'equal',  'Labels', labels),...
    mattest(rattus70f7b(:, 7:9), rattus70f7b(:, 10:12), 'VarType', 'equal',  'Labels', labels),...
    mattest(rattus70f7b(:, 4:6), rattus70f7a(:, 10:12), 'VarType', 'equal',  'Labels', labels)];
sumSorter = sum(pValueRattusa(:, [1 3]), 2) + sum(pValueRattusb(:, [1 3]), 2);
[~, p] = sort(sumSorter, 'ascend');
rattusTablea = rattusTablea(p, :);
rattusTableb = rattusTableb(p, :);
pValueRattusa = pValueRattusa(p, :);
pValueRattusb = pValueRattusb(p, :);
rattus70f7a = rattus70f7a(p, :);
rattus70f7b = rattus70f7b(p, :);
labels = labels(p);
figure(1)
t = tiledlayout(3, 5, 'TileSpacing', 'compact');
for i = 1:length(labels)
    nexttile
    bar(rattusTablea(i, 1:4)/rattusTablea(i, 1), 'FaceColor', [0.7 0.7 0.7]);
    hold on
    groups = {{1, 2}, {3, 4}, {2, 4}};
    H = sigstar(groups, [pValueRattusa(i, 1),pValueRattusa(i, 2), pValueRattusa(i, 3)]);
    if ((i == 11) || (i == 12)|| (i == 13)|| (i == 14) || (i == 15))
        set(gca,'XTickLabel',{'WT', 'WT+CDDO-Im', 'Nrf2 KO', 'Nrf2 KO+CDDO-Im'})
    else
        set(gca,'XTick',[])
    end  
    swarmchart(1:4, [rattus70f7a(i, 1:3);rattus70f7a(i, 4:6);rattus70f7a(i, 7:9);rattus70f7a(i, 10:12)]/rattusTablea(i, 1), 30,"k", 'LineWidth', 1)
    title(labels(i), 'FontSize', 14, 'FontAngle','italic')
    ylim([0 ((max(rattusTablea(i, 1:4)/rattusTablea(i, 1)))+max(rattusTablea(i, 5:8)/rattusTablea(i, 1)))*1.25])
    %er = errorbar(1:4, rattusTablea(i, 1:4)/rattusTablea(i, 1), zeros(1, 4), rattusTablea(i, 5:8)/rattusTablea(i, 1));
    %er.Color = [0 0 0]; er.LineStyle = 'none';
    set(gcf,'color','w');
    set(H(:,2),'FontSize',15)
    hold off
end
title(t, 'Nrf2 Status Rat CDDO-IM +1', 'FontSize', 20)
ylabel(t, 'Relative expression', 'FontSize', 14)
set(gcf,'color','w');
figure(2)
t = tiledlayout(3, 5, 'TileSpacing', 'compact');
for i = 1:length(labels)
    nexttile
    bar(rattusTableb(i, 1:4)/rattusTableb(i, 1), 'FaceColor', [0.7 0.7 0.7])
    hold on
    groups = {{1, 2}, {3, 4}, {2, 4}};
    H = sigstar(groups, [pValueRattusb(i, 1),pValueRattusb(i, 2), pValueRattusb(i, 3)]);
    if ((i == 11) || (i == 12)|| (i == 13)|| (i == 14) || (i == 15))
        set(gca,'XTickLabel',{'WT', 'WT+CDDO-Im', 'Nrf2 KO', 'Nrf2 KO+CDDO-Im'})
    else
        set(gca,'XTick',[])
    end
    swarmchart(1:4, [rattus70f7b(i, 1:3);rattus70f7b(i, 4:6);rattus70f7b(i, 7:9);rattus70f7b(i, 10:12)]/rattusTableb(i, 1), 30,"k", 'LineWidth', 1)
    title(labels(i), 'FontSize', 14, 'FontAngle','italic')
    ylim([0 ((max(rattusTableb(i, 1:4)/rattusTableb(i, 1)))+max(rattusTableb(i, 5:8)/rattusTableb(i, 1)))*1.3])
    %er = errorbar(1:4, rattusTableb(i, 1:4)/rattusTableb(i, 1), zeros(1, 4), rattusTableb(i, 5:8)/rattusTableb(i, 1));
    %er.Color = [0 0 0]; er.LineStyle = 'none';
    set(gcf,'color','w');
    set(H(:,2),'FontSize',15)
    hold off
end
title(t, 'Nrf2 Status Rat CDDO-IM +7', 'FontSize', 20)
ylabel(t, 'Relative expression', 'FontSize', 18)
set(gcf,'color','w');
clear  t er  i er xlab h p H
%% Rattus 5and6 of 7
rattus70f7a = [[-0.024309635	0.2139926	-0.22748995	0.034966946	-0.3429103	-0.30669165	0.07469416	-0.23674059	0.19831896	-0.07158995	0.2703576	0.024309635];...
[0.32384682	-0.15820503	-0.3271265	2.725583	-0.12853146	0.53150463	-0.24593735	0.049355507	0.00249958	-0.002498627	-0.25111008	0.10221386];...
[-0.2533245	-0.3211851	0.006337166	0.14132214	0.21152306	0.37306595	-0.10473919	-0.04147625	0.08232975	-0.09297848	-0.006337166	0.13135433];...
[-0.041243076	-0.1378665	-0.10901308	0.55721045	0.9105954	1.1078763	-0.052167416	-0.11745691	0.28719568	0.070491314	0.041243076	-0.18495417];...
[-0.17100334	-0.22534943	0.53462315	1.6206541	1.7538776	1.6562433	-0.32063293	-0.2346077	-0.14770794	-0.007837296	0.42085266	0.007837296];...
[0.40741682	0.5156908	0.95534897	2.1366396	2.6372423	2.7409563	-0.55621433	-0.40741682	-1.4950204	-0.75598145	-0.8887744	-0.76053476];...
[0.13582754	0.003647327	-0.019053936	0.8799386	0.9960437	1.1144814	-0.09474993	-0.26600695	-0.003647327	0.04838991	-0.043563366	-0.12425947];...
[0.31585932	0.053005695	0.29734945	2.6723666	3.1343646	3.1517959	-0.48055887	-0.615458	-0.67626524	-0.33456945	-0.11301851	-0.053005695];...
[0.15320301	0.03199196	0.044569016	0.333292	-0.12605476	-0.390563	0.39863014	-0.03199196	0.44451046	-0.03269005	-0.24123001	-0.3260765];...
[0.10362339	0.04309082	0.011972427	0.093447685	0.106770515	0.29817677	-0.25274754	-0.02114296	-0.08005333	-0.15245152	-0.011972427	-0.17060566];...
[0.066058636	-0.04005766	0.018002987	0.22977877	0.24251413	0.13926458	-0.24222231	-0.018002987	-0.16732359	-0.1011796	-0.034392834	0.13076067];...
[0.97262573	0.6354532	0.4768238	-0.059555054	0.036860466	-0.036860466	0.21544743	-0.62705564	0.0673008	-0.616004	-0.8527198	-0.5067706];...
[0.11948347	0.13800573	0.26773787	-0.011759281	-0.33079004	-0.01670599	0.027006626	-0.04985094	0.35677576	-0.08782244	-0.07817888	0.011759281]];
rattus70f7a = 2.^rattus70f7a;
rattus70f7b = [[-0.006602287	-0.25734043	0.4364109	-0.24585438	-0.9584875	-0.49352074	0.21082306	0.20637703	0.006602287	0.15765285	-0.10517979	0.3648739];...
[-0.012938499	-0.6316824	-0.01925373	-0.45675755	0.89490604	0.34914398	0.012937546	-0.7197647	0.18612003	0.78840256	0.17678928	-0.5072756];...
[-0.26650333	-0.02823162	0.10086918	0.21471405	0.34101677	0.27810287	-0.12873554	-0.046105385	-0.16163921	0.008061409	0.009460449	-0.008061409];...
[-0.27798986	0.007277012	-0.06777525	0.66845083	0.42570734	1.0772338	-0.52053213	-0.3116498	0.08707762	-0.007277012	-0.061067104	0.06589079];...
[-0.11428547	-0.07938576	0.2479763	1.2125607	1.280076	1.2519312	-0.19542885	-0.18046188	-0.5861559	0.16526985	0.07938576	-0.54963684];...
[0.116031885	0.54473376	0.56969523	2.7551239	2.1201484	2.2261436	-0.6465781	-0.92950416	-0.27708268	-0.69370675	-0.116031885	-0.9707172];...
[0.064359665	0.06017971	-0.19471359	0.87583447	0.9038048	1.1146641	-0.06017971	-0.5530758	0.14289284	-0.42973328	-0.21482182	-0.07439804];...
[0.042840958	0.5415058	-0.018321037	3.2256117	3.0724487	3.5776186	-0.04246998	-0.33377743	-0.22860336	-0.5075073	0.018321037	-0.37198257];...
[-0.17088795	-0.1979866	-0.04166031	0.015813828	-0.08634758	0.037903786	0.11102009	-0.05678463	0.097268105	0.031736374	-0.015813828	0.026112556];...
[-0.27986765	-0.10424662	-0.11945772	0.35583925	0.031401157	0.2606969	-0.28611517	-0.09249544	-0.008018017	0.008018017	0.13730669	0.05544901];...
[0.096223354	-0.012323856	0.03261423	0.43440294	0.108043194	0.23175478	-0.32266474	-0.14012575	-0.074615955	-0.103168964	-0.13313341	0.012323856];...
[-0.35656166	-0.14224005	-0.29296446	0.019867897	0.59594774	0.001530647	-0.001530647	-0.6786113	0.1618495	-0.31370878	0.15657663	0.09400129];...
[0.032081604	0.09673977	0.26266956	-0.28508282	-0.57311344	-0.21555233	0.022130013	0.07612133	-0.034269333	0.09518051	-0.1754036	-0.022130013]];
rattus70f7b = 2.^rattus70f7b;
labels = {'Akr1c1', 'Trim16', 'Abcb6', 'Slc48a1', 'Mafg', 'Panx2', 'Sqstm1', 'Txnrd1', 'Aifm2', 'Rit1', 'Prdx1', 'Slc6a6', 'Nqo2'};
rattusTablea = [mean(rattus70f7a(:, 1:3), 2), mean(rattus70f7a(:, 4:6), 2), mean(rattus70f7a(:, 7:9), 2), mean(rattus70f7a(:, 10:12), 2),...
    std(rattus70f7a(:, 1:3), 0, 2), std(rattus70f7a(:, 4:6), 0, 2), std(rattus70f7a(:, 7:9), 0, 2), std(rattus70f7a(:, 10:12),0,  2)];
rattusTableb = [mean(rattus70f7b(:, 1:3), 2), mean(rattus70f7b(:, 4:6), 2), mean(rattus70f7b(:, 7:9), 2), mean(rattus70f7b(:, 10:12), 2),...
    std(rattus70f7b(:, 1:3), 0, 2), std(rattus70f7b(:, 4:6), 0, 2), std(rattus70f7b(:, 7:9), 0, 2), std(rattus70f7b(:, 10:12),0,  2)];
pValueRattusa = [mattest(rattus70f7a(:, 1:3), rattus70f7a(:, 4:6), 'VarType', 'equal',  'Labels', labels),...
    mattest(rattus70f7a(:, 7:9), rattus70f7a(:, 10:12), 'VarType', 'equal',  'Labels', labels),...
        mattest(rattus70f7a(:, 4:6), rattus70f7a(:, 10:12), 'VarType', 'equal',  'Labels', labels)];
pValueRattusb = [mattest(rattus70f7b(:, 1:3), rattus70f7b(:, 4:6), 'VarType', 'equal',  'Labels', labels),...
    mattest(rattus70f7b(:, 7:9), rattus70f7b(:, 10:12), 'VarType', 'equal',  'Labels', labels),...
    mattest(rattus70f7b(:, 4:6), rattus70f7b(:, 10:12), 'VarType', 'equal',  'Labels', labels)];
figure(1)
t = tiledlayout(3, 5, 'TileSpacing', 'compact');
% for i = 1:length(labels)
%     nexttile
%     bar(rattusTablea(i, 1:4), 'FaceColor', [0.7 0.7 0.7])    
%     if ((i == 11) || (i == 12)|| (i == 13)|| (i == 14) || (i == 15))
%         set(gca,'XTickLabel',{'WT', 'WT+CDDO-Im', 'NRF2 KO', 'NRF2 KO+CDDO-Im'})
%     else
%         set(gca,'XTick',[])
%     end   
%     swarmchart(1:4, [rattus70f7a(i, 1:3);rattus70f7a(i, 4:6);rattus70f7a(i, 7:9);rattus70f7a(i, 10:12)]/rattusTablea(i, 1), 30,"k", 'LineWidth', 1)
%     groups = {{1, 2}, {3, 4}, {2, 4}};
%     H = sigstar(groups, [pValueRattusa(i, 1),pValueRattusa(i, 2), pValueRattusa(i, 3)]);
%     hold on
%     title(labels(i), 'FontSize', 16, 'FontAngle','italic')
%     ylim([0 ((max(rattusTablea(i, 1:4)))+max(rattusTablea(i, 5:8)))*1.25])
%     er = errorbar(1:4, rattusTablea(i, 1:4), rattusTablea(i, 5:8));
%     er.Color = [0 0 0]; er.LineStyle = 'none';
%     hold off
% end

for i = 1:length(labels)
    nexttile
    bar(rattusTablea(i, 1:4)/rattusTablea(i, 1), 'FaceColor', [0.7 0.7 0.7]);
    hold on
    groups = {{1, 2}, {3, 4}, {2, 4}};
    H = sigstar(groups, [pValueRattusa(i, 1),pValueRattusa(i, 2), pValueRattusa(i, 3)]);
    if ((i == 11) || (i == 12)|| (i == 13)|| (i == 14) || (i == 15))
        set(gca,'XTickLabel',{'WT', 'WT+CDDO-Im', 'Nrf2 KO', 'Nrf2 KO+CDDO-Im'})
    else
        set(gca,'XTick',[])
    end  
    swarmchart(1:4, [rattus70f7a(i, 1:3);rattus70f7a(i, 4:6);rattus70f7a(i, 7:9);rattus70f7a(i, 10:12)]/rattusTablea(i, 1), 30,"k", 'LineWidth', 1)
    title(labels(i), 'FontSize', 14, 'FontAngle','italic')
    ylim([0 ((max(rattusTablea(i, 1:4)/rattusTablea(i, 1)))+max(rattusTablea(i, 5:8)/rattusTablea(i, 1)))*1.25])
    %er = errorbar(1:4, rattusTablea(i, 1:4)/rattusTablea(i, 1), zeros(1, 4), rattusTablea(i, 5:8)/rattusTablea(i, 1));
    %er.Color = [0 0 0]; er.LineStyle = 'none';
    set(gcf,'color','w');
    set(H(:,2),'FontSize',15)
    hold off
end
title(t, 'Rat CDDO-IM Nrf2 +1', 'FontSize', 24)
ylabel(t, 'Relative expression', 'FontSize', 18)
set(gcf,'color','w');
figure(2)
t = tiledlayout(3,5, 'TileSpacing', 'compact');
% for i = 1:length(labels)
%     nexttile
%     bar(rattusTableb(i, 1:4), 'FaceColor', [0.7 0.7 0.7])
%     if ((i == 11) || (i == 12)|| (i == 13)|| (i == 14) || (i == 15))
%         set(gca,'XTickLabel',{'WT', 'WT+CDDO-Im', 'NRF2 KO', 'NRF2 KO+CDDO-Im'})
%     else
%         set(gca,'XTick',[])
%     end   
%     hold on
%     groups = {{1, 2}, {3, 4}, {2, 4}};
%     swarmchart(1:4, [rattus70f7b(i, 1:3);rattus70f7b(i, 4:6);rattus70f7b(i, 7:9);rattus70f7b(i, 10:12)]/rattusTableb(i, 1), 30,"k", 'LineWidth', 1)
%     H = sigstar(groups, [pValueRattusb(i, 1),pValueRattusb(i, 2),pValueRattusb(i, 3)]);
%     title(labels(i), 'FontSize', 16, 'FontAngle','italic')
%     ylim([0 ((max(rattusTableb(i, 1:4)))+max(rattusTableb(i, 5:8)))*1.25])
%     er = errorbar(1:4, rattusTableb(i, 1:4), rattusTableb(i, 5:8));
%     er.Color = [0 0 0]; er.LineStyle = 'none';
%     hold off
% end
for i = 1:length(labels)
    nexttile
    bar(rattusTableb(i, 1:4)/rattusTableb(i, 1), 'FaceColor', [0.7 0.7 0.7])
    hold on
    groups = {{1, 2}, {3, 4}, {2, 4}};
    H = sigstar(groups, [pValueRattusb(i, 1),pValueRattusb(i, 2), pValueRattusb(i, 3)]);
    if ((i == 11) || (i == 12)|| (i == 13)|| (i == 14) || (i == 15))
        set(gca,'XTickLabel',{'WT', 'WT+CDDO-Im', 'Nrf2 KO', 'Nrf2 KO+CDDO-Im'})
    else
        set(gca,'XTick',[])
    end
    swarmchart(1:4, [rattus70f7b(i, 1:3);rattus70f7b(i, 4:6);rattus70f7b(i, 7:9);rattus70f7b(i, 10:12)]/rattusTableb(i, 1), 30,"k", 'LineWidth', 1)
    title(labels(i), 'FontSize', 14, 'FontAngle','italic')
    ylim([0 ((max(rattusTableb(i, 1:4)/rattusTableb(i, 1)))+max(rattusTableb(i, 5:8)/rattusTableb(i, 1)))*1.3])
    %er = errorbar(1:4, rattusTableb(i, 1:4)/rattusTableb(i, 1), zeros(1, 4), rattusTableb(i, 5:8)/rattusTableb(i, 1));
    %er.Color = [0 0 0]; er.LineStyle = 'none';
    set(gcf,'color','w');
    set(H(:,2),'FontSize',15)
    hold off
end
title(t, 'Rat CDDO-IM Nrf2 +7', 'FontSize', 24)
ylabel(t, 'Relative expression', 'FontSize', 18)
set(gcf,'color','w');
clear rattusTablea rattusTableb t er  i er xlab
%% SNF 24Hrs
    excelData = '/Users/harshitakumar/Documents/Research/Letterio Research/Letterio Data/RNA Seq/GSE48812_24h_LNCAP_SFN.xlsx';
    [~, ~, dataFile1] = xlsread(excelData);
    dataLNSNF24 = cell2mat(dataFile1(2:end, [5 3]));
    labelLNSNF24 = cell2table(dataFile1(2:end, 1));
    upRegLNSNF24 = [];
    SNFLNVal = [];
    for i = 1:length(ylabelUniversal)
        gene = ylabelUniversal(i);
        location = find(string(labelLNSNF24.(1)) == string(gene));
        if (isempty(location) == 0)
            if ((dataLNSNF24(location, 1) > 1) && (dataLNSNF24(location, 2) < 0.05))
                upRegLNSNF24 = [upRegLNSNF24; labelLNSNF24.(1)(location), dataLNSNF24(location, 1), dataLNSNF24(location, 2)];
                if (dataLNSNF24(location, 1) >= 2)
                    SNFLNVal = [SNFLNVal; [1]];
                else
                    SNFLNVal = [SNFLNVal; [0]];
                end
            else
                SNFPC3Val = [SNFPC3Val; [-1]];
            end
        else
            SNFLNVal = [SNFLNVal; [-1]];
        end
    end
   
    excelData = '/Users/harshitakumar/Documents/Research/Letterio Research/Letterio Data/RNA Seq/GSE48812_24h_PC3_SFN.xlsx';
    [~, ~, dataFile1] = xlsread(excelData);
    dataPC3SNF24 = cell2mat(dataFile1(2:end, [5 3]));
    labelPC3SNF24 = cell2table(dataFile1(2:end, 1));
    upRegPC3SNF24 = [];
    SNFPC3Val = [];
    for i = 1:length(ylabelUniversal)
        gene = ylabelUniversal(i);
        location = find(string(labelPC3SNF24.(1)) == string(gene));
        if (isempty(location) == 0)
            if ((dataPC3SNF24(location, 1) > 1) && (dataPC3SNF24(location, 2) < 0.05))
                upRegPC3SNF24 = [upRegPC3SNF24; labelPC3SNF24.(1)(location), dataPC3SNF24(location, 1), dataPC3SNF24(location, 2)];
                if (dataPC3SNF24(location, 1) >= 2)
                    SNFPC3Val = [SNFPC3Val; [1]];
                else
                    SNFPC3Val = [SNFPC3Val; [0]];
                end
            else
                SNFPC3Val = [SNFPC3Val; [-1]];
            end
        else
            SNFPC3Val = [SNFPC3Val; [-1]];
        end
    end
    excelData = '/Users/harshitakumar/Documents/Research/Letterio Research/Letterio Data/RNA Seq/GSE48812_24h_PREC_SFN.xlsx';
    [~, ~, dataFile1] = xlsread(excelData);
    dataPRECSNF24 = cell2mat(dataFile1(2:end, [2 3]));
    dataPRECSNF24(:, 1) = 2.^dataPRECSNF24(:, 1);
    labelPRECSNF24 = cell2table(dataFile1(2:end, 1));
    upRegPRECSNF24 = [];
    SNFPRECVal = [];
    for i = 1:length(ylabelUniversal)
        gene = ylabelUniversal(i);
        location = find(string(labelPRECSNF24.(1)) == string(gene));
        if (isempty(location) == 0)
            if ((dataPRECSNF24(location, 1) > 1) && (dataPRECSNF24(location, 2) < 0.05))
                upRegPRECSNF24 = [upRegPRECSNF24; labelPRECSNF24.(1)(location), dataPRECSNF24(location, 1), dataPRECSNF24(location, 2)];
                if (dataPRECSNF24(location, 1) >= 2)
                    SNFPRECVal = [SNFPRECVal; [1]];
                else
                    SNFPRECVal = [SNFPRECVal; [0]];
                end
            end
        else
            SNFPRECVal = [SNFPRECVal; [-1]];
        end
    end
    %HEATMAP
        ylab = {'\it ABHD4'	'\it AKR1C3'	'\it ~CBR3'	'\it EPHX1'	'\it FTH1'	'\it FTL'	'\it GCLC'	'\it GCLM'	'\it GSR'	'\it ME1'	'\it NQO1'	'\it OSGIN1'	'\it PIR'	'\it SLC7A11'	'\it SRXN1'	'\it ABCB6'	'\it AIFM2'	'\it AKR1C1'	'\it GLA'	'\it MAFG'	'\it NQO2'	'\it PANX2'	'\it PGD'	'\it PRDX1'	'\it RIT1'	'\it SLC48A1'	'\it SLC6A6'	'\it SQSTM1'	'\it TRIM16'	'\it TRIM16L'	'\it TXNRD1'};
        figure(1)
        xlabels = char({'LNCaP', 'PC-3', 'PREC'});
        heatmap(xlabels, ylab, double([SNFLNVal, SNFPC3Val, SNFPRECVal]), 'Title', 'Validation of NRF2 Genes with Sulforaphane', 'FontSize', 20, 'CellLabelColor', 'none', 'ColorbarVisible', 'on');
        cmap = colormap([ 0.7, 0.7, 0.7; 1, 0.8, 0.8; 1, 0, 0]);
        set(gcf,'color','w');
        axs = struct(gca); %ignore warning that this should be avoided
        colorbar('Location','south')
        c = axs.Colorbar;
        c.Ticks = [-1; -0.9; 0; 1];
        c.TickLabelsMode = 'manual';
        c.TickLabels = {['Unchanged', 10, 'Not expressed/ '],'Upregulated', 'Highly Upregulated'};
        c.FontSize = 12;
        c.location = 'south';
        set(gcf,'color','w');
 clear excelData i location dataFile1 
%% Revised SNF 24 Hrs Validation
excelData = '/Users/harshitakumar/Documents/Research/Letterio Research/Letterio Data/RNA Seq/quant/LNCAP 24 SFN.xlsx';
[~, ~, dataFile1] = xlsread(excelData);
dataLNSNF24 = cell2mat(dataFile1(2:end, [4:9]));
labelLNSNF24 = cell2table(dataFile1(2:end, 10));
excelData = '/Users/harshitakumar/Documents/Research/Letterio Research/Letterio Data/RNA Seq/quant/PC3 24 SFN.xlsx';
[~, ~, dataFile1] = xlsread(excelData);
dataPC3SNF24 = cell2mat(dataFile1(2:end, [4:9]));
labelPC3SNF24 = cell2table(dataFile1(2:end, 10));
excelData = '/Users/harshitakumar/Documents/Research/Letterio Research/Letterio Data/RNA Seq/quant/Prec 24 SFN.xlsx';
[~, ~, dataFile1] = xlsread(excelData);
dataPRECSNF24 = cell2mat(dataFile1(2:end, [4:9]));
labelPRECSNF24 = cell2table(dataFile1(2:end, 10));

%% SNF 6Hrs
    excelData = '/Users/harshitakumar/Documents/Research/Letterio Research/Letterio Data/RNA Seq/GSE48812_6h_LNCAP_SFN.xlsx';
    [~, ~, dataFile1] = xlsread(excelData);
    dataLNSNF24 = cell2mat(dataFile1(2:end, [3 2]));
    labelLNSNF24 = cell2table(dataFile1(2:end, 1));
    upRegLNSNF24 = [];
    SNFLNVal = [];
    for i = 1:length(ylabelUniversal)
        gene = ylabelUniversal(i);
        location = find(string(labelLNSNF24.(1)) == string(gene));
        if (isempty(location) == 0)
            if ((dataLNSNF24(location, 1) > 1) && (dataLNSNF24(location, 2) < 0.05))
                upRegLNSNF24 = [upRegLNSNF24; labelLNSNF24.(1)(location), dataLNSNF24(location, 1), dataLNSNF24(location, 2)];
                if (dataLNSNF24(location, 1) >= 2)
                    SNFLNVal = [SNFLNVal; [1]];
                else
                    SNFLNVal = [SNFLNVal; [0]];
                end
            end
        else
            SNFLNVal = [SNFLNVal; [-1]];
        end
    end
   
    excelData = '/Users/harshitakumar/Documents/Research/Letterio Research/Letterio Data/RNA Seq/GSE48812_6h_PC3_SFN.xlsx';
    [~, ~, dataFile1] = xlsread(excelData);
    dataPC3SNF24 = cell2mat(dataFile1(2:end, [3 2]));
    labelPC3SNF24 = cell2table(dataFile1(2:end, 1));
    upRegPC3SNF24 = [];
    SNFPC3Val = [];
    for i = 1:length(ylabelUniversal)
        gene = ylabelUniversal(i);
        location = find(string(labelPC3SNF24.(1)) == string(gene));
        if (isempty(location) == 0)
            if ((dataPC3SNF24(location, 1) > 1) && (dataPC3SNF24(location, 2) < 0.05))
                upRegPC3SNF24 = [upRegPC3SNF24; labelPC3SNF24.(1)(location), dataPC3SNF24(location, 1), dataPC3SNF24(location, 2)];
                if (dataPC3SNF24(location, 1) >= 2)
                    SNFPC3Val = [SNFPC3Val; [1]];
                else
                    SNFPC3Val = [SNFPC3Val; [0]];
                end
            else
                SNFPC3Val = [SNFPC3Val; [-1]]
            end
        else
            SNFPC3Val = [SNFPC3Val; [-1]]
        end
    end
    excelData = '/Users/harshitakumar/Documents/Research/Letterio Research/Letterio Data/RNA Seq/GSE48812_6h_PREC_SFN.xlsx';
    [~, ~, dataFile1] = xlsread(excelData);
    dataPRECSNF24 = cell2mat(dataFile1(2:end, [3 2]));
    labelPRECSNF24 = cell2table(dataFile1(2:end, 1));
    upRegPRECSNF24 = [];
    SNFPRECVal = [];
    for i = 1:length(ylabelUniversal)
        gene = ylabelUniversal(i);
        location = find(string(labelPRECSNF24.(1)) == string(gene));
        if (isempty(location) == 0)
            if ((dataPRECSNF24(location, 1) > 1) && (dataPRECSNF24(location, 2) < 0.05))
                upRegPRECSNF24 = [upRegPRECSNF24; labelPRECSNF24.(1)(location), dataPRECSNF24(location, 1), dataPRECSNF24(location, 2)];
                if (dataPRECSNF24(location, 1) >= 2)
                    SNFPRECVal = [SNFPRECVal; [1]];
                else
                    SNFPRECVal = [SNFPRECVal; [0]];
                end
            else
                SNFPRECVal = [SNFPRECVal; [-1]];
            end
        else
            SNFPRECVal = [SNFPRECVal; [-1]];
        end
    end
    %HEATMAP
        figure(1)
        xlabels = char({'LNCaP', 'PC-3', 'PREC'});
        ylab = {'\it ABHD4'	'\it AKR1C3'	'\it ~CBR3'	'\it EPHX1'	'\it FTH1'	'\it FTL'	'\it GCLC'	'\it GCLM'	'\it GSR'	'\it ME1'	'\it NQO1'	'\it OSGIN1'	'\it PIR'	'\it SLC7A11'	'\it SRXN1'	'\it ABCB6'	'\it AIFM2'	'\it AKR1C1'	'\it GLA'	'\it MAFG'	'\it NQO2'	'\it PANX2'	'\it PGD'	'\it PRDX1'	'\it RIT1'	'\it SLC48A1'	'\it SLC6A6'	'\it SQSTM1'	'\it TRIM16'	'\it TRIM16L'	'\it TXNRD1'};
        heatmap(xlabels, ylab, double([SNFLNVal, SNFPC3Val, SNFPRECVal]), 'Title', 'Validation of NRF2 Genes with Sulforaphane at 6 hours', 'FontSize', 20, 'CellLabelColor', 'none', 'ColorbarVisible', 'on');
        cmap = colormap([ 0.7, 0.7, 0.7; 1, 0.8, 0.8; 1, 0, 0]);
        axs = struct(gca); %ignore warning that this should be avoided
        c = axs.Colorbar;
        c.Ticks = [-1; -0.9; 0; 1];
        c.TickLabelsMode = 'manual';
        c.TickLabels = {['Unchanged', 10, 'Not expressed/ '],'Upregulated', 'Highly Upregulated'};
        c.FontSize = 12;
        set(gcf,'color','w');
 clear excelData i location dataFile1 

%% Functions 
function varargout=sigstar(groups,stats,nosort)
    % sigstar - Add significance stars to bar charts, boxplots, line charts, etc,
    %
    % H = sigstar(groups,stats,nsort)
    %
    % Purpose
    % Add stars and lines highlighting significant differences between pairs of groups. 
    % The user specifies the groups and associated p-values. The function handles much of 
    % the placement and drawing of the highlighting. Stars are drawn according to:
    %   * represents p<=0.05
    %  ** represents p<=1E-2
    % *** represents p<=1E-3
    %
    %
    % Inputs
    % groups - a cell array defining the pairs of groups to compare. Groups defined 
    %          either as pairs of scalars indicating locations along the X axis or as 
    %          strings corresponding to X-tick labels. Groups can be a mixture of both 
    %          definition types.
    % stats -  a vector of p-values the same length as groups. If empty or missing it's 
    %          assumed to be a vector of 0.05s the same length as groups. Nans are treated
    %          as indicating non-significance.
    % nsort -  optional, 0 by default. If 1, then significance markers are plotted in 
    %          the order found in groups. If 0, then they're sorted by the length of the 
    %          bar.
    %
    % Outputs
    % H - optionally return handles for significance highlights. Each row is a different
    %     highlight bar. The first column is the line. The second column is the text (stars).
    %     
    %
    % Examples
    % 1. 
    % bar([5,2,1.5])
    % sigstar({[1,2], [1,3]})
    %
    % 2. 
    % bar([5,2,1.5])
    % sigstar({[2,3],[1,2], [1,3]},[nan,0.05,0.05])
    %
    % 3.  **DOESN'T WORK IN 2014b**
    % R=randn(30,2);
    % R(:,1)=R(:,1)+3;
    % boxplot(R)
    % set(gca,'XTick',1:2,'XTickLabel',{'A','B'})
    % H=sigstar({{'A','B'}},0.01);
    % ylim([-3,6.5])
    % set(H,'color','r')
    %
    % 4. Note the difference in the order with which we define the groups in the 
    %    following two cases. 
    % x=[1,2,3,2,1];
    % subplot(1,2,1)
    % bar(x)
    % sigstar({[1,2], [2,3], [4,5]})
    % subplot(1,2,2)
    % bar(x)
    % sigstar({[2,3],[1,2], [4,5]})
    %
    % ALSO SEE: demo_sigstar
    %
    % KNOWN ISSUES:
    % 1. Algorithm for identifying whether significance bar will overlap with 
    %    existing plot elements may not work in some cases (see line 277)
    % 2. Bars may not look good on exported graphics with small page sizes.
    %    Simply increasing the width and height of the graph with the 
    %    PaperPosition property of the current figure should fix things.
    %
    % Rob Campbell - CSHL 2013
    %Input argument error checking
    %If the user entered just one group pair and forgot to wrap it in a cell array 
    %then we'll go easy on them and wrap it here rather then generate an error
    if ~iscell(groups) & length(groups)==2
        groups={groups};
    end
    if nargin<2 
        stats=repmat(0.05,1,length(groups));
    end
    if isempty(stats)
        stats=repmat(0.05,1,length(groups));
    end
    if nargin<3
        nosort=0;
    end
    %Check the inputs are of the right sort
    if ~iscell(groups)
        error('groups must be a cell array')
    end
    if ~isvector(stats)
        error('stats must be a vector')
    end
    if length(stats)~=length(groups)
        error('groups and stats must be the same length')
    end
    %Each member of the cell array groups may be one of three things:
    %1. A pair of indices.
    %2. A pair of strings (in cell array) referring to X-Tick labels
    %3. A cell array containing one index and one string
    %
    % For our function to run, we will need to convert all of these into pairs of
    % indices. Here we loop through groups and do this. 
    xlocs=nan(length(groups),2); %matrix that will store the indices 
    xtl=get(gca,'XTickLabel');  
    for ii=1:length(groups)
        grp=groups{ii};
        if isnumeric(grp)
            xlocs(ii,:)=grp; %Just store the indices if they're the right format already
        elseif iscell(grp) %Handle string pairs or string/index pairs
            if isstr(grp{1})
                a=strmatch(grp{1},xtl);
            elseif isnumeric(grp{1})
                a=grp{1};
            end
            if isstr(grp{2})
                b=strmatch(grp{2},xtl);
            elseif isnumeric(grp{2})
                b=grp{2};
            end
            xlocs(ii,:)=[a,b];
        end
        %Ensure that the first column is always smaller number than the second
        xlocs(ii,:)=sort(xlocs(ii,:));
    end
    %If there are any NaNs we have messed up. 
    if any(isnan(xlocs(:)))
        error('Some groups were not found')
    end
    %Optionally sort sig bars from shortest to longest so we plot the shorter ones first
    %in the loop below. Usually this will result in the neatest plot. If we waned to 
    %optimise the order the sig bars are plotted to produce the neatest plot, then this 
    %is where we'd do it. Not really worth the effort, though, as few plots are complicated
    %enough to need this and the user can define the order very easily at the command line. 
    if ~nosort
        [~,ind]=sort(xlocs(:,2)-xlocs(:,1),'ascend');
        xlocs=xlocs(ind,:);groups=groups(ind);
        stats=stats(ind);
    end
    %-----------------------------------------------------
    %Add the sig bar lines and asterisks 
    holdstate=ishold;
    hold on
    H=ones(length(groups),2); %The handles will be stored here
    y=ylim;
    yd=myRange(y)*0.05; %separate sig bars vertically by 5% 
    for ii=1:length(groups)
        thisY=findMinY(xlocs(ii,:))+yd;
        H(ii,:)=makeSignificanceBar(xlocs(ii,:),thisY,stats(ii));
    end
    %-----------------------------------------------------
    %Now we can add the little downward ticks on the ends of each line. We are
    %being extra cautious and leaving this it to the end just in case the y limits
    %of the graph have changed as we add the highlights. The ticks are set as a
    %proportion of the y axis range and we want them all to be the same the same
    %for all bars.
    yd=myRange(ylim)*0.01; %Ticks are 1% of the y axis range
    for ii=1:length(groups)
        y=get(H(ii,1),'YData');
        y(1)=y(1)-yd;
        y(4)=y(4)-yd;   
        set(H(ii,1),'YData',y*1.2)
    end
    %Be neat and return hold state to whatever it was before we started
    if ~holdstate
        hold off
    elseif holdstate
        hold on
    end
    %Optionally return the handles to the plotted significance bars (first column of H)
    %and asterisks (second column of H).
    if nargout>0
        varargout{1}=H;
    end
end %close sigstar
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Internal functions
function H=makeSignificanceBar(x,y,p)
    %makeSignificanceBar produces the bar and defines how many asterisks we get for a 
    %given p-value
    if p<=1E-3
        stars='***'; 
    elseif p<=1E-2
        stars='**';
    elseif p<=0.05
        stars='*';
    elseif isnan(p)
        stars='n.s.';
    else
        stars='^{NS}';
    end
            
    x=repmat(x,2,1);
    y=repmat(y,4,1);
    H(1)=plot(x(:),y,'-k','LineWidth',1.5,'Tag','sigstar_bar');
    %Increase offset between line and text if we will print "n.s."
    %instead of a star. 
    if ~isnan(p)
        offset=0.005;
    else
        offset=0.02;
    end
    starY=mean(y)+myRange(ylim)*offset;
    H(2)=text(mean(x(:)),starY*1.2,stars,...
        'HorizontalAlignment','Center',...
        'BackGroundColor','none',...
        'Tag','sigstar_stars');
    Y=ylim;
    if Y(2)<starY
        ylim([Y(1),starY+myRange(Y)*0.05])
    end
end %close makeSignificanceBar
function Y=findMinY(x)
    % The significance bar needs to be plotted a reasonable distance above all the data points
    % found over a particular range of X values. So we need to find these data and calculat the 
    % the minimum y value needed to clear all the plotted data present over this given range of 
    % x values. 
    %
    % This version of the function is a fix from Evan Remington
    oldXLim = get(gca,'XLim');
    oldYLim = get(gca,'YLim');
    axis(gca,'tight')
    
    %increase range of x values by 0.1 to ensure correct y max is used
    x(1)=x(1)-0.1;
    x(2)=x(2)+0.1;
    
    set(gca,'xlim',x) %Matlab automatically re-tightens y-axis
    yLim = get(gca,'YLim'); %Now have max y value of all elements within range.
    Y = max(yLim);
    axis(gca,'normal')
    set(gca,'XLim',oldXLim,'YLim',oldYLim)
end %close findMinY
function rng=myRange(x)
    %replacement for stats toolbox range function
    rng = max(x) - min(x);
end %close myRange
%function to convert negative fold changes to ratio of fold changes between
%0 and 1
%     function y =  negToAbs(x)
%         y = x;
%         for i = 1:length(x)
%             if (x(i) < 0)
%             y(i) = 1/abs(x(i));
%             end 
%         end
%     end 
