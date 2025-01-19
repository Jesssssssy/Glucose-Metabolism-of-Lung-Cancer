%% AIM:  Build Dynamic networks based on Ki mean

addpath(genpath('/Volumes/Extreme Pro/糖代谢/code/circularGraph'));
addpath(genpath('/Volumes/Extreme Pro/糖代谢/code'));

patinfo = readtable('/Volumes/Extreme Pro/糖代谢/其他/info/patinfo.xlsx');

HCPlabels = logical(table2array(patinfo(1:110,8))); % Post-covid control 
ADClabels = logical(table2array(patinfo(1:110,5))); % adenocarcinoma
SCClabels = logical(table2array(patinfo(1:110,6))); % squamous cell carcinoma
% NoduleLabls = logical(table2array(patinfo(1:110,7))); % benign nodule
% SmallCellLabls = logical(table2array(patinfo(1:110,9))); % small cell 

% readin patnames_1 for patinfo table
in_dir = dir('/Volumes/Extreme Pro/糖代谢/3DKiDistr_body');
patnames_1 = {in_dir(1:end).name}';
patnames_1 = patnames_1(~startsWith(patnames_1, '.'));
destination_folder = '/Volumes/Extreme Pro/糖代谢/3DKiDistr_body_brain_merged';

%% Combine [body and brain] Ki distribution data
in_dir = dir('/Volumes/Extreme Pro/糖代谢/3DKiDistr_body');
patnames_1 = {in_dir(1:end).name}';
patnames_1 = patnames_1(~startsWith(patnames_1, '.'));
in_dir_2 = dir('/Volumes/Extreme Pro/糖代谢/3DKiDistr_brain');
patnames_2 = {in_dir_2(1:end).name}';
patnames_2 = patnames_2(~startsWith(patnames_2, '.'));

% Find matching identifiers
matching_files = {};
for i = 1:length(patnames_1)
    if ismember(patnames_1{i}, patnames_2)
        matching_files{end+1} = patnames_1{i}; 
    end
end

% Define the destination folder
destination_folder = '/Volumes/Extreme Pro/糖代谢/3DKiDistr_body_brain_merged';
if ~exist(destination_folder, 'dir')
    mkdir(destination_folder);
end

% Copy matched files to the destination folder
for i = 1:length(matching_files)
    source_path_body = fullfile(in_dir(3).folder, matching_files{i});
    source_path_brain = fullfile(in_dir_2(3).folder, matching_files{i});
    destination_path = fullfile(destination_folder, matching_files{i});
    body = load(source_path_body);
    brain = load(source_path_brain);
    featcell_comB = [body.featcell, brain.featcell];
    save(destination_path, 'featcell_comB');

end

%% Selection of ROIs
load('/Volumes/Extreme Pro/糖代谢/3DKiDistr_body_brain_merged/024_CHANG_LIANG.mat');
ROInames = featcell_comB(1,:);

ROInames{1} = 'L-Adrenal';
ROInames{2} = 'R-Adrenal';
ROInames{10} = 'Thyroid';

% ROInames{11} = 'Heart';

% ROInames{22} = 'Rib';
% ROInames{31} = 'Vertebra';
% ROInames{5} = 'Pelvics';

ROInames{41} = 'Muscle';    

ROInames{14} = 'L-Kidney';
ROInames{15} = 'R-Kidney';
ROInames{16} = 'Liver';

ROInames{17} = 'L-Lung';
ROInames{18} = 'R-Lung';
ROInames{19} = 'Pancreas';
ROInames{27} = 'Spleen';
ROInames{26} = 'SpinalC';

ROInames{44} = 'Brainstem';
ROInames{45} = 'L-Hipp';
ROInames{46} = "R-Hipp";
% ROInames{47} = 'L-Temporal gyrus (mi)'; % middle and inferior temporal gyrus
% ROInames{48} = 'R-Temporal gyrus (mi)'; 
ROInames{49} = 'L-PC'; %Cingulate gyrus cinguli posterior part
ROInames{50} = 'R-PC';
ROInames{51} = 'L-MF'; %middle frontal gyrus
ROInames{52} = 'R-MF';
ROInames{53} = 'L-MAT'; %Anterior temporal lobe medial part
ROInames{54} = 'R-MAT';
ROInames{55} = 'L-LAT'; %Anterior temporal lobe lateral part
ROInames{56} = 'R-LAT';
ROInames{57} = 'L-Cerebe';
ROInames{58} = 'R-Cerebe';

selectedIndices = [1,2,10,41,14,15,16,19,17,18,27,26,44,45,46,49,50,51,52,53,54,55,56,57,58];

ROInames_selected = ROInames(selectedIndices); 

num_roi = length(ROInames_selected); % 25 ROI selected
alpha = 0.05; % Original significance level
num_tests = nchoosek(num_roi, 2); % Number of pairwise correlations
corrected_alpha = alpha / num_tests;

%% SECTION 1.1: Calculate Group RefNet, ADC_Net, SCC_Net (using spearman)
close all

HCPnames = patnames_1(HCPlabels); % 新冠康复?
[refNet,ref_p] = getRmean(destination_folder, HCPnames, num_roi,selectedIndices);
% corrplotvis_jessy(refNet,ROInames_selected); 
% refNet_ZCC = zscore(refNet); % normalized network 
% idx_ = [find(refNet_ZCC > 1.96);find(refNet_ZCC < -1.96)];
idx_ = find(ref_p < corrected_alpha); 
adjMat_NonS = zeros(size(refNet));
adjMat_NonS(idx_) = 1;
corrMat_NonS = adjMat_NonS.*refNet;
A_covid_spear = corrMat_NonS;
symmetricMatrix_refnet = A_covid_spear + A_covid_spear.' - diag(diag(A_covid_spear)); % normalized & filtered network 

% R_upper = A_NonS;
% [row, col] = find(R_upper ~= 0);
% values = R_upper(R_upper ~= 0);
% figure;
% scatter3(row, col, values, 50, values, 'filled'); % 3D scatter plot
% colorbar;
% colormap(parula);
% set(gca, 'XTick', 1:length(ROInames_selected), 'XTickLabel', ROInames_selected, ...
%          'YTick', 1:length(ROInames_selected), 'YTickLabel', ROInames_selected);
% xtickangle(45);
% xlabel('Feature Names (Rows)');
% ylabel('Feature Names (Columns)');
% zlabel('Correlation Coefficients');
% title('Scatter Plot of Upper Triangle Correlation Coefficients');

% ADC tumor
ADCnames = patnames_1(ADClabels); % 共 ? 腺癌
[ADC_Net,ADC_p] = getRmean(destination_folder, ADCnames, num_roi,selectedIndices);
% corrplotvis_jessy(ADC_Net,ROInames_selected);
% ADC_Net_ZCC = zscore(ADC_Net); % normalized network 
% idx_adc = [find(ADC_Net_ZCC > 1.96);find(ADC_Net_ZCC < -1.96)];
idx_adc = find(ADC_p < corrected_alpha);
adjMat_NonS = zeros(size(ADC_Net));
adjMat_NonS(idx_adc) = 1;
corrMat_NonS = adjMat_NonS.*ADC_Net;
A_ADC_spear = corrMat_NonS;
symmetricMatrix_adcnet = A_ADC_spear + A_ADC_spear.' - diag(diag(A_ADC_spear)); % normalized & filtered network 


% SCC tumor
SCCnames = patnames_1(SCClabels); % 共11个鳞状细胞癌
[SCC_Net,SCC_p] = getRmean(destination_folder, SCCnames, num_roi,selectedIndices);
% corrplotvis_jessy(SCC_Net,ROInames_selected);
% SCC_Net_ZCC = zscore(SCC_Net); % normalized network 
% idx_scc = [find(refNet_ZCC > 1.96);find(refNet_ZCC < -1.96)];
idx_scc = find(SCC_p < corrected_alpha);
adjMat_NonS = zeros(size(SCC_Net));
adjMat_NonS(idx_scc) = 1;
corrMat_NonS = adjMat_NonS.*SCC_Net;
A_SCC_spear = corrMat_NonS;
symmetricMatrix_sccnet = A_SCC_spear + A_SCC_spear.' - diag(diag(A_SCC_spear)); % normalized & filtered network 

% % Benign Nodule % 共 8 个良性肺结节
% Nodulenames = patnames_1(NoduleLabls); % 共?个nodule
% [Nodule_Net,Nodule_p] = getRmean(destination_folder, Nodulenames, num_roi,selectedIndices);
% % corrplotvis_jessy(Nodule_Net,ROInames_selected);
% % Nodule_Net_ZCC = zscore(Nodule_Net); % normalized network 
% % idx_nol = [find(Nodule_Net_ZCC > 1.96);find(Nodule_Net_ZCC < -1.96)];
% idx_nol = find(Nodule_p < corrected_alpha);
% adjMat_NonS = zeros(size(Nodule_Net));
% adjMat_NonS(idx_nol) = 1;
% corrMat_NonS = adjMat_NonS.*Nodule_Net;
% A_Nodule_spear = corrMat_NonS;
% symmetricMatrix_nolnet = A_Nodule_spear + A_Nodule_spear.' - diag(diag(A_Nodule_spear)); % normalized & filtered network 
% 
% % small cell 
% SmallCellnames = patnames_1(SmallCellLabls); % 共11个small cell
% [SmallCell_Net, SmallCell_p] = getRmean(destination_folder, SmallCellnames, num_roi,selectedIndices);
% % corrplotvis_jessy(SmallCell_Net,ROInames_selected);
% % SmallCell_Net_ZCC = zscore(SmallCell_Net); % normalized network 
% % idx_SmallCell = [find(SmallCell_Net_ZCC > 1.96);find(SmallCell_Net_ZCC < -1.96)];
% idx_SmallCell = find(SmallCell_p < corrected_alpha);
% adjMat_NonS = zeros(size(SmallCell_Net));
% adjMat_NonS(idx_SmallCell) = 1;
% corrMat_NonS = adjMat_NonS.*SmallCell_Net;
% A_SCLC_spear = corrMat_NonS;
% symmetricMatrix_SmallCell = A_SCLC_spear + A_SCLC_spear.' - diag(diag(A_SCLC_spear)); % normalized & filtered network 

% Figure Plotting
figure,circularGraph(symmetricMatrix_refnet,symmetricMatrix_refnet,'Label',ROInames_selected); % covid net
figure,circularGraph(symmetricMatrix_adcnet,symmetricMatrix_adcnet,'Label',ROInames_selected); % adc net
figure,circularGraph(symmetricMatrix_sccnet,symmetricMatrix_sccnet,'Label',ROInames_selected); % scc net
% figure,circularGraph(symmetricMatrix_nolnet,symmetricMatrix_nolnet,'Label',ROInames_selected); % nodule net
% figure,circularGraph(symmetricMatrix_SmallCell,symmetricMatrix_SmallCell,'Label',ROInames_selected); % Small cell net

%% SECTION 1.2: Calculate Group RefNet, ADC_Net, SCC_Net (using mutual information)
close all

HCPnames = patnames_1(HCPlabels); % 新冠康复?
[refNet,ref_p] = getRmean_mutual(destination_folder, HCPnames, num_roi, selectedIndices);
% corrplotvis_jessy(refNet,ROInames_selected); 
% refNet_ZCC = zscore(refNet); % normalized network 
% idx_ = [find(refNet_ZCC > 1.96);find(refNet_ZCC < -1.96)];
idx_ = find(ref_p < corrected_alpha); 
adjMat_NonS = zeros(size(refNet));
adjMat_NonS(idx_) = 1;
corrMat_NonS = adjMat_NonS.*refNet;
A_Refnet = corrMat_NonS;
symmetricMatrix_refnet = A_Refnet + A_Refnet.' - diag(diag(A_Refnet)); % normalized & filtered network 

% R_upper = A_NonS;
% [row, col] = find(R_upper ~= 0);
% values = R_upper(R_upper ~= 0);
% figure;
% scatter3(row, col, values, 50, values, 'filled'); % 3D scatter plot
% colorbar;
% colormap(parula);
% set(gca, 'XTick', 1:length(ROInames_selected), 'XTickLabel', ROInames_selected, ...
%          'YTick', 1:length(ROInames_selected), 'YTickLabel', ROInames_selected);
% xtickangle(45);
% xlabel('Feature Names (Rows)');
% ylabel('Feature Names (Columns)');
% zlabel('Correlation Coefficients');
% title('Scatter Plot of Upper Triangle Correlation Coefficients');

% ADC tumor
ADCnames = patnames_1(ADClabels); % 共 ? 腺癌
[ADC_Net,ADC_p] = getRmean_mutual(destination_folder, ADCnames, num_roi, selectedIndices);
% corrplotvis_jessy(ADC_Net,ROInames_selected);
% ADC_Net_ZCC = zscore(ADC_Net); % normalized network 
% idx_adc = [find(ADC_Net_ZCC > 1.96);find(ADC_Net_ZCC < -1.96)];
idx_adc = find(ADC_p < corrected_alpha);
adjMat_NonS = zeros(size(ADC_Net));
adjMat_NonS(idx_adc) = 1;
corrMat_NonS = adjMat_NonS.*ADC_Net;
A_ADC = corrMat_NonS;
symmetricMatrix_adcnet = A_ADC + A_ADC.' - diag(diag(A_ADC)); % normalized & filtered network 


% SCC tumor
SCCnames = patnames_1(SCClabels); % 共11个鳞状细胞癌
[SCC_Net,SCC_p] = getRmean_mutual(destination_folder, SCCnames, num_roi,selectedIndices);
% corrplotvis_jessy(SCC_Net,ROInames_selected);
% SCC_Net_ZCC = zscore(SCC_Net); % normalized network 
% idx_scc = [find(refNet_ZCC > 1.96);find(refNet_ZCC < -1.96)];
idx_scc = find(SCC_p < corrected_alpha);
adjMat_NonS = zeros(size(SCC_Net));
adjMat_NonS(idx_scc) = 1;
corrMat_NonS = adjMat_NonS.*SCC_Net;
A_SCC = corrMat_NonS;
symmetricMatrix_sccnet = A_SCC + A_SCC.' - diag(diag(A_SCC)); % normalized & filtered network 

% % Benign Nodule % 共 8 个良性肺结节
% Nodulenames = patnames_1(NoduleLabls); % 共?个nodule
% [Nodule_Net,Nodule_p] = getRmean_mutual(destination_folder, Nodulenames, num_roi,selectedIndices);
% % corrplotvis_jessy(Nodule_Net,ROInames_selected);
% % Nodule_Net_ZCC = zscore(Nodule_Net); % normalized network 
% % idx_nol = [find(Nodule_Net_ZCC > 1.96);find(Nodule_Net_ZCC < -1.96)];
% idx_nol = find(Nodule_p < corrected_alpha);
% adjMat_NonS = zeros(size(Nodule_Net));
% adjMat_NonS(idx_nol) = 1;
% corrMat_NonS = adjMat_NonS.*Nodule_Net;
% A_LN = corrMat_NonS;
% symmetricMatrix_nolnet = A_LN + A_LN.' - diag(diag(A_LN)); % normalized & filtered network 
% 
% % small cell 
% SmallCellnames = patnames_1(SmallCellLabls); % 共11个small cell
% [SmallCell_Net, SmallCell_p] = getRmean_mutual(destination_folder, SmallCellnames, num_roi,selectedIndices);
% % corrplotvis_jessy(SmallCell_Net,ROInames_selected);
% % SmallCell_Net_ZCC = zscore(SmallCell_Net); % normalized network 
% % idx_SmallCell = [find(SmallCell_Net_ZCC > 1.96);find(SmallCell_Net_ZCC < -1.96)];
% idx_SmallCell = find(SmallCell_p < corrected_alpha);
% adjMat_NonS = zeros(size(SmallCell_Net));
% adjMat_NonS(idx_SmallCell) = 1;
% corrMat_NonS = adjMat_NonS.*SmallCell_Net;
% A_SCLC = corrMat_NonS;
% symmetricMatrix_SmallCell = A_SCLC + A_SCLC.' - diag(diag(A_SCLC)); % normalized & filtered network 

% Figure Plotting
figure,circularGraph(symmetricMatrix_refnet,symmetricMatrix_refnet,'Label',ROInames_selected); % covid net
figure,circularGraph(symmetricMatrix_adcnet,symmetricMatrix_adcnet,'Label',ROInames_selected); % adc net
figure,circularGraph(symmetricMatrix_sccnet,symmetricMatrix_sccnet,'Label',ROInames_selected); % scc net
% figure,circularGraph(symmetricMatrix_nolnet,symmetricMatrix_nolnet,'Label',ROInames_selected); % nodule net
% figure,circularGraph(symmetricMatrix_SmallCell,symmetricMatrix_SmallCell,'Label',ROInames_selected); % Small cell net

%% ADC-SqCC
close all
cancer_diff = A_ADC - A_SCC;
symmetricMatrix_cancer_diff = cancer_diff + cancer_diff.' - diag(diag(cancer_diff)); 
figure,circularGraph(symmetricMatrix_cancer_diff,symmetricMatrix_cancer_diff,'Label',ROInames_selected); 

%% SECTION 2.1: Calculate group difference and plot (spearman)
close all

ADC_Normal = A_ADC_spear - A_covid_spear;
SCC_Normal = A_SCC_spear - A_covid_spear;
% Nodule_Normal = A_Nodule_spear - A_covid_spear;
% SCLC_Normal = A_SCLC_spear - A_covid_spear;
% 
% ADC_SCC = A_ADC_spear - A_SCC_spear; 
% ADC_SCLC = A_ADC_spear - A_SCLC_spear;
% SCC_SCLC = A_SCC_spear - A_SCLC_spear;


% figure,circularGraph(ADC_Normal,ADC_Normal,'Label',ROInames_selected);
% title('Adenocarcinoma vs. Post-Covid')
% figure,circularGraph(SCC_Normal,SCC_Normal,'Label',ROInames_selected);
% title('Squamous cell carcinoma  vs. Post-Covid')
% figure,circularGraph(ADC_SCC,ADC_SCC,'Label',ROInames_selected);
% title('Adenocarcinoma vs. Squamous cell carcinoma')

% Compute ZCC 
% ADC_Normal_ZCC = zscore(ADC_Normal);
% idx_ADC = [find(ADC_Normal_ZCC > 1.96);find(ADC_Normal_ZCC < -1.96)];
% adjMat_NonS = zeros(size(ADC_Normal));
% adjMat_NonS(idx_ADC) = 1;
% corrMat_NonS = adjMat_NonS.*ADC_Normal;
% A_NonS = corrMat_NonS;
% symmetricMatrix_NonS = A_NonS + A_NonS.' - diag(diag(A_NonS));
% 
% 
% SCC_Normal_ZCC= zscore(SCC_Normal);
% idx_MidS = [find(SCC_Normal_ZCC> 1.96);find(SCC_Normal_ZCC < -1.96)];
% adjMat_MidS = zeros(size(SCC_Normal));
% adjMat_MidS(idx_MidS) = 1;
% corrMat_MidS=adjMat_MidS.*SCC_Normal;
% A_MidS = corrMat_MidS;
% symmetricMatrix_MidS = A_MidS + A_MidS.' - diag(diag(A_MidS));
% 
% 
% ADC_SCC_ZCC = zscore(ADC_SCC);
% idx_HigS = [find(ADC_SCC_ZCC > 1.96);find(ADC_SCC_ZCC < -1.96)];
% adjMat_HigS = zeros(size(ADC_SCC));
% adjMat_HigS(idx_HigS)=1;
% corrMat_HigS = adjMat_HigS.*ADC_SCC;
% A_HigS = corrMat_HigS;
% symmetricMatrix_HigS = A_HigS + A_HigS.' - diag(diag(A_HigS));

figure, h=circularGraph(ADC_Normal,ADC_Normal,'Label',ROInames_selected);
figure,circularGraph(SCC_Normal,SCC_Normal,'Label',ROInames_selected);
% figure,circularGraph(Nodule_Normal,Nodule_Normal,'Label',ROInames_selected);
% figure,circularGraph(SCLC_Normal,SCLC_Normal,'Label',ROInames_selected);
% 
% figure,circularGraph(ADC_SCC,ADC_SCC,'Label',ROInames_selected);
% figure,circularGraph(ADC_SCLC,ADC_SCLC,'Label',ROInames_selected);
% figure,circularGraph(SCC_SCLC,SCC_SCLC,'Label',ROInames_selected);

%% SECTION 2.2: Calculate group difference and plot (Mutual information)
close all

ADC_Normal_m = A_ADC - A_Refnet;
SCC_Normal_m = A_SCC - A_Refnet;
% Nodule_Normal_m = A_LN - A_Refnet;
% SCLC_Normal_m = A_SCLC - A_Refnet;

% ADC_SCC_m = A_ADC - A_SCC; 
% ADC_SCLC_m = A_ADC - A_SCLC;
% SCC_SCLC_m = A_SCC - A_SCLC;

figure,circularGraph(ADC_Normal_m,ADC_Normal_m,'Label',ROInames_selected);
figure,circularGraph(SCC_Normal_m,SCC_Normal_m,'Label',ROInames_selected);
% figure,circularGraph(Nodule_Normal_m,Nodule_Normal_m,'Label',ROInames_selected);
% figure,circularGraph(SCLC_Normal_m,SCLC_Normal_m,'Label',ROInames_selected);
% 
% figure,circularGraph(ADC_SCC_m,ADC_SCC_m,'Label',ROInames_selected);
% figure,circularGraph(ADC_SCLC_m,ADC_SCLC_m,'Label',ROInames_selected);
% figure,circularGraph(SCC_SCLC_m,SCC_SCLC_m,'Label',ROInames_selected);

%% SECTION 3.1: Calculate different systems (brain networks| Spearman)
close all

categories = {
    [12], ...spinal cord
    13, ...brainstem
    [14,15,16,17,18,19],  ... DMN Default mode network  
    [20,21,22,23], ...  Anterior temporal networks
    [24,25], ... Cerebellum network
    [9,10],        ...  呼吸系统
    [7,8],     ...消化系统
    [1,2,3],    ...内分泌系统
     11,           ...免疫系统    
    [5,6],      ...泌尿系统    
    [4]   ...运动系统
};

close all

Systemnames = [{'SpinalCord'}, {'Brainstem'},{'DMN'},{'ATN'},{'CN'},{'Respiratory'},{'Digestive'},{'Endocrine'},{'Lymphatic'},{'Urinary'},{'Muscular'}];
R_mean_HC_system = getRmean_system(A_covid_spear,categories);
refNet_sys = R_mean_HC_system;

% corrplotvis_jessy(refNet_sys,Systemnames);
% refNet_sys_ZCC = zscore(refNet_sys);
% idx_ = [find(refNet_sys_ZCC > 1.96);find(refNet_sys_ZCC < -1.96)];
% adjMat_NonS = zeros(size(refNet_sys));
% adjMat_NonS(idx_) = 1;
% corrMat_NonS = adjMat_NonS.*refNet_sys;
% A_NonS = corrMat_NonS;
symmetricMatrix_refnet = refNet_sys + refNet_sys.' - diag(diag(refNet_sys));
figure,circularGraph(refNet_sys,refNet_sys,'Label',Systemnames);

R_mean_HC_system = getRmean_system(A_ADC_spear,categories);
ADC_Net_sys = R_mean_HC_system; 
figure,circularGraph(ADC_Net_sys,ADC_Net_sys,'Label',Systemnames);


R_mean_HC_system = getRmean_system(A_SCC_spear,categories);
SCC_Net_sys = R_mean_HC_system;
figure,circularGraph(SCC_Net_sys,SCC_Net_sys,'Label',Systemnames);

% R_mean_HC_system = getRmean_system(A_Nodule_spear,categories);
% LN_Net_sys = R_mean_HC_system;
% figure,circularGraph(LN_Net_sys,LN_Net_sys,'Label',Systemnames);
% 
% R_mean_HC_system = getRmean_system(A_SCLC_spear,categories);
% SCLC_system = R_mean_HC_system;
% figure,circularGraph(SCLC_system,SCLC_system,'Label',Systemnames);

%% SECTION 3.2: Calculate different systems (brain networks| Mutual information)
close all

R_mean_HC_system = getRmean_system(A_Refnet,categories);
refNet_sys_m = R_mean_HC_system;
figure,circularGraph(refNet_sys_m,refNet_sys_m,'Label',Systemnames);

R_mean_HC_system = getRmean_system(A_ADC,categories);
ADC_Net_sys_m = R_mean_HC_system; 
figure,circularGraph(ADC_Net_sys_m,ADC_Net_sys_m,'Label',Systemnames);


R_mean_HC_system = getRmean_system(A_SCC,categories);
SCC_Net_sys_m = R_mean_HC_system;
figure,circularGraph(SCC_Net_sys_m,SCC_Net_sys_m,'Label',Systemnames);

% R_mean_HC_system = getRmean_system(A_LN,categories);
% LN_Net_sys_m = R_mean_HC_system;
% figure,circularGraph(LN_Net_sys_m,LN_Net_sys_m,'Label',Systemnames);
% 
% R_mean_HC_system = getRmean_system(A_SCLC,categories);
% SCLC_Net_sys_m = R_mean_HC_system;
% figure,circularGraph(SCLC_Net_sys_m,SCLC_Net_sys_m,'Label',Systemnames);
%% ADC-SqCC(system)
close all
cancer_diff = ADC_Net_sys_m - SCC_Net_sys_m;
symmetricMatrix_cancer_diff = cancer_diff + cancer_diff.' - diag(diag(cancer_diff)); 
figure,circularGraph(symmetricMatrix_cancer_diff,symmetricMatrix_cancer_diff,'Label',ROInames_selected); 

%% SECTION 4.1 SYSTEM DIFF (Spearman)
close all

ADC_Normal_m = ADC_Net_sys - refNet_sys;
SCC_Normal_m = SCC_Net_sys - refNet_sys;
Nodule_Normal_m = LN_Net_sys - refNet_sys;
SCLC_Normal_m = SCLC_system - refNet_sys;

ADC_SCC_m = ADC_Net_sys - SCC_Net_sys; 
ADC_SCLC_m = ADC_Net_sys - SCLC_system;
SCC_SCLC_m = SCC_Net_sys - SCLC_system;

figure,circularGraph(ADC_Normal_m,ADC_Normal_m,'Label',Systemnames);
figure,circularGraph(SCC_Normal_m,SCC_Normal_m,'Label',Systemnames);
figure,circularGraph(Nodule_Normal_m,Nodule_Normal_m,'Label',Systemnames);
figure,circularGraph(SCLC_Normal_m,SCLC_Normal_m,'Label',Systemnames);

figure,circularGraph(ADC_SCC_m,ADC_SCC_m,'Label',Systemnames);
figure,circularGraph(ADC_SCLC_m,ADC_SCLC_m,'Label',Systemnames);
figure,circularGraph(SCC_SCLC_m,SCC_SCLC_m,'Label',Systemnames);


%% SECTION 4.2 SYSTEM DIFF (Mutual)
close all

ADC_Normal_m = ADC_Net_sys_m - refNet_sys_m;
SCC_Normal_m = SCC_Net_sys_m - refNet_sys_m;
% Nodule_Normal_m = LN_Net_sys_m - refNet_sys_m;
% SCLC_Normal_m = SCLC_Net_sys_m - refNet_sys_m;
% 
% ADC_SCC_m = ADC_Net_sys_m - SCC_Net_sys_m; 
% ADC_SCLC_m = ADC_Net_sys_m - SCLC_Net_sys_m;
% SCC_SCLC_m = SCC_Net_sys_m - SCLC_Net_sys_m;

figure,circularGraph(ADC_Normal_m,ADC_Normal_m,'Label',Systemnames);
figure,circularGraph(SCC_Normal_m,SCC_Normal_m,'Label',Systemnames);
% figure,circularGraph(Nodule_Normal_m,Nodule_Normal_m,'Label',Systemnames);
% figure,circularGraph(SCLC_Normal_m,SCLC_Normal_m,'Label',Systemnames);
% 
% figure,circularGraph(ADC_SCC_m,ADC_SCC_m,'Label',Systemnames);
% figure,circularGraph(ADC_SCLC_m,ADC_SCLC_m,'Label',Systemnames);
% figure,circularGraph(SCC_SCLC_m,SCC_SCLC_m,'Label',Systemnames);

%% Calculate Individual difference Adenocarcinoma 
%% Section 5.0 Individual Network visulization only
% ADC (organs)
i = 3; 
[ADC_Net_diff, ADC_p_ind] = getRmean_mutual(destination_folder, [HCPnames; ADCnames(i)], num_roi, selectedIndices);
idx_ = find(ADC_p_ind < corrected_alpha); 
adjMat_NonS = zeros(size(ADC_Net_diff));
adjMat_NonS(idx_) = 1;
corrMat_NonS = adjMat_NonS.*ADC_Net_diff;
ADC_Net_ind = corrMat_NonS - A_Refnet;
symmetricMatrix_ADC_Net_ind  = ADC_Net_ind + ADC_Net_ind.' - diag(diag(ADC_Net_ind)); 
figure,circularGraph(symmetricMatrix_ADC_Net_ind,symmetricMatrix_ADC_Net_ind,'Label',ROInames_selected); 


% SCC (organs)
[SCC_Net_diff, SCC_p_ind] = getRmean_mutual(destination_folder, [HCPnames; SCCnames(i)], num_roi, selectedIndices);
idx_ = find(SCC_p_ind < corrected_alpha); 
adjMat_NonS = zeros(size(SCC_Net_diff));
adjMat_NonS(idx_) = 1;
corrMat_NonS = adjMat_NonS.*ADC_Net_diff;
SCC_Net_ind = corrMat_NonS - A_Refnet;
symmetricMatrix_SCC_Net_ind  = SCC_Net_ind + SCC_Net_ind.' - diag(diag(SCC_Net_ind)); 
figure,circularGraph(symmetricMatrix_SCC_Net_ind,symmetricMatrix_SCC_Net_ind,'Label',ROInames_selected); 
%% system
i=1;
R_mean_HC_system = getRmean_system(ADC_Net_ind,categories);
figure,circularGraph(R_mean_HC_system,R_mean_HC_system,'Label',Systemnames);


%% Function in use
function [R_mean, P_mean] = getRmean(dir,grpnames,num_roi,selectedIndices)
    %%% function to compute mean of entire Ki distribution and 
    % Computes the Pearson/spearman correlation coefficient matrix for the rows in the matrix
     
    %%% dir: Ki distribution path; grpnames: filename; num_roi: number of ROIs

    n = numel(grpnames);
    meanfeat = zeros([n num_roi]);
    for i = 1:n
        disp(['Working on patient ' grpnames{i} '...']);
        featcell = load(fullfile(dir,grpnames{i}));
        % add one row for feature selection 
        featcell_comb = featcell.featcell_comB(1:2, selectedIndices); 
        meanfeat(i,:) = cellfun(@mean,featcell_comb(2,:));
    end
    [R_mean,P_mean] = corr(meanfeat,'Rows','pairwise','type','Pearson');
R_mean = triu(R_mean,1);
P_mean = triu(P_mean,1);
end

function [MI_matrix, p_matrix] = getRmean_mutual(dir,grpnames,num_roi,selectedIndices)
    %%% function to compute mean of entire Ki distribution and 
    % Computes Mutual information using kernel density estimation (KDE)
    %%% dir: Ki distribution path; grpnames: filename; num_roi: number of ROIs

    n = numel(grpnames);
    meanfeat = zeros([n num_roi]);
    for i = 1:n
        disp(['Working on patient ' grpnames{i} '...']);
        featcell = load(fullfile(dir,grpnames{i}));
        % add one row for feature selection 
        featcell_comb = featcell.featcell_comB(1:2, selectedIndices); 
        meanfeat(i,:) = cellfun(@mean,featcell_comb(2,:));
    end
    [numSamples, numRegions] = size(meanfeat);
    MI_matrix = zeros(numRegions, numRegions);
    p_matrix = zeros(numRegions, numRegions); 
    % Permutation parameters
    numPermutations = 1000; % Number of permutations for null distribution

    % Compute pairwise MI and p-values
    for i = 1:numRegions
        for j = 1:numRegions
            if i ~= j
                % Compute observed MI
                MI_matrix(i, j) = computeMutualInformation(meanfeat(:, i), meanfeat(:, j));

                % Generate null distribution via permutations
                null_MI = zeros(numPermutations, 1);
                for k = 1:numPermutations
                    shuffled_y = meanfeat(randperm(numSamples), j); % Shuffle one variable
                    null_MI(k) = computeMutualInformation(meanfeat(:, i), shuffled_y);
                end

                % Compute p-value
                p_matrix(i, j) = sum(null_MI >= MI_matrix(i, j)) / numPermutations; % Right-tailed test
            end
        end
    end

    % Keep only upper triangle of matrices (for symmetry)
    MI_matrix = triu(MI_matrix, 1);
    p_matrix = triu(p_matrix, 1);
end

function MI = computeMutualInformation(x, y)
    % Mutual Information calculation for continuous variables using KDE 

    % Check for constant data
    if std(x) == 0 || std(y) == 0
        MI = 0; % No variability means no information
        return;
    end
    
    % Bandwidth selection for KDE
    bw_x = max(1.06 * std(x) * length(x)^(-1/5), eps); % Ensure positive bandwidth
    bw_y = max(1.06 * std(y) * length(y)^(-1/5), eps);
    
    % Estimate joint entropy
    joint_entropy = kde_entropy([x, y], [bw_x, bw_y]);
    
    % Estimate marginal entropies
    hx = kde_entropy(x, bw_x);
    hy = kde_entropy(y, bw_y);
    
    % Compute Mutual Information
    MI = hx + hy - joint_entropy;
end

function entropy = kde_entropy(data, bandwidth)
    % Kernel Density Estimation-based entropy calculation
    [~, n] = size(data);
    kde_pdf = mvksdensity(data, data, 'Bandwidth', bandwidth); % KDE PDF
    kde_pdf(kde_pdf <= 0) = eps; % Avoid log(0)
    entropy = -mean(log(kde_pdf)); % Shannon entropy
end



function R_mean_system = getRmean_system(R_mean_HC,categories)
% function to compute mean of different systems
% R_mean_HC is the network using mean value

R_symmetric = R_mean_HC + R_mean_HC.';


R_new = zeros(length(categories),length(categories));

for i = 1:length(categories)
    for j = i:length(categories) 
        if i == j
            % 类别内部的元素取均值
            tt = R_symmetric(categories{i}, categories{j});
            R_new(i,j) = mean(mean(tt));
        else
            % 不同类别间的元素直接从 R_mean_HC 中提取
            % 由于可能涉及多个元素，我们这里取它们的均值
            temp = R_symmetric(categories{i}, categories{j});
            R_new(i,j) = mean(temp(:));
        end
    end
end

R_mean_system = R_new;
end