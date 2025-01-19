%% AIM:  Build networks based on Static SUV mean

addpath(genpath('/Volumes/Extreme Pro/糖代谢/code/circularGraph'));
addpath(genpath('/Volumes/Extreme Pro/糖代谢/code'));

patinfo = readtable('/Volumes/Extreme Pro/糖代谢/其他/info/patinfo.xlsx');

HCPlabels = logical(table2array(patinfo(1:110,8))); % Post-covid control 
ADClabels = logical(table2array(patinfo(1:110,5))); % adenocarcinoma
SCClabels = logical(table2array(patinfo(1:110,6))); % squamous cell carcinoma

%% Selection of ROIs
load('/Volumes/Extreme Pro/糖代谢/Images/SUV/SUV_Distr_body_brain_merged/024_CHANG_LIANG.mat');
ROInames = featcell_comB(1,:);
ROInames{1} = 'L-Adrenal';
ROInames{2} = 'R-Adrenal';
ROInames{10} = 'Thyroid';
ROInames{41} = 'Muscle';    
ROInames{14} = 'L-Kidney';
ROInames{15} = 'R-Kidney';
ROInames{16} = 'Liver';
ROInames{17} = 'L-Lung';
ROInames{18} = 'R-Lung';
ROInames{19} = 'Pancreas';
ROInames{27} = 'Spleen';
ROInames{26} = 'SpinalC';% SpinalCord
ROInames{44} = 'Brainstem';
ROInames{49} = 'L-Hipp';
ROInames{55} = "R-Hipp";% 在文中代表 Hippocampus
ROInames{48} = 'L-PC'; % Cingulate gyrus cinguli posterior part在文中代表Posterior Cingulate
ROInames{54} = 'R-PC';
ROInames{50} = 'L-MF'; %middle frontal gyrus 在文中代表Medial Frontal
ROInames{56} = 'R-MF';
ROInames{46} = 'L-MAT'; %Anterior temporal lobe medial part. 在文中代表Medial Anterior Temporal
ROInames{52} = 'R-MAT';
ROInames{45} = 'L-LAT'; %Anterior temporal lobe lateral part在文中代表Lateral Anterior Temporal
ROInames{51} = 'R-LAT';
ROInames{47} = 'L-Cerebe'; % 在文中代表Cerebellum
ROInames{53} = 'R-Cerebe';

selectedIndices = [1,2,10,41,14,15,16,19,17,18,27,26,44,49,55,48,54,50,56,46,52,45,51,47,53];

ROInames_selected = ROInames(selectedIndices); 

% Beferoni test
num_roi = length(ROInames_selected); % 25 ROI selected
alpha = 0.05; % Original significance level
num_tests = nchoosek(num_roi, 2); % Number of pairwise correlations
corrected_alpha = alpha / num_tests;

%% section 1.1 Calculate Group RefNet, ADC_Net, SCC_Net (using spearman information)
close all
destination_folder = '/Volumes/Extreme Pro/糖代谢/Images/SUV/SUV_Distr_body_brain_merged';
% Get patient indices
in_dir = dir('/Volumes/Extreme Pro/糖代谢/Updated_110_Body_ROI&VOI_MASK43');
patnames_1 = {in_dir(1:end).name}';
patnames_1 = patnames_1(~startsWith(patnames_1, '.'));

HCPnames = patnames_1(HCPlabels); 
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
ADCnames = patnames_1(ADClabels); 
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
SCCnames = patnames_1(SCClabels); 
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

% Figure Plotting
figure,circularGraph(symmetricMatrix_refnet,symmetricMatrix_refnet,'Label',ROInames_selected); % covid net
figure,circularGraph(symmetricMatrix_adcnet,symmetricMatrix_adcnet,'Label',ROInames_selected); % adc net
figure,circularGraph(symmetricMatrix_sccnet,symmetricMatrix_sccnet,'Label',ROInames_selected); % scc net


%% SECTION 1.2: Calculate Group RefNet, ADC_Net, SCC_Net (using mutual information)
close all
destination_folder = '/Volumes/Extreme Pro/糖代谢/Images/SUV/SUV_Distr_body_brain_merged';
% Get patient indices
in_dir = dir('/Volumes/Extreme Pro/糖代谢/Updated_110_Body_ROI&VOI_MASK43');
patnames_1 = {in_dir(1:end).name}';
patnames_1 = patnames_1(~startsWith(patnames_1, '.'));

HCPnames = patnames_1(HCPlabels); % 新冠康复
[refNet,ref_p] = getRmean_mutual(destination_folder, HCPnames, num_roi, selectedIndices);
idx_ = find(ref_p < corrected_alpha); 
adjMat_NonS = zeros(size(refNet));
adjMat_NonS(idx_) = 1;
corrMat_NonS = adjMat_NonS.*refNet;
A_Refnet = corrMat_NonS;
symmetricMatrix_refnet = A_Refnet + A_Refnet.' - diag(diag(A_Refnet)); % normalized & filtered network 

% ADC tumor
ADCnames = patnames_1(ADClabels); % 腺癌
[ADC_Net,ADC_p] = getRmean_mutual(destination_folder, ADCnames, num_roi, selectedIndices);
idx_adc = find(ADC_p < corrected_alpha);
adjMat_NonS = zeros(size(ADC_Net));
adjMat_NonS(idx_adc) = 1;
corrMat_NonS = adjMat_NonS.*ADC_Net;
A_ADC = corrMat_NonS;
symmetricMatrix_adcnet = A_ADC + A_ADC.' - diag(diag(A_ADC)); % normalized & filtered network 

% SCC tumor
SCCnames = patnames_1(SCClabels); % 鳞状细胞癌
[SCC_Net,SCC_p] = getRmean_mutual(destination_folder, SCCnames, num_roi,selectedIndices);
idx_scc = find(SCC_p < corrected_alpha);
adjMat_NonS = zeros(size(SCC_Net));
adjMat_NonS(idx_scc) = 1;
corrMat_NonS = adjMat_NonS.*SCC_Net;
A_SCC = corrMat_NonS;
symmetricMatrix_sccnet = A_SCC + A_SCC.' - diag(diag(A_SCC)); % normalized & filtered network 

% Figure Plotting
figure,circularGraph(symmetricMatrix_refnet,symmetricMatrix_refnet,'Label',ROInames_selected); % covid net
figure,circularGraph(symmetricMatrix_adcnet,symmetricMatrix_adcnet,'Label',ROInames_selected); % adc net
figure,circularGraph(symmetricMatrix_sccnet,symmetricMatrix_sccnet,'Label',ROInames_selected); % scc net

%% ADC- SqCC
close all
cancer_diff = A_ADC - A_SCC;
symmetricMatrix_cancer_diff = cancer_diff + cancer_diff.' - diag(diag(cancer_diff)); 
figure,circularGraph(symmetricMatrix_cancer_diff,symmetricMatrix_cancer_diff,'Label',ROInames_selected); 

%% Section 2.1.1 Calculate the network statbility 
% Cacluate another matrix with # of HC
NumrandomSelect = 15; % Number of people to randomly select
numIterations = 20; % Number of iterations
rhoValues = zeros(numIterations, 1); % Preallocate array for rho values

for iter = 1:numIterations
    disp(['Iteration: ', num2str(iter)]);
    
    % Randomly select HCP names
    HCPnames = patnames_1(HCPlabels);
    selectedPeople = randperm(20, NumrandomSelect); % Randomly select indices
    selectedHCPnames = HCPnames(selectedPeople);
    
    % Compute mutual information matrix
    [matrix2, matrix_p] = getRmean_mutual(destination_folder, selectedHCPnames, num_roi, selectedIndices);
    
    % Extract upper triangular elements (excluding diagonal)
    upper_tri1 = refNet(triu(true(size(refNet)), 1)); % Upper triangle of reference network
    upper_tri2 = matrix2(triu(true(size(matrix2)), 1)); % Upper triangle of computed network
    
    % Compute Spearman correlation
    rho = corr(upper_tri1, upper_tri2, 'Type', 'Spearman'); % Spearman correlation
    rhoValues(iter) = rho; % Store rho value for this iteration
    
    % Optionally display results for each iteration
    disp(['Spearman correlation (rho): ', num2str(rho)]);
end

% Compute mean and standard deviation of rho values
meanRho = mean(rhoValues);
stdRho = std(rhoValues);

% Display final results
disp(['Mean Spearman correlation (rho): ', num2str(meanRho)]);
disp(['Standard deviation of rho: ', num2str(stdRho)]);


%% Section 2.2 Calculate difference (mutual)
close all

ADC_Normal_m = A_ADC - A_Refnet;
SCC_Normal_m = A_SCC - A_Refnet;

figure,circularGraph(ADC_Normal_m,ADC_Normal_m,'Label',ROInames_selected);
figure,circularGraph(SCC_Normal_m,SCC_Normal_m,'Label',ROInames_selected);

%% Section 3.1 Calculate system (spearman)
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
R_mean_HC_system = getRmean_system(A_covid_spear,categories);
refNet_sys = R_mean_HC_system;
Systemnames = [{'SpinalCord'}, {'Brainstem'},{'DMN'},{'ATN'},{'CN'},{'Respiratory'},{'Digestive'},{'Endocrine'},{'Lymphatic'},{'Urinary'},{'Muscular'}];
% Systemnames = [{'SpinalCord'}, {'Brainstem'},{'Default mode network'},{'Left anterior temporal network'},{'Right anterior temporal network'},{'Cerebellum network'},{'Respiratory'},{'Digestive'},{'Endocrine'},{'Lymphatic'},{'Urinary'},{'Muscular'}];


symmetricMatrix_refnet = refNet_sys + refNet_sys.' - diag(diag(refNet_sys));
figure,circularGraph(refNet_sys,refNet_sys,'Label',Systemnames);

R_mean_HC_system = getRmean_system(A_ADC_spear,categories);
ADC_Net_sys = R_mean_HC_system; 
figure,circularGraph(ADC_Net_sys,ADC_Net_sys,'Label',Systemnames);


R_mean_HC_system = getRmean_system(A_SCC_spear,categories);
SCC_Net_sys = R_mean_HC_system;
figure,circularGraph(SCC_Net_sys,SCC_Net_sys,'Label',Systemnames);

%% Section 3.2 Calculate System (mutual info)
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

%% SqCC-ADC (system)
close all
cancer_diff = ADC_Net_sys_m - SCC_Net_sys_m;
symmetricMatrix_cancer_diff = cancer_diff + cancer_diff.' - diag(diag(cancer_diff)); 
figure,circularGraph(symmetricMatrix_cancer_diff,symmetricMatrix_cancer_diff,'Label',ROInames_selected); 


%% Section 4.1 System difference (mutual information)
close all

ADC_Normal_m = ADC_Net_sys_m - refNet_sys_m;
SCC_Normal_m = SCC_Net_sys_m - refNet_sys_m;


figure,circularGraph(ADC_Normal_m,ADC_Normal_m,'Label',Systemnames);
figure,circularGraph(SCC_Normal_m,SCC_Normal_m,'Label',Systemnames);

%% Section 5.0 Individual Network visulization only
% ADC (organs)
close all
i = 3; 
[ADC_Net_diff, ADC_p_ind] = getRmean_mutual(destination_folder, [HCPnames; ADCnames(i)], num_roi, selectedIndices);
% idx_ = find(ADC_p_ind < corrected_alpha); 
% adjMat_NonS = zeros(size(ADC_Net_diff));
% adjMat_NonS(idx_) = 1;
% corrMat_NonS = adjMat_NonS.*ADC_Net_diff;
% ADC_Net_ind = corrMat_NonS - A_Refnet;
ADC_Net_ind = ADC_Net_diff - refNet;
symmetricMatrix_ADC_Net_ind  = ADC_Net_ind + ADC_Net_ind.' - diag(diag(ADC_Net_ind)); 
figure,circularGraph(symmetricMatrix_ADC_Net_ind,symmetricMatrix_ADC_Net_ind,'Label',ROInames_selected); 


% SCC (organs)
[SCC_Net_diff, SCC_p_ind] = getRmean_mutual(destination_folder, [HCPnames; SCCnames(i)], num_roi, selectedIndices);
% idx_ = find(SCC_p_ind < corrected_alpha); 
% adjMat_NonS = zeros(size(SCC_Net_diff));
% adjMat_NonS(idx_) = 1;
% corrMat_NonS = adjMat_NonS.*ADC_Net_diff;
% SCC_Net_ind = corrMat_NonS - A_Refnet;
SCC_Net_ind = SCC_Net_diff - refNet;
symmetricMatrix_SCC_Net_ind  = SCC_Net_ind + SCC_Net_ind.' - diag(diag(SCC_Net_ind)); 
figure,circularGraph(symmetricMatrix_SCC_Net_ind,symmetricMatrix_SCC_Net_ind,'Label',ROInames_selected); 

%% system 
close all
R_mean_HC_system_ADC = getRmean_system(ADC_Net_ind,categories);
figure,circularGraph(R_mean_HC_system_ADC,R_mean_HC_system_ADC,'Label',Systemnames);

R_mean_HC_system_SCC = getRmean_system(SCC_Net_ind,categories);
figure,circularGraph(R_mean_HC_system_SCC,R_mean_HC_system_SCC,'Label',Systemnames);


%% Section 5.1 Individual network for feature extraction (ADC)
close all

numFeatures = (25 * (25 - 1)) / 2;
featureMatrix_ADC = zeros(numel(ADCnames), numFeatures); % Preallocate feature matrix
columnNames = cell(1, numFeatures); % Preallocate column names
% Generate column names for features
idx = 1;
for row = 1:25
for col = row + 1:25 % Upper triangle only
columnNames{idx} = sprintf('R%d_C%d', row, col);
idx = idx + 1;
end
end
% Iterate through ADCnames and compute feature vectors
for i = 1:numel(ADCnames)
[ADC_Net_diff, ADC_p_ind] = getRmean_mutual(destination_folder, [HCPnames; ADCnames(i)], num_roi, selectedIndices);
idx_ = find(ADC_p_ind < corrected_alpha); 
adjMat_NonS = zeros(size(ADC_Net_diff));
adjMat_NonS(idx_) = 1;
corrMat_NonS = adjMat_NonS.*ADC_Net_diff;
ADC_Net_ind = corrMat_NonS - A_Refnet;
% Extract upper triangular elements (excluding diagonal)
upperTriElements = ADC_Net_ind(triu(true(size(ADC_Net_ind)), 1));
% Store the feature vector as a row in the feature matrix
featureMatrix_ADC(i, :) = upperTriElements;
end
% Convert results into a table for easier visualization and access
featureTable_ADC = array2table(featureMatrix_ADC, 'VariableNames', columnNames);
% Display the feature table
disp(featureTable_ADC);

% Section 5.2 Calculate individual network (SqCC)
numFeatures = (25 * (25 - 1)) / 2;
featureMatrix_SCC = zeros(numel(SCCnames), numFeatures); % Preallocate feature matrix
columnNames = cell(1, numFeatures); % Preallocate column names
% Generate column names for features
idx = 1;
for row = 1:25
for col = row + 1:25 % Upper triangle only
columnNames{idx} = sprintf('R%d_C%d', row, col);
idx = idx + 1;
end
end
% Iterate through ADCnames and compute feature vectors
for i = 1:numel(SCCnames)
[SCC_Net_diff, SCC_p_ind] = getRmean_mutual(destination_folder, [HCPnames; SCCnames(i)], num_roi, selectedIndices);
idx_ = find(SCC_p_ind < corrected_alpha); 
adjMat_NonS = zeros(size(SCC_Net_diff));
adjMat_NonS(idx_) = 1;
corrMat_NonS = adjMat_NonS.*SCC_Net_diff;
SCC_Net_ind = corrMat_NonS - A_Refnet;
% Extract upper triangular elements (excluding diagonal)
upperTriElements = SCC_Net_ind(triu(true(size(SCC_Net_ind)), 1));
% Store the feature vector as a row in the feature matrix
featureMatrix_SCC(i, :) = upperTriElements;
end
% Convert results into a table for easier visualization and access
featureTable_SCC = array2table(featureMatrix_SCC, 'VariableNames', columnNames);
% Display the feature table
disp(featureTable_SCC);

%% Section 5.3 Combine ADC and SCC table --> feature selection --> multivariate analysis
% STEP 1: COMBINE TWO TABLES
% Create the cancer type column for each table
cancerType_ADC = zeros(height(featureTable_ADC), 1); % 0 for ADC
cancerType_SCC = ones(height(featureTable_SCC), 1); % 1 for SCC
% Add the cancer type column to each table
featureTable_ADC.Cancer_Type = cancerType_ADC;
featureTable_SCC.Cancer_Type = cancerType_SCC;
% Combine the two tables
combinedFeatureTable = [featureTable_ADC; featureTable_SCC]; % Concatenate tables

%% Nested cross validation (not use this method in this study)
% Parameters
numOuterFolds = 10; % Outer loop folds
numInnerFolds = 3; % Inner loop folds
k = 10; % Number of top features to select

% Data
features = combinedFeatureTable{:, 1:end-1}; % Features (exclude 'Cancer_Type')
labels = combinedFeatureTable.Cancer_Type; % Labels
features = zscore(features); % Standardize features

% Outer Cross-Validation
cvOuter = cvpartition(labels, 'KFold', numOuterFolds); % Outer loop partition
oddsRatiosAll = cell(numOuterFolds, 1); % Store odds ratios for each fold

for outerFold = 1:numOuterFolds
    % Get training and testing indices for the outer fold
    trainIdx = training(cvOuter, outerFold);
    testIdx = test(cvOuter, outerFold);
    
    % Outer training and testing data
    trainFeatures = features(trainIdx, :);
    trainLabels = labels(trainIdx);
    testFeatures = features(testIdx, :);
    testLabels = labels(testIdx);
    
    % Inner Cross-Validation for Feature Selection
    cvInner = cvpartition(trainLabels, 'KFold', numInnerFolds);
    innerFeatureRanks = zeros(size(features, 2), numInnerFolds); % Rank features
    
    for innerFold = 1:numInnerFolds
        % Inner training and validation indices
        innerTrainIdx = training(cvInner, innerFold);
        innerValIdx = test(cvInner, innerFold);
        
        % Inner training data
        innerTrainFeatures = trainFeatures(innerTrainIdx, :);
        innerTrainLabels = trainLabels(innerTrainIdx);
        
        % Perform t-tests for feature selection
        pValues = zeros(size(innerTrainFeatures, 2), 1);
        for i = 1:size(innerTrainFeatures, 2)
            [~, pValues(i)] = ttest2(innerTrainFeatures(innerTrainLabels == 0, i), ...
                                     innerTrainFeatures(innerTrainLabels == 1, i));
        end
        [~, rankIndices] = sort(pValues); % Rank features by p-value
        innerFeatureRanks(:, innerFold) = rankIndices; % Store ranks
    end
    
    % Aggregate feature ranks across inner folds and select top k features
    meanRanks = mean(innerFeatureRanks, 2); % Average rank for each feature
    [~, topFeatureIndices] = mink(meanRanks, k); % Select top k features
    selectedFeatures = trainFeatures(:, topFeatureIndices); % Features for training
    
    % Train LASSO Logistic Regression on Outer Training Data
    [B, FitInfo] = lassoglm(selectedFeatures, trainLabels, 'binomial', 'Link', 'logit');
    [~, idxLambdaMinDeviance] = min(FitInfo.Deviance); % Best lambda
    bestCoefficients = B(:, idxLambdaMinDeviance); % Coefficients
    intercept = FitInfo.Intercept(idxLambdaMinDeviance); % Intercept
    
    % Compute Odds Ratios for Selected Features
    oddsRatios = exp([intercept; bestCoefficients]);
    selectedFeatureNames = combinedFeatureTable.Properties.VariableNames(topFeatureIndices);
    oddsRatiosTable = table(['Intercept'; selectedFeatureNames(:)], [intercept; bestCoefficients], oddsRatios, ...
        'VariableNames', {'Feature', 'Coefficient', 'OddsRatio'});
    
    % Store Results
    oddsRatiosAll{outerFold} = oddsRatiosTable;
    
    % Optionally: Evaluate on the test set (for overall model performance)
    testSelectedFeatures = testFeatures(:, topFeatureIndices);
    logOdds = intercept + testSelectedFeatures * bestCoefficients; % Linear combination
    predictedProbs = exp(logOdds) ./ (1 + exp(logOdds)); % Sigmoid function for probabilities
    predictedLabels = predictedProbs >= 0.5; % Threshold at 0.5
    accuracy = mean(predictedLabels == testLabels); % Calculate accuracy
    disp(['Test Set Accuracy: ', num2str(accuracy)]);
    % Evaluate predictions (e.g., AUC, accuracy)
end

% Summarize Results
disp('Odds Ratios for Each Fold:');
for fold = 1:numOuterFolds
    fprintf('Outer Fold %d:\n', fold);
    disp(oddsRatiosAll{fold});
end


%% Normal two-step logistic regression
% STEP 2: FEATURE SELECTION
% Step 1: Preprocessing and standardization
features = combinedFeatureTable{:, 1:end-1}; % All columns except 'Cancer_Type'
labels = combinedFeatureTable.Cancer_Type;
features = zscore(features); % standardization

% Step 3: Feature Selection
numFeatures = size(features, 2);
pValues = zeros(numFeatures, 1);
for i = 1:numFeatures
    [~, pValues(i)] = ttest2(features(labels == 0, i), features(labels == 1, i)); % Perform univariate statistical test (t-test) for each feature
end

% Select top k features (e.g., top 10 features)
k = 5;
[~, topFeatureIndices] = mink(pValues, k); % Indices of top k features
selectedFeatures = features(:, topFeatureIndices);
selectedFeatureNames = combinedFeatureTable.Properties.VariableNames(topFeatureIndices);

% Standardize the selected features (required for lasso)
features = zscore(selectedFeatures);

% Fit a logistic regression model with lasso regularization
[B, FitInfo] = lassoglm(features, labels, 'binomial', 'Link', 'logit');

% Find the lambda value that minimizes the deviance
[~, idxLambdaMinDeviance] = min(FitInfo.Deviance); % Index of the lambda with minimum deviance
bestCoefficients = B(:, idxLambdaMinDeviance); % Coefficients at the best lambda
intercept = FitInfo.Intercept(idxLambdaMinDeviance); % Intercept at the best lambda

% Compute standard errors for coefficients
% Standard errors are approximated using the Hessian matrix
X = [ones(size(selectedFeatures, 1), 1), selectedFeatures]; % Add intercept to design matrix
y = labels; % Binary outcome (0 or 1)
p = 1 ./ (1 + exp(-(X * [intercept; bestCoefficients]))); % Predicted probabilities
W = diag(p .* (1 - p)); % Weight matrix for logistic regression
Hessian = X' * W * X; % Hessian matrix
covarianceMatrix = inv(Hessian); % Covariance matrix
standardErrors = sqrt(diag(covarianceMatrix)); % Standard errors

% Compute Wald test statistics and p-values
zValues = [intercept; bestCoefficients] ./ standardErrors; % Z-statistics
pValues = 2 * (1 - normcdf(abs(zValues))); % Two-tailed p-values

% Compute odds ratios for the selected coefficients
oddsRatios = exp([intercept; bestCoefficients]);

% Create a results table including p-values
disp('Best Coefficients, Odds Ratios, and P-Values:');
resultsTable = table(['Intercept'; selectedFeatureNames(:)], ...
                     [intercept; bestCoefficients], ...
                     oddsRatios, ...
                     pValues, ...
                     'VariableNames', {'Feature', 'Coefficient', 'OddsRatio', 'PValue'});
disp(resultsTable);

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