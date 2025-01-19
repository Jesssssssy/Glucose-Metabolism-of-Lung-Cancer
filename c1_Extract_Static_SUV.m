addpath(genpath('/Volumes/Extreme Pro/糖代谢/code/circularGraph'));
addpath(genpath('/Volumes/Extreme Pro/糖代谢/code'));
addpath(genpath('/Volumes/Extreme Pro/糖代谢/code姜然'));
addpath(genpath('/Volumes/Extreme Pro/糖代谢/code姜然/ToolBox'));

%% extract static SUV distribution within given ROIs 
PET_dir_name = '/Volumes/Extreme Pro/糖代谢/Images/PET'; 

% Filter patnames 
patinfo = readtable('/Volumes/Extreme Pro/糖代谢/其他/info/patinfo.xlsx');
HCPlabels = logical(table2array(patinfo(1:110,8))); % Post-covid control 
ADClabels = logical(table2array(patinfo(1:110,5))); % adenocarcinoma
SCClabels = logical(table2array(patinfo(1:110,6))); % squamous cell carcinoma
% Combine the labels using a logical OR operation
combinedLabels = HCPlabels | ADClabels | SCClabels;
% Select PatientName where any of the labels are 1
selectedPatientNames = patinfo.PatientName(combinedLabels);


out_dir_name = '/Volumes/Extreme Pro/糖代谢/Images/SUV/SUV_Distr_body';

mask_dir_name = '/Volumes/Extreme Pro/糖代谢/Updated_110_Body_ROI&VOI_MASK43';

% Loop through each patient:
for i = 1:numel(selectedPatientNames) 

    disp(['Working on patient ' selectedPatientNames{i} '...']);

    pat_dir_ = dir(fullfile(PET_dir_name,selectedPatientNames{i},"PET_Body"));
    pat_dir = pat_dir_(~startsWith({pat_dir_.name}, '.'));
    mask_dir_ = dir(fullfile(mask_dir_name,selectedPatientNames{i}));
    mask_dir = mask_dir_(~startsWith({mask_dir_.name}, '.'));
    
    vol_SUV = zeros([192 192 numel(pat_dir)]);
    for ii = 1:numel(pat_dir) 
        raw = dicomread(fullfile(pat_dir(ii).folder,pat_dir(ii).name));
        info = dicominfo(fullfile(pat_dir(ii).folder,pat_dir(ii).name));
        
        if all(raw(:)== 0) % to handle slices zero-padded for spatial alignment
            vol_SUV(:,:,ii) = double(raw);
        else
            vol_SUV(:,:,ii) = double(raw) * info.RescaleSlope + ...
                    info.RescaleIntercept;
        end
    end

    % Calculate SUV distribution per Body ROI
    featcell = cell([2 numel(mask_dir)]); 
    featcell(1,1:end) = strrep(erase({mask_dir(1:end).name}','.nii'),'_',' ');

    for j = 1:numel(mask_dir)
        disp(['Extracting from ROI: ' mask_dir(j).name '...']);
        mask = niftiread(fullfile(mask_dir(j).folder,mask_dir(j).name));

        SUV_distr = vol_SUV(mask==1);
        if size(mask,3) == size(vol_SUV,3)% check vol and mask length
            featcell{2,j} = SUV_distr;
        else
            featcell{2,j} = [];
            disp('Vol and mask not aligned!');
        end
    end
    save(fullfile(out_dir_name,selectedPatientNames{i}),"featcell");
end

%% Calculate SUV distribution for Brain ROI
% PET_dir_name = '/Volumes/Extreme Pro/糖代谢/Images/PET'; 
% % Filter patnames 
% patinfo = readtable('/Volumes/Extreme Pro/糖代谢/其他/info/patinfo.xlsx');
% HCPlabels = logical(table2array(patinfo(1:110,8))); % Post-covid control 
% ADClabels = logical(table2array(patinfo(1:110,5))); % adenocarcinoma
% SCClabels = logical(table2array(patinfo(1:110,6))); % squamous cell carcinoma
% % Combine the labels using a logical OR operation
% combinedLabels = HCPlabels | ADClabels | SCClabels;
% % Select PatientName where any of the labels are 1
% selectedPatientNames = patinfo.PatientName(combinedLabels);

out_dir_name = '/Volumes/Extreme Pro/糖代谢/Images/SUV/SUV_Distr_brain';
mask_dir_name = '/Volumes/Extreme Pro/糖代谢/Selected_Brain_ROI_mask';

% Loop through each patient:
for i = 1:numel(selectedPatientNames) 

    disp(['Working on patient ' selectedPatientNames{i} '...']);

    pat_dir_ = dir(fullfile(PET_dir_name,selectedPatientNames{i},"PET_Body"));
    pat_dir = pat_dir_(~startsWith({pat_dir_.name}, '.'));
    mask_dir_ = dir(fullfile(mask_dir_name,selectedPatientNames{i}));
    mask_dir = mask_dir_(~startsWith({mask_dir_.name}, '.'));
    vol_SUV = zeros([192 192 numel(pat_dir)]);
    for ii = 1:numel(pat_dir) 
        raw = dicomread(fullfile(pat_dir(ii).folder,pat_dir(ii).name));
        info = dicominfo(fullfile(pat_dir(ii).folder,pat_dir(ii).name));
        
        if all(raw(:)== 0) % to handle slices zero-padded for spatial alignment
            vol_SUV(:,:,ii) = double(raw);
        else
            vol_SUV(:,:,ii) = double(raw) * info.RescaleSlope + ...
                    info.RescaleIntercept;
        end
    end

    % Calculate SUV distribution per Body ROI
    featcell = cell([2 numel(mask_dir)]); 
    featcell(1,1:end) = strrep(erase({mask_dir(1:end).name}','.nii'),'_',' ');

    for j = 1:numel(mask_dir)
        disp(['Extracting from ROI: ' mask_dir(j).name '...']);
        mask = niftiread(fullfile(mask_dir(j).folder,mask_dir(j).name));
        % for brain ROI, rotate and flip
        mask_rotated = imrotate(mask, 90);
        mask_rotated_flipped = flip(mask_rotated, 3);

        SUV_distr = vol_SUV(mask_rotated_flipped==1);
        if size(mask,3) == size(vol_SUV,3) % check vol and mask length
            featcell{2,j} = SUV_distr;
        else
            featcell{2,j} = []; 
            disp('Vol and mask not aligned!');
        end
    end
    save(fullfile(out_dir_name,selectedPatientNames{i}),"featcell");
end




%% Merge brain and body
in_dir = dir('/Volumes/Extreme Pro/糖代谢/Images/SUV/SUV_Distr_body');
patnames_1 = {in_dir(1:end).name}';
patnames_1 = patnames_1(~startsWith(patnames_1, '.'));
in_dir_2 = dir('/Volumes/Extreme Pro/糖代谢/Images/SUV/SUV_Distr_brain');
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
destination_folder = '/Volumes/Extreme Pro/糖代谢/Images/SUV/SUV_Distr_body_brain_merged';
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

