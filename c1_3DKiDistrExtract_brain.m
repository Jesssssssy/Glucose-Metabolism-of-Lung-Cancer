%% Ki brain region feature extraction

Total_folder_dir_name = '/Volumes/Extreme Pro/糖代谢/HC_LC_Smoking_60_PET_Brain_SEG_83ROI_MASK';
Four_categories_folder_dir = dir(Total_folder_dir_name);
categories = {Four_categories_folder_dir(3:end).name};
Ki_directory = '/Volumes/Extreme Pro/糖代谢/Ki';
out_dir_name = '/Volumes/Extreme Pro/糖代谢/3DKiDistr_brain';


for c = 1:numel(categories)
    disp(['Working on category ' categories{c} '...']);

    % make directory if not present
    if exist([out_dir_name,'/',categories{c}],"dir")==0 
        mkdir(out_dir_name,categories{c});
    end

   
    mask_dir = dir(fullfile(Total_folder_dir_name,categories{c})); 
    patnames = {mask_dir(3:end).name}';
    patnames = patnames(~startsWith(patnames, '.'));
    for i = 1:numel(patnames)
        disp(['Working on patient ' patnames{i} '...']);
        
        % get the Ki map 
        pat_dir_ = dir(fullfile(Ki_directory,patnames{i})); 
        pat_dir = pat_dir_(~startsWith({pat_dir_.name}, '.'));
    
        vol = zeros([192 192 numel(pat_dir)]); % 192 x 192 x n 
        for ii=1:numel(pat_dir)
            raw = dicomread(fullfile(pat_dir(ii).folder,pat_dir(ii).name));
            info = dicominfo(fullfile(pat_dir(ii).folder,pat_dir(ii).name));
            if all(raw(:)== 0) % to handle slices zero-padded for spatial alignment
                vol(:,:,ii) = double(raw);
            else
                vol(:,:,ii) = double(raw) * info.RescaleSlope + ...
                    info.RescaleIntercept;
            end
        end

        % get the mask file
        mask_file_d_ = dir(fullfile(Total_folder_dir_name, categories{c}, patnames{i}));
        mask_file_d = mask_file_d_(endsWith({mask_file_d_.name}, '.nii.gz') & ...
                             ~startsWith({mask_file_d_.name}, '.'));
        mask_file = mask_file_d.name;

        % create feature cell
        related_brain_structures = { "Brainstem"; %19
                                     % DMN 
                                     "L-Hippocampus";  %2
                                     "R-Hippocampus";  %1
                                     "L-Middle and inferior temporal gyrus"; %14
                                     "R-Middle and inferior temporal gyrus"; %13
                                     "L-Cingulate gyrus gyrus cinguli posterior part"; %26
                                     "R-Cingulate gyrus gyrus cinguli posterior part"; %27
                                     "L-Middle frontal gyrus"; %28
                                     "R-Middle frontal gyrus"; %29
                                     % Anterior temporal networks
                                    "L-Anterior-temporal-lobe-medial-part"; %6 
                                    "R-Anterior-temporal-lobe-medial-part"; %5
                                    "L-Anterior-temporal-lobe-lateral-part"; %8
                                    "R-Anterior-temporal-lobe-lateral-part"; %7
                                    % Cerebellem network
                                    "L-Cerebellum"; %18
                                    "R-Cerebellum"; %17
                                   };
        featcell = cell([2 numel(related_brain_structures)]);
        featcell(1,1:end) = related_brain_structures;
        mask = niftiread(fullfile(Total_folder_dir_name, categories{c}, patnames{i},mask_file));
        mask_rotated = imrotate(mask, 90);
        mask_rotated_flipped = flip(mask_rotated, 3);
        
        % n_slices = numel(pat_dir);
        % figure(1);
        % imshow(vol(:,:,77));
        % colormap('jet'); 
        % clim([0 6]);
        % figure(2);
        % imshow(mask_rotated_flipped(:,:,77));
        % colormap('jet'); 
        % clim([0 50]);
        % figure(3);
        % imshow(rotated_mask(:,:,n_slices-77+1));
        % colormap('jet'); 
        % clim([0 50]);
        
        related_brain_structures_number = [19,2,1,14,13,26,27,28,29,6,5,8,7,18,17];
        for j=1:numel(related_brain_structures_number)
            numj = related_brain_structures_number(j);
            Ki_distr = vol(mask_rotated_flipped(:,:,1:200)==numj); % 1.check voi&mask读取顺序; 2.取mask前200个slice的roi.
            featcell{2,j} = Ki_distr;
        end
        save(fullfile(out_dir_name,categories{c},patnames{i}),"featcell");
    end
end
