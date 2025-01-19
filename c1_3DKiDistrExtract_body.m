%% extract 3D Ki distribution within BODY ROIs
mask_dir_name = '/Volumes/Extreme Pro/糖代谢/Updated_110_Body_ROI&VOI_MASK43';
Ki_dir_name = '/Volumes/Extreme Pro/糖代谢/Ki';
Ki_dir_ = dir(Ki_dir_name);
Ki_dir = Ki_dir_([Ki_dir_.isdir] & ~ismember({Ki_dir_.name}, {'.', '..'}));
patnames = {Ki_dir.name}';
out_dir_name = '/Volumes/Extreme Pro/糖代谢/3DKiDistr';


%% Loop circulation 
for i = 1:numel(patnames)

    disp(['Working on patient ' patnames{i} '...']);

    pat_dir_ = dir(fullfile(Ki_dir_name,patnames{i}));
    pat_dir = pat_dir_(~startsWith({pat_dir_.name}, '.'));
    mask_dir_ = dir(fullfile(mask_dir_name,patnames{i}));
    mask_dir = mask_dir_(~startsWith({mask_dir_.name}, '.'));


    vol = zeros([192 192 numel(pat_dir)]); % 192 x 192

    for ii = 1:numel(pat_dir) 
        raw = dicomread(fullfile(pat_dir(ii).folder,pat_dir(ii).name));
        info = dicominfo(fullfile(pat_dir(ii).folder,pat_dir(ii).name));
        
        if all(raw(:) == 0) % to handle slices zero-padded for spatial alignment
            vol(:,:,ii) = double(raw);
        else
            vol(:,:,ii) = double(raw) * info.RescaleSlope + ...
                info.RescaleIntercept;
        end
    end

    featcell = cell([2 numel(mask_dir)]);
    featcell(1,1:end) = strrep(erase({mask_dir(1:end).name}','.nii'),'_',' ');

    for ii = 1:numel(mask_dir)
        disp(['Extracting from ROI: ' mask_dir(ii).name '...']);

        mask = niftiread(fullfile(mask_dir(ii).folder,mask_dir(ii).name));
        Ki_distr = vol(mask==1); % Return a 1D array containing the extracted values from vol.
        if size(mask,3) == size(vol,3) % check vol and mask length
            featcell{2,ii} = Ki_distr;
        else
            featcell{2,ii} = [];
            disp('Vol and mask not aligned!');
        end
    end

    save(fullfile(out_dir_name,patnames{i}),"featcell");

end
