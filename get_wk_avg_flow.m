clear;
close all;clc;

%% working path includes patient data
%% patient folder
folder_name_orig = 'D:\Ying\Anatomical_Model\';

dirOutput = dir(fullfile(folder_name_orig,'17*'));
fileNames = cell(length(dirOutput),1);
Name = {dirOutput.name}';
for i = 1: length(dirOutput)
    fileNames{i,1} = [folder_name_orig Name{i,1}];
end

%% altas and atlas mask path, mask is the outline of the patient
average_name = 'D:\Ying\Patients\average.nii.gz';
average_mask = 'D:\Ying\Patients\average_mask.nii.gz';

all = 1 : length(Name);
%% get deformable image registration fields for each patient
for p = 1 : length(Name)
    %% in patient folder, wk*.nii is the nifty data of the weekly CT
    list_nii = dir(fullfile(fileNames{p,1},'wk*.nii'));
    list_nii_name = {list_nii.name};
    tmp = regexpi(list_nii_name,'wk..nii','match');
    ind = ~cellfun(@isempty, tmp);
    list_nii_name = list_nii_name(ind);
    %% in patient folder, wk*_mask.nii is the mask of the body of the weekly CTs
    list_mask = dir(fullfile(fileNames{p,1},'wk*_mask.nii'));
    list_mask_name = {list_mask.name};
    tmp = regexpi(list_mask_name,'wk._mask.nii','match');
    ind = ~cellfun(@isempty, tmp);
    list_mask_name = list_mask_name(ind);
    
    flo_name = list_nii_name{1,1};
    flo_mask = list_mask_name{1,1};
    
    % rigid registration between pCT and atlas
    
    rig_pCT_name = 'rig_pCT_avg.nii';
    rig_pCT_mask = 'rig_mask_pCT_avg.nii';
    
    rig_command = [ 'cd ' fileNames{p,1} ' && reg_aladin  -ref ' average_name ...
        ' -rmask ' average_mask ' -flo ' flo_name ' -fmask ' flo_mask ' -rigOnly ' ...
            ' -aff aff_pCT_average.txt  -res ' rig_pCT_name];
    [status_rig_avg,cmdout_rig_avg] = system(rig_command);
    
    rig_mask = [ 'cd ' fileNames{p,1} ' && reg_resample  -ref ' average_name ' -flo ' flo_mask ...
        ' -trans aff_pCT_average.txt ' ' -res '  rig_pCT_mask];   
    [status_rig_mk,cmdout_rig_mk] = system(rig_mask);
    
    % deformalbe registration between pCT and atlas

    dir_command = [ 'cd '  fileNames{p,1} ' && reg_f3d  -ref ' average_name ' -rmask ' average_mask ' -flo ' rig_pCT_name  ...
        ' -fmask ' rig_pCT_mask ...
        ' -be 0.01 -vel --nmi -sx 7 -sy 7 -sz 7 -res dir_pCT_avg.nii -cpp cpp_pCT_avg.nii'];    
    [status_dir_avg,cmdout_dir_avg] = system(dir_command);

    ref_name = rig_pCT_name;
    ref_mask = rig_pCT_mask;
    
    for i = 2 : length(list_nii_name)
        wk = num2str(i)+1;
        flo_name = list_nii_name{1,i};
        flo_mask = list_mask_name{1,i};
        aff_txt = ['aff_wk' wk '_pCT.txt'];
        rig_res_name = ['rig_wk' wk '_pCT.nii'];
        rig_mask_res_name = ['rig_wk' wk '_mask.nii'];
        dir_res_name = ['dir_wk' wk '_pCT.nii'];
        cpp_name =['cpp_wk' wk '_pCT.nii'];
        vel2avg_res_name =[folder_name_orig 'wk' wk '\' 'vel_wk_avg_' Name{p,1} '.nii'];
        
        %% command
        % rigid regisration between wCT and pCT
        rig_command = [ 'cd ' fileNames{p,1} ' && reg_aladin  -ref ' ref_name ...
            ' -rmask ' ref_mask ' -flo ' flo_name ' -fmask ' flo_mask ' -cog '...
            ' -affDirect  -aff ' aff_txt ' -res ' rig_res_name];
        
        [status_rig_wk,cmdout_rig_wk] = system(rig_command);
        
        rig_mask = [ 'cd ' fileNames{p,1} ' && reg_resample  -ref ' ref_name ' -flo ' flo_mask ...
            ' -trans ' aff_txt ' -res ' rig_mask_res_name];
        
        [status_rig_mk,cmdout_rig_mk] = system(rig_mask);
        
        % deformable regisration between wCT and pCT
        dir_command = [ 'cd '  fileNames{p,1} ' && reg_f3d -ref ' ref_name ' -rmask ' ref_mask ' -flo ' rig_res_name ...
            ' -fmask ' rig_mask_res_name ...
            ' -be 0.1 -vel --nmi -sx 5 -sy 5 -sz 5 -res ' dir_res_name ' -cpp ' cpp_name];
        [status_dir_wk,cmdout_dir_wk] = system(dir_command); 
   
%%      convert wk to average anatomy
        trans_command = ['cd ' fileNames{p,1} ' && reg_transform -trans_vel cpp_pCT_avg.nii ' cpp_name ' cpp_pCT_avg_backward.nii ' vel2avg_res_name ...
            ' -ref ' average_name   ' -ref2 ' ref_name];
        [status_wk_avg_trans,cmdout_wk_avg_trans] = system(trans_command);        
        
    end
      
end




