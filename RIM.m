clear;
close all;clc;

folder_name_orig = 'D:\Ying\Anatomical_Model\';
result_path = 'D:\Ying\YingModel\RIM_2\';
dirOutput = dir(fullfile(folder_name_orig,'17*'));
fileNames = cell(length(dirOutput),1);
Name = {dirOutput.name}';
for i = 1: length(dirOutput)
    fileNames{i,1} = [result_path Name{i,1}];
end

dirOutput = dir(fullfile(folder_name_orig,'wk*'));
wkfolder = cell(length(dirOutput),1);
wkName = {dirOutput.name}';
for i = 1: length(dirOutput)
    wkfolder{i,1} = [folder_name_orig wkName{i,1}];
end

average_name = 'D:\Ying\Anatomical_Model_New\average.nii.gz';
average_mask = 'D:\Ying\Anatomical_Model_New\average_mask.nii.gz';
mask = niftiread(average_mask);
len = length(mask(mask>0));

ref_name = 'rig_pCT_avg.nii';
ref_mask = 'rig_wk0_organ_masks_avg.nii';

record = zeros(length(wkfolder)-1,1);
week_vector= 1:6;

for p = 1 : 5
    if p == 13
        week_vector = [1 3:6];
    elseif p == 17
        week_vector = 1:5;
    else
        week_vector = 1:6 ;
    end
    for kk = 2 : length(week_vector)
        
        wk = week_vector(kk);
        
        list_nii = dir(fullfile(wkfolder{wk,1},'vel*.nii'));
        list_nii_name = {list_nii.name};
        id = list_nii_name(contains(list_nii_name,Name{p,1})>0);
        train = list_nii_name(6:end)';
        
        info = niftiinfo([wkfolder{wk,1} '\' train{1,1}]);
        pca_matrix = zeros(length(train),len*3);
        
        N = length(train);
        for i = 1 : length(train)
            aa = niftiread([wkfolder{wk,1} '\' train{i,1}]);
            x = aa(:,:,:,1); %% displcement along the axis x
            y = aa(:,:,:,2); %% displcement along the axis y
            z = aa(:,:,:,3); %% displcement along the axis z
            pca_matrix(i,:) = [x(mask>0);y(mask>0);z(mask>0)]';
        end
        
        mu = mean(pca_matrix,1);
        

        if kk == 2
            %% get patient-specific progressive information
            aa = niftiread([wkfolder{wk,1} '\' id{1,1}]);
            x = aa(:,:,:,1); %% displcement along the axis x
            y = aa(:,:,:,2); %% displcement along the axis y
            z = aa(:,:,:,3); %% displcement along the axis z
            tt = [x(mask>0);y(mask>0);z(mask>0)]';
            dmean = zeros(size(niftiread([wkfolder{wk,1} '\' train{1,1}])));
            rand_wk = tt - mu;
        else
            recons = rand_wk + mu; % refine the prediction for the following weeks
            
            dx = dmean(:,:,:,1);
            dx(mask>0) = recons(1:len);
            
            dy = dmean(:,:,:,2);
            dy(mask>0) = recons(len+1:2*len);
            
            dz = dmean(:,:,:,3);
            dz(mask>0) = recons(2*len+1:end);
            
            dmean(:,:,:,1) = dx;
            dmean(:,:,:,2) = dy;
            dmean(:,:,:,3) = dz;
            
            dmean_name = [result_path Name{p,1} '\wk' num2str(wk+2) '_pca.nii'];
            niftiwrite(single(dmean),dmean_name,info);
            
            vel_info = niftiinfo(dmean_name);
            vel = niftiread(vel_info);
            vel = - vel;
            inv_vel_name = [result_path Name{p,1} '\' 'inv_pca_wk' num2str(wk+2) '_avg.nii'];
            niftiwrite(vel,inv_vel_name,vel_info);
            
            comp_name = [fileNames{p,1} '\wk' num2str(wk+2) '_comp_transformed.nii'];
            trans_command = ['cd ' fileNames{p,1} ' && reg_transform -trans_vel cpp_pCT_avg_backward.nii ' inv_vel_name ' cpp_pCT_avg.nii ' comp_name ...
                ' -ref ' ref_name ' -ref2 ' average_name];
            [status_wk_avg_trans,cmdout_wk_avg_trans] = system(trans_command);
            
            flow_name = [fileNames{p,1} '\wk' num2str(wk+2) '_defflow_reversed.nii'];
            flow_command = ['cd ' fileNames{p,1} ' && reg_transform -flow ' comp_name ' ' flow_name];
            [status_flow,cmdout_flow] = system(flow_command);
            
            disp_name = [fileNames{p,1} '\wk' num2str(wk+2) '_dispflow_reversed.nii'];
            disp_command = ['reg_transform -disp ' flow_name ...
                ' ' disp_name];
            [status_disp,cmdout_disp] = system(disp_command);
            
            aff_txt = 'aff_pCT_average.txt';
            save_name = [fileNames{p,1} '\inv_' aff_txt];
            trans_command = [ 'cd ' fileNames{p,1} ' && reg_transform -invAff ' aff_txt ...
                ' ' save_name];
            [status_trans,cmdout_trans] = system(trans_command);
            
            
            resample_name = [fileNames{p,1} '\wk' num2str(wk+2) '_pred_trans.nii'];
            resample_command = ['cd ' result_path Name{p,1} ' && reg_resample -ref ' ref_name ' -flo ' ref_name ...
                ' -inter 0 -trans ' flow_name ' -res ' resample_name];
            [status_resample,cmdout_resample] = system(resample_command);
            
            
            ref_wk ='wk0.nii';
            inv_trans = [fileNames{p,1} '\' 'predicted_wk' num2str(wk+2) '.nii'];
            trans_command2 = ['cd ' fileNames{p,1} ' && reg_resample -ref ' ref_wk ...
                ' -flo ' resample_name ' -inter 0 -trans ' save_name  ' -res ' inv_trans ];
            [status_wk_inv,cmdout_wk_inv] = system(trans_command2);
            
            %% mask
            resample_mask_name = [fileNames{p,1} '\wk' num2str(wk+2) '_pred_mask_avg.nii'];
            resample_command2 = ['cd ' result_path Name{p,1} ' && reg_resample -ref ' ref_mask ' -flo ' ref_mask ...
                ' -inter 0 -trans ' flow_name ' -res ' resample_mask_name];
            [status_resample1,cmdout_resample1] = system(resample_command2);
            
            
            inv_trans_mask = [fileNames{p,1} '\' 'predicted_wk' num2str(wk+2) '_mask.nii'];
            trans_command3 = ['cd ' fileNames{p,1} ' && reg_resample -ref ' ref_wk ...
                ' -flo ' resample_mask_name ' -inter 0 -trans ' save_name  ' -res ' inv_trans_mask ];
            [status_wk_inv1,cmdout_wk_inv1] = system(trans_command3);

        end
        
    end
    
end




