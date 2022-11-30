clear;
close all;clc;

folder_name_orig = 'D:\Ying\Anatomical_Model_New\';
result_path = 'D:\Ying\YingModelNew\ProbabilityModel\';
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

average_name = 'D:\Ying\Patients\average.nii.gz';
average_mask = 'D:\Ying\Patients\average_mask.nii.gz';

mask_orig = niftiread(average_mask);
id_s = 23;
id_e = 112;
mask = mask_orig(:,:,id_s:id_e);
ref_mask = 'rig_wk0_organ_masks_avg.nii';

average_img = niftiread(average_name);
bone = average_img > 300;
se = strel('disk',1);
bone_mask = imdilate(bone,se);

len = length(mask(mask>0));

ref_name = 'rig_pCT_avg.nii';
sample_N = 100;
id = 1 : 20;
prob_para_all = struct;

for p = 1 : length(fileNames)
    
    for wk = 1 : 1%length(wkfolder)
        
        list_nii = dir(fullfile(wkfolder{wk,1},'vel*.nii'));
        list_nii_name = {list_nii.name};
        
        id = list_nii_name(contains(list_nii_name,Name{p,1})>0);
        train = setdiff(list_nii_name,id)';
        dmean = single(zeros(size(niftiread([wkfolder{wk,1} '\' list_nii_name{1,1}]))));
        
        info = niftiinfo([wkfolder{wk,1} '\' list_nii_name{1,1}]);
        pca_matrix = zeros(length(train),len*3);
        
        N = length(train);
        for i = 1 : length(train)
            aa = niftiread([wkfolder{wk,1} '\' train{i,1}]);
            x = aa(:,:,id_s:id_e,1); %% displcement along the axis x
            y = aa(:,:,id_s:id_e,2); %% displcement along the axis y
            z = aa(:,:,id_s:id_e,3); %% displcement along the axis z
            pca_matrix(i,:) = [x(mask>0);y(mask>0);z(mask>0)]';
        end
        
        [coeff,score,latent,tsqured,explained,mu] = pca(pca_matrix);
        
        num = select_N(latent,0.9);
        
        min_range = min(score);
        inter = floor(log10(abs(min_range)))-1;
        max_range = max(score);
        sigma_all = 1.06/N^(0.2)*std(score);
        
        joint_N = 15;
        
        edge_all = zeros(joint_N+1,num);
        joint_prob = zeros(joint_N,num);
        joint_prob2 = zeros(sample_N,num);
        para = zeros(sample_N,num);
        
        for i = 1 : num
            range_template = min_range(i) - 10*10^(inter(i)+1) : 10^(inter(i)) : max_range(i) + 10*10^(inter(i)+1);
            score_l = score(:,i);
            range = repmat(range_template,length(score_l),1);
            sigma = sigma_all(i);
            
            P = sum(1/(N*sigma*sqrt(2*pi))*exp(-(range-score_l).^2/(2*sigma^2)));   
            [a,b] = find_ab(range_template,P);
            
            theta = 1/sum(P);
            
            sample = sampleDist(max(P),sample_N,[a,b],score_l,sigma,N);

            sample_ext = repmat(sample,1,length(score_l))';
            
            prob = theta*sum(1/(N*sigma*sqrt(2*pi))*exp(-(sample_ext-score_l).^2/(2*sigma^2)));
            
            joint_prob2(:,i) = prob;
            [jn,edges] = histcounts(sample,joint_N);
            edge_all(:,i) = edges;
            joint_prob(:,i) = jn/sum(jn);
            
            para(:,i) = sample;
        end
        
        [prob_para2,prob_orig] = cal_joint(joint_prob2);
        prob_para = select_prob_para(joint_prob,edge_all,para);
        para_com3 = ['save ' fileNames{p,1} '\sample_para_' wkName{wk,1} '.mat' ' para'];
        para_com2 = ['save ' fileNames{p,1} '\prob2_' wkName{wk,1} '.mat' ' prob_para2'];
        para_com1 = ['save ' fileNames{p,1} '\prob_' wkName{wk,1} '.mat' ' prob_para'];
        para_com4 = ['save ' fileNames{p,1} '\prob_orig_' wkName{wk,1} '.mat' ' prob_orig'];
        eval(para_com2)
        eval(para_com1)
        eval(para_com3)
        eval(para_com4)
        
        for iter = 1 : sample_N
            recons = para(iter,:)*coeff(:,1:num)' + mu;
            dx = dmean(:,:,id_s:id_e,1);
            dx(mask>0) = recons(1:len);
            
            dy = dmean(:,:,id_s:id_e,2);
            dy(mask>0) = recons(len+1:2*len);
            
            dz = dmean(:,:,id_s:id_e,3);
            dz(mask>0) = recons(2*len+1:end);
            
            dmean(:,:,id_s:id_e,1) = dx;
            dmean(:,:,id_s:id_e,2) = dy;
            dmean(:,:,id_s:id_e,3) = dz;
            
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
            inv_trans = [fileNames{p,1} '\' 'predicted_wk' num2str(wk+2) '_' num2str(iter) '.nii'];
            trans_command2 = ['cd ' fileNames{p,1} ' && reg_resample -ref ' ref_wk ...
                ' -flo ' resample_name ' -inter 0 -trans ' save_name  ' -res ' inv_trans ];
            [status_wk_inv,cmdout_wk_inv] = system(trans_command2);
            
            %% mask
            resample_mask_name = [fileNames{p,1} '\wk' num2str(wk+2) '_pred_mask_avg.nii'];
            resample_command2 = ['cd ' result_path Name{p,1} ' && reg_resample -ref ' ref_mask ' -flo ' ref_mask ...
                ' -inter 0 -trans ' flow_name ' -res ' resample_mask_name];
            [status_resample1,cmdout_resample1] = system(resample_command2);
            
            
            inv_trans_mask = [fileNames{p,1} '\' 'predicted_wk' num2str(wk+2) '_' num2str(iter) '_mask.nii'];
            trans_command3 = ['cd ' fileNames{p,1} ' && reg_resample -ref ' ref_wk ...
                ' -flo ' resample_mask_name ' -inter 0 -trans ' save_name  ' -res ' inv_trans_mask ];
            [status_wk_inv1,cmdout_wk_inv1] = system(trans_command3);
            
        end
        
    end

end



