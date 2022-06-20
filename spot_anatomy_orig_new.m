clear;
clc;
close all;
addpath D:\Ying\Proton_Robust\code_WEPL_wk_prob

folder_name = 'D:\Ying\ProtonDataFromMeg\';
mask_path = 'D:\Ying\Anatomical_Model\';
%workpath = 'D:\Ying\code_WEPL_new\';
RD_folder = 'D:\Ying\ProtonDataFromMeg\beam_dose\';
dirOutput = dir(fullfile(folder_name,'17*'));
fileNames = cell(length(dirOutput),1);
RDName = cell(length(dirOutput),1);
maskName = cell(length(dirOutput),1);
Name = {dirOutput.name}';
for i = 1: length(dirOutput)
    fileNames{i,1} = [folder_name Name{i,1}];
    RDName{i,1} = [RD_folder Name{i,1}];
    maskName{i,1} = [mask_path Name{i,1}];
end
% this mat includes the stopping power curve data for different energies. 
load protons_Generic.mat
%%  save wepl for 7 weeks
week_vector = [0 3 4 5 6 7 8];
total_data = struct;
roi_center = zeros(length(week_vector),3);
mvct = true;
spatial_distribution = struct;
%%
color = ['r';'g';'y'];
for kk = 3 : 3%length(fileNames)
    %% get file name
    record_dis = struct;
    %% path of the plan on planning CT and weekly CTs
    wk_folder = dir(fullfile([fileNames{kk,1},'\DICOM\'],'3B*'));
    wk_name = {wk_folder.name}';
    
    wk_path =  cell(length(wk_name),1);
    for tt = 1: length(wk_name)
        wk_path{tt,1} = [[fileNames{kk,1},'\DICOM\'] wk_name{tt,1} ];
    end
    
    masks = dir(fullfile(maskName{i,1},'wk*mask.nii'));
    maskfile = {masks.name}';
    
    %% path of the dose on planning CT and weekly CT
    RD_folder = dir(fullfile(RDName{kk,1},'W*'));
    RD_name = {RD_folder.name}';
    
    for wk = 1 : length(wk_name)
        
        %% here we get the RN of each weeks
        RN = dir(fullfile(wk_path{wk,1},'RN*'));
        RN_name = {RN.name}';
        full_RN_name = strcat(wk_path{wk,1},'\',RN_name);
        full_RN_name = full_RN_name{1,1};
        
        %% here we get the nifty of the each week
        nii = dir(fullfile(wk_path{wk,1},'*.nii*'));
        nii_name = {nii.name}';
        full_nii_name = strcat(wk_path{wk,1},'\',nii_name); %% 
        nii_name_wk = full_nii_name{1,1};
        
        %% here we get the coresponding dose matric. 
        RD = dir(fullfile([RDName{kk,1} '\' RD_name{wk,1}],'RD*'));
        RD = {RD.name}';
        
        mask = niftiread([maskName{kk,1} '\' maskfile{wk,1}]);

        %% info
        nii_info = niftiinfo(nii_name_wk);
        ct_img = niftiread(nii_info);
        ct_img = int16(ct_img);
        ct_orig = ct_img*nii_info.raw.scl_slope + nii_info.raw.scl_inter;

        ct_orig(mask ==0) = min(ct_orig(:)); 
        ct_orig = rot90(ct_orig);
        
        %% this block is used to interpolate the relationship between HU and RSP
        HU_CT = [-200 -120 -20 35 100 140];
        Y = [0.793 0.957 1.013 1.025 1.09 1.106];
        hu = -200:1:140;
        yy = interp1(HU_CT,Y,hu);
        
        %% convert to RSP map
        %% convert HU is a function used to convert HU to RSP for whole images
        RSP = convertHU(double(ct_orig),mvct,hu,yy);
        
        RN_info = dicominfo(full_RN_name);
        dx=nii_info.raw.pixdim(2); % length between each planes in x axis (mm)
        dy=nii_info.raw.pixdim(3);
        dz=nii_info.raw.pixdim(4);
        
        %% beam info
        beamFieldNames = fieldnames(RN_info.IonBeamSequence);
        nb_fields = length(beamFieldNames);
        myBeamData = cell(nb_fields,1);
        beam_spots = zeros(size(ct_orig));
        field_spatial = struct;
        ind = 1;
        record_distribution = zeros(size(ct_orig,1)*size(ct_orig,2)*size(ct_orig,3),13);
        
        for i = 1:nb_fields
            img_beam = zeros(size(ct_orig));
            beamSequence = RN_info.IonBeamSequence.(beamFieldNames{i}).IonControlPointSequence;
            finalCumulativeMetersetWeight = RN_info.IonBeamSequence.(beamFieldNames{i}).FinalCumulativeMetersetWeight;
            layerFieldNames = fieldnames(RN_info.IonBeamSequence.(beamFieldNames{i}).IonControlPointSequence);
            nb_layers = length(layerFieldNames)/2;
            BeamMeterset = RN_info.FractionGroupSequence.Item_1.ReferencedBeamSequence.(beamFieldNames{i}).BeamMeterset;
            %% dose_registration is used to map the dose to the CT images. 
            [cube_d_new,xVec,yVec,zVec,dose_info] = dose_registration(wk_path{wk,1},[RDName{kk,1} '\' RD_name{wk,1}],RD{i,1});
            for j = 1:nb_layers
                
                if j == 1
                    gantry_angle = beamSequence.(layerFieldNames{2*j-1}).GantryAngle;
                    rangeShifterWaterEquivalentThickness = beamSequence.(layerFieldNames{2*j-1}).RangeShifterSettingsSequence.Item_1.RangeShifterWaterEquivalentThickness;
                    isocenterToRangeShifterDistance  = beamSequence.(layerFieldNames{2*j-1}).RangeShifterSettingsSequence.Item_1.IsocenterToRangeShifterDistance;
                    isocenter_orig = beamSequence.(layerFieldNames{2*j-1}).IsocenterPosition;
                    isocenter(1:2) = -isocenter_orig(1:2);
                    isocenter(3) = isocenter_orig(3);
                    center = iso_center_change(isocenter,nii_info);
                    iso_center_im(2) = center(1);
                    iso_center_im(1) = nii_info.ImageSize(1) - center(2);
                    iso_center_im(3) = center(3);
                    couch_angle = 0;
                    Length = 1000;
                end
                cumulativeMetersetWeight = beamSequence.(layerFieldNames{2*j-1}).CumulativeMetersetWeight;
                scanSpotMetersetWeights = beamSequence.(layerFieldNames{2*j-1}).ScanSpotMetersetWeights;
                spotMeterset = scanSpotMetersetWeights;
                
                xy = reshape(beamSequence.(layerFieldNames{2*j-1}).ScanSpotPositionMap,2,beamSequence.(layerFieldNames{2*j-1}).NumberOfScanSpotPositions).';
                energy = beamSequence.(layerFieldNames{2*j-1}).NominalBeamEnergy;
                energy_z = energy_to_peakPos(energy,machine.data,'range') - rangeShifterWaterEquivalentThickness;
                
                xy_virtual = zeros(size(xy,1),3);
                theta = gantry_angle/180*pi;
                
                spatial_info = zeros(size(xy,1),6);
                
                xy_virtual(:,1) = xy(:,1)*sin(theta)/dx + double(iso_center_im(1));
                xy_virtual(:,2) = xy(:,1)*cos(theta)/dy + double(iso_center_im(2));
                xy_virtual(:,3) = double(xy(:,2)/dz) + double(iso_center_im(3));
                xy_virtual = round(xy_virtual);
                
                for spot_num = 1 : size(xy,1)
                    spot_loc = double(xy_virtual(spot_num,:));
                    [i_label_t,j_label_t,~] = cal_wepl_z(nii_info,spot_loc,gantry_angle,couch_angle,Length,ct_orig,RSP,energy_z);
                    if ~isempty(i_label_t) && ~isempty(j_label_t)
                        v_index =[i_label_t,j_label_t,ceil(xy_virtual(spot_num,3))];
                        row = v_index(1);
                        col = v_index(2);
                        page = v_index(3);
                        img_beam(row,col,page) = img_beam(row,col,page) + spotMeterset(spot_num);
                        Eud = sqrt(((row - iso_center_im(1))*dx).^2 + ((col-iso_center_im(2))*dy).^2 + ((page-iso_center_im(3))*dz).^2);
                        record_distribution(ind,1) = row;
                        record_distribution(ind,2) = col;
                        record_distribution(ind,3) = page;
                        record_distribution(ind,4) = (row - iso_center_im(1))*dx;
                        record_distribution(ind,5) = (col - iso_center_im(2))*dy;
                        record_distribution(ind,6) = (page - iso_center_im(3))*dz;
                        record_distribution(ind,7) = Eud;
                        record_distribution(ind,8) = gantry_angle;
                        record_distribution(ind,9) = xy(spot_num,1);
                        record_distribution(ind,10) = xy(spot_num,2);
                        record_distribution(ind,11) = energy;
                        record_distribution(ind,12) = cube_d_new(row,col,page);
%                         record_distribution(ind,13) = energy_z;
                        record_distribution(ind,13) = img_beam(row,col,page);
                        ind = ind +1;
                    end
                end
            end  
            beam_spots = beam_spots + img_beam;
        end
       
%         range_end = find(sum(sum(beam_spots))>0,1,'last');
%         range_start = find(sum(sum(beam_spots))>0,1,'first');
            
%         for slice = range_start : range_end
%             figure(1);
%             imshow(ct_orig(:,:,slice),[]);
%             show_structure(0.2,beam_spots(:,:,slice),'r');
%             pause;
%         end 
%         eval(['save ' Name{kk,1} '_' wk_name{wk,1} '.mat' ' beam_spots']);

        record_distribution = record_distribution(record_distribution(:,end)>0,:);
        record_dis.(['week' num2str(wk)]) = record_distribution;
    end
    
    eval(['save ref\' Name{kk,1} '_pos_orig_all_new2.mat record_dis'])
end