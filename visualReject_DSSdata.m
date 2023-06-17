

clear all;
input_folder = 'D:\Results\Trigger_analysis_PRE_HP0_LP2\DSS\DSStransformed';
output_folder = 'D:\Results\Trigger_analysis_PRE_HP0_LP2\DSS\DSStransformed\visualRejection';
% mkdir(output_folder);
addpath('D:\fieldtrip-20220707'); 

n_components = 3; 
trigger_list            =   [10, 20]; 
subject_list = [2 3 4 5 6 7 8 9 10 11 12 13 15 16 17 18 19 20 21 22 23 24]; 

for subject_ind = 1:length(subject_list)
   
    % load channels
    load(fullfile('D:\MEGGAP\Channels_DSS',...
    sprintf('Channels-SUBJ_%d',subject_list(subject_ind))),'channels', 'channels_num');
        
    for trigger_ind = 1:length(trigger_list)

     load(fullfile(input_folder,sprintf('DSSdata_subject-TRIG_%d-SUBJ_%d-COMP_%d.mat',...
     trigger_list(trigger_ind), subject_list(subject_ind), n_components)),'dss_data_subject');

     %% Visually reject bad trials
     cfg.channel = 'all';
     dss_data_subject = ft_rejectvisual(cfg,dss_data_subject);
     
     save(fullfile(output_folder,sprintf('DSSdata_subject-TRIG_%d-SUBJ_%d-COMP_%d.mat',...
     trigger_list(trigger_ind), subject_list(subject_ind), n_components)),'dss_data_subject');
     
         

    end

end