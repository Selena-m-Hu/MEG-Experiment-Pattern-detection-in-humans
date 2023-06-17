%% Fix the missing channels
% This script is for fixing the dataset with dead channels
% After the fixation we combine all the DSSed data into one matrix for purpose of tone
% analysis 
% adapted by Mingyue Hu, Sep, 2022

clear all;
clc;
config.hpfreq           =   0;                     
config.lpfreq           =   2; 
config.out_folder       =   sprintf('Trigger_analysis_PRE_HP%d_LP%d', config.hpfreq, config.lpfreq);  % Output data folder. High passed data (2Hz). Detrended.
config.channels_path    =   fullfile('..','Results','Trigger_analysis_PRE_HP0_LP30','Channels'); 
trigger_list            =   [10, 20]; 
subject_list = [2 3 4 5 6 7 8 9 10 11 12 13 15 16 17 18 19 20 21 22 23 24];
dssLoad = 1;       %1: load DSSed data 0:load raw data
out_folder = config.out_folder;
n_components = 3; 
store_data = 0;  %1: save the data matrix of all subjects, useful for tone analysis 


for subject_ind = 1:length(subject_list)
    
    for trigger_ind = 1:length(trigger_list)
        
     load(fullfile('..','Results',out_folder,'DSS_components','normalizedDSS_cfg',...
     sprintf('Xdss-TRIG_%d-SUBJ_%d-SINGLE-Clean-COMP_%d.mat',trigger_list(trigger_ind), subject_list(subject_ind), n_components)),'dss_data_subject'); 

     cfg = [];
     cfg.layout = 'CTF275';
     cfg.method = 'template';
     neighbours = ft_prepare_neighbours(cfg);

     cfg = [];
     cfg.neighbours = neighbours;
     cfg.method = 'spline'; % spline can handle missing grad struct
     cfg.layout = 'CTF275'; % use this rather than grad structure to find sens info

     data = dss_data_subject;
     % find missing channels in this dataset
     cfg.missingchannel = {neighbours((~ismember({neighbours(:).label}, data.label))).label};
     data = ft_channelrepair(cfg,data);
     assert(length(data.label) == 275, 'CTF275 not completely interpolated!');

     dss_data_subject = data; 
     
     %% We save the data with corrected channel number
     save(fullfile('..','Results',out_folder,'DSS_components','normalizedDSS_cfg',...
     sprintf('ChannFixedXdss-TRIG_%d-SUBJ_%d-SINGLE-Clean-COMP_%d.mat',trigger_list(trigger_ind), subject_list(subject_ind), n_components)),'dss_data_subject');
     
     x_dss2 = dss_data_subject.trial; 
     x_dss = cat(3,x_dss2{:}); %convert the dss data into 3 dimension and save it into one cell
     
     %% Put data of all subjects into one mat file, for further analysis purpose
     dss_comp(:,:, subject_ind, trigger_ind) = mean(x_dss,3);  %calculate the average of trials and save it into one cell              
         
    end
      
end 

if store_data
 save(fullfile('..','Results',out_folder,'DSS_components','Transformed',...
 sprintf('ALLchann_allSubj_Xdss-TRIG_%d-%d-SINGLE-COMP_%d.mat',trigger_list(1), trigger_list(2), n_components)),'dss_comp','-v7.3');
end 



