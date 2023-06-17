
%% Visual rejection, Channel fixation and data saving

%___updated 03/10/2022, Mingyue Hu
    
    clear all;
    clc;
    %----Parameters required-------------
    single                  = 0;   % 1: if you wanna apply the DSS components based on conditions
    config.hpfreq           =   0;                     
    config.lpfreq           =   30; 
    config.out_folder       =   sprintf('Trigger_analysis_PRE_HP%d_LP%d', config.hpfreq, config.lpfreq);  % Output data folder. High passed data (2Hz). Detrended.
    config.channels_path    =   fullfile('..','Results','Trigger_analysis_PRE_HP0_LP30','Channels'); 
    trigger_list            =   [5 15];
    oldSub  = 1;  % load the subjects data from Antonio's folder
%      subject_list = [24];
     subject_list = [2 3 4 5 6 7 8 9 10 11 12 13 15 16 17 18 19 20 21 22 23 24];   %all subjects
    dssAdd = 1;       %1: load DSSed data 0:load raw data
    out_folder = config.out_folder;
    n_components = 3;
    store_data = 1; % 1: save the data matrix for further analysis
    %-----add path---------------
    fieldtrip_path      =   'C:\Users\Chait Lab\Documents\MATLAB\Toolboxes\fieldtrip-20170831';
    addpath(genpath(fieldtrip_path)); 
   

  for subject_ind = 1:length(subject_list)
    for trigger_ind = 1:length(trigger_list)
      
        % If the file exist, we load it
     if (single == 0 & exist(fullfile('..','Results',out_folder,'DSS_components','normalizedDSS_cfg',...
        sprintf('Xdss-TRIG_%d-SUBJ_%d-SINGLE-Clean-COMP_%d.mat',trigger_list(trigger_ind), subject_list(subject_ind), n_components))) == 2) | ...
        (single == 1 &  exist(fullfile('..','Results',out_folder,'DSS_components','normalizedDSS_cfg',...
        sprintf('Xdss-TRIG_%d-SUBJ_%d-Clean-COMP_%d.mat',trigger_list(trigger_ind), subject_list(subject_ind), n_components))) == 2)
    
        if single
        load(fullfile('..','Results',out_folder,'DSS_components','normalizedDSS_cfg',...
        sprintf('Xdss-TRIG_%d-SUBJ_%d-SINGLE-Clean-COMP_%d.mat',trigger_list(trigger_ind), subject_list(subject_ind), n_components)),'dss_data_subject');  
        else 
        load(fullfile('..','Results',out_folder,'DSS_components','normalizedDSS_cfg',...
        sprintf('Xdss-TRIG_%d-SUBJ_%d-Clean-COMP_%d.mat',trigger_list(trigger_ind), subject_list(subject_ind), n_components)),'dss_data_subject');  
        end
       
     else  % If the file does not exist, we generate it      
        %% Load short raw data
        load(fullfile('..','Results',out_folder, 'Preprocessed_data_AllChannels','short_timelock',...
        sprintf('Short_timelock-TRIG_%d-SUBJ_%d.mat',trigger_list(trigger_ind), subject_list(subject_ind)))); 
            
        %% Load DSS data
        if single            
        load(fullfile('..','Results',out_folder,'DSS_components','New_Transformed',...
        sprintf('Xdss-TRIG_%d-SUBJ_%d-SINGLE-COMP_%d.mat',trigger_list(trigger_ind), subject_list(subject_ind), n_components)));    
        else 
        load(fullfile('..','Results',out_folder,'DSS_components','New_Transformed',...
        sprintf('Xdss-TRIG_%d-SUBJ_%d-COMP_%d.mat',trigger_list(trigger_ind), subject_list(subject_ind), n_components)));        
        end 
        data_subject = short_data_subject;

        if dssAdd
        %% replace the data of preprocessing data with the DSSed data
        data_subject.trial = x_dss2;
        else 
        end  

        %% Visual rejection
         cfg = [];
         cfg.channel = 'all';
         data_subject = ft_rejectvisual(cfg,data_subject);
         dss_data_subject = data_subject;
         
        %% Fix channels
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
    
        if single
        save(fullfile('..','Results',out_folder,'DSS_components','normalizedDSS_cfg',...
        sprintf('Xdss-TRIG_%d-SUBJ_%d-SINGLE-Clean-COMP_%d.mat',trigger_list(trigger_ind), subject_list(subject_ind), n_components)),'dss_data_subject');  
        else 
        save(fullfile('..','Results',out_folder,'DSS_components','normalizedDSS_cfg',...
        sprintf('Xdss-TRIG_%d-SUBJ_%d-Clean-COMP_%d.mat',trigger_list(trigger_ind), subject_list(subject_ind), n_components)),'dss_data_subject');  
        end
        
    end
         %% Timelock analysis
         
         dataDSS_timelock = ft_timelockanalysis([],dss_data_subject);       
         timelock_all(:,:, trigger_ind, subject_ind) = dataDSS_timelock.avg;
        
        %% save the out put DSSed data in the fieldTrip structure    
        if exist(fullfile('..','Results',out_folder,'DSS_components','normalizedDSS_cfg')) == 0
        mkdir(fullfile('..','Results',out_folder,'DSS_components','normalizedDSS_cfg'));
        else 
        end 
               
     end
  end 
  
    if store_data
      if single
     save(fullfile('..','Results',out_folder,'DSS_components','normalizedDSS_cfg',...
     sprintf('ALLchann_allSubj_timelock-TRIG_%d-%d-SINGLE-COMP_%d.mat',trigger_list(1), trigger_list(2), n_components)),'timelock_all','-v7.3');
      else
     save(fullfile('..','Results',out_folder,'DSS_components','normalizedDSS_cfg',...
     sprintf('ALLchann_allSubj_timelock-TRIG_%d-%d-COMP_%d.mat',trigger_list(1), trigger_list(2), n_components)),'timelock_all','-v7.3');
      end
      
    end
    
    clear timelock_all
