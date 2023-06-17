% File to be used in order to compute temporal analysis of MEG data
% including the following tasks:
% * Data preprocessing (Smoothing, detrending, trial extraction, 
%   baselining, LP filtering).
% * Timelock extraction
% * Plotting of timelock data
% 
% Visitor: 
% Antonio Rodriguez Hidalgo 
% Dept. of Signal Theory and Communications
% Universidad Carlos III de Madrid
% arodh91@gmail.com
%
% Principal Investigator:
% Maria Chait 
% Ear Institute
% University College London
% m.chait@ucl.ac.uk
%
% Last update: 18/June/2018
% Last update: July 2021 (Roberta Bianco)

clc
clear all

subject_list            = [2 3 4 5 6 7 8 9 10 11 12 13 15 16 17 18 19 20 21 22 23 24]; 
%subject_list            = [2:13, 15]; % old subjects
% subject_list = [16:24];   % new subjects
config.hpfreq           =   0;                     
config.lpfreq           =   2; 
% config.baseline         =   'tone';           % Options: 'tone', 'silence' (original prestimuli data), 'activity' (2 to 2.5 seconds)
config.n_components     =   3;                      % DSS components to keep. We usually keep 3 components
config.single           =   0;                      % Set to 1 to compute a DSS matrix for each one of the conditions.
config.out_folder       =   sprintf('Trigger_analysis_PRE_HP%d_LP%d', config.hpfreq, config.lpfreq);  % Output data folder. High passed data (2Hz). Detrended.
config.channels_path    =   'D:\MEGGAP\Channels_DSS\';  % channels selected based on the data after proprocessing 
config.fs               =   600;                      % Sampling frequency for MEG is 600, for EEG is 2048
% config.load_channels    =   0;                      % Sampling frequency
 trigger_list            =   [10, 20];                 % old triggers[5 10 15 20], 5=RAND SHORT, 10=RAND LONG, 15=REG SHORT, 20=REG LONG. Put them in couples (5,15) and (10,20). Indicate the trigger list.
%  trigger_list            =   [5 15]; 
%  trigger_list            =   [10, 20];                  % we use the old triggers in our following analysis for simplicity
%  trigger_list2           =   [40, 80];                % new triggers[40 60 80 100] 40 = RAND SHORT, 60 RAND LONG, 80 = REG SHORT, 100 = REG LONG
trigger_list2           =   [60, 100];  
% addpath('D:\fieldtrip-20220707'); 
addpath('C:\Users\i7 System\Desktop\fieldtrip-20190819');
addpath('D:\NoiseTools');

% We define which parts of the code should be (re)computed.
setup.preprocessingDetrend      =   0;   % Detrending and smoothing
setup.preprocessingMEM          =   1;   % Filtering and more (epoching etc)
% setup.channel_selection         =   0; % Set to 1 to select automatically the optimal channels.
setup.compute_DSScomputation    =   0;
% setup.compute_ToneButterfly     =   0; % Must be performed for HP=2.
% setup.compute_PSD               =   0;

% setup.to_SPM                    =   0;
nicolas_procedure               =   0; % 1: Used to compute DSS maximizing differences between conditions. 0: Computes DSS after grouping the trials of both conditions.




%% Preprocessing
if setup.preprocessingMEM == 1
    % We have considered the standard configuration:
    % * Allows to reject trials visually.
    % * It stores the output structure that contains the merged info of all the
    %   blocks of the subject.
    config.reject_visual    = 1; % Set to 1 to reject trials visually. Use 0 otherwise.
    config.store_data       = 1; % Set to 1 to store the data in a .mat file. Use 0 otherwise (useful for testing purposes).
    config.load_channels    = 0;
    config.pre_filter       = 1; %filter before epoching (critical for the edge artefacts)
    
   for subject = subject_list
       
     if subject < 16 %analyse old subjects
      
       [data_subject] = pre_TempBlockMEM(trigger_list, subject, config);
      
     else  %analyse new subjects

      [data_subject] = pre_TempBlockMEMnewsubj(trigger_list2, subject, config, trigger_list);
      
     end 

   end
    
 end


% %% Channel selection (temporal channels).
% % Channel selection runs independently for each subject, but is obtained
% % from the average signals of all the trigger conditions.
% 
% if setup.channel_selection == 1
%     config.channel_modality = 'temporal';   % 'temporal', 'occipital', 'auto'
%     config.plot_channels = 1;               % In this case, it would plot topography using temporal information.
%     config.store_data = 1;
%     config.DSS = 1;
%     averageTrig = 1;     %1: we average the information of all triggers for channel selection
%                          %0: we average the information of each pair of
%                          %triggers
%     if averageTrig
%        trigger_list = [5,10,15,20];
%     else
%     end 
%     
%     for subject = subject_list
% %       pre_ChannelSelection(trigger_list, subject, config);  
%         pre_ChannelSelection_newSubject(trigger_list, subject, config);   %adapted for newly aquired subjects
%     end
% end


% DSS analysis
 if setup.compute_DSScomputation == 1
   
%      oldSubjects = 1; % 1: we analyse the old subjects 0: we analyse the new subjects
     config.reject_visual    = 0; % Set to 1 to reject trials visually. Use 0 otherwise.
     config.store_data       = 1; % Set to 1 to store the data in a .mat file. Use 0 otherwise (useful for testing purposes).
     DSS_project             = 1; % project the DSS components to the sensor space
   % DSS computation
    if nicolas_procedure == 1
       pre_DSSNicola(trigger_list, subject_list, config); 

    elseif config.single == 1   % single = 1; compute the dss matrix within condition
        for trigger = trigger_list         
            pre_DSScomputation(trigger, subject_list, config);  % We compute the DSS matrix for each trigger.
        end
    else      
       pre_DSScomputation(trigger_list, subject_list, config);  % We compute the DSS matrix.
    end
    
    %% DSS projection
    if DSS_project 

      if config.single == 0  %we only project the components based on two conditions
      DSSprojection_allSubject(trigger_list, subject_list, config);
      else 
      end 

    end 


 end 


% %% Butter fly plotting for tone 
% if setup.compute_ToneButterfly == 1
% 
%     config.reject_visual        = 0; % Set to 1 to reject trials visually. Use 0 otherwise.
%     config.store_data           = 1; % Set to 1 to store the data in a .mat file. Use 0 otherwise (useful for testing purposes).
%     config.load_channels        = 1;
%     config.DSS                  = 0; %0: analyse the raw data; 1: analyse the DSS data
%     config.nicolas_procedure    = nicolas_procedure;
% %     newSubject                  = 0; % set to 1 we analyse new subjects
% %     oldSubject                  = 0; % set to 1 we analyse old subjects
%     config.temporal_frame   = [8,14]; % We determine the temporal frame in seconds that we want to process. dafult [8,14]
% %     config.temporal_frame   = [0, 2.5]; % We determine the temporal frame in seconds that we want to process.
%     config.baseline         = 'tone'; % 'tone', 'silence', 'activity'
%     config.firstCycle       = 1; % Analyse the response of first cycle
%     
% %     if newSubject
% %         
% %     post_ToneButterfly_newSubject(trigger_list, subject_list, config);
% %     
% %     elseif oldSubject
% %         
% %     config.channels_path  = fullfile('..','Results_Antonio_S1_S15','Channels_DSS');  % channels selected based on the data after proprocessing     
% %     subject_list = [2 3 4 5 6 7 8 9 10 11 12 13 15];       
% %     post_ToneButterfly_oldSubject(trigger_list, subject_list, config);
% %     
% %     else
%         
% %     subject_list = [2 3 4 5 6 7 8 9 10 11 12 13 15 16 17 18 19 20 21 22 23 24];  
%     post_ToneButterfly_allSubject(trigger_list, subject_list, config);
% %   post_ToneButterfly(trigger_list, subject_list, config);
% %     end     
% end


%% PSD estimation for and individual subject and condition (trigger).
% In this section we make use of the preprocessed (and stored) information
% and the precomputed channels for a set of subjects.
% if setup.compute_PSD == 1
%     config.store_data           = 1;            % Set to 1 to store the data in a .mat file. Use 0 otherwise (useful for testing purposes).
%     config.channel_modality     = 'temporal';   % 'temporal', 'occipital'
% 
%     [psd_data] = proc_DSS_PSD(trigger_list, subject_list, config);
% 
% end


