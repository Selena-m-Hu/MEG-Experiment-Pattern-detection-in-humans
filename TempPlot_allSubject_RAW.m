    %% This script put all separate raw ddata files into one file
    % trigger_list: formed by the elements {5, 10, 15, 20]. 
    % * 5 : 3 second RAND sequences.
    % * 10: 15 second RAND sequences.
    % * 15: 3 second REG sequences.
    % * 20: 15 second REG sequences.
    % Updated by Mingyue Hu, 21/11/2022
    
    %--------notes---------------
    % Look at subject 9, whose performance is the outlier in long trial REG
    % condition
    
    clear all;
    clc
    fieldtrip_path      =   'C:\Users\Chait Lab\Documents\MATLAB\Toolboxes\fieldtrip-20170831';
    addpath(genpath(fieldtrip_path)); 
%   in_folder = 'Trigger_analysis_PRE_HP0_LP30';  %channels are saved in the 0-30Hz folder
%   out_folder = 'Trigger_analysis_PRE_HP0_LP2';
%   store_output = config.store_data;    
    config.hpfreq           =   0;                     
    config.lpfreq           =   2; 
    config.out_folder       =   sprintf('Trigger_analysis_PRE_HP%d_LP%d', config.hpfreq, config.lpfreq);  % Output data folder. High passed data (2Hz). Detrended.
    % config.channels_path    =   fullfile('..','Results','Trigger_analysis_PRE_HP0_LP30','Channels');
    trigger_list            =   [10 20]; 
    subject_list = [2 3 4 5 6 7 8 9 10 11 12 13 15 16 17 18 19 20 21 22 23 24];
    out_folder = config.out_folder;
 
    
    allTimelock = []; % save all data into this matrix
 for trigger_ind = 1:length(trigger_list) 
     
     for subject_ind = 1:length(subject_list)
        
            %load channels
%           load(fullfile('..','Results',in_folder,'Channels_DSS',...
%           sprintf('Channels-SUBJ_%d',subject_list(subject_ind))),...
%           'channels', 'channels_num');
        
            load(fullfile('..','Results',out_folder,'Preprocessed_data_AllChannels','short_timelock',...
            sprintf('Short_timelock-TRIG_%d-SUBJ_%d.mat',trigger_list(trigger_ind), subject_list(subject_ind))),'short_data_subject');
            timelock = ft_timelockanalysis([],short_data_subject);
%             rmsData = rms(timelock.avg(channels_num,:),1);
%             all_subjects_timelock(:,trigger_ind,subject_ind) = rmsData;      
      
      %Put all subjects' timelock data from each trigger into one file 
      allTimelock=[allTimelock timelock];
     end     
      
      switch trigger_list(trigger_ind) 
      
          case 5
             save(fullfile('..','Results',out_folder,'Preprocessed_data_AllChannels','short_timelock',...
             sprintf('allRAWDATA_timelock-TRIG_%d.mat',5)),'allTimelock', '-v7.3'); 
          case 10
             save(fullfile('..','Results',out_folder,'Preprocessed_data_AllChannels','short_timelock',...
             sprintf('allRAWDATA_timelock-TRIG_%d.mat',10)),'allTimelock', '-v7.3'); 
          case 15
             save(fullfile('..','Results',out_folder,'Preprocessed_data_AllChannels','short_timelock',...
             sprintf('allRAWDATA_timelock-TRIG_%d.mat',15)),'allTimelock', '-v7.3'); 
          case 20
             save(fullfile('..','Results',out_folder,'Preprocessed_data_AllChannels','short_timelock',...
             sprintf('allRAWDATA_timelock-TRIG_%d.mat',20)),'allTimelock', '-v7.3'); 
      end    
       
  end

  disp('..............finish writing the files................') 

 %Load raw data
%        load();
%         if trigger == 5 | trigger == 15
%             cfg.toilim = [-0.2, 4];   % fast trial
%         else
%            cfg.toilim = [-0.2, 16];  % slow trial
% %         end
%      short_data_subject = ft_redefinetrial(cfg,data_subject);
%     save(fullfile('..','Results',out_folder, 'Preprocessed_data_AllChannels','short_timelock',...
%       sprintf('Short_timelock-TRIG_%d-SUBJ_%d.mat',10, 18)),'short_data_subject'); 
%   
