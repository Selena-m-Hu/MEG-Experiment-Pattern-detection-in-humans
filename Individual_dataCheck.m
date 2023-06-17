
%% Sanity check for individual brain response 
% --- Mingyue Hu, 05/10/2022

 clear all;
 clc;
 %parameter
%  subject_list = [2 3 4 5 6 7 8 9 10 11 12 13 15];
%  subject_list = [16:24];
 subject_list = [2 3 4 5 6 7 8 9 10 11 12 13 15 16 17 18 19 20 21 22 23 24]; 
 out_folder = 'Trigger_analysis_PRE_HP2_LP30'; 
 trigger_list = [10 20];
 n_components = 3; 
 channels_path = fullfile('..','Results','Trigger_analysis_PRE_HP0_LP30','Channels_DSS');
 single = 0;
 dss = 0;  %1:load the DSS data; 0: load the raw data

 %% Load data matrix
 if dss
     if single
     load(fullfile('..','Results',out_folder,'DSS_components','normalizedDSS_cfg',...
     sprintf('ALLchann_allSubj_timelock-TRIG_%d-%d-SINGLE-COMP_%d.mat',trigger_list(1), trigger_list(2), n_components)),'timelock_all');
     else 
     load(fullfile('..','Results',out_folder,'DSS_components','normalizedDSS_cfg',...
     sprintf('ALLchann_allSubj_timelock-TRIG_%d-%d-COMP_%d.mat',trigger_list(1), trigger_list(2), n_components)),'timelock_all');
     end
     
     timelock = ft_timelockanalysis([],short_data_subject); 
 else
     for subject_ind = 1:length(subject_list)
         
         %load selected channels
         load(fullfile(channels_path, sprintf('Channels-SUBJ_%d', subject_list(subject_ind))));
         
         for trigger_ind = 1:length(trigger_list)
              load(fullfile('..','Results',out_folder,'Preprocessed_data_AllChannels','short_timelock',...
              sprintf('Short_timelock-TRIG_%d-SUBJ_%d.mat',trigger_list(trigger_ind),subject_list(subject_ind))),'short_data_subject');
              timelock = ft_timelockanalysis([],short_data_subject);
          
              if trigger_list(trigger_ind) == 10
              RAN_all(:,subject_ind) = rms(timelock.avg(channels_num,:),1)*1e15;
              else 
              REG_all(:,subject_ind) = rms(timelock.avg(channels_num,:),1)*1e15;
              end        
         end
     end
          
 end
 timelock_all(:,:,:,1) = RAN_all;
 timelock_all(:,:,:,2) = REG_all;

 
 save(fullfile('..','Results',out_folder, 'Preprocessed_data_AllChannels','short_timelock',...
 sprintf('All_Subjects_Short_timelock-TRIG-10-20.mat')), 'timelock_all');
 

 load(fullfile('..','Results',out_folder, 'Preprocessed_data_AllChannels','short_timelock',...
 sprintf('Short_timelock-TRIG_%d-SUBJ_%d.mat',10, 23)));
  
 

 %% Plot the data of individuals 
    
    for i = 1:length(subject_list)               
       figure(i)
        
       load(fullfile(channels_path, sprintf('Channels-SUBJ_%d', subject_list(i))));
       randData =  rms(timelock_all(channels_num,:, 1, i),1)*1e15;
       regData  =  rms(timelock_all(channels_num,:, 2, i),1)*1e15;
       plot(timelock.time, randData, 'color', 'k', 'Linewidth',1);   % RAND
       hold on
       plot(timelock.time, regData, 'color', 'r', 'Linewidth',1);   % REG
       
       hold on
%        plot(timelock.time, k1*abs(s),'Linewidth', 3);
%        hold on
%        plot(timelock.time, k2*abs(s2),'Linewidth', 3);
       xlabel('Time (ms)')
       ylabel('RMS magnitude (fT)')
       grid on
%        ylim([0,15])
       RANall(:,i) = randData;
       REGall(:,i) = regData; 
    end 
     
    figure;
    plot(timelock.time, mean(RANall(:,:),2), 'color', 'k', 'Linewidth',1);   % RAND
    hold on
    plot(timelock.time, mean(REGall(:,:),2), 'color', 'r', 'Linewidth',1);   % REG
       
    Diff=RANall-REGall; 

    % We use bootstrap and compute if there is a significant difference
    % between conditions.
    perc = 0.05;
    perc2 = 0.01;
    dataB=bootstrap(Diff'); 
    s=findSigDiff(dataB, perc);
    s2=findSigDiff(dataB, perc2);
    
    s_out = proc_DiffPruning(s, trigger_list);
    s_out2 = proc_DiffPruning(s2, trigger_list);
    
%     k1 = 1;
%     k2 = 2;
    k3 = 3;
    k4 = 5;    
%      hold on
%     plot(timelock.time, k1*abs(s),'Linewidth', 3);
%     hold on
%     plot(timelock.time, k2*abs(s2),'Linewidth', 3);
%     hold on
    plot(timelock.time, k3*abs(s_out),'Linewidth', 12);
    hold on
    plot(timelock.time, k4*abs(s_out2),'Linewidth', 12);
    hold on
       