    %% Temporal plot for all subjects
    % Scripts designed to plot the timelock information for all the
    % triggers (5, 10, 15, 20) for a specific subject (or group of 
    % subjects). Esentially, it plots two separate graphs for LONG and 
    % SHORT, considering the modalities REG and RAND for each one of the
    % configurations. NO TIMELOCK DATA IS COMPUTED IN THIS SCRIPT.
    % 
    % trigger_list: formed by the elements {5, 10, 15, 20]. 
    % * 5 : 3 second RAND sequences.
    % * 10: 15 second RAND sequences.
    % * 15: 3 second REG sequences.
    % * 20: 15 second REG sequences.
    % {Subject 9, the 8th subject is the outlier, whose performance (RMS) in long
    % REG is extremly high during the last half time of the trial.}
    % Last update: 07/June/2018
    % Adapted by Mingyue Hu, Aug, 2022
    
    %--------notes---------------
    % Look at subject 9, whose performance is the outlier in long trial REG
    % condition
    
    clear all;
    clc
    in_folder = 'Trigger_analysis_PRE_HP0_LP30';  %channels are saved in the 0-30Hz folder
    out_folder = 'Trigger_analysis_PRE_HP0_LP2';
%   store_output = config.store_data;    
    config.hpfreq           =   0;                     
    config.lpfreq           =   2; 
    config.out_folder       =   sprintf('Trigger_analysis_PRE_HP%d_LP%d', config.hpfreq, config.lpfreq);  % Output data folder. High passed data (2Hz). Detrended.
    % config.channels_path    =   fullfile('..','Results','Trigger_analysis_PRE_HP0_LP30','Channels');
    trigger_list            =   [10, 20]; 
    subject_list = [2 3 4 5 6 7 8 9 10 11 12 13 15 16 17 18 19 20 21 22 23 24];
%     out_folder = config.out_folder;
    n_components = 3; 
    short = 0;  %1 if we plot the short condition; 0: we plot the long condition
    single = 0; 
   %% We create some matrices with the timelock information.
   %load the timelock data from all new subjects
   
   if single       
   load(fullfile('..','Results',out_folder,'DSS_components','normalizedDSS_cfg',...
   sprintf('ALLchann_allSubj_timelock-TRIG_%d-%d-SINGLE-COMP_%d.mat',trigger_list(1), trigger_list(2), n_components)),'timelock_all');    
   SUBJECT_GLOBAL = [2 3 4 5 6 7 8 9 10 11 12 13 15 16 17 18 19 20 21 22 23 24];    
   
   else
%    load(fullfile('..','Results',out_folder,'DSS_components','normalizedDSS_cfg',...
%    sprintf('ALLchann_allSubj_timelock-TRIG_%d-%d-COMP_%d.mat',trigger_list(1), trigger_list(2), n_components)),'timelock_all');    
%    SUBJECT_GLOBAL = [2 3 4 5 6 7 8 9 10 11 12 13 15 16 17 18 19 20 21 22 23 24];
  load(fullfile('D:\Results',out_folder,'DSS','DSStransformed',...
  sprintf('DSSdata_subject-TRIG_%d-SUBJ_%d-COMP_%d.mat',trigger_list(trigger_ind), subject_list(subject_ind), n_components)),'dss_data_subject');
   SUBJECT_GLOBAL = [2 3 4 5 6 7 8 9 10 11 12 13 15 16 17 18 19 20 21 22 23 24];
   
   end 
   %we load the cfg information from subject 2 (any subject should be ok), as the time information is
   % needed for further plotting 
   load(fullfile('..','Results',out_folder,'DSS_components','normalizedDSS_cfg',...
   sprintf('Xdss-TRIG_%d-SUBJ_%d-SINGLE-Clean-COMP_%d.mat',trigger_list(1), 2, n_components)),'dss_data_subject');
   timelock = ft_timelockanalysis([],dss_data_subject); % we need the time information from the timelock

  for subject_ind = 1:length(subject_list)         
          if subject_list(subject_ind) == SUBJECT_GLOBAL(subject_ind)  
            %load channels
            load(fullfile('..','Results',in_folder,'Channels_DSS',...
            sprintf('Channels-SUBJ_%d',subject_list(subject_ind))),...
            'channels', 'channels_num');
        
            all_subjects_rms(:,:,subject_ind) = rms(timelock_all(channels_num,:,:, subject_ind));
%             all_subjects_rms(:,:,subject_ind) = rms(timelock_all(:,:,:, subject_ind),1);

          end 
   end 
     
   for trigger_ind = 1:length(trigger_list)    
     switch trigger_list(trigger_ind)
        case 5
            out_5 = all_subjects_rms(:,1,:); %rms for selected channels for each subject
            time_5 = timelock.time;
            out_5 = squeeze(out_5);
        case 10
            out_10 = all_subjects_rms(:,1,:); %rms for selected channels for each subject
            time_10 = timelock.time;
            out_10 = squeeze(out_10);

        case 15
            out_15 = all_subjects_rms(:,2,:); %rms for selected channels for each subject
            time_15 = timelock.time;
            out_15 = squeeze(out_15);

        case 20
            out_20 = all_subjects_rms(:,2,:); %rms for selected channels for each subject
            time_20 = timelock.time;
            out_20 = squeeze(out_20);
     end
  end  
       
 
    %% We plot data.   
    if short 
        %% SHORT condition (3 seconds long).        
        figure;
        shadedErrorBar(time_5(1,:),mean(out_5,2)*1e15,(std(out_5')/sqrt(size(out_5,2)))*1e15,'lineProps','k');
        hold on 
        shadedErrorBar(time_15(1,:),mean(out_15,2)*1e15,(std(out_15')/sqrt(size(out_15,2)))*1e15,'lineProps','r');
%         f1 = figure;
%         plot(time_10(1,:),mean(out_10,1)*1e15,'color','k','Linewidth', 1.5); hold on;
%         plot(time_20(1,:),mean(out_20,1)*1e15,'color','r','Linewidth', 1.5)
        xlabel('Time (s)');
        ylabel('Magnitude (fT)');  
        xlim([min(time_5(1,:)), max(time_5(1,:))])
%         legend('RAND','REG')
        title(sprintf('Timelock. LONG. GLOBAL, N=22'))
                
        hold on
        %Bootstrap analysis
          diffREGRAN= out_5 - out_15; 
          dataB=bootstrap(diffREGRAN'); 
          s=findSigDiff(dataB, 0.01);
          s1=findSigDiff(dataB, 0.05);
          
         s_out = proc_DiffPruning(s, trigger_list);
         s_out1 = proc_DiffPruning(s1, trigger_list);
          plot(timelock.time,10*abs(s_out),'color',[0.9290 0.6940 0.1250],'Linewidth', 12);
          plot(timelock.time,20*abs(s_out1),'color',[0.8500 0.3250 0.0980],'Linewidth', 12);  
   else 
        %% LONG condition (15 seconds long). 
        figure;
        shadedErrorBar(time_10(1,:),mean(out_10,2)*1e15,(std(out_10')/sqrt(size(out_10,2)))*1e15,'lineProps','k');
        hold on 
        shadedErrorBar(time_20(1,:),mean(out_20,2)*1e15,(std(out_20')/sqrt(size(out_20,2)))*1e15,'lineProps','r');
%         f1 = figure;
%         plot(time_10(1,:),mean(out_10,1)*1e15,'color','k','Linewidth', 1.5); hold on;
%         plot(time_20(1,:),mean(out_20,1)*1e15,'color','r','Linewidth', 1.5)
        xlabel('Time (s)');
        ylabel('Magnitude (fT)');  
        xlim([min(time_10(1,:)), max(time_10(1,:))])
%         legend('RAND','REG')
        title(sprintf('Timelock. LONG. GLOBAL, N=22'))
                
        hold on
        %Bootstrap analysis
          diffREGRAN= out_10 - out_20; 
          dataB=bootstrap(diffREGRAN'); 
          s=findSigDiff(dataB, 0.01);
         s1=findSigDiff(dataB, 0.05);
         
         s_out = proc_DiffPruning(s, trigger_list);
         s_out1 = proc_DiffPruning(s1, trigger_list);
         
          plot(timelock.time,2*abs(s_out),'color',[0.9290 0.6940 0.1250],'Linewidth', 12);
          plot(timelock.time,4*abs(s_out1),'color',[0.8500 0.3250 0.0980],'Linewidth', 12);   
        
    end
  
  %% Topograph
  
%     load(fullfile('..','Results',out_folder,'Preprocessed_data_AllChannels','short_timelock',...
%                 sprintf('Short_timelock-TRIG_%d-SUBJ_%d.mat',trigger_list(trigger_ind),subject_list(subject_ind))),'short_data_subject')
%                                  
%     timelock = ft_timelockanalysis([],short_data_subject); 
% 
%     timelock.avg = mean(mean(timelock_all(:,:,:,:),4),3);
%     
%     low_t = 0.08;
%     high_t = 0.12;
%     figure;
%     cfg = [];
%     cfg.parameter = 'avg';
%     cfg.layout = 'CTF275.lay';
%     cfg.xlim            = [low_t, high_t]';
%     cfg.marker          = 'labels';
%     cfg.interactive     = 'yes';
%     cfg.comment         = 'no';
%     cfg.colorbar        = 'yes';
%     cfg.highlight       = 'off';
%     cfg.highlightcolor  = 'k';
%     cfg.highlightsymbol = '.';
%     cfg.highlightsize   = 10;
% %     cfg.highlightchannel = [mI; nI];
%     cfg.markersymbol    = '.';
%     cfg.interpolatenan  = 'no';
%     cfg.markersize      = 8;
%     ft_topoplotER(cfg,timelock);

  
  %% Save figure
%   mkdir(fullfile('..','Results',out_folder,'DSS timelock','allRMSplot'));    
%    plotFolder = fullfile('..','Results',out_folder,'DSS timelock','allRMSplot'...
%     ,sprintf('Power-TRIG_%d_%d-allSUBJ',trigger_list(1),trigger_list(2)));
%     savefig(f1,plotFolder)

  
