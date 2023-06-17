function post_ToneButterfly_allSubject(trigger_list, subject_list, config)
    % trigger_list: 
    % * [5, 15]     : SHORT sequences, 3 seconds
    % * [10, 20]    : long sequences, 15 seconds
    %
    % subject_list:
    % * 2-15, the index of the subject aquired by Antonio.
    % * 16-14, subjects aquired by RB and MYH
    % config: allows for certain configurations.
    %   .out_folder: name of the folder where data will be stored.
    %   .store_data: set to 1 to store in a .mat file the whole data of the
    %       subject once that it has been processed.
    %   .temporal_frame: temporal region where the tones of interest are
    %       located. It gets divided into T_init and T_end, which indicate
    %       the beginning and the end of the data of interest. Is a vector
    %       with [T_init, T_end] in seconds.
    %   .hpfreq: in this case, useful to store the results in the proper
    %       folder.
    %   .lpfreq: in this case, useful to store the results in the proper 
    %       folder.
    %   .baseline: baseline modality that we want to use. We have three
    %       different ones:
    %       * 'tone': baselined independently for each tone using its
    %           initial 50ms (from the total of 250ms per tone).
    %       * 'silence': baseline data is computed using the prestimuli 
    %           time. In our case, 500ms.
    %       * 'activity': baseline data is acquired from 2 to 2.5 seconds
    %           poststimuli, just 500ms before the second repetition of the
    %           sequence occurs (in case of REG).    
    %   .channels_path: path indicating where to find the channel
    %       information file.
    %
    % OUTPUT FOLDER (one of the following):
    %   * ../Results/***/Global_trigger,
    %   * ../Results/***/Global_trigger_silence
    %   * ../Results/***/Global_trigger_activity
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
    % Last update: 16/July/2018
    % Last update: 07/07/2022 By Mingyue Hu, Ear Institute
    
    
    out_folder = config.out_folder;
    firstCycle      = config.firstCycle;
    store_output    = config.store_data;
    T_init          = config.temporal_frame(1);
    T_end           = config.temporal_frame(2);
    hpfreq          = config.hpfreq;
    lpfreq          = config.lpfreq;
    baseline        = config.baseline;
    fs              = 600; % Sampling frequency
    pre_estim       = 0.2; % Prestimuli time. In our case, 200ms.
    window_size     = 0.250*fs; % 250 ms, the length of the stimuli (50ms signal + 200ms silence).
    channels_path   = config.channels_path;
    dss_flag        = config.DSS;
    n_components    = config.n_components;
    single = config.single;    
%    trigger_length = 0.25; % 250ms.
   
    %% If you analyse the data before transition
    if firstCycle == 1
    T_init          = 0;
    T_end           = 2.5;
    end 
    


for trigger_ind = 1:length(trigger_list)
   for subject_ind = 1:length(subject_list)

     
    %% We load the timelock information that we computed in a previous stage.
            
            
%             load(fullfile('..','Results',out_folder,'Preprocessed_data_AllChannels','short_timelock',...
%                 sprintf('Short_timelock-TRIG_%d-SUBJ_%d.mat',trigger_list(trigger_ind),subject_list(subject_ind))),'short_data_subject')
      if dss_flag == 1
           if single    
            load(fullfile('..','Results',out_folder,'DSS_components','normalizedDSS_cfg',...
             sprintf('Xdss-TRIG_%d-SUBJ_%d-SINGLE-Clean-COMP_%d',trigger_list(trigger_ind),subject_list(subject_ind), n_components)), 'dss_data_subject') 
           else
            load(fullfile('..','Results',out_folder,'DSS_components','normalizedDSS_cfg',...
             sprintf('Xdss-TRIG_%d-SUBJ_%d-Clean-COMP_%d',trigger_list(trigger_ind),subject_list(subject_ind), n_components)), 'dss_data_subject')    
           end
           timelock = ft_timelockanalysis([],dss_data_subject);
      else
           load(fullfile('..','Results',out_folder,'Preprocessed_data_AllChannels','short_timelock',...
               sprintf('Short_timelock-TRIG_%d-SUBJ_%d.mat',trigger_list(trigger_ind),subject_list(subject_ind))),'short_data_subject');
           timelock = ft_timelockanalysis([],short_data_subject);
      end
                                       
            timelock_subject{trigger_ind} = timelock;
            timelock275 = timelock; 
            

            
%           We load the channels that we obtained in channel selection.
            load(fullfile(channels_path, sprintf('Channels-SUBJ_%d', subject_list(subject_ind))));
%           channels_num = 40;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            counter = 1;
            
            for t_ind = 1:window_size:length(timelock.time)
                % We fill the first one (which will not be used) with
                % zeros, since it will be incomplete. This way, the
                % structure of the data is easier and we consider a 50ms
                % to baseline and 200ms of signal. (so the trial structure is pre-stim(50 ms) + tone(50ms) + silence(150ms) =
                % 150 ms). Here we use the pre-stim to conduct 'tone'
                % baseline 
                %-----------------------------------------------------------
                % Notbly, the new subjects (16-24) only include 273 channels,
                % particularly, subject 16 only has 272 channels, all
                % channels have been fixed as 275
                %-----------------------------------------------------------
                if t_ind == 1                 
%                     trigger_shape275(:,:,counter) = zeros(275,150);
                    
                    trigger_shape(:,:,counter) = zeros(length(channels_num),150); %first counter include no information
                else  %[-0.2,0] include no information,the tone onset corresponding to time point of 120, 90 is -0.5 sec(no signal), followed with 200 ms signal
                    trigger_shape(:,:,counter) = timelock.avg(channels_num,t_ind-.1*fs:t_ind+0.15*fs-1); %second counter(first tone), information should start from -0.05
%                     trigger_shape275(:,:,counter) = timelock275.avg(:,t_ind-.1*fs:t_ind+0.15*fs-1);
                    % Security check (at t_ind == 151). % Should be equal to [-0.05, 0.2] 
%                     timelock.time(t_ind-.1*fs:t_ind+0.15*fs-1) 
                end
                counter = counter + 1;
            end
            
%            Multiply by a constant in order to get the units into
            % femtoTeslas.
            trigger_shape = trigger_shape*1e15;
%             trigger_shape275 = trigger_shape275*1e15;

            % We get the value of the initial time (0 seconds).
            t0 = 0.25/(window_size/fs);
            
            % We get the temporal index of the beginning and the end of our
            % data from the global matrix.
            t_stable = max(t0,1) + ((T_init/(window_size/fs)):(T_end/(window_size/fs)));
            stable_average_shape(:,:,subject_ind, trigger_ind)     = mean(trigger_shape(:,:,t_stable),3);
%             stable_average_shape275(:,:,subject_ind, trigger_ind)  = mean(trigger_shape275(:,:,t_stable),3);
            
            % We baseline the data according to the criteria we want to
            % use.
            % 'tone'    :  each tone is baselined using its initial 50ms.
            % 'silence' :  each tone is baselined using the 500ms 
            %              prestimuli data.  
            % 'activity':  each tone is baselined using the data from 2 to
            %              2.5 seconds poststimuli, just before the tone
            %              sequence starts repeating.
            switch baseline
                
                case 'tone'
%                     bl_window = 11:30;  % define the time window you want to use for tone baseline
%                     bl_tone = '11-30'; %useful for saving file name 
%                       bl_window = 21:30;  % define the time window you want to use for tone baseline
%                       bl_tone = '21-30'; %useful for saving file name 
                        bl_window = 30;  % define the time window you want to use for tone baseline
                        bl_tone = '-30'; %useful for saving file name 
                    baseline_data = mean(stable_average_shape(:,bl_window,subject_ind, trigger_ind), 2); % we only baseline from -30ms from the onset of the tone
%                     baseline_data275 = mean(stable_average_shape275(:,1:30,subject_ind, trigger_ind), 2);
                    output_appendix = '_tone';
                case 'silence'
                    baseline_data = mean(timelock.avg(channels_num, 1:fs*pre_estim),2)*1e15;
%                     baseline_data275 = mean(timelock275.avg(:, 1:fs*pre_estim),2)*1e15;
                    bl_tone = '-';
                    output_appendix = '_silence';
                case 'activity'
                    baseline_data = mean(timelock.avg(channels_num, fs*(.5+2.0):fs*(.5+2.5)),2)*1e15;
%                     baseline_data275 = mean(timelock275.avg(:, fs*(.5+2.0):fs*(.5+2.5)),2)*1e15;
                    output_appendix = '_activity';
                    bl_tone = '-';

            end
          
                 
           % stable_shape(:,:,:,subject_ind, trigger_ind) = trigger_shape - repmat(mean(stable_average_shape(:,1:30,subject_ind, trigger_ind) ,2), 1, 150, size(trigger_shape,3));
            stable_average_shape(:,:,subject_ind, trigger_ind)  =  stable_average_shape(:,:,subject_ind, trigger_ind) - repmat(baseline_data, 1, 150);
%            stable_average_shape275(:,:,subject_ind, trigger_ind)  =  stable_average_shape275(:,:,subject_ind, trigger_ind) - repmat(baseline_data275, 1, 150);
            
            clear trigger_shape
  
    end    
end
    
  %% We save the matrix into memory 
if dss_flag == 1
 if single 
%   output_appendix = 'no_baseline';
  if firstCycle
    save(fullfile('..','Results',out_folder,'tone_analysis',...
    sprintf('Tone_40Channels_allSub-FirstCycle-SINGLE-COMP_%d_BL_%s_BLwindow_%s', n_components, output_appendix, bl_tone)), 'stable_average_shape');

%     save(fullfile('..','Results',out_folder,'tone_analysis',...
%     sprintf('Tone_275Channels_allSub-FirstCycle-SINGLE-COMP_%d_BL_%s_BLwindow_%s', n_components, output_appendix, bl_tone)), 'stable_average_shape275');
         
  else 
    save(fullfile('..','Results',out_folder,'tone_analysis',...
    sprintf('Tone_40Channels_allSub-SINGLE-COMP_%d_BL_%s_BLwindow_%s', n_components, output_appendix, bl_tone)), 'stable_average_shape');

%     save(fullfile('..','Results',out_folder,'tone_analysis',...
%     sprintf('Tone_275Channels_allSub-SINGLE-COMP_%d_BL_%s_BLwindow_%s', n_components, output_appendix, bl_tone)), 'stable_average_shape275');
  end 
  
 else 
  if firstCycle
    save(fullfile('..','Results',out_folder,'tone_analysis',...
    sprintf('Tone_40Channels_allSub-FirstCycle-COMP_%d_BL_%s_BLwindow_%s', n_components, output_appendix, bl_tone)), 'stable_average_shape');

%     save(fullfile('..','Results',out_folder,'tone_analysis',...
%     sprintf('Tone_275Channels_allSub-FirstCycle-COMP_%d_BL_%s_BLwindow_%s', n_components, output_appendix, bl_tone)), 'stable_average_shape275');
         
  else 
    save(fullfile('..','Results',out_folder,'tone_analysis',...
    sprintf('Tone_40Channels_allSub-COMP_%d_BL_%s_BLwindow_%s', n_components, output_appendix, bl_tone)), 'stable_average_shape');

%     save(fullfile('..','Results',out_folder,'tone_analysis',...
%     sprintf('Tone_275Channels_allSub-COMP_%d_BL_%s_BLwindow_%s', n_components, output_appendix, bl_tone)), 'stable_average_shape275');
  end
 end 
 
else
    
    if firstCycle
    save(fullfile('..','Results',out_folder,'tone_analysis',...
    sprintf('Tone_40Channels_allSub-FirstCycle-RAWdata_BL_%s_BLwindow_%s', output_appendix, bl_tone)), 'stable_average_shape');

%     save(fullfile('..','Results',out_folder,'tone_analysis',...
%     sprintf('Tone_275Channels_allSub-FirstCycle-RAWdata_BL_%s_BLwindow_%s', output_appendix, bl_tone)), 'stable_average_shape275');
         
  else 
    save(fullfile('..','Results',out_folder,'tone_analysis',...
    sprintf('Tone_40Channels_allSub-RAWdata_BL_%s_BLwindow_%s', output_appendix, bl_tone)), 'stable_average_shape');

%     save(fullfile('..','Results',out_folder,'tone_analysis',...
%     sprintf('Tone_275Channels_allSub-RAWdata_BL_%s_BLwindow_%s', output_appendix, bl_tone)), 'stable_average_shape275');
    end
end   
    %% Subject curve-topography graphs
%     figure('units','normalized','outerposition',[0 0 1 1])  
    for subject_ind = 1:length(subject_list)
        subject_data.avg = mean(stable_average_shape(:,:, subject_ind,:),4);
        
        subject_data.time = (1:150)/150*250-50;
        subject_data.label = timelock.label;
        subject_data.dimord = 'chan_time';
        
        
        % We load and plot the temporal data.
        channels_num = 1:length(channels_num);
        mean_subjects_data(:, subject_ind,:) = squeeze(rms(stable_average_shape(channels_num,:, subject_ind, :),1));
    end
    
%     close all
    %% Global results (average of all the subjects)
    subject_list_short = 1:length(subject_list);
%     figure('units','normalized','outerposition',[0 0 1 1])
%     subplot(2,5,6:10)
    Diff=mean_subjects_data(:,subject_list_short,1)- mean_subjects_data(:,subject_list_short,2); 

    % We use bootstrap and compute if there is a significant difference
    % between conditions.
    perc = 0.05;
    perc2 = 0.01;
    dataB=bootstrap(Diff'); 
    s=findSigDiff(dataB, perc);
    s2=findSigDiff(dataB, perc2);
    k1 = 1;
    k2 = 2;
    
    figure;
    plot(subject_data.time, squeeze(mean(mean_subjects_data(:,subject_list_short,:),2)), 'Linewidth',3)
    hold on
    plot(subject_data.time, k1*abs(s),'Linewidth', 3);
    hold on
    plot(subject_data.time, k2*abs(s2),'Linewidth', 3);
    xlabel('Time (ms)')
    ylabel('RMS magnitude (fT)')
    title(sprintf('All subjects. Baseline: %s', baseline));
    legend('RAND','REG','p=0.05','p=0.01')
%     legend('RAND','REG','p=0.01')
    grid on
   ylim([0,20])
%      ylim([0,35])
%     subject_data.label = timelock275.label;
%     subject_data.avg = mean(mean(stable_average_shape274(:,:,subject_list_short,:),4),3);
%     for ind_T = 1:5
%         switch ind_T
%             case 1
%                 T_plot = [-50,0];
%             case 2
%                 T_plot = [0,50];
%             case 3
%                 T_plot = [50,100];
%             case 4
%                 T_plot = [100,150];
%             case 5
%                 T_plot = [150,200];
%         end
%         subplot(2,5,ind_T)
%         cfg = [];
%         cfg.parameter = 'avg';
%         cfg.layout='CTF275.lay';
%         cfg.xlim=[T_plot(1), T_plot(2)]';
%         cfg.marker = 'none';
%         cfg.interactive = 'yes';
%         cfg.colorbar = 'no';
%         ft_topoplotER(cfg, subject_data); title (sprintf('%dms -- %dms', T_plot(1), T_plot(2)));
% 
%     end
%     
%     %% save the data
%     if store_output == 1
%         mkdir(fullfile('..','Results',out_folder,'ToneTopography'))
%         savefig(fullfile('..','Results',out_folder,'ToneTopography',sprintf('GLOBAL_%s_Topography_%.2f-T%.2f_HP%d_LP%d_BLwindow_%s.fig', baseline, T_init, T_end,hpfreq,lpfreq,bl_tone)));
%     end
%     
%     if store_output == 1
%         mkdir(fullfile('..','Results',sprintf('Global_trigger%s', output_appendix)))
%         save(fullfile('..','Results',sprintf('Global_trigger%s', output_appendix),sprintf('GLOBAL_T%.2f-T%.2f_HP%d_LP%d.mat', T_init, T_end, hpfreq, lpfreq)), 'aux', 'tone_aux');
%       
%     end
    
end
