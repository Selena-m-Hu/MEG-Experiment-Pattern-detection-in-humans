%% Tone epoching
% DSS analysis on tone epoch
% Subject list: [2 3 4 5 6 7 8 9 10 11 12 13 15(old) 16 17 18 19 20 21 22 23 24(new)];
% Frequency band: [2-30 Hz]
% Time window of interests [8-14 s], saved as name Time 8-14 as part of the file name
% Time window for control [0-2.5 s], saved as name Time 0-3 as part of the file name
%_________Mingyue Hu, 03/11/2022


%% Cut the signal sequence into individual epoch of 250 ms/150 samples 
%parameter
    clear all;
    clc;
    trigger_list    = [10 20];
%     subject_list    = [2 3 4 5 6 7 8 9 10 11 12 13 15 16 17 18 19 20 21 22 23 24];
    subject_list    = [2];
    %     time_frame      = [0.5 2.5]; % the time window of interest; in seconds
    time_frame      = [8 14]; % the time window of interest; in seconds

    T_init          = time_frame(1);
    T_end           = time_frame(2);
    hpfreq          = 2;
    lpfreq          = 30;
    out_folder      = sprintf('Trigger_analysis_PRE_HP%d_LP%d',hpfreq,lpfreq);
%   baseline        = config.baseline;
    fs              = 600; % Sampling frequency
    pre_estim       = 0.2; % Prestimuli time. In our case, 200ms.
    window_size     = 0.250*fs; % 250 ms, the length of the stimuli (50ms signal + 200ms silence).
    channels_path   = fullfile('..','Results','Trigger_analysis_PRE_HP0_LP30','Channels_DSS');
    in_folder       = 'Trigger_analysis_PRE_HP0_LP30';  %channels are saved in the 0-30Hz folder

    n_components    = 5;  
    addpath('D:\NoiseTools\');
    addpath('D:\fieldtrip-20220707'); 

    dss_computation = 1;  %1: compute the DSS components
    dss_projection = 0; %project the DSS components into the sensor space, and baseline
    dss_plot  = 0;
%% tone DSS components computation

if dss_computation 
    dsscomp = 1;
    cut_epoch = 0;
    toneDSS_computation(trigger_list, subject_list,time_frame,dsscomp,cut_epoch);
end 

%% tone DSS Projection, baseline

if dss_projection    
    
  for trigger_ind = 1:length(trigger_list)               
       
    for subject_ind = 1:length(subject_list)
        
            % Load raw data
            load(fullfile('..','Results',out_folder,'ToneDSS','DSS_components',...
            sprintf('toneData-TRIG_%d_%d_Time_%d_%d-SUBJ_%d.mat',trigger_list(1), trigger_list(2), round(T_init), round(T_end),...
            subject_list(subject_ind))),'x_orig'); 
     
            % Load DSS components
%             load(fullfile('..','Results',out_folder,'ToneDSS','DSS_components',...
%             sprintf('toneDSS-TRIG_%d_%d_Time_%d_%d-SUBJ_%d-COMP_%d.mat',trigger_list(1), trigger_list(2), round(T_init), round(T_end),...
%             subject_list(subject_ind), 274)));     
            load(fullfile('..','Results',out_folder,'ToneDSS','DSS_components',...
            sprintf('toneDSS-TRIG_%d_%d_Time_%d_%d-SUBJ_%d-COMP_%d.mat',trigger_list(1), trigger_list(2), 0, 3,...
            subject_list(subject_ind), 274)));  

            t_cond{1} = 1:z_timelock.samples_cond1;
            t_cond{2} = z_timelock.samples_cond1+(1:z_timelock.samples_cond2);                

            %% Here we back-project the data into the original 274
            % channels, using the number of components that we want to
            % keep.

            raw_data = x_orig{trigger_ind};
            trans_raw_data = permute(raw_data(:,:,:), [2,1,3]);
            c = nt_xcov(z_timelock.avg(:,:,t_cond{trigger_ind}),trans_raw_data); % c is cross-covariance between z(raw data) and x(DSS components)

            tone_dss = nt_mmat(z_timelock.avg(:,1:n_components,t_cond{trigger_ind}),c(1:n_components,:)); % project from component to sensor space, only using the KEEP components
            tone_dss = nt_mat2trial(tone_dss);  %convert the dss projection into the configuration that fieldtrip expects


            %% We save the output matrix into folder

            mkdir(fullfile('..','Results',out_folder,'ToneDSS','DSS_trials')) % data structure is saved as 3D 
            save(fullfile('..','Results',out_folder,'ToneDSS','DSS_trials',...
            sprintf('toneDSS_Trials-TRIG_%d_Time_%d_%d-SUBJ_%d-COMP_%d.mat',trigger_list(trigger_ind),round(T_init),round(T_end),...
            subject_list(subject_ind), n_components)),'tone_dss');


            %convert the dss data into 3 dimension and save it into one cell
            tone_dss_transform = cat(3,tone_dss{:}); 

            %% Put data from all all triggers into one mat file, for further analysis purpose
            tone_dss_all(:, :, subject_ind, trigger_ind) = mean(tone_dss_transform(:,:,:),3); % mean of all trials from raw data after preprocessing

            clear tone_dss c z
            clear z_timelock  


            %% baseline: method(pre-tone activity)

            %load channels
            load(fullfile('..','Results',in_folder,'Channels_DSS',sprintf('Channels-SUBJ_%d',subject_list(subject_ind))),...
            'channels', 'channels_num');

            bl_window = 21:30;  %define the time window you want to use for tone baseline
            bl_tone = '21-30';  %useful for saving file name 
            baseline_data = mean(tone_dss_all(:,bl_window,subject_ind,trigger_ind), 2);
            output_appendix = '_tone';

            BL40_tone_dss_allsubj(:,:,subject_ind, trigger_ind)  =  tone_dss_all(channels_num,:,subject_ind, trigger_ind) - repmat(baseline_data(channels_num), 1, 150);

            clear tone_dss_all
                               
    end
                        
  end
  
  %% Save the baselined DSS data
 save(fullfile('..','Results',out_folder,'ToneDSS','DSS_trials',...
 sprintf('TI_%d-%d_40Channels_allsubj_toneData-TRIG_%d_%d_BL_%s.mat',round(T_init),round(T_end),trigger_list(1), trigger_list(2),output_appendix)),'BL40_tone_dss_allsubj');  
end 

%% Plot the results

if dss_plot
 
 output_appendix = '_tone';
 load(fullfile('..','Results',out_folder,'ToneDSS','DSS_trials',...
 sprintf('TI_%d-%d_40Channels_allsubj_toneData-TRIG_%d_%d_BL_%s.mat',T_init,round(T_end),trigger_list(1), trigger_list(2),output_appendix)),'BL40_tone_dss_allsubj');  




    for subject_ind = 1:length(subject_list)

            subject_data.time = (1:150)/150*250-50;
    %       subject_data.label = timelock.label;
            subject_data.dimord = 'chan_time';       
            %We load and plot the temporal data.
    %       channels_num = 40;
            mean_subjects_data(:, subject_ind,:) = squeeze(rms(BL40_tone_dss_allsubj(:,:, subject_ind, :),1));
    end
    
    mean_subjects_data=mean_subjects_data*1e15; % convert to femotesla
    %% Global results (average of all the subjects)
    subject_list_short = 1:length(subject_list);

    Diff=mean_subjects_data(:,subject_list_short,1) - mean_subjects_data(:,subject_list_short,2); 

    % We use bootstrap and compute if there is a significant difference
    % between conditions
    perc = 0.05;
    perc2 = 0.01;
    dataB=bootstrap(Diff'); 
    s=findSigDiff(dataB, perc);
    s2=findSigDiff(dataB, perc2);
    k1 = 1;
    k2 = 3;
    
    condi = squeeze(mean(mean_subjects_data(:,subject_list_short,:),2));
    REG = condi(:,2);
    RAND = condi(:,1);
    REGall = mean_subjects_data(:,:,2);
    RANall = mean_subjects_data(:,:,1);
    REGstd = std(REGall')/sqrt(size(REGall,2));
    RANstd = std(RANall')/sqrt(size(RANall,2));
    time = 16:150;
    figure;
    shadedErrorBar(subject_data.time(time),RAND(time),RANstd(time),'lineProps','k');
    hold on 
    shadedErrorBar(subject_data.time(time),REG(time),REGstd(time),'lineProps','r');
    xlim([-25,200])
    
    hold on
    plot(subject_data.time(time), k1*abs(s(time)),'Linewidth', 12);
    hold on
    plot(subject_data.time(time), k2*abs(s2(time)),'Linewidth', 12);
    xlabel('Time (ms)')
    ylabel('Magnitude (fT)')
    title ('0.5-2.5 s')
%     title(sprintf('All subjects. Baseline: %s', baseline));
    
        
end










