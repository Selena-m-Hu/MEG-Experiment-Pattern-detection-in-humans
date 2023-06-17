function post_ToneButterfly(trigger_list, subject_list, config)
    % This script is modified by Mingyue Hu, August, 2022    
    % trigger_list: 
    % * [5, 15]     : SHORT sequences, 3 seconds
    % * [10, 20]    : long sequences, 15 seconds
   
    

    out_folder = config.out_folder;
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
%     trigger_length = 0.25; % 250ms.
   

    
    %% We create some matrices with the timelock information.
    SUBJECT_GLOBAL = [2:13,15:24];   %Subjects 2:15, aquired by Antonio %Subjects 16:24 aquired by MYH, RB

    if dss_flag == 1     
           load(fullfile('..','Results',out_folder,'DSS_components','Transformed',...
           sprintf('Xdss-TRIG_%d-%d-COMP_%d',trigger_list(1), trigger_list(2), n_components)));
    else  
    end
    
    for trigger_ind = 1:length(trigger_list)
        for subject_ind = 1:length(subject_list)
            if subject_ind == 1
               subject_aux = 1;  % only one subject
            else
            subject_aux = find(subject_list(subject_ind) == SUBJECT_GLOBAL);% original input
            end 
            % We load the timelock information that we computed in a
            % previous stage.
            
            % load the timelock data
            load(fullfile('..','Results',out_folder,'Preprocessed_data_AllChannels',...
                sprintf('Long_timelock-TRIG_%d-SUBJ_%d.mat',trigger_list(trigger_ind),subject_list(subject_ind))),'timelock')
            
            preTimelock = timelock;
            OriginalTimelock = timelock; 
            if dss_flag == 1   % use DSSed data
                timelock.avg = dss_comp(:, :, subject_aux, trigger_ind);
                timelock.time = (1:size(dss_comp,2))/fs-pre_estim;
            end
            timelock_subject{trigger_ind} = timelock;
            timelock274 = timelock; 
            
            clear timelock
%              load(fullfile('..','Results',out_folder,'Preprocessed_data_AllChannels',...
%              sprintf('Long_timelock-TRIG_%d-SUBJ_%d.mat',trigger_list(trigger_ind),subject_list(subject_ind))), 'data_subject')
%             
            timelock = preTimelock;  %timelock data from the non-DSSed data
            timelock = ft_timelockanalysis([], data_subject); %timelock data from the non-DSSed data
            
            
            %%%%%%%%%%%%what is this part for??? This is Antonio's note%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%             We load the channels that we obtained during the PSD analysis.
%             load(fullfile(channels_path, sprintf('Channels-SUBJ_%d', subject_list(subject_ind))));
%              channels_num = 1:274; % this is for data aquired by Antonio
             
%              channels_num = 1:273; % this is for data aquired by Rb,MYH
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % We load the channels that we obtained during the channel
            % selection based on M100 response
            
                      
            load(fullfile(channels_path, sprintf('Channels-SUBJ_%d', subject_list(subject_ind))));
           figure
           plot(timelock274.time,rms(timelock274.avg(channels_num, :)),'Color','k','Linewidth', 3);
           hold on; plot(timelock.time,rms(timelock.avg(channels_num, :)), 'Color', 'r', 'Linewidth', 3);
            
                 
            counter = 1;
            
            for t_ind = 1:window_size:length(timelock.time)
                % We fill the first one (which will not be used) with
                % zeros, since it will be incomplete. This way, the
                % structure of the data is easier and we consider a 50ms
                % to baseline and 200ms of signal.
                if t_ind == 1
                    trigger_shape274(:,:,counter) = zeros(273,150);
                    trigger_shape(:,:,counter) = zeros(length(channels_num),150);
                else
                    trigger_shape(:,:,counter) = timelock.avg(channels_num,t_ind-.1*fs:t_ind+0.15*fs-1);
                    trigger_shape274(:,:,counter) = timelock274.avg(:,t_ind-.1*fs:t_ind+0.15*fs-1);
                    % Security check (at t_ind == 151). % Should be equal to [-0.05, 0.2] 
%                     timelock.time(t_ind-.1*fs:t_ind+0.15*fs-1) 
                end
                counter = counter + 1;
            end
            
            % Multiply by a constant in order to get the units into
            % femtoTeslas.
            trigger_shape = trigger_shape*1e15;
            trigger_shape274 = trigger_shape274*1e15;

            % We get the value of the initial time (0 seconds).
            t0 = 0.25/(window_size/fs);
            
            % We get the temporal index of the beginning and the end of our
            % data from the global matrix.
            t_stable = max(t0,1) + ((T_init/(window_size/fs)):(T_end/(window_size/fs)));
            stable_average_shape(:,:,subject_ind, trigger_ind)      = mean(trigger_shape(:,:,t_stable),3);
            stable_average_shape274(:,:,subject_ind, trigger_ind)   = mean(trigger_shape274(:,:,t_stable),3);
            
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
                    baseline_data = mean(stable_average_shape(:,1:30,subject_ind, trigger_ind) ,2);
                    baseline_data274 = mean(stable_average_shape274(:,1:30,subject_ind, trigger_ind) ,2);
                    output_appendix = '';
                case 'silence'
                    baseline_data = mean(timelock.avg(channels_num, 1:fs*pre_estim),2)*1e15;
                    baseline_data274 = mean(timelock274.avg(:, 1:fs*pre_estim),2)*1e15;
                    output_appendix = '_silence';
                case 'activity'
                    baseline_data = mean(timelock.avg(channels_num, fs*(.5+2.0):fs*(.5+2.5)),2)*1e15;
                    baseline_data274 = mean(timelock274.avg(:, fs*(.5+2.0):fs*(.5+2.5)),2)*1e15;
                    output_appendix = '_activity';
            end
                    
%             stable_shape(:,:,:,subject_ind, trigger_ind) = trigger_shape - repmat(mean(stable_average_shape(:,1:30,subject_ind, trigger_ind) ,2), 1, 150, size(trigger_shape,3));
            stable_average_shape(:,:,subject_ind, trigger_ind)  =  stable_average_shape(:,:,subject_ind, trigger_ind) - repmat(baseline_data, 1, 150);
            stable_average_shape274(:,:,subject_ind, trigger_ind)  =  stable_average_shape274(:,:,subject_ind, trigger_ind) - repmat(baseline_data274, 1, 150);
            
            clear trigger_shape
  
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
    figure('units','normalized','outerposition',[0 0 1 1])
    subplot(2,5,6:10)
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
%     ylim([0,12])
%     ylim([0,25])
    subject_data.label = timelock274.label;
    subject_data.avg = mean(mean(stable_average_shape274(:,:,subject_list_short,:),4),3);
    for ind_T = 1:5
        switch ind_T
            case 1
                T_plot = [-50,0];
            case 2
                T_plot = [0,50];
            case 3
                T_plot = [50,100];
            case 4
                T_plot = [100,150];
            case 5
                T_plot = [150,200];
        end
        subplot(2,5,ind_T)
        cfg = [];
        cfg.parameter = 'avg';
        cfg.layout='CTF275.lay';
        cfg.xlim=[T_plot(1), T_plot(2)]';
        cfg.marker = 'none';
        cfg.interactive = 'yes';
        cfg.colorbar = 'no';
        ft_topoplotER(cfg, subject_data); title (sprintf('%dms -- %dms', T_plot(1), T_plot(2)));

    end
    
    %%
    if store_output == 1
        mkdir(fullfile('..','Results',out_folder,'ToneTopography'))
        savefig(fullfile('..','Results',out_folder,'ToneTopography',sprintf('GLOBAL_%s_Topography_%.2f-T%.2f_HP%d_LP%d.fig', baseline, T_init, T_end,hpfreq,lpfreq)));
    end
%     
%     if store_output == 1
%         mkdir(fullfile('..','Results',sprintf('Global_trigger%s', output_appendix)))
%         save(fullfile('..','Results',sprintf('Global_trigger%s', output_appendix),sprintf('GLOBAL_T%.2f-T%.2f_HP%d_LP%d.mat', T_init, T_end, hpfreq, lpfreq)), 'aux', 'tone_aux');
%       
%     end


    
end
