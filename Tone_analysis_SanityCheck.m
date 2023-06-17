
   %% sanity check for tone analysis 
    config.hpfreq           =   2;                     
    config.lpfreq           =   30; 
    config.channels_path    =   fullfile('..','Results','Trigger_analysis_PRE_HP0_LP30','Channels'); 
    trigger_list            =   [10, 20]; 
    subject_list = [2 3 4 5 6 7 8 9 10 11 12 13 15 16 17 18 19 20 21 22 23 24];

      
    load('D:\Antonio\Results\Trigger_analysis_PRE_HP2_LP30\tone_analysis\Tone_40Channels_allSub-SINGLE-COMP_3_BL__tone_50ms');

    load('D:\Antonio\Results\Trigger_analysis_PRE_HP2_LP30\tone_analysis\Tone_275Channels_allSub-SINGLE-COMP_3_BL__tone_50ms');
    
    load('D:\Antonio\Results\Trigger_analysis_PRE_HP2_LP30\Preprocessed_data_AllChannels\short_timelock\Short_timelock-TRIG_10-SUBJ_2');
    timelock = ft_timelockanalysis([],short_data_subject);             

    %% Subject curve-topography graphs
    for subject_ind = 1:length(subject_list)
        subject_data.avg = mean(stable_average_shape(:,:, subject_ind,:),4);
        
        subject_data.time = (1:150)/150*250-50;
        subject_data.label = timelock.label;
        subject_data.dimord = 'chan_time';
        
        
        % We load and plot the temporal data.
        channels_num = 1:40;
        mean_subjects_data(:, subject_ind,:) = squeeze(rms(stable_average_shape(channels_num,:, subject_ind, :),1));
    end
    
    
    %% Global results (average of all the subjects)
    subject_list_short = 1:length(subject_list);
    subplot(2,5,6:10)
    Diff=mean_subjects_data(:,subject_list_short,1)- mean_subjects_data(:,subject_list_short,2); 
    timewindow = (1:30);
    Diff =Diff(timewindow,:);
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
    plot(subject_data.time(timewindow), squeeze(mean(mean_subjects_data(timewindow,subject_list_short,:),2)*100), 'Linewidth',3)
    hold on
    plot(subject_data.time(timewindow), k1*abs(s),'Linewidth', 3);
    hold on
    plot(subject_data.time(timewindow), k2*abs(s2),'Linewidth', 3);
    xlabel('Time (ms)')
    ylabel('RMS magnitude (fT)')
    title(sprintf('All subjects. Baseline: %s', 'tone'));
    legend('RAND','REG','p=0.05','p=0.01')
%     legend('RAND','REG','p=0.01')
    grid on
%     ylim([0,12])
     ylim([0,35])
    subject_data.label = timelock275.label;
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
        ft_topoplotER(cfg, subject_data); title(sprintf('%dms -- %dms', T_plot(1), T_plot(2)));

    end
   
