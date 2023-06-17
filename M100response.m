
%% Ploting the time series data of all channels for M100 response

    % This function plots a topography map showing the data for all the MEG
    % sensors using a time interval specifically for M100 response, this
    % time interval has to be manipulated based on its topography and
    % time-series response
    % The data during the tigger interval are averaged accross all conditions[5,15;10,20]
    % For subject 16 to subject 24
    %_____________________Adapted by Mingyue Hu, Aug,2022
    
    clear all;
    clc;
    config.hpfreq           =   0;                     
    config.lpfreq           =   30; 
    config.out_folder       =   sprintf('Trigger_analysis_PRE_HP%d_LP%d', config.hpfreq, config.lpfreq);  % Output data folder. High passed data (2Hz). Detrended.
    config.channels_path    =   fullfile('..','Results','Trigger_analysis_PRE_HP0_LP30','Channels'); 
    trigger_list            =   [5,15,10,20]; % we analyse data from all trigger conditions
    subject = 23;
    dssLoad = 1;       %1:load DSSed data 0:load raw data
    out_folder = config.out_folder;
    n_components = 3; 
    fieldtrip_path      =   'C:\Users\Chait Lab\Documents\MATLAB\Toolboxes\fieldtrip-20170831';
    addpath(genpath(fieldtrip_path)); 
 
    %% Trigger data averaging
    % In order to select the proper channels, we average the information from all
    % the triggers.
for subject_ind = 1:length(subject_ind)
       
    for trigger_ind = 1:length(trigger_list)
      
      if dssLoad
      %Load DSSed clean data
       load(fullfile('..','Results',out_folder,'DSS_components','normalizedDSS_cfg',...
       sprintf('Xdss-TRIG_%d-SUBJ_%d-SINGLE-Clean-COMP_%d.mat',trigger_list(trigger_ind), subject(subject_ind), n_components)),'dss_data_subject');     
       data_subject = dss_data_subject;
      else    
      %Load raw short timelock data
       load(fullfile('..','Results',out_folder, 'Preprocessed_data_AllChannels','short_timelock',...
       sprintf('Short_timelock-TRIG_%d-SUBJ_%d.mat',trigger_list(trigger_ind), subject(subject_ind))),'short_data_subject');
       data_subject = short_data_subject;
      end  
        
        cfg=[];
        cfg.method = 'summary';
        cfg.channel = {'MEG'};
        cfg.vartrllength       = 2;
        timelock = ft_timelockanalysis(cfg, data_subject);
     %% We temporarily keep the information from t = [pre_stim_time, 0.3]
        data_subject.fsample = 600;
        
        t0 = 0.2*round(data_subject.fsample); % This represent the position of t = 0. The trial was epoched from -0.2 sec
        % We will temporarily keep the information from t = [pre_stim_time, 0.3]
        % The assumption of M100 response occurs around 80ms to 120ms, so
        % we keep the time interval from -0.2 to 0.3 sec      
        timelock_trigger(:,:,trigger_ind)=timelock.avg(:,1:t0+0.3*round(data_subject.fsample));
 
    end

    timelock.avg = mean(timelock_trigger,3); % We average the information from all the triggers.
    timelock.time = timelock.time(1:t0+0.3*round(data_subject.fsample));
   
   %% Plot all channels and identify the time window of M100
    figure(1);
    plot(timelock.time,timelock.avg*1E15)
    
    title(['subject' num2str(subject)])
    
    lowT = 0.095;
    highT= 0.121;
    %% Plot topography 
    figure(2)
    cfg = [];
    cfg.parameter = 'avg';
    cfg.layout = 'CTF275.lay';
    cfg.xlim            = [lowT, highT]';
    cfg.marker          = 'labels';
    cfg.interactive     = 'yes';
    cfg.comment         = 'no';
    cfg.colorbar        = 'yes';
    cfg.highlight       = 'off';
    cfg.highlightcolor  = 'k';
    cfg.highlightsymbol = '.';
    cfg.highlightsize   = 10;
%   cfg.highlightchannel = [mI; nI];
    cfg.markersymbol    = '.';
    cfg.interpolatenan  = 'no';
    cfg.markersize      = 10;
    ft_topoplotER(cfg,timelock);   
    
    dim = [.2 .5 .3 .3];
    str = ['time interval(in sec): ' num2str(lowT) ' - ' num2str(highT)];
    annotation('textbox',dim,'String',str,'FitBoxToText','on');

end 

load('D:\Antonio\Results\Localiser_analysis_PRE_HP0_LP\30Localiser_timelock-SUBJ_20')
    
    
    
    
