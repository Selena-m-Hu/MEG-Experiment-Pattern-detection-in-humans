    
    %% Localiser analysis
    % localiser data preprocessing
    % M100 response plot
    % Channel selection
    % _______Written by Mingyue Hu, 18/08/2022
    clc
    clear all
   %% Global Parameters -- we use the same parameters we used for the Main_temporal analysis
%     subject_list            =   [18:24];                   % [2:13, 15]. Indicate the subjects to analyze.
    subject_list            =   17; %trying just one subject for now
    config.hpfreq           =   0;                     
    config.lpfreq           =   30; 
    config.baseline         =   'silence';              % Options: 'tone', 'silence' (original prestimuli data), 'activity' (2 to 2.5 seconds)
    config.single           =   1;                      % Set to 1 to compute a DSS matrix for each one of the conditions.
    config.out_folder       =   sprintf('Localiser_analysis_PRE_HP%d_LP%d', config.hpfreq, config.lpfreq);  % Output data folder. High passed data (2Hz). Detrended.
    % config.channels_path    =   fullfile('..','Results','Channels_DSS');
    % %channels selected after DSS
    config.fs               =   600;                      % Sampling frequency for MEG is 600, for EEG is 2048
    % fieldtrip_path      =   fullfile('..','fieldtrip'); % Fieldtrip path.  
    fieldtrip_path      =   'C:\Users\Chait Lab\Documents\MATLAB\Toolboxes\fieldtrip-20170831';
    addpath(genpath(fieldtrip_path)); 
    addpath(genpath('..\NoiseTools'))
    save_Localiser = 1; 

   % Local paramters -- specifically for localiser analysis
    out_folder      =           config.out_folder;          % Output data folder
    hpfreq          =           config.hpfreq;
    lpfreq          =           config.lpfreq;
    
    mkdir(fullfile('..','Results'));
    mkdir(fullfile('..','Results',out_folder));
  
    %load the data
    
 for i = 1:length(subject_list)
     
    data_folder = fullfile('..','data_longshort', sprintf('Subj%d/',subject_list(i)));
    filelist_bad = dir([data_folder,'*.ds']);
    counter = 1;
    for ind = 1:length(filelist_bad) 
        if filelist_bad(ind).name(1)~= '.'
            filelist{counter} = filelist_bad(ind).name;     % This variable contains the names of the MEG data files.
            counter = counter+1;
        end
    end
    clear counter filelist_bad
    
   
    %% Define trial and epoching
    loc = 1;  % the first file(block) is the localiser data
    str = filelist{loc};
    newStr = split(str, '.');
    name = newStr{1};
    
    cfg = [];
    cfg.headerformat = 'ctf_res4';
    cfg.trigchan = 'UADC013';
    cfg.conditionlabels = {'loc'};
    if subject_list(i) == 17
    name = 'MG06250_Roberta_20210505_01';
    else 
    end 
    cfg.dataset =  {[data_folder filelist{loc}, '/' name '.res4' ]}; % your filename with file extension; MUST be the res4 file for MEG not the .ds folder
    hdr = ft_read_header(cfg.dataset);
    cfg.header = hdr;
    cfg.channel = 'MEG';  % We capture all the MEG channels.
    cfg.trialdef.eventtype  = 'trigger_up'; % does the onset of the trigger go down(negative ) or up (positive)?
    cfg.trialdef.eventvalue = 6; % your event value just for the localiser tone bip
    cfg.trialdef.conditionlabel = {'loc'};
    cfg.trialdef.prestim    = 0.2;  % before stimulation (sec), positive value means trial begins before trigger
    cfg.trialdef.poststim   = 0.5; % after stimulation (sec) , only use positive value
    [cfg] = ft_definetrial_filMEG(cfg, 1); % our triggers are audio pulses sent through the sound card

     
    % Preprocessing
     cfg = ft_definetrial(cfg);
     data = ft_preprocessing(cfg);
     
     cfg=[];  
     cfg.lpfilter = 'yes';
     cfg.lpfreq = lpfreq;
     if hpfreq > 0
        cfg.hpfilter = 'yes';
        cfg.hpfreq = hpfreq;
     end   
     
     data = ft_preprocessing(cfg,data);
    
    % Visual rejection
    cfg = [];
    clean = ft_rejectvisual(cfg, data);
    
    data=clean; 
    % baseline correction
    cfg = [];
    cfg.baseline =[-0.1 0]; %in seconds
    data = ft_timelockbaseline(cfg,data);
    
    % average across trials
    cfg = [];
    cfg.trials = 'all';
    cfg.covariance         = 'no';
    cfg.covariancewindow   = 'all';
    cfg.keeptrials         = 'no';
    cfg.removemean         = 'no';
    cfg.vartrllength       = 0;
    timelock = ft_timelockanalysis(cfg,data);
    
       
    %% Plot the localiser data
    
    figure(i);
    plot(timelock.time,timelock.avg);
    
    
    if save_Localiser == 1
                save(fullfile('..','Results',out_folder,sprintf('Localiser_timelock-SUBJ_%d.mat',subject_list(i))),'timelock')
    end
    
    
     %% Topography
    %load biosemi64_AB_layout;
% 
%     hold on
%     cfg = [];
%     cfg.parameter = 'avg';
%     %cfg.layout = layout_AB;
%     cfg.layout = 'Biosemi64.lay';
%     cfg.xlim            = [.08 .12];
%     %cfg.zlim            = [-.5 .5];
%     cfg.marker          = 'on';
%     cfg.interactive     = 'yes';
%     cfg.comment         = 'no';
%     cfg.colorbar        = 'no';
%     cfg.highlight       = 'on';
%     cfg.highlightcolor  = 'k';
%     cfg.highlightsymbol = '*';
%     cfg.highlightsize   = 10;
%     cfg.highlightchannel = [mI; nI];
%     cfg.markersymbol    = '*';
%     cfg.interpolatenan  = 'no';
%     cfg.markersize      = 3;
%     subplot(1,3,3)
%     ft_topoplotER(cfg,timelock);   
    
 end
   
   
   
   
   
   
   
   
   
   
   
   