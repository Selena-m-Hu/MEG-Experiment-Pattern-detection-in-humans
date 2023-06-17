function psd_data = proc_PSD(trigger, subject, config)
    % *********************************************************************
    % WARNING: THIS FUNCTION NEEDS TO BE TESTED ONCE AGAIN TO VERIFY THAT
    % IT STORES THE OUTPUT DATA PROPERLY.
    % *********************************************************************
    % Function in charge of computing and storing the PSD data for a
    % specific subject and trigger. It does not apply any additional
    % computation, and it can be useful to compute conditions whose files
    % are missing or wrong (for some reason).
    % 
    % trigger: value of the trigger representing the modality whose PSD we
    %       want to compute.
    %
    % subject: index representing the subject.
    %
    % config 
    %   .out_folder : path of a folder to store data.
    %   .channel_modality : in our case, it can be 'temporal' or
    %       'occipital', and is related with the field config.channels 
    %       below.
    %   .channels  : list of selected channels, whose structure is the one
    %       provided by FT_CHANNELSELECTION.
    %   .store_data : set to 1 to store data in a .mat file and figures in
    %       a .fig file.
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
    % Last update: 07/June/2018
    out_folder          = config.out_folder;  %'Trigger_analysis_Post16_Pre500';
    channel_modality    = config.channel_modality;
    channels            = config.channels;
    store_data          = config.store_data;
        
    %% Data loading and reduction
    % Trial data of the selected subject.
    load(fullfile('..','Results',out_folder, 'Preprocessed_data',sprintf('Long_timelock-TRIG_%d-SUBJ_%d.mat',trigger, subject)))
    aux = data_subject;
        
    % We keep only the temporal information when the data is stable. That
    % happens approximately 20 cycles after stimuli.
    cfg = [];
    if (trigger == 5) | (trigger == 15)
        cfg.toilim      = [1 3];        % In this case, the initial time is 20*50ms = 1 second.        
    else
        cfg.toilim      = [5 15];       % In this case, 20*250ms = 5 seconds.        
    end
    data_chunk = ft_redefinetrial(cfg, aux); % This is the data structure that we will use to compute the PSD.
    
    %% Data concatenation.

    f_res = 3072;
    if trigger == 5 | trigger == 15
        % The following code is used for the 3 second-long signals,
        % in order to concatenate a certain number of them
        % (average_num) and while we keep so many trials as posible.
        data_chunk.trial = {horzcat(data_chunk.trial{:})};
        data_chunk.time = {(1/aux.fsample:1/aux.fsample:length(data_chunk.trial{1})/aux.fsample)};
    else
        % In this case (LONG), we are not supposed to concatenate data.
        % However, it may happen that its temporal length is not long
        % enough to get a proper amount of nfft elements. Consequently, we 
        % will concatenate some of the trials into a kind of 'super-trial' 
        % structure.
        t_trial = size(data_chunk.trial{1},2);
        if f_res*4 > t_trial
            n_concatenate = ceil(f_res*32/t_trial);
            n_trials = length(data_chunk.trial) - mod(length(data_chunk.trial),n_concatenate);
            counter = 1;
            clear data_aux
            for i = 1:n_concatenate:n_trials
                data_aux{counter} = horzcat(data_chunk.trial{i:i+n_concatenate-1});
                counter = counter + 1;
            end
        end
        data_chunk.trial = data_aux;
        clear data_aux
        [data_aux{1:length(data_chunk.trial)}]  = deal(1/aux.fsample:1/aux.fsample:length(data_chunk.trial{1})/aux.fsample);
        data_chunk.time = data_aux;
    end
                
    %% PSD using fieldtrip
    cfg             = [];            
    cfg.channel     = 'MEG';
    cfg.method      = 'mtmfft'; % mtmfft
    cfg.output      = 'pow';
    cfg.taper       = 'hanning';
    cfg.channel     = channels; % Preselected channels.
%     cfg.foi = 0:.1:30;%linspace(0,30,1024);
    cfg.foi      = linspace(0,30,f_res);%[0,30]; % Bandwith of interest.
    cfg.keeptrials  = 'no';

    % PSD for all the individual trials, and then averaged.
    psd_data = ft_freqanalysis(cfg, data_chunk);
    
    %% Results storage
    if store_data == 1
        mkdir(fullfile('..','Results',out_folder,'PSD'));
        save(fullfile('..','Results',out_folder,'PSD', sprintf('PSD-MOD_%s-TRIG_%d-SUBJ_%d', channel_modality, trigger, subject)),'psd_data','channels');
    end

end