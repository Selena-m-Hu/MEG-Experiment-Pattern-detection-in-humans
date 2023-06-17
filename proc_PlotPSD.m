function psd_data = proc_PlotPSD(trigger, subject, config)
    % Computes the PSD for a subject using a list of triggers. Then, it 
    % plots each one of the PSD's and stores them in individual files.
    % 
    % trigger: list of triggers to use for a certain subject.
    %
    % subject: index representing the subject.
    %
    % config 
    %   .out_folder : path of a folder to store data.
    %   .channel_modality : in our case, it can be 'temporal' or
    %       'occipital', and is related with the field config.channels 
    %       below.
    %   .channels : list of selected channels, whose structure is the one
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
    out_folder          = config.out_folder; 
    channel_modality    = config.channel_modality;
    channels            = config.channels;
    store_data          = config.store_data;
    
    %% Data loading and reduction
    % Trial data of the selected subject and all the triggers.
    for trigger_ind = 1:length(trigger)
        load(fullfile('..','Results',out_folder, 'Preprocessed_data',sprintf('Long_timelock-TRIG_%d-SUBJ_%d.mat',trigger(trigger_ind), subject)))
        aux = data_subject;

        % We keep only the temporal information when the data is stable. That
        % happens approximately 20 cycles after stimuli.
        cfg = [];
        if (trigger(trigger_ind) == 5) | (trigger(trigger_ind) == 15)
            cfg.toilim      = [1 3];        % In this case, the initial time is 20*50ms = 1 second.        
        else
            cfg.toilim      = [5 15];       % In this case, 20*250ms = 5 seconds.        
        end
        data_chunk = ft_redefinetrial(cfg, aux); % This is the data structure that we will use to compute the PSD.

        %% Data concatenation.
        % The following code is used for the 3 second-long signals,
        % in order to concatenate a certain number of them
        % (average_num) and while we keep so many trials as posible.
        f_res = 3072;
        if trigger(trigger_ind) == 5 | trigger(trigger_ind) == 15
            data_chunk.trial = {horzcat(data_chunk.trial{:})};
            data_chunk.time = {(1/aux.fsample:1/aux.fsample:length(data_chunk.trial{1})/aux.fsample)};
        else
            % In this case, we are not supposed to concatenate data.
            % However, it may happen that its temporal length is not long
            % enough. Consequently, we will concatenate some of the trials
            % into a kind of 'super-trial' structure.
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
        % We cut some temporal information from the signal, so it fits the
        % PSD function properly.
        time_cut = round(mod(length(data_chunk.time{1}), aux.fsample));
        for ind = 1:length(data_chunk.trial)
            data_chunk.trial{ind} = data_chunk.trial{ind}(:,floor(time_cut/2):end-ceil(time_cut/2)-1);
            data_chunk.time{ind} = data_chunk.time{ind}(:,1:end-time_cut);
        end
        
        
        cfg             = [];            
        cfg.channel     = 'MEG';
        cfg.method      = 'mtmfft'; % mtmfft
        cfg.output      = 'pow';
        cfg.taper       = 'hanning';
        cfg.channel     = channels; % Preselected channels.
        cfg.pad         = 200;
    %     cfg.foi = 0:.1:30;%linspace(0,30,1024);
        % We get the padding.
%         padding = round(length(data_chunk.trial{1})/aux.fsample);
        
        cfg.foi      = 1/100:1/100:30;%[0,30]; % Bandwith of interest.
%         cfg.foilim = [0,30];
        cfg.keeptrials  = 'no';

        
        % PSD for all the individual trials, and then averaged.
        psd_data{trigger_ind} = ft_freqanalysis(cfg, data_chunk);
        

        
        
    end
    %%
    figure('units','normalized','outerposition',[0 0 1 1])        
    subplot(221)
    ind = 1;
    semilogy(psd_data{ind}.freq, sqrt(rms(psd_data{ind}.powspctrm))*1e15);
    legend('TRIG:5. SHORT - RAND');
    xlabel('Frequency (Hz)');
    ylabel('Magnitude (fT)');
    title(sprintf('PSD. Subject: %d. TRIG: %d.',subject, trigger(ind)))
    xlim([0,30]);
    hold on;
    subplot(223)
    ind = 3;
    semilogy(psd_data{ind}.freq, sqrt(rms(psd_data{ind}.powspctrm))*1e15);
    legend('TRIG:15. SHORT - REG');
    xlim([0,30]);
    title(sprintf('PSD. Subject: %d. TRIG: %d.',subject, trigger(ind)))
    xlabel('Frequency (Hz)');
    ylabel('Magnitude (fT)');


    subplot(222)
    ind = 2;
    semilogy(psd_data{ind}.freq, sqrt(rms(psd_data{ind}.powspctrm))*1e15,'r');
    legend('TRIG:10. LONG - RAND');
    xlabel('Frequency (Hz)');
    ylabel('Magnitude (fT)');
    xlim([0,30]);
    title(sprintf('PSD. Subject: %d. TRIG: %d.',subject, trigger(ind)))
    hold on;
    
    subplot(224)
    ind = 4;
    semilogy(psd_data{ind}.freq, sqrt(rms(psd_data{ind}.powspctrm))*1e15,'r');
    legend('TRIG:20. LONG - REG');
    xlim([0,30]);
    title(sprintf('PSD. Subject: %d. TRIG: %d.',subject, trigger(ind)))
    xlabel('Frequency (Hz)');
    ylabel('Magnitude (fT)');
    
    %% Results storage
    if store_data == 1
        mkdir(fullfile('..','Results',out_folder,'PSD'));
        % We store each psd in a separate file.
        psd_aux = psd_data;
        for ind = 1:length(psd_data)
            psd_data = psd_aux{ind};
            save(fullfile('..','Results',out_folder,'PSD', sprintf('PSD-MOD_%s-TRIG_%d-SUBJ_%d', channel_modality, trigger(ind), subject)),'psd_data','channels');
        end
        
        savefig(fullfile('..','Results',out_folder,'PSD', sprintf('GRAPH-MOD_%s-SUBJ_%d', channel_modality, subject)));
        close all
        pause(1)
    end

end