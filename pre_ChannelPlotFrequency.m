function [] = pre_ChannelPlotFrequency(temporal_data, config)
    % This function plots a topography map showing the data for all the MEG
    % sensors using a certain FREQUENCY BANDWIDTH, whose limit values can 
    % be set below.
    % 
    % temporal_data: temporal information in the structure format from 
    %   FT_PREPROCESSING. 
    %
    % config 
    %   .low_f  : lower bound of our bandwidth of interest.
    %   .high_f : upper bound of our bandwidth of interest.
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
    
    low_freq        = config.low_f;
    high_freq       = config.high_f;
    
    %% Frequency analysis computation
    cfg             = [];
    cfg.channel     = 'MEG';
    cfg.method      = 'mtmfft';     
    cfg.output      = 'pow';
    cfg.taper       = 'hanning';    
    cfg.foi         = 0:.1:30;      % Frequency bandwitdh to compute the PSD, set to [0,30] Hz.
    
    if length(temporal_data) == 4 % We check if we are using the four conditions of our experiment.
        % We concatenate our two SHORT condition data, that come after
        % averaging the trial information. For a better resolution, the
        % information stored in temporal_data should be formed by
        % concatenated trials.
        data_short.trial    = {horzcat(temporal_data{1}.trial, temporal_data{3}.trial)};
        data_short.time     = {1/temporal_data{1}.fsample:1/temporal_data{1}.fsample:length(data_short.trial{1})/temporal_data{1}.fsample};
        data_short.label    = temporal_data{1}.label;
        
        % Here we concatenate the two LONG condition data. Notice that we
        % will plot the topographies of LONG and SHORT in different graphs.
        data_long.trial     = {horzcat(temporal_data{2}.trial, temporal_data{4}.trial)};
        data_long.time      = {1/temporal_data{2}.fsample:1/temporal_data{2}.fsample:length(data_long.trial{1})/temporal_data{2}.fsample};
        data_long.label     = temporal_data{2}.label;
        
        % We perform the frequency analysis (PSD).
        freq_short          = ft_freqanalysis(cfg, data_short);
        freq_long           = ft_freqanalysis(cfg, data_long);
        
       
    end

    %% Freq analysis plot
    cfg                     = [];     
    cfg.xlim                = [low_freq,high_freq];  % Frequency range that we want to analyze.
    cfg.parameter           = 'powspctrm';
    cfg.layout              = 'CTF275.lay';
    cfg.interactive         = 'yes';
    cfg.colorbar            = 'yes';
    cfg.highlightfontsize   = 20;
    
    figure('units','normalized','outerposition',[0 0 1 1])        
    subplot(121)
    ft_topoplotER(cfg, freq_short); title ('Alpha response (SHORT)');
    subplot(122)
    ft_topoplotER(cfg, freq_long); title ('Alpha response (LONG)');
end