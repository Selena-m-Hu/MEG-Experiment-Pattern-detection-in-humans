function pre_ICAcomputation(trigger, subject, config)
    % Function used to apply ft_timelockanalysis to the data of a subject
    % for a specific trigger. It is useful to reduce computational time,
    % since we precompute all the information for later plots. In order to
    % do that, it uses the pre-processed data. It is computed only over the
    % channels of interest in order to reduce the computational time.
    %
    % trigger: 
    % * 5 : 3 second RAND sequences.
    % * 10: 15 second RAND sequences.
    % * 15: 3 second REG sequences.
    % * 20: 15 second REG sequences.
    %
    % subject:
    % * 2-15, the index of the subject.
    %
    % config: allows for certain configurations.
    %   .out_folder: name of the folder where data will be stored.
    %   .store_data: set to 1 to store in a .mat file the whole data of the
    %       subject once that it has been processed.
    %   .channels_path: path indicating where to find the channel
    %       information file.
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
    % Last update: 18/June/2018
    

    out_folder = config.out_folder;
    store_output = config.store_data;
    channels_path = config.channels_path;
    %%   
    load(fullfile('..','Results',out_folder,'Preprocessed_data', sprintf('Long_timelock-TRIG_%d-SUBJ_%d', trigger ,subject)));

    % We perform artifact removal
    % perform the independent component analysis (i.e., decompose the data)
    cfg        = [];
    cfg.method = 'runica'; % this is the default and uses the implementation from EEGLAB

    components = ft_componentanalysis(cfg, data_subject);
%     %%
%     % plot the components for visual inspection
%     figure
%     cfg = [];
%     cfg.component = 1:20;       % specify the component(s) that should be plotted
%     cfg.layout    = 'CTF275.lay'; % specify the layout file that should be used for plotting
%     cfg.comment   = 'no';
%     ft_topoplotIC(cfg, comp)
%     %%
%     cfg          = [];
% cfg.channel  = [2,10]; % components to be plotted
% cfg.viewmode = 'component';
% cfg.layout   = 'CTF275.lay'; % specify the layout file that should be used for plotting
% ft_databrowser(cfg, comp)
    %% We can store the block information in a .mat file.
    if store_output == 1
        mkdir(fullfile('..','Results',out_folder,'ICA_components'))
        save(fullfile('..','Results',out_folder,'ICA_components',sprintf('ICA-TRIG_%d-SUBJ_%d.mat',trigger,subject)),'components')
    end
    

    
end
