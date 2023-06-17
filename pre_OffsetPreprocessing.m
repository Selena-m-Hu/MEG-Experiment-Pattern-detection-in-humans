function data_subject = pre_OffsetPreprocessing(trigger, subject, config)
    % This function reads the epoched (and preprocessed) data and changes
    % its temporal reference. Instead of 0 seconds, we set it to be placed
    % in the offset temporal area (15 seconds).
    %    
    % trigger_list:
    % * 5 : 3 second RAND sequences.
    % * 10: 15 second RAND sequences.
    % * 15: 3 second REG sequences.
    % * 20: 15 second REG sequences.
    %
    % subject: we can choose a single subject or a list of them.
    % * 2-15, the index of the subject.
    %
    % config: allows for certain configurations.
    %   .out_folder: name of the folder where data will be stored.
    %   .store_data: set to 1 to store in a .mat file the whole data of the
    %       subject once that it has been processed.
    %   .channels_path: path where we can find the channels.
    %
    % OUTPUT FOLDER:
    % ../Results/***/Offset
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
    

    out_folder      = config.out_folder;
    store_output    = config.store_data;
    channels_path   = config.channels_path;

    %%   
    load(fullfile('..','Results',out_folder,'Preprocessed_data_AllChannels_OLD', sprintf('Long_timelock-TRIG_%d-SUBJ_%d', trigger ,subject)));

    % We need to select the channels for each configuration
    % automatically, since it proved before that a preselection does
    % not work properly.
    load(fullfile(channels_path, sprintf('Channels-SUBJ_%d', subject)));
    
    % Here we redefine the temporal area of interest. We must be careful
    % not to get out of the epoching area (in our case, from -500ms to 18
    % seconds originally). We choose the region from 14 to 16 seconds.
    cfg = [];
    cfg.toilim = [-0.5,16];
    data_new = ft_redefinetrial(cfg,data_subject );
    
    % We baseline the data considering a new reference.
    cfg = [];
    cfg.demean = 'yes'; % Necessary to baseline.
    % The last tone should have been heard at 15s minus 200ms.
    % Consequently, our prestim window will begin just before the tone was
    % played. (15.00-0.25) = 14.75
    cfg.baselinewindow = [14.5 15];% in seconds
    data_short = ft_preprocessing(cfg, data_new);
   
    cfg = [];
    cfg.channel = channels;

    % We compute the average and other statistics from the trials using
    % ft_timelockanalysis.
    timelock = ft_timelockanalysis(cfg, data_short);
  
    %% We can store the block information in a .mat file.
    if store_output == 1
        mkdir(fullfile('..','Results',out_folder,'Offset'))
        save(fullfile('..','Results',out_folder,'Offset',sprintf('Timelock-TRIG_%d-SUBJ_%d.mat',trigger,subject)),'timelock', 'data_new')
    end
  
        
end


