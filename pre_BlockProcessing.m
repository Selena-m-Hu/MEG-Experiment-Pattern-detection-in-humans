function data_subject = pre_BlockProcessing(trigger, subject, config)
    % Function in charge of merging the information available in all the
    % blocks of a single subject, considering only the trigger that we are
    % interested at. 
    %
    % We also preprocess the information:
    % * Pre-stim and post-stim values depend on the modality that we are
    %       using.
    % * Data baselined using the pre-stim average value.
    % * Low-pass filtering: 30 Hz.
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
    %   .reject_visual: set to 1 if we want to reject trials visually once
    %       that we have grouped them.
    %   .store_data: set to 1 to store in a .mat file the whole data of the
    %       subject once that it has been processed.
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
    
    
    % Some internal paramters.
    out_folder      =           config.out_folder;          % Output data folder
    reject_visual   =           config.reject_visual;       % Set to 1 if we want to reject trials visually.
    store_output    =           config.store_data;          % Set to 1 if we want to store a .mat file with the output block data.
    
    
    mkdir(fullfile('..','Results'));
    mkdir(fullfile('..','Results',out_folder));
    mkdir(fullfile('..','Results',out_folder,'Preprocessed_data'));
    
    
    % Input data information.
    % We get the names of the files, considering that are some files with a
    % '.' in the beginning of their name. (Hidden files)
    data_folder = fullfile('..','data_longshort', sprintf('Subj%d/',subject));
    filelist_bad = dir([data_folder,'*.ds']);
    counter = 1;
    for ind = 1:length(filelist_bad)
        if filelist_bad(ind).name(1)~= '.'
            filelist{counter} = filelist_bad(ind).name; % This variable contains the names of the MEG data files.
            counter = counter+1;
        end
    end
    clear counter filelist_bad
        
    %% Data preprocessing and merging
    % We load the information of each one of the blocks using a loop.
    % ** WARNING ** 
    % This loop can be integrated with the previous one, but we leave it
    % here for the sake of clarity.
    % *************
    data_subject = []; % We use this variable to concatenate the blocks.
    for ind = 2:length(filelist) % We start in the second file, since the first one is the LOC one.
        cfg = [];
        cfg.feedback = 'no';
        cfg.channel = 'MEG';  % We capture all the MEG channels.
        cfg.trialdef.eventtype  = 'UPPT001'; 
        cfg.trialdef.eventvalue = trigger;
        cfg.dataset = fullfile(data_folder, filelist{ind});  % File to read.
        % Depending on the trigger we consider a different post-stimuli
        % time.
        if (trigger == 5) | (trigger == 15)
            cfg.trialdef.prestim    = 0.5; 
            cfg.trialdef.poststim   = 3.5;
        else
            cfg.trialdef.prestim    = 0.5;
            cfg.trialdef.poststim   = 16;
        end
        try
            cfg = ft_definetrial(cfg);
        catch
            continue
        end
        % We baseline using the average value of the pre-stimuli
        % information.
        cfg.demean = 'yes'; % Necessary to baseline.
        cfg.baselinewindow = [-cfg.trialdef.prestim 0];% in seconds
        data = ft_preprocessing(cfg);
        
        % We low-pass filter the information. In our example, the cutoff
        % freq is on 30 Hz.
        cfg=[];
        cfg.lpfilter = 'yes';
        cfg.lpfreq = 30;
        data = ft_preprocessing(cfg, data);
        
        % Once that we have all the data of the block, we concatenate it
        % with the previous ones.
        if isempty(data_subject) 
            data_subject = data;
        else
            data_subject = ft_appenddata(cfg, data_subject, data);
        end
        clear data

        
    end

    % We can reject some of the trials visually.
    if reject_visual == 1
        data_subject = ft_rejectvisual(cfg,data_subject);   
    end
    
    % We can store the block information in a .mat file.
    if store_output == 1
        save(fullfile('..','Results',out_folder,'Preprocessed_data',sprintf('Long_timelock-TRIG_%d-SUBJ_%d.mat',trigger,subject)),'data_subject')
    end
        
end


