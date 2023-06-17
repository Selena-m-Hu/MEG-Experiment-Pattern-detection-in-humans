function pre_ComputeTemporal(trigger, subject, config)
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
    % OUTPUT FOLDER: 
    % * ../Results/**/Timelock
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

    % We need to select the channels for each configuration
    % automatically, since it proved before that a preselection does
    % not work properly.
    load(fullfile(channels_path, sprintf('Channels-SUBJ_%d', subject)));
    cfg = [];
    cfg.channel = channels;

    % We compute the average and other statistics from the trials using
    % ft_timelockanalysis.
    timelock = ft_timelockanalysis(cfg, data_subject);
       
    % This function detects if some of the data coming from the
    % pre-processing stage has some NaN. In that case, it pauses the
    % execution.
    for ind = 1:length(data_subject.trial)
        aux = data_subject.trial{ind};
        ind_nan(ind) = sum(isnan(aux(:)));
    end        
    if sum(ind_nan) > 0
        fprintf('NaN IS COMING FOR YOU!\n')
        keyboard
    end
    clear ind_nan      
    
    %% We can store the block information in a .mat file.
    if store_output == 1
        mkdir(fullfile('..','Results',out_folder,'Timelock'))
        save(fullfile('..','Results',out_folder,'Timelock',sprintf('Timelock-TRIG_%d-SUBJ_%d.mat',trigger,subject)),'timelock')
    end
    

    
end
