function [] = pre_TriggerPlot(subject, config)
    % Function used to plot trigger signals. With the config parameter we
    % can choose whether we want to plot the LOC signals or the triggers of
    % any other block. Useful to analyze LOC, but also to determine what
    % trigger values have been used in the experimental setup.
    %
    % subject: index of the subject.
    % * 2-15
    % 
    % config: 
    %   .config.localizer. If its set to 1, it plots the localizer signal.
    %       If its value is 0, it plots the triggers for another selected block 
    %       (can be modified in the code below).
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
    % Last update: 05/June/2018
    % Some internal paramters.
    out_folder      =           config.out_folder;  % Output data folder
    localizer       =           config.localizer;   % Set to 1 to use plot the localizer signal. 0 to plot the triggers of other block.
    block           =           2;                  % Block to choose if we don't want to plot LOC.
    
    
    
    % Input data information.
    % We get the names of the files, considering that are some files with a
    % '.' in the beginning of their name. (Hidden files)
    data_folder = fullfile('..','data_longshort', sprintf('Subj%d/',subject));
    filelist_bad = dir([data_folder,'*.ds']);
    counter = 1;
    for ind = 1:length(filelist_bad) 
        if filelist_bad(ind).name(1)~= '.'
            filelist{counter} = filelist_bad(ind).name;     % This variable contains the names of the MEG data files.
            counter = counter+1;
        end
    end
    clear counter filelist_bad
    
    %% Trigger extraction
    % We choose the triggers that we want to plot.
    if localizer == 1
        username = filelist{1};                     % LOC information is stored in the first block (in our case).
    else
        username = filelist{block};                     % We choose another arbitrary block.
    end
    cfg = [];
    cfg.feedback = 'no';
    cfg.channel = 'MEG';                            % Le decimos que queremos capturar todos los canales de MEG.
    cfg.dataset = fullfile(data_folder,username);   % Fichero a leer.
    hdr = ft_read_header(cfg.dataset);
    cfg.header = hdr;
    cfg.headerformat = 'ctf_res4';
    cfg.trialdef.eventtype  = 'UADC013';
    % We obtain the trigger information from the events of the dataset.
    hdr   = ft_read_header(cfg.dataset);
    event = ft_read_event(cfg.dataset);
    
    Fs = hdr.Fs;                                    % Sampling frequency.    
%     tr_channel = 'UPPT001';                         % Label of the trigger channel.
    tr_channel = 'UADC013';   %new trigger channel
    % We create a vector with the triggers and their latency values.
    counter = 1;
    for i = 1:length(event)
        if strcmp(event(i).type, tr_channel)
            sample(counter) = event(i).sample;
            value(counter) = event(i).value;
            counter = counter+1;
        end

    end
    clear event hdr

%% Trigger plotting
    sample_time = sample/Fs;
    stem(sample_time, value)
    xlabel('Time (s)');
    ylabel('Trigger val');
    if localizer == 1
        title(sprintf('Localizer signal'));
    else
        title(sprintf('Triggers for block %d',block))        
    end
    
    
end
