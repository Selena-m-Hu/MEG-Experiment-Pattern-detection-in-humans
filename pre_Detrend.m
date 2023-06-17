function data_store = pre_Detrend(subject, config)
    % This function is designed to DETREND over the whole data file under
    % consideration (block-wise). 
    %
    % subject:
    % * [2:13, 15], the index of the subject. Only one at a time.
    %
    % config: allows for certain configurations.
    %   .out_folder: name of the folder where data will be stored.
    %   .store_data: set to 1 to store in a .mat file the whole data of the
    %       subject once that it has been processed.
    %   .load_channels: set to 1 to load the channels from a file (in our
    %       case, computed in a previous analysis). Set to 0  to use all
    %       the channels available.
    %   .channels_path: path indicating where to find the channel
    %       information file.
    % 
    % OUTPUT FOLDER (depending on .load_channels):
    % * load_channels = 1 : ../Results/***/Preprocessed_data
    % * load_channels = 0 : ../Results/***/Preprocessed_data_AllChannels
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
    % Last update: 06/August/2018
    
    
    % Some internal paramters.
    out_folder      =           config.out_folder;          % Output data folder
    store_output    =           config.store_data;          % Set to 1 if we want to store a .mat file with the output block data.
    load_channels   =           config.load_channels;
    channels_path   =           config.channels_path;
     
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
    data_store = []; % We use this variable to concatenate the blocks.
    for ind = 2:length(filelist) % We start in the second file, since the first one is the LOC one.
        % We check if the detrended files exist. If they don't, we compute
        % them.
        if isdir(fullfile('..','Results','Detrended',sprintf('SUBJ_%d', subject))) == 0 | exist(fullfile('..','Results','Detrended',sprintf('SUBJ_%d', subject), sprintf('BLOCK_%d.mat',ind))) == 0
      
            % We load the whole MEG file, which means that we consider the
            % whole block instead of its trials.
            cfg = [];
            cfg.feedback = 'no';
%           cfg.headerformat = 'ctf_res4';
%           cfg.trigchan = 'UADC013';
            %cfg.trialdef.eventtype  = 'UADC013'; 
%           cfg.trialdef.eventtype  = 'trigger_up'; % does the onset of the trigger go down(negative ) or up (positive)?
%             cfg.headerformat = 'ctf_res4';
%             cfg.trigchan = 'UADC013';
            cfg.channel = 'MEG';  % We capture all the MEG channels.
           
            cfg.trialdef.triallength = Inf;
            cfg.trialdef.ntrials = 1;
            cfg.dataset = fullfile(data_folder, filelist{ind});  % File to read.
            cfg = ft_definetrial(cfg);
            data = ft_preprocessing(cfg);
            hdr = data.hdr;
            cfg.header = hdr;
%            cfg = ft_definetrial(cfg);
           % cfg = ft_definetrial_filMEG_RB(cfg);   %reading triggers
%             data = ft_preprocessing(cfg);
            % If the block is broken for some reason, we don't process it.
            % CASE FOR SUBJECT NUMBER 9.
            if isempty(data.trial{1}) == 1 
                continue
            end


            % We load the channels that we obtained during the PSD analysis.
            if load_channels == 1
                load(fullfile(channels_path, sprintf('Channels-SUBJ_%d', subject)));
            else
                channels_num = 1:size(data.trial{1},1);
            end

            % Since there exist a peak at the end of the block, we won't
            % consider its information during the detrending process. We simply
            % keep the temporal information from the beggining until just
            % before it occurs. This peak probably has to do with the
            % acquisition process, and it does not interfere with the data
            % under analysis.
            end_block = cfg.event(end).sample + 20*data.fsample;


            fs = data.fsample; % sampling rate
        
        
        
              mkdir(fullfile('..','Results','Detrended',sprintf('SUBJ_%d', subject)));
            % We use only the channels of interest.
            x = data.trial{1}(channels_num,1:min(end_block+1, length(data.trial{1}))); % extract data from fieldtrip structure
            x = x';

            % We use Noisetools toolbox to smoth and detrend the whole block 
            % information.       
            x = nt_smooth(x,fs/50,3,1);   % smooth with 1/50 Hz kernel (similar to Low-Pass 50Hz)

            w = ones(size(x,1),1);
            x_dt = nt_detrend(x,10,w);     % remove polynomial fit over entire data

            % We organize the processed data according to the Fieldtrip
            % structure.
            data.trial = {x_dt'};
            data.time =  {data.time{1}(1:min(end_block+1, length(data.trial{1})))};
            data.sampleinfo = [1 min(end_block+1, length(data.trial{1}))];
            data.label = data.label(channels_num) ;
            clear x x_dt w

            if store_output == 1
                save(fullfile('..','Results','Detrended',sprintf('SUBJ_%d', subject), sprintf('BLOCK_%d',ind)), 'data');
            end
        end
    end
end