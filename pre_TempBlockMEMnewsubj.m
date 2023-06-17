 function data_store = pre_TempBlockMEMnewsubj(trigger_list_new, subject, config, trigger_list_old)
    % Function designed to group the trial information from blocks for a
    % subject/trigger combo. In particular, this function applies our
    % standard preprocessing (baselining and filtering), but it also can
    % apply DETRENDING over the whole data file under consideration
    % (block-wise). In order to accelerate the process, we can use the
    % information from the CHANNELS that we obtained in previous analysis
    % (btw for PSD analysis).
    %
    % We also preprocess the information:
    % * Block-wise detrending (if it was not done before).
    % * Filtering (before OR after epoching).
    % * Pre-stim and post-stim epoching depending on the modality that we 
    %       are using.
    % * Data baselined using the pre-stim average value.
    % * Append the trials of all the blocks concatenation.
    % * Visual rejection.
    % **************************************************************
    % trigger: in couples, (5,15) or (10, 20)
    % (New triggers were used (40,80) or (60,100), for simplicity, we
    % replaced the label of triggers with the old triggers after the
    % epoch)
    % Fast sequence
    % * 5 : 3 second RAND sequences.
    % * 15: 3 second REG sequences.
    % Slow sequence
    % * 10: 15 second RAND sequences.  
    % * 20: 15 second REG sequences.
    % **************************************************************
    % subject:
    % * 2-15, the index of the subject, aquired by Antonio.
    % * 16-24, newly aquired by RB and MYH in 2021, this script adapted to this group of
    % subjects
    
    % config: allows for certain configurations.
    %   .out_folder: name of the folder where data will be stored.
    %   .reject_visual: set to 1 if we want to reject trials visually once
    %       that we have grouped them.
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
    % Last update: 25/November/2021: Roberta Bianco, UCL Ear Institute
    % Last update: 07/07/2022: Mingyue Hu, UCL Ear Institute
    
    
    % Some internal paramters.
    out_folder      =           config.out_folder;          % Output data folder
    reject_visual   =           config.reject_visual;       % Set to 1 if we want to reject trials visually.
    store_output    =           config.store_data;          % Set to 1 if we want to store a .mat file with the output block data.
    load_channels   =           config.load_channels;
    channels_path   =           config.channels_path;
    hpfreq          =           config.hpfreq;
    lpfreq          =           config.lpfreq;
    pre_filter      =           config.pre_filter;          % Set to 1 to filter the whole data. Set to 0 to filter on the Epoched data.
    
%     mkdir(fullfile('D:\Results'));
%     mkdir(fullfile('D:\Results',out_folder));
  
    % Input data information.
    % We get the names of the files, considering that are some files with a
    % '.' in the beginning of their name. (Hidden files)
    data_folder = fullfile('D:\MEGGAP','data_longshort', sprintf('Subj%d/',subject));
    filelist_bad = dir([data_folder,'*.ds']);
    counter = 1;
    for ind = 1:length(filelist_bad)
        if filelist_bad(ind).name(1)~= '.'
            filelist{counter} = filelist_bad(ind).name; % This variable contains the names of the MEG data files.
            counter = counter+1;
        end
    end
    clear counter filelist_bad
        
    %% Data reading into memory
    % We load the information of each one of the blocks using a loop.
    % ** WARNING ** 
    % This loop can be integrated with the previous one (pre_Detrend function), but we leave it
    % here for the sake of clarity.
    % *************
    
    data_store = []; % We use this variable to concatenate the blocks. 
     for ind = 2:length(filelist) % We start in the second file, since the first one is the LOC one.
        % We check if the detrended files exist. If they don't, we compute
        % them.
        if isdir(fullfile('D:\Results','Detrended',sprintf('SUBJ_%d', subject))) == 0 |...
                exist(fullfile('D:\Results','Detrended',sprintf('SUBJ_%d', subject), sprintf('BLOCK_%d.mat',ind))) == 0
      
            % We load and read the whole MEG file into the work space, no
            % filtering/referencing or epoching are required in this step
            cfg = [];
            cfg.feedback = 'no';
            cfg.channel = 'MEG';  % We capture all the MEG channels.
            cfg.trialdef.triallength = Inf;
            cfg.trialdef.ntrials = 1;
            cfg.dataset = fullfile(data_folder, filelist{ind});  % File to read.
            
            cfg = ft_definetrial(cfg);
            data = ft_preprocessing(cfg);
    
            % Here, we find the events(Trials) which were triggered by our audio trigger channel,
            % Then, we define the time point of the onset of the last event
            % We don't do epoching in this step as we need to apply
            % detrending on continuous data first
 
            %% Find the time point of triggered events and identify the last event
            str = filelist{ind};
            newStr = split(str, '.');
            name = newStr{1};
            
            cfg = []; 
            cfg.headerformat = 'ctf_res4';
            cfg.trigchan = 'UADC013'; 
            cfg.trialdef.eventtype  = 'trigger_up'; % does the onset of the trigger go down(negative ) or up (positive)?
            cfg.dataset =  [data_folder filelist{ind}, '/' name '.res4' ]; % your filename with file extension; MUST be the res4 file for MEG not the .ds folder
            cfg.trialdef.conditionlabel = {'currtrg'}; 
            
            %cfg.dataset = fullfile(data_folder, filelist{ind});  % File to read.
            cfg.channel = 'MEG';  % We capture all the MEG channels.
            % Depending on the trigger we consider a different post-stimuli
            % time.
            events = ft_definetrial_events_MYH(cfg);  % This function was adapted to capture the all events based on the audio channel triggers
         
           
            % If the block is broken for some reason, we don't process it.
            if isempty(data.trial{1}) == 1 
                continue
            end


            % We load the channels that we obtained during the PSD analysis.
            if load_channels == 1
                load(fullfile(channels_path, sprintf('Channels-SUBJ_%d', subject)));
            else
                channels_num = 1:size(data.trial{1},1);
            end

            %% Apply the Detrending on the continuous data 
            
            % Since there exist a peak at the end of the block, we won't
            % consider its information during the detrending process. We simply
            % keep the temporal information from the beggining until just
            % before it occurs. This peak probably has to do with the
            % acquisition process, and it does not interfere with the data
            % under analysis.
            cfg.event = events; % The captured events/trials relative to its time point 
            end_block = cfg.event(end) + 20*data.fsample;

            fs = data.fsample; % sampling rate
            
            mkdir(fullfile('D:\Results','Detrended',sprintf('SUBJ_%d', subject)));
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
            
            save(fullfile('D:\Results','Detrended',sprintf('SUBJ_%d', subject), sprintf('BLOCK_%d',ind)), 'data');
            
        else 
            load(fullfile('D:\Results','Detrended',sprintf('SUBJ_%d', subject), sprintf('BLOCK_%d',ind)), 'data');
           
        end
        
        
 %% Filtering, preprocessing, epoching based on triggers, merge the blocks
 % After detrending, we do filtering, and epoching based on the triggers,
 % and then we merge the data from 6 blocks into one file
 
        if pre_filter == 1
            cfg=[];
        
      
            cfg.lpfilter = 'yes';
            cfg.lpfreq = lpfreq;
            if hpfreq > 0
                cfg.hpfilter = 'yes';
                cfg.hpfreq = hpfreq;
            end
            data = ft_preprocessing(cfg, data);
        end


        for trigger_ind = 1:length(trigger_list_new)
            % Here, we epoch the continuous data in each block into trials after it has been
            % detrended and filtered.
            
            str = filelist{ind};
            newStr = split(str, '.');
            name = newStr{1};
            
            hdr = data.hdr;
            cfg = []; 
            cfg.header = hdr;
            cfg.headerformat = 'ctf_res4';
            cfg.trigchan = 'UADC013';
            %cfg.trialdef.eventtype  = 'UADC013'; 
            cfg.trialdef.eventtype  = 'trigger_up'; % does the onset of the trigger go down(negative ) or up (positive)?
            cfg.dataset =  [data_folder filelist{ind}, '/' name '.res4' ]; % your filename with file extension; MUST be the res4 file for MEG not the .ds folder
            cfg.trialdef.eventvalue = trigger_list_new(trigger_ind);
            cfg.trialdef.conditionlabel = {'currtrg'}; 
            
            %cfg.dataset = fullfile(data_folder, filelist{ind});  % File to read.
            cfg.channel = 'MEG';  % We capture all the MEG channels.
            % Depending on the trigger we consider a different post-stimuli
            % time.
            prestim =  0.2;    %200ms latency before sound onset
            if (trigger_list_new(trigger_ind) == 40) || (trigger_list_new(trigger_ind) == 80)
                cfg.trialdef.prestim    = prestim;  
                cfg.trialdef.poststim   = 4;   % in seconds, sample = 2520
            else
                cfg.trialdef.prestim    = prestim;
                cfg.trialdef.poststim   = 16;  % sample = 9721
            end
            try
                cfg = ft_definetrial_filMEG_RB(cfg, 1);  % second argument 1 if plot triggers, 0 to avoid plot
            catch
               continue
            end
            
            data_short = ft_redefinetrial(cfg, data);
            
            % find Nan trials and remove from analysis
            good_trials = [];
            for tr = 1:length(data_short.trial)
                trial = data_short.trial{tr};
                idx = find(isnan(trial(:,:)));
                if isempty(idx)
                    good_trials = [good_trials tr];
                end
            end
            
            cfg = [];
            cfg.trials = good_trials;
            data_short = ft_selectdata(cfg, data_short); % amend data to only keep good trials
            
            % remap newtriggers with old triggers
            data_short.trialinfo = repmat(trigger_list_old(trigger_ind),length(data_short.trialinfo),1);
            
            % We baseline using the average value of the pre-stimuli
            % information.
            cfg = [];
            cfg.preproc.demean = 'yes'; % Necessary to baseline.
            cfg.preproc.baselinewindow = [-prestim 0];% in seconds
            data_short = ft_preprocessing(cfg, data_short);
    %       aux = ft_timelockanalysis([], data_short);
    
            % We low-pass filter the information. In our example, the cutoff
            % freq is on 30 Hz.
            if pre_filter == 0
                cfg=[];
                cfg.lpfilter = 'yes';
                cfg.lpfreq = lpfreq;
                if hpfreq > 0
                    cfg.hpfilter = 'yes';
                    cfg.hpfreq = hpfreq;
                end
                data_short = ft_preprocessing(cfg, data_short);    
            end

    %           cfg=[];
    %         timelock = ft_timelockanalysis(cfg, data_short);
    %         figure;plot(timelock.time, rms(timelock.avg))
            % Once that we have all the data of the block, we concatenate it
            % with the previous ones.
            if isempty(data_store) || length(data_store) < length(trigger_list_new)
                data_store{trigger_ind} = data_short;
            else
                data_store{trigger_ind} = ft_appenddata(cfg, data_store{trigger_ind}, data_short);
            end
            clear data_short
         end
        
    end
    
    
    
    for trigger_ind = 1:length(trigger_list_new)
        %% We can reject some of the trials visually.
        if reject_visual == 1
            cfg.channel = 'all';
            data_store{trigger_ind} = ft_rejectvisual(cfg,data_store{trigger_ind});   
        end

        data_subject = data_store{trigger_ind};
        
        %% We can store the block information in a .mat file.
        if store_output == 1
            if load_channels == 1
                mkdir(fullfile('D:\Results',out_folder,'Preprocessed_data'));
                save(fullfile('D:\Results',out_folder,'Preprocessed_data',sprintf('data_subject-TRIG_%d-SUBJ_%d.mat',trigger_list_old(trigger_ind),subject)),'data_subject')
            else
                cfg = [];


                % We compute the average and other statistics from the trials using
                % ft_timelockanalysis.
                timelock = ft_timelockanalysis(cfg, data_subject);

                mkdir(fullfile('D:\Results',out_folder,'Preprocessed_data_AllChannels'))
                save(fullfile('D:\Results',out_folder,'Preprocessed_data_AllChannels',sprintf('data_subject-TRIG_%d-SUBJ_%d.mat',trigger_list_old(trigger_ind),subject)),'timelock','data_subject','-v7.3')
            end
        end
    
    end
        
end


