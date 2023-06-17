function [channels, channels_num] = pre_ChannelSelectionORIGINALCODE(trigger_list, subject, config)
    % Function created to determine and store which channels are more relevant
    % depending on their magnitude. It also allows to depict a topography
    % map using temporal information (for 'temporal') or frequency
    % information (for 'occipital'). 
    %
    % trigger: list of triggers to be processed
    % * 5 : 3 second RAND sequences.
    % * 10: 15 second RAND sequences.
    % * 15: 3 second REG sequences.
    % * 20: 15 second REG sequences.
    %
    % subject: index of the subject.
    % * 2-15
    %
    % config: allows for certain configurations.
    %   .out_folder: name of the folder where data will be stored.
    %   .channel_modality: set to 'temporal' to obtain automatically the
    %       optimal channels for the temporal information, or 'occipital' to
    %       get the channels for the occipital region. The 'auto' mode is
    %       used does not set any particular limitation for the subjects
    %       low_t and high_t.
    %   .plot_channels: set to 1 to plot a topography map for the chosen
    %       modality.
    %   .store_data: set to 1 to store in a .mat file the selected channel
    %       names and indexes.
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
    plot_channels       = config.plot_channels;
    store_data          = config.store_data;
    DSS_signals         = config.DSS;
    %% Trigger data averaging
    % In order to select the proper channels, we average the information from all
    % the triggers.
    for trigger_ind = 1:length(trigger_list)
        % We load the data that we stored previously in the block
        % processing stage.
        load(fullfile('..','Results',out_folder, 'Preprocessed_data_AllChannels',sprintf('Long_timelock-TRIG_%d-SUBJ_%d.mat',trigger_list(trigger_ind), subject)))
        
        % We can load the DSS data instead.
        if DSS_signals == 1
            n_components = 3;
            load(fullfile('..','Results',out_folder,'DSS_components','Transformed',sprintf('Xdss-TRIG_%d-SUBJ_%d-COMP_%d.mat',trigger_list(trigger_ind), subject, n_components)),'x_dss');
            data_subject.trial = {};
            data_subject.time = {};
            time_data = (1:size(x_dss,2))/600-0.2;
            for trial_ind = 1:size(x_dss,3)
                data_subject.trial{trial_ind} = x_dss(:, :, trial_ind);
                data_subject.time{trial_ind} = time_data;
                
            end
        end
        aux = data_subject;

        cfg=[];
        cfg.method = 'summary';
        cfg.channel = {'MEG'};
        timelock = ft_timelockanalysis(cfg, aux);
        t0 = 0.2*round(aux.fsample); % This represent the position of t = 0. In our case, it happens after a 0.5 seconds pre-stimuli interval.
        % We will temporarily keep the information from t = [pre_stim_time, 0.3]
        % seconds to average trigger data.
        timelock_trigger(:,:,trigger_ind)=timelock.avg(:,1:t0+0.3*round(aux.fsample));
        
        % If we are processing the occipital information, we store the
        % temporal data once that it gets stable. This only works if we
        % intend to ask for the graphical plot of this information.
        if (strcmp(channel_modality, 'occipital') == 1) & (plot_channels == 1)
            if trigger_list(trigger_ind) == 5 | trigger_list(trigger_ind) == 15
                trial_tail = timelock;
                trial_tail.trial = timelock.avg(:,t0+1*round(aux.fsample):t0+3*round(aux.fsample));
                trial_tail.time = timelock.time(t0+1*round(aux.fsample):t0+3*round(aux.fsample));
                trial_tail.fsample = aux.fsample;
                timelock_data_tail{trigger_ind} = trial_tail;
            else
                trial_tail = timelock;
                trial_tail.trial = timelock.avg(:,t0+5*round(aux.fsample):t0+15*round(aux.fsample));
                trial_tail.time = timelock.time(t0+5*round(aux.fsample):t0+15*round(aux.fsample));
                trial_tail.fsample = aux.fsample;
                timelock_data_tail{trigger_ind} = trial_tail;
            end
        end

    end

    timelock.avg = mean(timelock_trigger,3); % We average the information from all the triggers.
    timelock.time = timelock.time(1:t0+0.3*round(aux.fsample));
    
    %% Channel selection 
    % We setup the temporal intervals independently for each subject.
    % Consequently, we will get a different set of channels for each
    % subject.
    if strcmp(channel_modality,'temporal')
        %% Temporal interval
        switch subject
            case 2 
                low_t       = .05;
                high_t      = .076;
            case 3
%                 low_t = .14;
%                 high_t = .15;
                low_t       = .075;
                high_t      = .093;
            case 4
                low_t       = .06;
                high_t      = .085;
            case 5 
                low_t       = .105;
                high_t      = .12
            case 6
                low_t       = .05;
                high_t      = .06;
            case 7
                low_t       = .09;
                high_t      = .11;
            case 8 
%                 low_t = .21;
%                 high_t = .23;
                low_t       = .048;
                high_t      = .073;
            case 9 
                low_t       = .12;
                high_t      = .16;
            case 10 
                low_t       = .07;
                high_t      = .09;
            case 11
                low_t       = .05;
                high_t      = .06;
            case 12
                low_t       = .09;
                high_t      = .11;  
            case 13
                low_t       = .065;
                high_t      = .085;  
            case 14 
                low_t       = .048;
                high_t      = .055;
            case 15 
                low_t       = .11;
                high_t      = .12;  
            otherwise
                low_t = .09;
                high_t = .11;
        end

        % Data extraction
        M100dat = timelock.avg';
        amps=mean(M100dat((t0+low_t*aux.fsample):(t0+high_t*aux.fsample), :),1);
        [ampsSorted,idx]= sort(amps,2,'descend');
        chnsSorted = timelock.label(idx);

        % De-facto channel selection proccess
        % Here we specify which group of channels a subject is bounded to.
    
        chns_selectedLpos=[];
        chns_selectedRpos=[];
        chns_selectedLneg=[];
        chns_selectedRneg=[];
        leftChansCountPos=0;
        rightChansCountPos=0;
        leftChansCountNeg=0;
        rightChansCountNeg=0;

        chns_occipital = [];
        chansCountOcci = 0;

        for count=1:length(ampsSorted)
            strPos=chnsSorted(count);
            strNeg=chnsSorted(end-count+1);

            % Lookup in Left Hemisphere Positive
            if  subject == 11
                if  ~isempty(strfind(strPos{1},'MLT'))
                    if(leftChansCountPos<10)
                        leftChansCountPos=leftChansCountPos+1;
                        chns_selectedLpos=[chns_selectedLpos strPos];
                    end
                end       
            elseif subject == 2
                if  ~isempty(strfind(strPos{1},'MLO')) | ~isempty(strfind(strPos{1},'MLT'))
                    if(leftChansCountPos<10)
                        leftChansCountPos=leftChansCountPos+1;
                        chns_selectedLpos=[chns_selectedLpos strPos];
                    end
                end
            else
                if  ~isempty(strfind(strPos{1},'ML'))
                    if(leftChansCountPos<10)
                        leftChansCountPos=leftChansCountPos+1;
                        chns_selectedLpos=[chns_selectedLpos strPos];
                    end
                end
            end
            
            % Lookup in Left Hemisphere Negative
            if  subject == 4 | subject == 5 | subject == 8  | subject == 9 |  subject == 11 | subject == 14
                if ~isempty(strfind(strNeg{1},'ML'))
                    if(leftChansCountNeg<10)
                            leftChansCountNeg=leftChansCountNeg+1;
                            chns_selectedLneg=[chns_selectedLneg strNeg];
                    end
                end
            elseif subject == 6 |  subject == 10| subject == 13
                if ~isempty(strfind(strNeg{1},'MLT')) | ~isempty(strfind(strNeg{1},'MLO'))
                    if(leftChansCountNeg<10)
                            leftChansCountNeg=leftChansCountNeg+1;
                            chns_selectedLneg=[chns_selectedLneg strNeg];
                    end
                end
            elseif 0%subject == 3
                if ~isempty(strfind(strNeg{1},'MLP')) 
                    if(leftChansCountNeg<10)
                            leftChansCountNeg=leftChansCountNeg+1;
                            chns_selectedLneg=[chns_selectedLneg strNeg];
                    end
                end
            else
                if ~isempty(strfind(strNeg{1},'MLT'))
                    if(leftChansCountNeg<10)
                        leftChansCountNeg=leftChansCountNeg+1;
                        chns_selectedLneg=[chns_selectedLneg strNeg];
                    end

                end
            end

            % Lookup in Right Hemisphere Positive
            if subject == 3 | subject == 4 | subject == 5 | subject == 6 | subject == 8 | subject == 9 | subject == 11 | subject == 14
                if ~isempty(strfind(strPos{1},'MR')) 
                    if(rightChansCountPos<10)
                        rightChansCountPos=rightChansCountPos+1;
                        chns_selectedRpos=[chns_selectedRpos strPos];
                    end
                end

            else
                if  ~isempty(strfind(strPos{1},'MRT'))
                    if(rightChansCountPos<10)
                        rightChansCountPos=rightChansCountPos+1;
                        chns_selectedRpos=[chns_selectedRpos strPos];
                    end
                end
            end

            % Lookup in Right Hemisphere Negative
            if  subject == 6 | subject == 8 | subject == 9  |  subject == 10
                if ~isempty(strfind(strNeg{1},'MRT'))
                    if(rightChansCountNeg<10)
                        rightChansCountNeg=rightChansCountNeg+1;
                        chns_selectedRneg=[chns_selectedRneg strNeg];
                    end
                end
            elseif subject == 3 | subject == 4
                if ~isempty(strfind(strNeg{1},'MRT')) | ~isempty(strfind(strNeg{1},'MRF'))
                    if(rightChansCountNeg<10)
                        rightChansCountNeg=rightChansCountNeg+1;
                        chns_selectedRneg=[chns_selectedRneg strNeg];
                    end
                end
            else
                if ~isempty(strfind(strNeg{1},'MR'))
                    if(rightChansCountNeg<10)
                        rightChansCountNeg=rightChansCountNeg+1;
                        chns_selectedRneg=[chns_selectedRneg strNeg];
                    end
                end

            end
        end

        chns_selectedL = [chns_selectedLpos chns_selectedLneg];
        chns_selectedR = [chns_selectedRpos chns_selectedRneg];

        chnsL_num=[];
        for count1=1:length(timelock.label)
            for count2=1:length(chns_selectedL)
                if (strcmp(timelock.label{count1},chns_selectedL{count2}) ~= 0)
                    chnsL_num=[chnsL_num count1];
                end
            end
        end

        chnsR_num=[];
        for count1=1:length(timelock.label)
            for count2=1:length(chns_selectedR)
                if (strcmp(timelock.label{count1},chns_selectedR{count2}) ~= 0)
                    chnsR_num=[chnsR_num count1];
                end
            end
        end
        timelock.low_t = low_t;
        timelock.high_t = high_t;
        
        channels_num = [chnsL_num,chnsR_num];
        channels = timelock.label(channels_num); 
        
        % We can plot the resultant channels in a topography map.
        if plot_channels == 1
            pre_ChannelPlotTime(timelock, channels_num)
        end
%         keyboard
    elseif strcmp(channel_modality,'occipital')        
        %% Occipital channel extraction
        chansCountOcci = 1;
        chns_occipital = [];
        chnsOcci_num = [];
        % Here we get the indexes of the occipital channels.
        for count=1:length(timelock.label)
            if  ~isempty(strfind(timelock.label{count},'MLO')) | ~isempty(strfind(timelock.label{count},'MRO')) | ~isempty(strfind(timelock.label{count},'MZO'))
                chnsOcci_num = [chnsOcci_num count];
            end     
        end

        channels = timelock.label([chnsOcci_num]);
        channels_num = chnsOcci_num;
        
        
        config.low_f = 8;
        config.high_f = 12;
        
        % Since we are working with the occipital channels, it might be
        % interesting to analyze the inhibition that should happen in the
        % alpha frequency bands [8, 12] Hz. In order to plot this
        % information over a layout, we can use the condition below.
        
        if plot_channels == 1
            pre_ChannelPlotFrequency(timelock_data_tail, config)
        end
        
        
        
    elseif strcmp(channel_modality,'auto')
        %% Temporal interval
        
        low_t = .09;
        high_t = .11;

        % Data extraction
        M100dat = timelock.avg';
        amps=mean(M100dat((t0+low_t*aux.fsample):(t0+high_t*aux.fsample), :),1);
        [ampsSorted,idx]= sort(amps,2,'descend');
        chnsSorted = timelock.label(idx);

        % De-facto channel selection proccess
        % Here we specify which group of channels a subject is bounded to.
    
        chns_selectedLpos=[];
        chns_selectedRpos=[];
        chns_selectedLneg=[];
        chns_selectedRneg=[];
        leftChansCountPos=0;
        rightChansCountPos=0;
        leftChansCountNeg=0;
        rightChansCountNeg=0;

        chns_occipital = [];
        chansCountOcci = 0;

        for count=1:length(ampsSorted)
            strPos=chnsSorted(count);
                strNeg=chnsSorted(end-count+1);
                % Lookup in Left Hemisphere
                if  ~isempty(strfind(strPos{1},'ML'))
                    if(leftChansCountPos<10)
                        leftChansCountPos=leftChansCountPos+1;
                        chns_selectedLpos=[chns_selectedLpos strPos];
                    end
                end
                if ~isempty(strfind(strNeg{1},'MLT'))
                    if(leftChansCountNeg<10)
                        leftChansCountNeg=leftChansCountNeg+1;
                        chns_selectedLneg=[chns_selectedLneg strNeg];
                    end

                end

                % Lookup in Right Hemisphere
                if  ~isempty(strfind(strPos{1},'MRT'))
                    if(rightChansCountPos<10)
                        rightChansCountPos=rightChansCountPos+1;
                        chns_selectedRpos=[chns_selectedRpos strPos];
                    end
                end
                if ~isempty(strfind(strNeg{1},'MR'))
                    if(rightChansCountNeg<10)
                        rightChansCountNeg=rightChansCountNeg+1;
                        chns_selectedRneg=[chns_selectedRneg strNeg];
                    end
                end

            
        end

        chns_selectedL = [chns_selectedLpos chns_selectedLneg];
        chns_selectedR = [chns_selectedRpos chns_selectedRneg];

        chnsL_num=[];
        for count1=1:length(timelock.label)
            for count2=1:length(chns_selectedL)
                if (strcmp(timelock.label{count1},chns_selectedL{count2}) ~= 0)
                    chnsL_num=[chnsL_num count1];
                end
            end
        end

        chnsR_num=[];
        for count1=1:length(timelock.label)
            for count2=1:length(chns_selectedR)
                if (strcmp(timelock.label{count1},chns_selectedR{count2}) ~= 0)
                    chnsR_num=[chnsR_num count1];
                end
            end
        end
        timelock.low_t = low_t;
        timelock.high_t = high_t;
        
        channels_num = [chnsL_num,chnsR_num];
        channels = timelock.label(channels_num); 
        
        % We can plot the resultant channels in a topography map.
        if plot_channels == 1
            pre_ChannelPlotTime(timelock, channels_num)
        end
%         keyboard
        
        
        
        
    end
    %% Here we store the output channels data.
    channels =ft_channelselection(channels, timelock.label);
    if store_data == 1
        if strcmp(channel_modality,'temporal') == 1 % Temporal case  
            mkdir(fullfile('..','Results',out_folder,'Channels'));
            save(fullfile('..','Results',out_folder,'Channels',sprintf('Channels-SUBJ_%d', subject)),'channels', 'channels_num');
            savefig(fullfile('..','Results',out_folder,'Channels',sprintf('Topography-SUBJ_%d', subject)));
        elseif strcmp(channel_modality,'auto') == 1
            mkdir(fullfile('..','Results',out_folder,'Channels'));
            save(fullfile('..','Results',out_folder,'Channels',sprintf('Channels-auto-SUBJ_%d', subject)),'channels', 'channels_num');
        else % Occipital case
            mkdir(fullfile('..','Results',out_folder,'Channels'));
            save(fullfile('..','Results',out_folder,'Channels','Channels-Occipital'),'channels', 'channels_num');
        end
    end
          
         
    close all
end
