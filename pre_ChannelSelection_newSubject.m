function [channels, channels_num] = pre_ChannelSelection_newSubject(trigger_list, subject, config)
    % Function created to determine and store which channels are more relevant
    % depending on their magnitude. It also allows to depict a topography
    % map using temporal information (for 'temporal') or frequency
    % information (for 'occipital'). 
    % This script is adapted for new subjects aquired by RB and MYH
    % trigger: list of triggers to be processed
    % * 5 : 3 second RAND sequences.
    % * 10: 15 second RAND sequences.
    % * 15: 3 second REG sequences.
    % * 20: 15 second REG sequences.
    %--------------------------------------------------------------------------
    % The time window identified for M100 is based on the average M100 response
    % of all triggers  
    % subject: index of the subject.
    % * 16-24  New subjects
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
    % Last update: 07/07/2022 by Mingyue Hu, Ear Insitute 

    out_folder          = config.out_folder;  
    channel_modality    = config.channel_modality;
    plot_channels       = config.plot_channels;
    store_data          = config.store_data;
    dssLoad             = config.DSS;
    n_components        =config.n_components;
    %% Trigger data averaging
    % In order to select the proper channels, we average the information from all
    % the triggers (Four triggers in this experiment).
    
             
       
      for trigger_ind = 1:length(trigger_list)

           if dssLoad
          %Load DSSed clean data
           load(fullfile('..','Results',out_folder,'DSS_components','normalizedDSS_cfg',...
           sprintf('Xdss-TRIG_%d-SUBJ_%d-SINGLE-Clean-COMP_%d.mat',trigger_list(trigger_ind), subject, n_components)),'dss_data_subject');     
           data_subject = dss_data_subject;
          else    
          %Load raw short timelock data
           load(fullfile('..','Results',out_folder, 'Preprocessed_data_AllChannels','short_timelock',...
           sprintf('Short_timelock-TRIG_%d-SUBJ_%d.mat',trigger_list(trigger_ind), subject)),'short_data_subject');
           data_subject = short_data_subject;
          end  

            cfg=[];
            cfg.method = 'summary';
            cfg.channel = {'MEG'};
            cfg.vartrllength       = 2;
            timelock = ft_timelockanalysis(cfg, data_subject);

         %% We temporarily keep the information from t = [pre_stim_time, 0.3]
            data_subject.fsample = 600;
            aux = data_subject; 
            t0 = 0.2*round(data_subject.fsample); % This represent the position of t = 0. The trial was epoched from -0.2 sec
            % We will temporarily keep the information from t = [pre_stim_time, 0.3]
            % The assumption of M100 response occurs around 80ms to 120ms, so
            % we keep the time interval from -0.2 to 0.3 sec      
            timelock_trigger(:,:,trigger_ind)=timelock.avg(:,1:t0+0.3*round(data_subject.fsample));


            %% If we are processing the occipital information, we store the
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
        
    %% We average the information from all the triggers.
    timelock.avg = mean(timelock_trigger,3); 
    timelock.time = timelock.time(1:t0+0.3*round(data_subject.fsample));
    
    if subject == 20   % we analyse subject 20 with localiser data
     load('D:\Antonio\Results\Localiser_analysis_PRE_HP0_LP30\Localiser_timelock-SUBJ_20.mat','timelock');
     timelock.time = timelock.time(1:300); %only want -0.2 to 0.3 sec
    end 
   
    %% Channel selection 
    % The channel selection is based on the activity of M100 response
    % The M100 response usually occurs at 80 - 120 ms after the sound onset
    % We setup the temporal intervals independently for each subject.
    % Consequently, we will get a different set of channels for each
    % subject.
    
    if strcmp(channel_modality,'temporal')
        %% Temporal interval
        % The time window for this section should reference to the M100
        % response, which needs to be manually analysed before this step
        % Run script M100 response to analyse M100
        switch subject
            % Subect 16 to 24 are newly aquired and the time interval was
            % identified by Mingyue Hu
            case 16
                low_t       = .073;
                high_t      = .115;
            case 17
                low_t       = .025;
                high_t      = .063;
            case 18
                low_t       = .061;
                high_t      = .106;
            case 19 
                low_t       = .063;
                high_t      = .098;
            case 20    %we look at the M100 response of localiser data for subject 20
                low_t       = .02;
                high_t      = .058;
            case 21
                low_t       = .08;
                high_t      = .108;
            case 22
                low_t       = .06;
                high_t      = .1;
            case 23
                low_t       = .095;
                high_t      = .121;
            case 24 
                low_t       = .07;
                high_t      = .1;
            otherwise       
                low_t = .09;
                high_t = .11;
        end

        % Data extraction
        M100dat = timelock.avg';
        amps=mean(M100dat((t0+low_t*aux.fsample):(t0+high_t*aux.fsample), :),1); % The value of average response of M100 for each channel
        [ampsSorted,idx]= sort(amps,2,'descend'); % Ranked the value based on the descending order
        chnsSorted = timelock.label(idx); 

        % De-facto channel selection proccess
        % Here we specify which group of channels a subject is bounded to.
    
        chns_selectedLpos=[]; % Left positive channels
        chns_selectedRpos=[]; % Right positive channels
        chns_selectedLneg=[]; % Left negative channels
        chns_selectedRneg=[]; % Right negative channels
        leftChansCountPos=0;   %start from 0
        rightChansCountPos=0; 
        leftChansCountNeg=0;
        rightChansCountNeg=0;

        chns_occipital = [];   % occipital channels
        chansCountOcci = 0;

        for count=1:length(ampsSorted)
            strPos=chnsSorted(count); % start from the most postive channel
            strNeg=chnsSorted(end-count+1);

            %% Lookup in Left Hemisphere Positive
            if  subject == 16 | subject == 19 | subject == 22 | subject == 21 | subject == 23
                if ~isempty(strfind(strPos{1},'ML'))
                    if(leftChansCountPos<10)
                            leftChansCountPos=leftChansCountPos+1;
                            chns_selectedLpos=[chns_selectedLpos strPos];
                    end
                end
           elseif subject == 17 | subject == 20 
                if ~isempty(strfind(strPos{1},'MLT'))      
                    if(leftChansCountPos<10)
                            leftChansCountPos=leftChansCountPos+1;
                            chns_selectedLpos=[chns_selectedLpos strPos];
                    end
                end
            elseif subject == 18 
                if ~isempty(strfind(strPos{1},'MLT'))|~isempty(strfind(strPos{1},'MLF')) 
                    if(leftChansCountPos<10)
                            leftChansCountPos=leftChansCountPos+1;
                            chns_selectedLpos=[chns_selectedLpos strPos];
                    end
                end
            elseif subject == 24
                if ~isempty(strfind(strPos{1},'MLT'))|~isempty(strfind(strPos{1},'MLO'))...
                        |~isempty(strfind(strPos{1},'MLO'))|~isempty(strfind(strPos{1},'MLP'))
                    if(leftChansCountPos<10)
                            leftChansCountPos=leftChansCountPos+1;
                            chns_selectedLpos=[chns_selectedLpos strPos];
                    end

                end
            else 
                if ~isempty(strfind(strPos{1},'MLT'))
                    if(leftChansCountPos<10)
                            leftChansCountPos=leftChansCountPos+1;
                            chns_selectedLpos=[chns_selectedLpos strPos];
                    end
                end
                
            end
            
            %% Lookup in Left Hemisphere Negative
            if  subject == 16 | subject == 19 | subject == 23 
                if ~isempty(strfind(strNeg{1},'MLT'))
                    if(leftChansCountNeg<10)
                            leftChansCountNeg=leftChansCountNeg+1;
                            chns_selectedLneg=[chns_selectedLneg strNeg];
                    end
                end
           elseif subject == 17 | subject == 20 | subject == 21
                if ~isempty(strfind(strNeg{1},'MLT')) | ~isempty(strfind(strNeg{1},'MLO'))...
                  | ~isempty(strfind(strNeg{1},'MLP'))      
                    if(leftChansCountNeg<10)
                            leftChansCountNeg=leftChansCountNeg+1;
                            chns_selectedLneg=[chns_selectedLneg strNeg];
                    end
                end
            elseif subject == 18 | subject == 24
                if ~isempty(strfind(strNeg{1},'MLT'))|~isempty(strfind(strNeg{1},'MLF')) 
                    if(leftChansCountNeg<10)
                            leftChansCountNeg=leftChansCountNeg+1;
                            chns_selectedLneg=[chns_selectedLneg strNeg];
                    end
                end
            elseif subject == 22 
                if ~isempty(strfind(strNeg{1},'MLT'))|~isempty(strfind(strNeg{1},'MLO')) 
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

            %% Lookup in Right Hemisphere Positive
            if subject == 16 | subject == 18 | subject == 19 | subject == 24
                if ~isempty(strfind(strPos{1},'MRT')) | ~isempty(strfind(strPos{1},'MRF'))
                    if(rightChansCountPos<10)
                        rightChansCountPos=rightChansCountPos+1;
                        chns_selectedRpos=[chns_selectedRpos strPos];
                    end
                end
                
             elseif subject == 20 | subject == 22 
                if  ~isempty(strfind(strPos{1},'MRT')) | ~isempty(strfind(strPos{1},'MRO')) 
                    if(rightChansCountPos<10)
                        rightChansCountPos=rightChansCountPos+1;
                        chns_selectedRpos=[chns_selectedRpos strPos];
                    end
                end
                
            elseif subject == 21
                if  ~isempty(strfind(strPos{1},'MRT'))
                    if(rightChansCountPos<10)
                        rightChansCountPos=rightChansCountPos+1;
                        chns_selectedRpos=[chns_selectedRpos strPos];
                    end
                end                  
            elseif subject == 17
                if  ~isempty(strfind(strPos{1},'MR'))
                    if(rightChansCountPos<10)
                        rightChansCountPos=rightChansCountPos+1;
                        chns_selectedRpos=[chns_selectedRpos strPos];
                    end
                end
                elseif subject == 23
                if  ~isempty(strfind(strPos{1},'MRC')) | ~isempty(strfind(strPos{1},'MRT'))...
                        | ~isempty(strfind(strPos{1},'MRF'))
                    if(rightChansCountPos<10)
                        rightChansCountPos=rightChansCountPos+1;
                        chns_selectedRpos=[chns_selectedRpos strPos];
                    end
                end
            end

            %% Lookup in Right Hemisphere Negative
            if  subject == 16 | subject == 19
                if ~isempty(strfind(strNeg{1},'MR'))
                    if(rightChansCountNeg<10)
                        rightChansCountNeg=rightChansCountNeg+1;
                        chns_selectedRneg=[chns_selectedRneg strNeg];
                    end
                end
            elseif subject == 17 | subject == 20
                if ~isempty(strfind(strNeg{1},'MRT')) | ~isempty(strfind(strNeg{1},'MRF'))
                    if(rightChansCountNeg<10)
                        rightChansCountNeg=rightChansCountNeg+1;
                        chns_selectedRneg=[chns_selectedRneg strNeg];
                    end
                end
            elseif subject == 18
                if ~isempty(strfind(strNeg{1},'MRO')) | ~isempty(strfind(strNeg{1},'MRT'))
                    if(rightChansCountNeg<10)
                        rightChansCountNeg=rightChansCountNeg+1;
                        chns_selectedRneg=[chns_selectedRneg strNeg];
                    end
                end
             elseif subject == 21
                if ~isempty(strfind(strNeg{1},'MRO')) | ~isempty(strfind(strNeg{1},'MRT'))...
                 ~isempty(strfind(strNeg{1},'MRP'))        
                   if(rightChansCountNeg<10)
                    rightChansCountNeg=rightChansCountNeg+1;
                    chns_selectedRneg=[chns_selectedRneg strNeg];
                   end
                end
             elseif subject == 22 | subject == 23 | subject == 24
                 if ~isempty(strfind(strNeg{1},'MRO')) | ~isempty(strfind(strNeg{1},'MRT'))|...
                 ~isempty(strfind(strNeg{1},'MRP'))| ~isempty(strfind(strNeg{1},'MRC'))         
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
%         %% Temporal interval
        
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
        
        if dssLoad == 1
            if strcmp(channel_modality,'temporal') == 1 % Temporal case  
                
                mkdir(fullfile('..','Results',out_folder,'Channels_DSS'));
                save(fullfile('..','Results',out_folder,'Channels_DSS',...
                    sprintf('Channels-SUBJ_%d', subject)),'channels', 'channels_num');
                
                savefig(fullfile('..','Results',out_folder,'Channels_DSS',...
                    sprintf('Topography-SUBJ_%d', subject)));
                
            elseif strcmp(channel_modality,'auto') == 1
                mkdir(fullfile('..','Results',out_folder,'Channels_DSS'));
                save(fullfile('..','Results',out_folder,'Channels_DSS',...
                    sprintf('Channels-auto-SUBJ_%d', subject)),'channels', 'channels_num');
            end
        else
            if strcmp(channel_modality,'temporal') == 1 % Temporal case  
                mkdir(fullfile('..','Results',out_folder,'Channels'));
                save(fullfile('..','Results',out_folder,'Channels',...
                 sprintf('Channels-SUBJ_%d', subject)),'channels', 'channels_num');
             
                savefig(fullfile('..','Results',out_folder,'Channels',...
                 sprintf('Topography-SUBJ_%d', subject)));
            elseif strcmp(channel_modality,'auto') == 1
                mkdir(fullfile('..','Results',out_folder,'Channels'));
                save(fullfile('..','Results',out_folder,'Channels',...
                    sprintf('Channels-auto-SUBJ_%d', subject)),'channels', 'channels_num');
            else % Occipital case
                mkdir(fullfile('..','Results',out_folder,'Channels'));
                save(fullfile('..','Results',out_folder,'Channels','Channels-Occipital'),'channels', 'channels_num');
            end
        end
    end
          
         
    close all
end
