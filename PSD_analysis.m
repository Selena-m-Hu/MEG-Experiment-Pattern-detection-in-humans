    function [] = PSD_analysis(trigger, store, filenum)
    %{
        trigger: 5 and 15 for 3 second files.
        trigger: 10 and 20 for 15 second files.
    %}

    dataset = 2; % 2 - Long dataset ; 1 - Short dataset;

    out_folder = 'Trigger_analysis_Post16_Pre500';
    erroneous_block = [];
    if (abs(store)== 1) | (abs(store) == 3)
        %{ 
            Section used to read the files and select the trials to keep. We store them in
            an individual file for each subject.            
        %}
        cfg = [];
        cfg.feedback = 'no';
        cfg.channel = 'MEG';


        error_counter = 1;
    %     trigger = [5]; % Trigger signal to read: (5,15) for short events - (10,20) for long events.
        subjects = [];
        blocks = [];
        for ind_subjects = 1:length(filenum)
           switch filenum(ind_subjects)
                case 2
                    n_blocks = 6;
                case 9 
                    n_blocks = 8;
               case  10
                    n_blocks = 6; 
                otherwise
                    n_blocks = 7;
            end
            [subjects_aux, blocks_aux ] = ndgrid( filenum(ind_subjects), 2:n_blocks);
            subjects = [subjects; subjects_aux(:)];
            blocks = [blocks; blocks_aux(:)];
        end
                 
        for trigger_ind = 1:length(trigger)   
            data_subject = [];
            
            if store == -1
                store_counter = 0;
                data_subject = [];
            end

            counter = 1;
            data_out = [];
            for ind_subjects = 1:length(subjects)
                if ind_subjects == 1
                    sub_num = subjects(ind_subjects);
                end
                 
                if sub_num ~= subjects(ind_subjects)
                    aux = ft_rejectvisual(cfg,data_subject);
                    aux_sum = aux.trial{1};
                    for ind_aux = 2:length(aux.trial)
                        aux_sum = aux_sum + aux.trial{ind_aux};
                    end
                    len_sum = length(aux.trial);
                    mkdir(fullfile('..','Results',out_folder));
                    save(fullfile('..','Results',out_folder,sprintf('Long_timelock-DATA_%d-TRIG_%d-SUBJ_%d.mat',dataset,trigger,sub_num)),'aux_sum','len_sum')
                     
                    clear aux len_sum aux_sum
                    data_subject = [];

                    sub_num = subjects(ind_subjects);
                end
                
                data_block = Trigger_reader(dataset, subjects(ind_subjects), blocks(ind_subjects), trigger(trigger_ind));
%                 counter = counter +1;
                if iscell(data_block.trial) == 0
                    error(error_counter) = ind_subjects;
                    error_counter = error_counter +1;
                else                    
                    if isempty(data_subject) 
                        data_subject = data_block;
                    else
                        data_subject = ft_appenddata(cfg, data_subject, data_block);
                    end
                end
                
                if  ind_subjects == length(subjects) % For the last iteration...
                                      
                    aux = ft_rejectvisual(cfg,data_subject);
                    aux_sum = aux.trial{1};
                    for ind_aux = 2:length(aux.trial)
                        aux_sum = aux_sum + aux.trial{ind_aux};
                    end
                    len_sum = length(aux.trial);
                    if filenum == 8
                        save(fullfile('..','Results',out_folder,sprintf('Long_timelock-DATA_%d-TRIG_%d-SUBJ_%d.mat',dataset,trigger,sub_num)),'aux_sum','len_sum','aux')
                    else
                        save(fullfile('..','Results',out_folder,sprintf('Long_timelock-DATA_%d-TRIG_%d-SUBJ_%d.mat',dataset,trigger,sub_num)),'aux_sum','len_sum')
                    end
                    %{
                        We store a file with the structure, that will be used
                        during the visual rejection of the channels (not the
                        trials, since that rejection has already been
                        performed).
                    %}                    
                    save(fullfile('..','Results',out_folder,sprintf('BAIT-DATA_%d-TRIG_%d.mat',dataset,trigger)),'aux')
                     
                    clear aux len_sum aux_sum
                    data_subject = [];

                    sub_num = subjects(ind_subjects);
                end
                
                clear data_block 
            end 
            clear data_subject
         


        end
    
        % Extraemos algunos estadísticos de los trials. Por ejemplo, el promedio de
        % la magnitud.    
        mkdir(fullfile('..','Results',out_folder));
       
                    

    elseif store== 4
        
        %{ 
            In this section of the code we determine which channels are
            more relevant for our task, and store them for future uses.
            We also show a first plot of the PSD for the subject
            considering the average of all the trials.
        %}
%         addpath(genpath('C:\Users\Chait Lab\Documents\MATLAB\Toolboxes\fieldtrip-20170831'))
        addpath(genpath(fullfile('..','fieldtrip')))
        dataset = 2;
        for trigger_ind = 1:length(trigger)
            load(fullfile('..','Results',out_folder,sprintf('Long_timelock-DATA_%d-TRIG_%d-SUBJ_%d',dataset,trigger(trigger_ind), filenum)))
            
            cfg = [];
            cfg.resamplefs = 120;
            aux = ft_resampledata(cfg, aux);
            
            
            cfg=[];
            cfg.method = 'summary';
            cfg.channel = {'MEG'};
            timelock = ft_timelockanalysis(cfg, aux);
     %% Channel selection for the auditory system
            % Ordenamos las amplitudes, para lo cual escogemos un rango de instantes de
            timelock.fsample = aux.fsample;
            M100dat=timelock.avg'; % Obsérvese que aquí se transpone la media, por lo que tenemos (tiempo x canal).
            t0 = 0.5*timelock.fsample;
            
            amps=mean(M100dat((t0+0.09*timelock.fsample):(t0+0.11*timelock.fsample), :),1);
            switch filenum
                case 2 
                    low_t = .05;
                    high_t = .07;
                case 3
                    low_t = .15;
                    high_t = .165;
                case 4
                    low_t = .07;
                    high_t = .09;
                case 5 
                    low_t = .105;
                    high_t = .12
                case 6
                    low_t = .05;
                    high_t = .06;
                case 7
                    low_t = .09;
                    high_t = .11;
                case 8 
                    low_t = .21;
                    high_t = .23;
                case 9 
                    low_t = .12;
                    high_t = .16;
                case 10 
                    low_t = .07;
                    high_t = .09;
                case 11
                    low_t = .12;
                    high_t = .15;
                case 12
                    low_t = .09;
                    high_t = .11;  
                case 13
                    low_t = .065;
                    high_t = .085;  
                case 14 
                    low_t = .15;
                    high_t = .17;  
                case 15 
                    low_t = .11;
                    high_t = .12;  
                otherwise
                    low_t = .09;
                    high_t = .11;
            end
            amps=mean(M100dat((t0+low_t*timelock.fsample):(t0+high_t*timelock.fsample), :),1);
            
            [ampsSorted,idx]= sort(amps,2,'descend');

            chnsSorted = timelock.label(idx);

            % selecting channels:

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
                % Occipital channel
                 if  ~isempty(strfind(strPos{1},'MLO')) | ~isempty(strfind(strPos{1},'MRO')) | ~isempty(strfind(strPos{1},'MZO'))
                    if(chansCountOcci<length(ampsSorted))
                        chansCountOcci=chansCountOcci+1;
                        chns_occipital=[chns_occipital strPos];
                    end
                 end
                
                 
                 
                
                % Lookup in Left Hemisphere
                if  filenum == 11
                    if  ~isempty(strfind(strPos{1},'MLT'))
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
                
                if  filenum == 4 | filenum == 5 | filenum == 8  | filenum == 9 |  filenum == 11 | filenum == 14
%                     if ~isempty(strfind(strNeg{1},'MLT')) | ~isempty(strfind(strNeg{1},'MLP')) | ~isempty(strfind(strNeg{1},'MLO')) 
                    if ~isempty(strfind(strNeg{1},'ML'))
                        if(leftChansCountNeg<10)
                                leftChansCountNeg=leftChansCountNeg+1;
                                chns_selectedLneg=[chns_selectedLneg strNeg];
                        end
                    end
                elseif filenum == 6 | filenum == 13
                    if ~isempty(strfind(strNeg{1},'MLT')) | ~isempty(strfind(strNeg{1},'MLO'))
                        if(leftChansCountNeg<10)
                                leftChansCountNeg=leftChansCountNeg+1;
                                chns_selectedLneg=[chns_selectedLneg strNeg];
                        end
                    end
                elseif filenum == 3
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

                % Lookup in Right Hemisphere
                if filenum == 3 | filenum == 4 | filenum == 5 | filenum == 6 | filenum == 8 | filenum == 9 | filenum == 11 | filenum == 14
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
                
                
                if filenum == 6 | filenum == 8 | filenum == 9
                    if ~isempty(strfind(strNeg{1},'MRT'))
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

            chns_selectedL = timelock.label(chnsL_num); 
            chns_selectedR = timelock.label(chnsR_num);

            %% Occipital channel selection, used to test visual inhibition.
            chnsOcci_num=[];
            for count1=1:length(timelock.label)
                for count2=1:length(chns_occipital)
                    if (strcmp(timelock.label{count1},chns_occipital{count2}) ~= 0)
                        chnsOcci_num=[chnsOcci_num count1];
                    end
                end
            end
            channel_occipital = ft_channelselection(chns_occipital, timelock.label);

            
            channel = ft_channelselection(unique([chns_selectedL, chns_selectedR]), timelock.label);

            verbose = 0;
            if verbose == 1
                % Timelock plot
                figure(3);
                cfg = [];
                cfg.parameter = 'avg';
                cfg.layout='CTF275.lay';
                cfg.xlim=[low_t, high_t]';
                cfg.marker = 'labels';
                cfg.interactive = 'yes';
                cfg.colorbar = 'yes';

                cfg.markerfontsize = 8;
                cfg.highlight='on';
                cfg.highlightchannel=channel;
                cfg.highlightfontsize=20;
                subplot(121)
                ft_topoplotER(cfg, timelock); title ('M100 response');
                subplot(122)
                t_chunk = timelock.time;
                t_chunk = t_chunk(t_chunk < .5);
%                 plot(t_chunk,timelock.avg([chnsL_num, chnsR_num],t_chunk < .5))
                plot(t_chunk,timelock.avg(:,t_chunk < .5))
                
                figure(2)                 
                ft_topoplotER(cfg, timelock); title ('M100 response');
            end

            if trigger(trigger_ind) == 10
                save(fullfile('..','Results',out_folder,sprintf('Channels_SUBJ-%d.mat', filenum)),...
                            'chns_selectedL','chnsL_num','chns_selectedR','chnsR_num')
            else
                load(fullfile('..','Results',out_folder,sprintf('Channels_SUBJ-%d.mat', filenum)));
            end


            cfg = [];
            if (trigger(trigger_ind) == 5) | (trigger(trigger_ind) == 15)
                cfg.toilim = [1 3];
            else
                cfg.toilim = [5 15]; 
            end
            data_chunk = ft_redefinetrial(cfg, aux);
            channel{trigger_ind} = ft_channelselection(unique([chns_selectedL, chns_selectedR]), data_chunk.label);
            
            %% Data concatenation.
            average_num = 1; % Number of chunks to concatenate.
            counter = 1;
            clear buff_out
            % The following code is used for the 3 second-long signals,
            % in order to concatenate a certain number of them
            % (average_num) and while we keep so many trials as posible.
            if average_num == 1
                % This code can be simplified and removed, but since its
                % already working and being tested it will stay here
                % forever.
                buff_ave = zeros(size(data_chunk.trial{1}));
                for i = 1:length(data_chunk.trial)
                    if (mod(i, average_num) == 1) 
                        buff_ave = data_chunk.trial{i}/average_num;
                    else
                        buff_ave = data_chunk.trial{i}/average_num;
                        if mod(i, average_num) == 0
                            buff_out{counter} = buff_ave;
                            counter = counter + 1;

                        end                    
                    end
                end
                data_aux.trial = {horzcat(buff_out{:})};
%                 data_aux.time = {horzcat(data_chunk.time{1:length(buff_out)})};
                data_aux.time = {0.0033:1/aux.fsample:length(data_aux.trial{1})/aux.fsample};
                
            else               
                % Here we concatenate (temporal axis) a certain number of
                % chunks.
                trial_aux = {};
                counter = 1;
                for i = 1:length(data_chunk.trial)    
                    % We limit the number of chunks concatenated
                    % considering that all of them should have at least
                    % 1200 temporal samples. Consequently, we might remove
                    % the last cluster if its temporal dim. is smaller.
                    if counter <= floor(length(data_chunk.trial)/average_num)
                        if (mod(i, average_num) == 1) 
                            trial_aux{counter} = data_chunk.trial{i};                     
                        else
                            trial_aux{counter} = horzcat([trial_aux{counter},data_chunk.trial{i}]);
                        end

                        if mod(i,average_num) == 0
                            counter = counter + 1;
                        end
                    else
                        % Nothing at all. This occurs when the last cluster
                        % of chunks is way too small.
                    end
                    
                end
                data_aux.trial = trial_aux;
                [data_aux.time{1:length(trial_aux)}] = deal(0.0033:1/aux.fsample:length(data_aux.trial{1})/aux.fsample);
            end

            
            data_aux.label = data_chunk.label;
            if trigger(trigger_ind) == 5 | trigger(trigger_ind) == 15
                data_aux.fsample = aux.fsample;
            end


           
            fs = aux.fsample;
            % PSD using fieldtrip
            cfg = [];            
            cfg.channel = 'MEG';
            cfg.method = 'mtmfft'; % mtmfft
            cfg.output = 'pow';
            cfg.taper = 'hanning'
            cfg.channel = channel{trigger_ind};
%             cfg.foi = 0:.1:30;%linspace(0,30,1024);
            cfg.foi = linspace(0,30,4096);
%             cfg.foilim = [0,30];
            cfg.keeptrials = 'yes';
            
            % PSD for all the individual trials, and then averaged.
            [freq{trigger_ind}] = ft_freqanalysis(cfg, data_chunk)
%             plot(freq_chunk{1}.freq,reshape(rms(freq_chunk{1}.powspctrm,2),4761,1))

            % PSD for a huge signal formed by merged chunks.            
            [freq_chunk{trigger_ind}] = ft_freqanalysis(cfg, data_aux)
            verbose = 0;
            if verbose == 1
                % Freq analysis plot
                figure(1);       
                cfg = [];
                cfg.channel = 'MEG';
                cfg.method = 'mtmfft'; % mtmfft
                cfg.output = 'pow';
                cfg.taper = 'hanning'
                cfg.foi = 0:.1:30;
    
                if (trigger(trigger_ind) == 5) | (trigger(trigger_ind) == 15)
                     % PSD for all the individual trials, and then averaged.
                    [aux_freq] = ft_freqanalysis(cfg, data_aux)
                    cfg = [];    
                    cfg.xlim=[19.9,20.1]';  % Bandwidth of interest for the SHORT sequences.

                    
                else                    
                    % PSD for a huge signal formed by merged chunks.            
                    [aux_freq] = ft_freqanalysis(cfg, data_chunk)
                    cfg = [];    
                    cfg.xlim=[3.9,4.1]'; % Bandwidth of interest for the LONG sequences.
%                     aux_freq = freq{trigger_ind};
                end
                cfg.parameter = 'powspctrm';
                cfg.layout='CTF275.lay';
                                  

                cfg.foi = 0:.1:30;%linspace(0,30,1024);
                
                cfg.interactive = 'yes';
                cfg.colorbar = 'yes';
                cfg.highlightfontsize=20;
                subplot(121)
%                 figure('units','normalized','outerposition',[0 0 1 1])        
                ft_topoplotER(cfg, aux_freq); title ('Frequency response');
                subplot(122)
                t_chunk = aux_freq.freq;
                semilogy(t_chunk,rms(aux_freq.powspctrm([chnsL_num,chnsR_num],:)))
                
                figure(2)    
                cfg.xlim=[8.5,12]';
                ft_topoplotER(cfg, aux_freq); title ('Alpha response');
            end
            
            
            
            %% PSD for occipital channels
            % PSD using fieldtrip
            cfg = [];            
            cfg.channel = 'MEG';
            cfg.method = 'mtmfft'; % mtmfft
            cfg.output = 'pow';
            cfg.taper = 'hanning'
            cfg.channel = channel_occipital;
            cfg.foi =  0:.1:30;
            cfg.keeptrials = 'yes';
    
            
            % PSD for all the individual trials, and then averaged.
            [freq_occipital{trigger_ind}] = ft_freqanalysis(cfg, data_chunk)
            % PSD for a huge signal formed by merged chunks.            
            [freq_chunk_occipital{trigger_ind}] = ft_freqanalysis(cfg, data_aux)
            
            
        end
        

         %% We store the output results. 

        for trigger_ind = 1:length(trigger)            
            output.full_freq = freq;
            output.chunk_freq = freq_chunk;
            output.occipital_full_freq = freq_occipital;
            output.occipital_chunk_freq = freq_chunk_occipital;
            
            
        end
%         close all
        if (trigger(trigger_ind) == 5) | (trigger(trigger_ind) == 15)
            if length(trigger) > 1
                savefig(fullfile('..','Results',out_folder,sprintf('Short_PSD_TOP-channels_SUBJ-%d.fig', filenum)))
            else
                close all
                save(fullfile('..','Results',out_folder,sprintf('PSD_TOP-channels_TRIG-%d_SUBJ-%d.mat',trigger(trigger_ind), filenum)),'output')
            end
        else
            if length(trigger) > 1
                savefig(fullfile('..','Results',out_folder,sprintf('Long_PSD_TOP-channels_SUBJ-%d.fig', filenum)))
            else
                close all
                save(fullfile('..','Results',out_folder,sprintf('PSD_TOP-channels_TRIG-%d_SUBJ-%d.mat',trigger(trigger_ind), filenum)),'output')
            end
        end
    close all
       
    else
%         trigger = [15];
        dataset = 2;
        for trigger_ind = 1:length(trigger)
            counter = 1;
            
            counter_trial = 1;
            for filenum = [2,3,5:13]%2:15
            
                load(fullfile('..','Results',out_folder,sprintf('PSD_TOP-channels_TRIG-%d_SUBJ-%d',trigger(trigger_ind), filenum)))
                
                % Fieltrip for TRIAL-AVG signal.
                
                if trigger(trigger_ind) == 10 | trigger(trigger_ind) == 20
                    [n_trial, n_chann, n_freq] = size(output.full_freq{1}.powspctrm);
                    global_psd(counter,:) = sqrt(rms(reshape(mean(output.full_freq{1}.powspctrm,1), n_chann, n_freq)));
                    global_f(counter,:) = output.full_freq{1}.freq;
                    global_occipital_psd(counter,:) = sqrt(rms(reshape(mean(output.occipital_full_freq{1}.powspctrm,1), n_chann, n_freq)));
                    global_occipital_f(counter,:) = output.occipital_full_freq{1}.freq;
                    
                    fieldtrip_psd(counter_trial:counter_trial+n_trial-1,:) = permute(rms(output.full_freq{1}.powspctrm,2),[1,3,2]);
                    fieldtrip_f(counter_trial:counter_trial+n_trial-1,:) = repmat(output.full_freq{1}.freq, n_trial, 1);
                    occipital_psd(counter_trial:counter_trial+n_trial-1,:) =  permute(rms(output.occipital_full_freq{1}.powspctrm,2),[1,3,2]);
                    occipital_f(counter_trial:counter_trial+n_trial-1,:) = repmat(output.occipital_full_freq{1}.freq, n_trial, 1);
                else
                    [n_trial, n_chann, n_freq] = size(output.chunk_freq{1}.powspctrm);
                    global_psd(counter,:) = sqrt(rms(reshape(mean(output.chunk_freq{1}.powspctrm,1), n_chann, n_freq)));
                    global_f(counter,:) = output.chunk_freq{1}.freq;
                    global_occipital_psd(counter,:) = sqrt(rms(reshape(mean(output.occipital_chunk_freq{1}.powspctrm,1), n_chann, n_freq)));
                    global_occipital_f(counter,:) = output.occipital_chunk_freq{1}.freq;
                    
                    fieldtrip_psd(counter_trial:counter_trial+n_trial-1,:) = permute(rms(output.chunk_freq{1}.powspctrm,2),[1,3,2]);
                    fieldtrip_f(counter_trial:counter_trial+n_trial-1,:) = repmat(output.chunk_freq{1}.freq, n_trial, 1);                
                    occipital_psd(counter_trial:counter_trial+n_trial-1,:) = permute(rms(output.occipital_chunk_freq{1}.powspctrm,2),[1,3,2]);
                    occipital_f(counter_trial:counter_trial+n_trial-1,:) = repmat(output.occipital_chunk_freq{1}.freq, n_trial,1);
                    
                end
%                 semilogy(output.full_freq{1}.freq,reshape(output.full_freq{1}.powspctrm(4,:,:),40,301)')
                subject_vec(counter_trial:counter_trial+n_trial-1) = filenum;
                counter = counter + 1;
                counter_trial = counter_trial + n_trial;
            end
            
            psd_group{trigger_ind} = fieldtrip_psd;
            freq_group{trigger_ind} = fieldtrip_f;
            
            global_psd_group{trigger_ind} = global_psd;
            global_freq_group{trigger_ind} = global_f;
            
            occipital_psd_group{trigger_ind} = occipital_psd;
            occipital_freq_group{trigger_ind} = occipital_f;
            
            global_occipital_psd_group{trigger_ind} = global_occipital_psd;
            global_occipital_freq_group{trigger_ind} = global_occipital_f;
            
            subject_group{trigger_ind} = subject_vec;
            clear global_psd global_f global_occipital_psd global_occipital_f fieldtrip_psd fieldtrip_f occipital_psd occipital_f subject_vec
        end
        

        % Subject plot
        figure(1)
        plot(global_freq_group{1}',global_psd_group{1}')
        xlabel('Frequency (Hz)');
        
        % GLOBAL SHORT
        figure(2)
        ind_aux = 1;
        semilogy(global_freq_group{ind_aux}',mean(global_psd_group{ind_aux})','b','Linewidth',2)
        hold on
        ind_aux = 3;
        semilogy(global_freq_group{ind_aux}',mean(global_psd_group{ind_aux})','r','Linewidth',2)
        xlabel('Frequency (Hz)');
        legend('RAND','REG')
        
        
         % GLOBAL SHORT
        figure(2)
        ind_aux = 1;
        plot(global_occipital_freq_group{ind_aux}',mean(global_occipital_psd_group{ind_aux})','b','Linewidth',2)
        hold on
        ind_aux = 3;
        plot(global_occipital_freq_group{ind_aux}',mean(global_occipital_psd_group{ind_aux})','r','Linewidth',2)
        xlabel('Frequency (Hz)');
        legend('RAND','REG')
        
        
        % OCCIPITAL VS REST
        figure(3)
        ind_aux = 1;
        plot(global_freq_group{ind_aux}',mean(global_psd_group{ind_aux})','b','Linewidth',2)
        hold on
        ind_aux = 3;
        plot(global_occipital_freq_group{ind_aux}',mean(global_occipital_psd_group{ind_aux})','r','Linewidth',2)
        legend('Auditory', 'Occipital')
        

        for trigger_ind = 1:length(trigger)
            if (trigger(trigger_ind) == 10) | (trigger(trigger_ind) == 20)
                ind_signal = [39:40,41,42:43];
                ind_noise = [31:min(ind_signal)-1, max(ind_signal)+1:52];
            else
                ind_signal = [199:200,201,202:203];
                ind_noise = [190:min(ind_signal)-1, max(ind_signal)+1:212];
            end
            
            signal = abs(global_psd_group{trigger_ind}(:,ind_signal));
            noise = abs(global_psd_group{trigger_ind}(:,ind_noise));
            
            signal_trial_like = abs(psd_group{trigger_ind}(:,ind_signal));
            noise_trial_like = abs(psd_group{trigger_ind}(:,ind_noise));

            snr(trigger_ind, :) = (mean(signal,2)./mean(noise,2));   
            snr_trial_like{trigger_ind} = (mean(signal_trial_like,2)./mean(noise_trial_like,2));   

            
        end

        
        
        %% Alpha computations
        ind_signal = [91:122]; % Information from 9 Hz to 12 Hz.
        for trigger_ind = 1:length(trigger)
            signal = abs(global_occipital_psd_group{trigger_ind}(:,ind_signal));
            signal_trial_like = abs(occipital_psd_group{trigger_ind}(:,ind_signal));
            
            alpha_magnitude(trigger_ind,:) = mean(signal,2)*1E15;
            alpha_magnitude_trial_like{trigger_ind} = mean(signal_trial_like,2)*1E15;
            
        end 


        %% One way anova.
        
        anova_snr_short_RANDvsREG = anova_rm( snr([1,3],:)','off');
        anova_snr_long_RANDvsREG = anova_rm(snr([2,4],:)','off');
        
        anova_alpha_short_RANDvsREG = anova_rm(alpha_magnitude([1,3],:)','off');
        anova_alpha_long_RANDvsREG = anova_rm(alpha_magnitude([2,4],:)','off');
        
        anova_snr_RAND_LONGvsSHORT = anova_rm( snr([1,2],:)','off');
        anova_snr_REG_LONGvsSHORT = anova_rm(snr([3,4],:)','off');
        
        anova_alpha_RAND_LONGvsSHORT = anova_rm(alpha_magnitude([1,2],:)','off');
        anova_alpha_REG_LONGvsSHORT = anova_rm(alpha_magnitude([3,4],:)','off');
        
        
        
        %% Two way anova (average)
        snr_vector = [];
        REG_vector = [];
        LONG_vector = [];
        subj_vector = [];
        alpha_vector = [];
        % TRIGGER 5
        snr_vector = snr(1,:)';
        alpha_vector = alpha_magnitude(1,:)';
        REG_vector = zeros(size(snr(1,:)'));
        LONG_vector = zeros(size(snr(1,:)'));
        subj_vector = [1:11]';
        
        % TRIGGER 10
        snr_vector = [snr_vector; snr(2,:)'];
        alpha_vector = [alpha_vector; alpha_magnitude(2,:)'];
        REG_vector = [REG_vector; zeros(size(snr(2,:)'))];
        LONG_vector = [LONG_vector; ones(size(snr(2,:)'))];
        subj_vector = [subj_vector; [1:11]'];
        
        % TRIGGER 15
        snr_vector = [snr_vector; snr(3,:)'];
        alpha_vector = [alpha_vector; alpha_magnitude(3,:)'];
        REG_vector = [REG_vector; ones(size(snr(3,:)'))];
        LONG_vector = [LONG_vector; zeros(size(snr(3,:)'))];
        subj_vector = [subj_vector; [1:11]'];
        
        % TRIGGER 20
        snr_vector = [snr_vector; snr(4,:)'];
        alpha_vector = [alpha_vector; alpha_magnitude(4,:)'];
        REG_vector = [REG_vector; ones(size(snr(4,:)'))];
        LONG_vector = [LONG_vector; ones(size(snr(4,:)'))];
        subj_vector = [subj_vector; [1:11]'];
        
        snr_stats = rm_anova2(snr_vector,subj_vector,REG_vector,LONG_vector,{'REG','LONG'});
        alpha_stats = rm_anova2(alpha_vector,subj_vector,REG_vector,LONG_vector,{'REG','LONG'});

        
       
        %% Two way anova (trial-like)
        snr_vector = [];
        REG_vector = [];
        LONG_vector = [];
        subj_vector = [];
        alpha_vector = [];
        % TRIGGER 5
        for ind = 1:4
            [n_fil, n_col] = size(snr_trial_like{ind});
            switch ind 
                case 1
                    snr_vector = snr_trial_like{ind};
                    alpha_vector = alpha_magnitude_trial_like{ind};
                    REG_vector = zeros(n_fil, n_col);
                    LONG_vector = zeros(n_fil, n_col);
                    subj_vector = subject_group{ind}';
                case 2
                    snr_vector = [snr_vector; snr_trial_like{ind}];
                    alpha_vector = [alpha_vector; alpha_magnitude_trial_like{ind}];
                    REG_vector = [REG_vector; zeros(n_fil, n_col)];
                    LONG_vector = [LONG_vector; ones(n_fil, n_col)];
                    subj_vector = [subj_vector; subject_group{ind}'];
                case 3
                    snr_vector = [snr_vector; snr_trial_like{ind}];
                    alpha_vector = [alpha_vector; alpha_magnitude_trial_like{ind}];
                    REG_vector = [REG_vector; ones(n_fil, n_col)];
                    LONG_vector = [LONG_vector; zeros(n_fil, n_col)];
                    subj_vector = [subj_vector; subject_group{ind}'];
                case 4
                    snr_vector = [snr_vector; snr_trial_like{ind}];
                    alpha_vector = [alpha_vector; alpha_magnitude_trial_like{ind}];
                    REG_vector = [REG_vector; ones(n_fil, n_col)];
                    LONG_vector = [LONG_vector; ones(n_fil, n_col)];
                    subj_vector = [subj_vector; subject_group{ind}'];
            end
        end
        
        snr_stats_trial = rm_anova2(snr_vector,subj_vector,REG_vector,LONG_vector,{'REG','LONG'});
        alpha_stats_trial = rm_anova2(sqrt(alpha_vector),subj_vector,REG_vector,LONG_vector,{'REG','LONG'});
       
        
        
           
        
        
    end
    
    
end



function [data_block] = Trigger_reader(dataset, filenum, block, trigger)
% clear all;
% close all;
addpath(genpath(fullfile('..','fieldtrip')));
if nargin == 2
    block = [];
elseif nargin == 1
    data_folder = dataset.data_folder;
    username = dataset.username;
    dataset = -1;
    filenum = 1;
end
sampling = 600;
switch dataset 
    case 1
        data_folder = fullfile('..','data');
        switch filenum
            case 1
                username = 'ac040981_MChait2_20120417_01.ds';
            case 2  
                username = 'cp071189_MChait2_20120416_01.ds'; %% FILE NOT WORKING
            case 3
                username = 'lt290589_MChait2_20120417_01.ds';
            case 4
                username = 'ns070382_MChait2_20120423_01.ds';
            case 5
                username = 'rs120685_MChait2_20120417_01.ds';
            otherwise
                if filenum < 0
                    username = fullfile('Loc',...
                                        sprintf('Subj%d', abs(filenum)),...
                                        sprintf('Subj%d_loc', abs(filenum)));
                end
        end
    case 2
        data_folder = fullfile('..','data_longshort', sprintf('Subj%d/',filenum));
        filelist_bad = dir([data_folder,'*.ds']);
        contador = 1;
        for ind = 1:length(filelist_bad)
            if filelist_bad(ind).name(1)~= '.'
                filelist{contador} = filelist_bad(ind).name;
                contador = contador+1;
            end
        end
        fprintf('** Number of blocks: %d **\n', length(filelist));
        if block <= length(filelist)
            fprintf('** Reading block %d **\n', block);
            username = filelist{block};
        else
            fprintf('** BLOCK NUMBER OUT OF BOUNDS **\n');
            keyboard
        end
    
end

cfg = [];
cfg.feedback = 'no';
cfg.channel = 'MEG'; % Le decimos que queremos capturar todos los canales de MEG.
% filenum = 4;
switch dataset
    case 1
        signature = 'Analysis_Wave4'; % 'params'
    case 2
%         signature = 'Dataset2-Long_first';
        signature = 'Dataset2-Long_first-DETREND';
    case -1
        signature = 'null';
end

if filenum > 0
    cfg.dataset = fullfile(data_folder,username);  % Fichero a leer.
else
    load(fullfile('..','data',username));
end

% We create the results directory if it doesn't exist.
if isdir(fullfile('..','Results')) == 0
    mkdir(fullfile('..','Results'))    
end

if isdir(fullfile('..','Results', signature)) == 0
    mkdir(fullfile('..','Results', signature))    
end

if isdir(fullfile('..','Results', signature, username)) == 0
    mkdir(fullfile('..','Results', signature, username))    
end

% cfg.dataset = fullfile(dataset,username);
% data = ft_preprocessing(cfg);
% plot((1:length(data.trial{1}))/600,data.trial{1}(250:255,:)')



hdr   = ft_read_header(cfg.dataset);
event = ft_read_event(cfg.dataset);
tr_channel = 'UPPT001';

counter = 1;
for i = 1:length(event)
    if strcmp(event(i).type, tr_channel)% & event(i).value <= 40
        sample(counter) = event(i).sample;
        value(counter) = event(i).value;
        counter = counter+1;
    end
        
end
clear hdr event

%% Trigger signals
verbose = 0;
if verbose == 1
    subplot(211)
    sample_time = sample/sampling;
    stem(sample_time, value)
    xlabel('Time (s)');
    ylabel('Trigger val');
    mean((sample(1:end-1)-sample(2:end))/600)
    title(sprintf('Time diff.(s) Mean: %.2f, Std: %.2f', mean(diff(sample_time)), std(diff(sample_time))));
    diff(sample_time)
    subplot(212)
    channel_val =5;
    sample_time = sample(value == channel_val)/sampling;
    stem( sample_time, value(value == channel_val)) 
    xlabel('Time (s)');
    ylabel('Trigger val')
    title(sprintf('Time diff.(s) Mean: %.2f, Std: %.2f', mean(diff(sample_time)), std(diff(sample_time))));
    % fprintf('Mean delay: %.2f\n', mean(diff(sample_time)))
    diff(sample_time)
end
%% We read data according to each trigger value.
% We begin with number 5.
% data = ft_preprocessing(cfg);

% cfg.dataset = fullfile('am061187_MChait2_20120510_01.ds');
cfg.trialdef.eventtype  = 'UPPT001'; % Aquí indicamos el itpo de evento que nos ayudará a separar en epochs.
cfg.trialdef.eventvalue = trigger;

if (trigger == 5) | (trigger == 15)
    cfg.trialdef.prestim    = 0.5; % Espacio de preestímulo (antes del estímulo). Puede usarse para extraer la media.
    cfg.trialdef.poststim   = 3.5%.5; % Espacio de señal que nos interesa, después del evento.
else
    cfg.trialdef.prestim    = 0.5; % Espacio de preestímulo (antes del estímulo). Puede usarse para extraer la media.
    cfg.trialdef.poststim = 16;
end
try
    cfg = ft_definetrial(cfg);
    % cfg = [];
%     cfg.demean      = 'yes'; %%%% UNDO
    % cfg.detrend = 'yes';
    cfg.baselinewindow = [-cfg.trialdef.prestim 0];% in seconds
    data = ft_preprocessing(cfg);
    cfg=[];
    cfg.lpfilter ='yes';
    cfg.lpfreq = 30;
    data_block=ft_preprocessing(cfg, data);
catch
%     data_block = struct('hdr',{},'fsample',{},'sampleinfo',{},'trialinfo',{},'grad',{},'trial',-1,...
%                         'time',{},'label',{},'cfg',{});   
%                     data_block.trial = -1:;
    data_block.hdr = -1;
    data_block.fsample = -1;
    data_block.sampleinfo = -1;
    data_block.trialinfo = -1;
    data_block.grad = -1;
    data_block.trial = -1;    
    data_block.time = -1;
    data_block.label = -1;
    data_block.cfg = -1;    
    
    
    
    
    
    
    
end
end
% We append data to the final storage structure.
% ft_appenddata(cfg, data_pre, data)





