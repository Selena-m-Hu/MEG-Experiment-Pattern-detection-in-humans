    function [] = Temporal_analysis(trigger, store, filenum)
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
        
%         addpath(genpath('C:\Users\Chait Lab\Documents\MATLAB\Toolboxes\fieldtrip-20170831'))
        addpath(genpath(fullfile('..','fieldtrip')))
        dataset = 2;
        generate_cell = 0;
        if generate_cell == 1
            for subject_ind = 1:length(filenum)
                for trigger_ind = 1:length(trigger)
                    load(fullfile('..','Results',out_folder,sprintf('Long_timelock-DATA_%d-TRIG_%d-SUBJ_%d',dataset,trigger(trigger_ind), filenum(subject_ind))))


                    cfg = [];
                    cfg.resamplefs = 120;
                    aux = ft_resampledata(cfg, aux);

                    load(fullfile('..','Results',out_folder,sprintf('Channels_SUBJ-%d.mat', filenum(subject_ind))));


    
                    cfg=[];
                    cfg.method = 'summary';
    %                 cfg.channel = channel{1};
                    timelock{trigger_ind,filenum(subject_ind)} = ft_timelockanalysis(cfg, aux);

                    % CHANNEL SELECTION (AROUND 100MS)
                    timelock{trigger_ind,filenum(subject_ind)}.fsample = aux.fsample;
                    M100dat=timelock{trigger_ind,filenum(subject_ind)}.avg'; % Obsérvese que aquí se transpone la media, por lo que tenemos (tiempo x canal).
                    t0 = 0.5*timelock{trigger_ind,filenum(subject_ind)}.fsample;
                    % tiempo, y todos los canales.

                    amps=mean(M100dat((t0+0.09*timelock{trigger_ind,filenum(subject_ind)}.fsample):(t0+0.11*timelock{trigger_ind,filenum(subject_ind)}.fsample), :),1);
                    [ampsSorted,idx]= sort(amps,2,'descend');

                    chnsSorted = timelock{trigger_ind,filenum(subject_ind)}.label(idx);

                    %selecting channels:

                    chns_selectedLpos=[];
                    chns_selectedRpos=[];
                    chns_selectedLneg=[];
                    chns_selectedRneg=[];
                    leftChansCountPos=0;
                    rightChansCountPos=0;
                    leftChansCountNeg=0;
                    rightChansCountNeg=0;

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
                    for count1=1:length(timelock{trigger_ind,filenum(subject_ind)}.label)
                        for count2=1:length(chns_selectedL)
                            if (strcmp(timelock{trigger_ind,filenum(subject_ind)}.label{count1},chns_selectedL{count2}) ~= 0)
                                chnsL_num=[chnsL_num count1];
                            end
                        end
                    end

                    chnsR_num=[];
                    for count1=1:length(timelock{trigger_ind,filenum(subject_ind)}.label)
                        for count2=1:length(chns_selectedR)
                            if (strcmp(timelock{trigger_ind,filenum(subject_ind)}.label{count1},chns_selectedR{count2}) ~= 0)
                                chnsR_num=[chnsR_num count1];
                            end
                        end
                    end

                    chns_selectedL = timelock{trigger_ind,filenum(subject_ind)}.label(chnsL_num); 
                    chns_selectedR = timelock{trigger_ind,filenum(subject_ind)}.label(chnsR_num);

                    channel{trigger_ind} = unique([chns_selectedL, chns_selectedR]);

                    clear area trial_ind
                    int_area{trigger_ind,filenum(subject_ind)} = [];                    
                    int_mean{trigger_ind,filenum(subject_ind)} = [];
                    for trial_ind = 1:length(aux.trial)
%                         keyboard
                        int_area{trigger_ind,filenum(subject_ind)}= [int_area{trigger_ind,filenum(subject_ind)} trapz(rms(aux.trial{trial_ind}([chnsL_num, chnsR_num],:)))];
                        int_mean{trigger_ind,filenum(subject_ind)}= [int_mean{trigger_ind,filenum(subject_ind)} mean(rms(aux.trial{trial_ind}([chnsL_num, chnsR_num],:)))];
                                       
                      
                      
                      
                      
                    end

                end

            end
            % We store the output results. 
            keyboard
            save timelocks timelock  chnsL_num chnsR_num int_area int_mean
        else         
            load timelocks
        end
        %% RM ANOVA for the integrals
        area_vector = [];
        REG_vector = [];
        LONG_vector = [];
        subj_vector = [];
        mean_vector = [];
        % TRIGGER 5
        for subject_ind = [2:3,5:13]
            for trigger_ind = 1:4
                area_len = length(int_area{trigger_ind, subject_ind});
                area_vector = vertcat(area_vector, int_area{trigger_ind, subject_ind}');
                subj_vector = vertcat(subj_vector, subject_ind*ones(area_len,1));
                mean_vector = vertcat(mean_vector, int_mean{trigger_ind, subject_ind}');
                
                switch trigger_ind
                    case 1
                        REG_vector = vertcat(REG_vector, zeros(area_len,1));
                        LONG_vector = vertcat(LONG_vector, zeros(area_len,1));
                    case 2
                        REG_vector = vertcat(REG_vector, zeros(area_len,1));
                        LONG_vector = vertcat(LONG_vector, ones(area_len,1));
                    case 3
                        REG_vector = vertcat(REG_vector, ones(area_len,1));
                        LONG_vector = vertcat(LONG_vector, zeros(area_len,1));
                    case 4
                        REG_vector = vertcat(REG_vector, ones(area_len,1));
                        LONG_vector = vertcat(LONG_vector, ones(area_len,1));
                end
                    
            end
        end
      
        
        area_stats_trial = rm_anova2(area_vector,subj_vector,REG_vector,LONG_vector,{'REG','LONG'});
        mean_stats_trial = rm_anova2(mean_vector,subj_vector,REG_vector,LONG_vector,{'REG','LONG'});
        
        
        
        
        
        
        
        
        
        %%
        for trigger_ind = [2,4] 
            for subject_ind = 1:length(filenum)
                if subject_ind == 1
                    aux_avg{trigger_ind} = timelock{trigger_ind,filenum(subject_ind)}.avg([chnsL_num,chnsR_num],:)/length(filenum);
                else
                    aux_avg{trigger_ind} = aux_avg{trigger_ind} + timelock{trigger_ind,filenum(subject_ind)}.avg([chnsL_num,chnsR_num],:)/length(filenum);
                end
            end
%             plot(timelock{trigger_ind}.time,rms(timelock{trigger_ind}.avg([chnsL_num,chnsR_num],:)));
            figure(1)
            
            subplot(211)
            plot(timelock{trigger_ind,filenum(subject_ind)}.time,rms(aux_avg{trigger_ind})*1E15);
            xlabel('Time (s)');
            ylabel('Magnitude (fT)');
            title('RMS signal from selected channels');
            legend('RAND','REG');
            
            hold on;
            subplot(212)
%             poly = polyfit(timelock{trigger_ind,subject_ind}.time,abs(hilbert(rms(aux_avg{trigger_ind}))),15);
            [envelope, residual] = emd(abs(hilbert(rms(aux_avg{trigger_ind}))));
            envelope_sum = sum(envelope(:,3:end),2)+residual
            plot(timelock{trigger_ind,filenum(subject_ind)}.time,envelope_sum*1E15, 'Linewidth',3)
            xlabel('Time (s)');
            ylabel('Magnitude (fT)');
            title('EMD envelope');
            legend('RAND','REG');
            hold on
            
            figure
           
            plot(timelock{trigger_ind,filenum(subject_ind)}.time,rms(aux_avg{trigger_ind})*1E15);
            hold on;
            plot(timelock{trigger_ind,filenum(subject_ind)}.time,envelope_sum*1E15, 'Linewidth',3)
            xlabel('Time (s)');
            ylabel('Magnitude (fT)');
            title('Original vs Envelope');
            legend('RMS signal','EMD envelope');
            
            
        end
%        
% %           
% %         close all
%         if (trigger(trigger_ind) == 5) | (trigger(trigger_ind) == 15)
%             if length(trigger) > 1
%                 savefig(fullfile('..','Results',out_folder,sprintf('Short_PSD_TOP-channels_SUBJ-%d.fig', filenum)))
%             else
%                 close all
%                 save(fullfile('..','Results',out_folder,sprintf('PSD_TOP-channels_TRIG-%d_SUBJ-%d.mat',trigger(trigger_ind), filenum)),'output')
%             end
%         else
%             if length(trigger) > 1
%                 savefig(fullfile('..','Results',out_folder,sprintf('Long_PSD_TOP-channels_SUBJ-%d.fig', filenum)))
%             else
%                 close all
%                 save(fullfile('..','Results',out_folder,sprintf('PSD_TOP-channels_TRIG-%d_SUBJ-%d.mat',trigger(trigger_ind), filenum)),'output')
%             end
%         end
%     close all
%        
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
    cfg.demean      = 'yes';
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





