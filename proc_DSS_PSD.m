function psd_data = proc_DSS_PSD(trigger_list, subject_list, config)
    % Function in charge of computing and storing the PSD data for a
    % specific subject and trigger. It does not apply any additional
    % computation, and it can be useful to compute conditions whose files
    % are missing or wrong (for some reason).
    % 
    % trigger: value of the trigger representing the modality whose PSD we
    %       want to compute.
    %
    % subject: index representing the subject.
    %
    % config 
    %   .out_folder : path of a folder to store data.
    %   .channel_modality : in our case, it can be 'temporal' or
    %       'occipital', and is related with the field config.channels 
    %       below.
    %   .channels  : list of selected channels, whose structure is the one
    %       provided by FT_CHANNELSELECTION.
    %   .store_data : set to 1 to store data in a .mat file and figures in
    %       a .fig file.
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
    out_folder          = config.out_folder;  %'Trigger_analysis_Post16_Pre500';
    channel_modality    = config.channel_modality;
%     channels            = config.channels;
    store_data          = config.store_data;
    n_components        = config.n_components;
    fs                  = config.fs;
    %% Data loading and reduction
    % Trial data of the selected subject.
%     load(fullfile('..','Results',out_folder, 'Preprocessed_data',sprintf('Long_timelock-TRIG_%d-SUBJ_%d.mat',trigger, subject)))
%     load(fullfile('..','Results',out_folder,'DSS_components','Transformed',sprintf('Xdss-TRIG_%d-%d-COMP_%d.mat',trigger_list(1), trigger_list(2), n_components)),'dss_comp','x_comp');
    
    % We apply channel selection for each one of the subjects.
    for subject_ind = 1:length(subject_list)
        switch config.channel_modality
            case 'temporal' 
                load(fullfile(config.channels_path, sprintf('Channels-SUBJ_%d.mat', subject_list(subject_ind))))
            case 'occipital'
                load(fullfile(config.channels_path, 'Channels','Channels-Occipital.mat'))
        end
        
        for trigger_ind = 1:length(trigger_list)
            if isdir(fullfile('..','Results',out_folder,'DSS_components','PSD')) == 0 | exist(fullfile('..','Results',out_folder,'DSS_components','PSD',sprintf('PSD-TRIG_%d-SUBJ_%d-COMP_%d.mat',trigger_list(trigger_ind), subject_list(subject_ind), n_components))) == 0
                load(fullfile('..','Results',out_folder,'DSS_components','Transformed',sprintf('Xdss-TRIG_%d-SUBJ_%d-COMP_%d.mat',trigger_list(trigger_ind), subject_list(subject_ind), n_components)),'x_dss');

                % We load this file in order to have the structure valid for
                % the PSD computation. We use it as a reference for the DSS
                % data.
                load(fullfile('..','Results',out_folder, 'Preprocessed_data_AllChannels',sprintf('Long_timelock-TRIG_%d-SUBJ_%d.mat',trigger_list(trigger_ind), subject_list(subject_ind))))

                % Depending on the trigger value, we keep data from certain
                % temporal area.
                t0 = 0.2*fs;
                if (trigger_list(trigger_ind) == 5) | (trigger_list(trigger_ind) == 15)
                    T_init = 1;
                    T_end  = 3;
                else
                    T_init = 5;
                    T_end  = 15;
                end
                t_index = t0+(T_init*fs:T_end*fs-1);

                % We keep only the data that we are interested at.
                for ind_trials = 1:size(x_dss,3)
                    x_out{ind_trials} = x_dss(:, t_index,ind_trials);
                    t_out{ind_trials} = (1:size(x_out{ind_trials},2))/fs;
                end


                data_chunk = data_subject;
                x_out = {reshape(x_dss(channels_num, t_index, :), [length(channels_num), length(t_index)*size(x_dss,3)])};
                t_out = {(1:size(x_out{1},2))/fs};

                data_chunk.trial = x_out;
                data_chunk.time = t_out; 
                data_chunk.fsample = 600;
                data_chunk.label =  data_subject.label(channels_num);




                f_res = 3072;
                %% PSD using fieldtrip
                cfg             = [];            
                cfg.channel     = 'MEG';
                cfg.method      = 'mtmfft'; % mtmfft
                cfg.output      = 'pow';
                cfg.taper       = 'hanning';
                cfg.channel     = channels; % Preselected channels.
            %     cfg.foi = 0:.1:30;%linspace(0,30,1024);
                cfg.foi      = linspace(0,30,f_res);%[0,30]; % Bandwith of interest.
                cfg.keeptrials  = 'no';

                % PSD for all the individual trials, and then averaged.
                psd = ft_freqanalysis(cfg, data_chunk);
                mkdir(fullfile('..','Results',out_folder,'DSS_components','PSD'));
                save(fullfile('..','Results',out_folder,'DSS_components','PSD',sprintf('PSD-TRIG_%d-SUBJ_%d-COMP_%d.mat',trigger_list(trigger_ind), subject_list(subject_ind), n_components)),'psd');
            else
                load(fullfile('..','Results',out_folder,'DSS_components','PSD',sprintf('PSD-TRIG_%d-SUBJ_%d-COMP_%d.mat',trigger_list(trigger_ind), subject_list(subject_ind), n_components)),'psd');
            end
            
            psd_data(:,:,subject_ind, trigger_ind) = psd.powspctrm;
            psd_freq(:,subject_ind, trigger_ind) = psd.freq;
        end
        clear channels_num
        
    end
               
%     save(fullfile('..','Results',out_folder,'DSS_components','PSD',sprintf('PSD-TRIG_%d_%d-COMP_%d.mat',trigger_list(1), trigger_list(2), n_components)),'psd');

   %%
    subj_ind_short  = 1:length(subject_list);
    figure('units','normalized','outerposition',[0 0 1 1])        
    subplot(211)
    ind = 1;
    semilogy(mean(psd_freq(:,subj_ind_short,ind),2), squeeze(mean(rms(psd_data(:,:,subj_ind_short,ind),1),3)));
%     legend('TRIG:5. SHORT - RAND');
    xlabel('Frequency (Hz)');
    ylabel('Magnitude (fT)');
    title(sprintf('PSD. Subject: GLOBAL. MOD: %s. TRIG: %d.', channel_modality, trigger_list(ind)));
    xlim([0,30]);
    grid on
    
    subplot(212)
    ind = 2;
    semilogy(mean(psd_freq(:,subj_ind_short,ind),2), squeeze(mean(rms(psd_data(:,:,subj_ind_short,ind),1),3)));
%     legend('TRIG:5. SHORT - RAND');
    xlabel('Frequency (Hz)');
    ylabel('Magnitude (fT)');
    title(sprintf('PSD. Subject: GLOBAL. MOD: %s. TRIG: %d.', channel_modality, trigger_list(ind)));
    xlim([0,30]);
    grid on
    %% Results storage
    if store_data == 1
        mkdir(fullfile('..','Results',out_folder,'DSS_components','PSD_plots'));
        savefig(fullfile('..','Results',out_folder,'DSS_components','PSD_plots',sprintf('PSD-TRIG_%d_%d-COMP_%d.fig',trigger_list(1), trigger_list(2), n_components)));
    end
    
    
    %% Analysis of some individual bands.
    if trigger_list(1) == 10
        f_interest = [4, 8];     
    elseif trigger_list(1) == 5
        f_interest = [2, 20];        
    end
    
    for f_interest_ind = 1:length(f_interest)
        switch f_interest(f_interest_ind)
            case 2
                s_ind = 206
            case 4
                s_ind = 411;

            case 8
                s_ind = 820;  % 8.001 Hz.
            case 20
                s_ind = 2049;
        end
        n_ind = [(s_ind)-(50:-1:1), (s_ind)+(1:50)];
        for subject_ind = 1:length(subject_list)
            subplot(ceil(length(subject_list)/2), 2, subject_ind)
            semilogy(squeeze(mean(psd_freq(:,subject_ind,:),3)),squeeze(rms(psd_data(:,:,subject_ind,:),1)))
            title(sprintf('PSD. Subject: %d. MOD: %s.', subject_list(subject_ind), channel_modality));    


            s_magnitude = squeeze(rms(psd_data(:,s_ind, subject_ind, :),1));
            n_magnitude = mean(squeeze(rms(psd_data(:,n_ind, subject_ind, :),1)),1);

            snr(subject_ind, :, f_interest_ind) = s_magnitude'./n_magnitude;
            x_magnitude(subject_ind, :) = s_magnitude;
            signal_out(:, subject_ind,:) = (squeeze(rms(psd_data(:,[(s_ind)-(50:-1:1), s_ind, (s_ind)+(1:50)], subject_ind, :),1)))./repmat(n_magnitude,length(n_ind)+1,1);
        end
        subplot(ceil(length(subject_list)/2), 2, subject_ind+1)
           
        legend('RAND','REG')
        semilogy(mean(psd_freq(:,subj_ind_short,ind),2), squeeze(mean(rms(psd_data(:,:,subj_ind_short,:),1),3)));
%     legend('TRIG:5. SHORT - RAND');
        xlabel('Frequency (Hz)');
        ylabel('Magnitude (fT)');
        title(sprintf('PSD. Subject: GLOBAL. MOD: %s. TRIG: %d.', channel_modality, trigger_list(ind)));
        mkdir(fullfile('..','Results',out_folder,'DSS_components','PSD_plots'));
        savefig(fullfile('..','Results',out_folder,'DSS_components','PSD_plots',sprintf('PSD_SUBJECTS-F_%d-TRIG_%d_%d-COMP_%d.fig', f_interest(f_interest_ind), trigger_list(1), trigger_list(2), n_components)));
        close all
        
        % We plot the tone of interest mean PSD
        figure('units','normalized','outerposition',[0 0 1 .75])
        freq =  [(s_ind)-(50:-1:1), s_ind, (s_ind)+(1:50)];       
        plot(squeeze(mean(psd_freq(freq,subj_ind_short,:),2)),squeeze(mean(signal_out(:, subj_ind_short,:),2)), 'Linewidth',3)
        xlabel('Frequency (Hz)');
        ylabel('AVG SNR');
        switch trigger_list(1) 
            case 5
                condition = 'SHORT';
            otherwise
                condition = 'LONG';
        end
        title(sprintf('Seq: %s. Tone: %d Hz.', condition, f_interest(f_interest_ind)))        
        grid on
        hold on
        perc = 0.05;
        Diff = signal_out(:,:,1) - signal_out(:,:,2);
        dataB=bootstrap(Diff'); 
        s=findSigDiff(dataB, perc);
        plot(squeeze(mean(psd_freq(freq,subj_ind_short,1),2)),abs(s)*0.2, 'Linewidth',3)
        legend('RAND','REG',sprintf('p=%.2f', perc))
%         ylim([-1,35])
        
        mkdir(fullfile('..','Results',out_folder,'DSS_components','PSD_plots'));
        savefig(fullfile('..','Results',out_folder,'DSS_components','PSD_plots',sprintf('Tone-F_%d-TRIG_%d_%d-COMP_%d.fig', f_interest(f_interest_ind), trigger_list(1), trigger_list(2), n_components)));
        
        snr_out = snr(subj_ind_short,:,f_interest_ind);
        save(fullfile('..','Results',out_folder,'DSS_components','PSD_plots',sprintf('SNR-F_%d-TRIG_%d_%d-COMP_%d.mat', f_interest(f_interest_ind), trigger_list(1), trigger_list(2), n_components)), 'snr_out', 'x_magnitude', 'signal_out');
        
        close all
        
        
    end
end