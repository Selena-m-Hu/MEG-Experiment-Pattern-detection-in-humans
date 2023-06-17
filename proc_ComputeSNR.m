function snr_results = proc_ComputeSNR(trigger_list, subject_list, config)
    % This function is used to compute the SNR from the different subjects
    % and modalities. If the modality is SHORT the frequency of interest is
    % around 20 Hz. If it is LONG, it is placed around 4 Hz. 
    % The function not only computes an SNR table, but also plots a set of
    % graphs representing the magnitude of the frequency of interest 
    % normalized to the region noise for all our modalities.
    %
    % trigger_list: list of triggers whose average PSD we want to compute.
    %
    % subject_list: list of subjects whose PSDs we want to use to compute
    % an average PSD graph.
    %
    % config 
    %   .out_folder : path of a folder to store data.
    %   .channel_modality : in our case, it can be 'temporal' or
    %       'occipital', and is related with the field config.channels 
    %       below.
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

    out_folder          = config.out_folder; 
    channel_modality    = config.channel_modality;
    store_data          = config.store_data;

    %% Data loading and reduction
    % Trial data of the selected subject and all the triggers.
    clear psd_out frequency_out
    for trigger_ind = 1:length(trigger_list)
        for subject_ind = 1:length(subject_list)
            load(fullfile('..','Results',out_folder,'PSD', sprintf('PSD-MOD_%s-TRIG_%d-SUBJ_%d', channel_modality, trigger_list(trigger_ind), subject_list(subject_ind))),'psd_data','channels');
            
            if trigger_list(trigger_ind) == 5 | trigger_list(trigger_ind) == 15
                signal_ind = 1997:1999; % f = [19.97, 19.99] Hz %% [19.99, 20.01] Hz
                noise_ind = [1947:1996,2000:2047]; % f = [19.47, 19.96] U [20.00, 20.47] Hz
            else
                signal_ind = 399:401; % f = [3.99, 4.01] Hz
                noise_ind = [350:398,402:450]; % f = [3.5, 3.98] U [4.02, 4.5] Hz
            end
            
            
            signal_mean(trigger_ind, subject_ind) = mean(rms(psd_data.powspctrm(:, signal_ind)));
            noise_mean(trigger_ind, subject_ind) = mean(rms(psd_data.powspctrm(:, noise_ind)));
            snr_graph(trigger_ind, subject_ind, :) = rms(psd_data.powspctrm(:, min(noise_ind):max(noise_ind)))/noise_mean(trigger_ind, subject_ind);
            frequency_graph(trigger_ind, subject_ind, :) = psd_data.freq(:, min(noise_ind):max(noise_ind));
        end
        
    end
    %% SNR graphs for each subject, considering LONG and SHORT conditions.
    
    for subject_ind = 1:length(subject_list)
        
        figure('units','normalized','outerposition',[0 0 1 1])  
        % SHORT condition
        subplot(211)
        trigger_ind = 1;
        snr_data = reshape(snr_graph(trigger_ind, subject_ind, :), 1, 101);
        freq_data = reshape(frequency_graph(trigger_ind, subject_ind, :), 1, 101);
        plot(freq_data, snr_data, 'Linewidth',3)
        hold on
        trigger_ind = 3;
        snr_data = reshape(snr_graph(trigger_ind, subject_ind, :), 1, 101);
        freq_data = reshape(frequency_graph(trigger_ind, subject_ind, :), 1, 101);
        plot(freq_data, snr_data, 'Linewidth',3)
        legend('TRIG:5 - RAND', 'TRIG:15 - REG');       
        xlabel('Frequency (Hz)');
        ylabel('SNR')
        title(sprintf('SNR graph - Subject: %d',subject_list(subject_ind)))
        xlim([min(freq_data), max(freq_data)])
        
        % LONG condition
        subplot(212)
        trigger_ind = 2;
        snr_data = reshape(snr_graph(trigger_ind, subject_ind, :), 1, 101);
        freq_data = reshape(frequency_graph(trigger_ind, subject_ind, :), 1, 101);
        plot(freq_data, snr_data, 'Linewidth',3)
        hold on
        trigger_ind = 4;
        snr_data = reshape(snr_graph(trigger_ind, subject_ind, :), 1, 101);
        freq_data = reshape(frequency_graph(trigger_ind, subject_ind, :), 1, 101);
        plot(freq_data, snr_data, 'Linewidth',3)
        legend('TRIG:10 - RAND', 'TRIG:20 - REG');       
        xlabel('Frequency (Hz)');
        ylabel('SNR')
        xlim([min(freq_data), max(freq_data)])
        if store_data == 1
            mkdir(fullfile('..','Results',out_folder,'PSD_graphs'));
            % We store each subject graph in an individula file.
            savefig(fullfile('..','Results',out_folder,'PSD_graphs', sprintf('GRAPH-MOD_%s-SUBJ_%d', channel_modality, subject_list(subject_ind))));
            close all
    
        end
    end
    
    %%
    snr_results = (signal_mean./noise_mean);
    if store_data == 1
        mkdir(fullfile('..','Results',out_folder,'Numerical_results'));
        save(fullfile('..','Results',out_folder,'Numerical_results', sprintf('SNR')), 'snr_results');
        close all
        
    end

end