function psd_data = proc_PlotGroupPSD(trigger_list, subject_list, config)
    % This code does not compute any PSD. It simply reads it from the 
    % stored files that were obtained using proc_PSD.m or proc_PlotPSD.m,
    % and then plots four subplots (one per condition) using the average
    % subject PSD.
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
            frequency_out(trigger_ind, subject_ind,:) = psd_data.freq;
            psd_out(trigger_ind, subject_ind, :) = rms(psd_data.powspctrm);
        end
        
    end
    %%
    psd_global = permute(mean(sqrt(psd_out),2),[1,3,2])*1E15;
    freq_global = permute(mean(frequency_out,2),[1,3,2]);
    
    figure('units','normalized','outerposition',[0 0 1 1])        
    subplot(221)
    ind = 1;
    semilogy(freq_global(ind,:), psd_global(ind,:));
    legend('TRIG:5. SHORT - RAND');
    xlabel('Frequency (Hz)');
    ylabel('Magnitude (fT)');
    title(sprintf('PSD. Subject: GLOBAL. MOD: %s. TRIG: %d.', channel_modality, trigger_list(ind)));
    xlim([0,30]);
    
    
    subplot(223)
    ind = 3;
    semilogy(freq_global(ind,:), psd_global(ind,:));
    legend('TRIG:15. SHORT - REG');
    xlim([0,30]);
    title(sprintf('PSD. Subject: GLOBAL. MOD: %s. TRIG: %d.', channel_modality, trigger_list(ind)));
    xlabel('Frequency (Hz)');
    ylabel('Magnitude (fT)');

    subplot(222)
    ind = 2;
    semilogy(freq_global(ind,:), psd_global(ind,:));
    legend('TRIG:10. LONG - RAND');
    xlabel('Frequency (Hz)');
    ylabel('Magnitude (fT)');
    xlim([0,30]);
    title(sprintf('PSD. Subject: GLOBAL. MOD: %s. TRIG: %d.', channel_modality, trigger_list(ind)));
    hold on;
    
    subplot(224)
    ind = 4;
    semilogy(freq_global(ind,:), psd_global(ind,:));
    legend('TRIG:20. LONG - REG');
    xlim([0,30]);
    title(sprintf('PSD. Subject: GLOBAL. MOD: %s. TRIG: %d.', channel_modality, trigger_list(ind)));
    xlabel('Frequency (Hz)');
    ylabel('Magnitude (fT)');
    
    %% Results storage
    if store_data == 1
        mkdir(fullfile('..','Results',out_folder,'PSD_graphs'));
        % We store each psd in a separate file.
        savefig(fullfile('..','Results',out_folder,'PSD_graphs', sprintf('GRAPH-MOD_%s', channel_modality)));
%         close all
%         pause(1)
    end

end