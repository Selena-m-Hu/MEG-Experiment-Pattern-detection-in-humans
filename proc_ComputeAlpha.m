function alpha_results = proc_ComputeAlpha(trigger_list, subject_list, config)
    % This function is used to compute the average magnitude of the ALPHA 
    % WAVES for the whole set of subjects and modalities. In this case, the
    % region of interest is the bandwidth f = [8.5, 12] Hz, excluding the
    % area of the 8 Hz since we expect the aparition of certain armonics
    % caused in the LONG condition.
    % For this function no plot is depicted.
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
                        
            signal_ind = 850:1200; % f = [8.50, 12.00] Hz 

%             signal(trigger_ind, subject_ind, :) = (rms(psd_data.powspctrm(:, signal_ind))*1e15);
            alpha_results(trigger_ind, subject_ind) = mean(sqrt(rms(psd_data.powspctrm(:, signal_ind)))*1e15);
%             frequency_graph(trigger_ind, subject_ind, :) = psd_data.freq(:, signal_ind);
        end
        
    end    
 
    %%
   
    if store_data == 1
        mkdir(fullfile('..','Results',out_folder,'Numerical_results'));
        save(fullfile('..','Results',out_folder,'Numerical_results', sprintf('Alpha_magnitude')), 'alpha_results');
        close all
        
    end

end