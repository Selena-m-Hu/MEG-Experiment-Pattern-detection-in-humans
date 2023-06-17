% File to be used in order to compute frequency analysis from MEG data, 
% including the following tasks:
% * Data preprocessing (trial extraction, baselining, LP filtering).
% * Channel selection (from Temporal and Occipital sensors).
% * Computation of the Power Spectral Density (PSD).
% * Plotting of PSD results.
% * Computation of the SNR for frequencies of interest (4 Hz and 20 Hz).
% * Computation of the AVG magnitude of Alpha waves (from Occipital sensors
%   and considering the bandwidth from [8.5, 12] Hz.
% * Repeated measures anova for the results.
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
clc
trigger_list        =   [5, 10, 15, 20];               % Indicate the trigger list.
% subject_list        =   [2:15];                       % Indicate the subjects to analyze.
subject_list        =   20;  
config.out_folder   =   'Trigger_analysis_Module';  % Output data folder
% fieldtrip_path      =   fullfile('..','fieldtrip'); % Fieldtrip path.  
fieldtrip_path      =   'C:\Users\Chait Lab\Documents\MATLAB\Toolboxes\fieldtrip-20170831';
addpath(genpath(fieldtrip_path)); 

% We define which parts of the code should be (re)computed.
setup.trigger_plot          =   0; % Set to 1 to plot trigger information. Set to 0 otherwise.
setup.preprocessing         =   0; % Set to 1 to preprocess data and merge block information. Set to 0 otherwise.
setup.channel_selection     =   0; % Set to 1 to select channels from temporal or occipital regions (select modality below). Set to 0 otherwise.
setup.channel_occipital     =   0; % Set to 1 to select channels from temporal or occipital regions (select modality below). Set to 0 otherwise.
setup.compute_PSD           =   0; % Set to 1 if you have a PSD missing for a specific (subject, trigger) pair and you need to get it. Set to 0 otherwise.
setup.compute_PlotPSD       =   0; % Set to 1 to compute all the PSDs for a subject considering all the triggers. Set to 0 otherwise.
setup.compute_PlotGroupPSD  =   0; % Set to 1 to plot the AVERAGE PSD (through subjects) for the specified triggers. Set to 0 otherwise.
setup.compute_ComputeSNR    =   0; % Set to 1 to compute an SNR table for the combo SUBJxMODALITY, and depict their graphics. Set to 0 otherwise.
setup.compute_ComputeAlpha  =   0; % Same than the previous one, but for ALPHA WAVES.
setup.compute_AnovaSNR      =   0; % Set to 1 to compute Repeated Measures ANOVA for the SNR table previously computed. Set to 0 otherwise.
setup.compute_AnovaAlpha    =   0; % Same than the previous one, but computed instead for the ALPHA WAVES table.
%% Trigger plotting (LOC and regular experimental data).
if setup.trigger_plot == 1
    subject             = 3; % Indicate the index of the subject under analysis.
    config.localizer    = 1; % Set to 1 to plot LOC triggers. Use 0 if you want to check the triggers for one of the regular blocks of the experiment.
    pre_TriggerPlot(subject, config);
    
end

%% Preprocessing block by block, and merging the data in a single file.
if setup.preprocessing == 1
    % We have considered the standard configuration:
    % * Allows to reject trials visually.
    % * It stores the output structure that contains the merged info of all the
    %   blacks of the subject.
    config.reject_visual    = 1; % Set to 1 to reject trials visually. Use 0 otherwise.
    config.store_data       = 1; % Set to 1 to store the data in a .mat file. Use 0 otherwise (useful for testing purposes).
    
    for subject = subject_list
        for trigger = trigger_list
            [data_subject] = pre_BlockProcessing(trigger, subject, config);
        end
    end
    
end

%% Channel selection (temporal channels).
% Channel selection runs independently for each subject, but is obtained
% from the average signals of all the trigger conditions.
if setup.channel_selection == 1
    config.channel_modality = 'temporal';   % 'temporal', 'occipital'
    config.plot_channels = 1;               % In this case, it would plot topography using temporal information.
    config.store_data = 1;
    
    for subject = subject_list
        pre_ChannelSelection(trigger_list, subject, config);
    end
end

%% Occipital channel selection
% This case is similar to the previous one, although we only focus on
% getting the occipital channels. Consequently, we only need to perform
% this operation using a single subject.
if setup.channel_occipital == 1
    config.channel_modality = 'occipital';  % 'temporal', 'occipital'
    config.plot_channels = 0;               % In this case, it would plot topography using PSD information.
    config.store_data = 1;
    
    for subject = subject_list
        pre_ChannelSelection(trigger_list, subject, config);
    end
end

%% PSD estimation for and individual subject and condition (trigger).
% In this section we make use of the preprocessed (and stored) information
% and the precomputed channels for a set of subjects.
if setup.compute_PSD == 1
    config.store_data           = 1;            % Set to 1 to store the data in a .mat file. Use 0 otherwise (useful for testing purposes).
    config.channel_modality     = 'temporal';   % 'temporal', 'occipital'

    for subject = subject_list
        % Channels of the selected subject.
        switch config.channel_modality
            case 'temporal' 
                load(fullfile('..','Results',config.out_folder, 'Channels',sprintf('Channels-SUBJ_%d.mat', subject)))
            case 'occipital'
                load(fullfile('..','Results',config.out_folder, 'Channels','Channels-Occipital.mat'))
        end
        config.channels     = channels;
        
        for trigger = trigger_list
             [psd_data] = proc_PSD(trigger, subject, config);
        end
    end
end

%% Compute and plot PSD results for a specific subject considering all the triggers at the same time.
% Similar to the previous section, although it focuses on every individual
% subject and computes the PSD for the four triggers in a single run.
if setup.compute_PlotPSD == 1
    config.store_data           = 1;            % Set to 1 to store the data in a .mat file. Use 0 otherwise (useful for testing purposes).
    config.channel_modality     = 'temporal';   % 'temporal', 'occipital'
    for subject = subject_list
        % Channels of the selected subject.
        switch config.channel_modality
            case 'temporal' 
                load(fullfile('..','Results',config.out_folder, 'Channels',sprintf('Channels-SUBJ_%d.mat', subject)))
            case 'occipital'
                load(fullfile('..','Results',config.out_folder, 'Channels','Channels-Occipital.mat'))
        end
        config.channels         = channels;        

        [dummy] = proc_PlotPSD(trigger_list, subject, config);

    end
end

%% Plot PSD results for ALL the subjects and triggers. 
% This function is used to plot the global results for each modality
% (trigger). It uses the PSD files previously computed and average the data
% through all the subjects.
if setup.compute_PlotGroupPSD == 1
    config.store_data           = 1;            % Set to 1 to store the data in a .mat file. Use 0 otherwise (useful for testing purposes).
    config.channel_modality     = 'temporal';   % 'temporal', 'occipital'
  
    [dummy] = proc_PlotGroupPSD(trigger_list, subject_list, config);
end

%% Compute the SNR graphs and results.
% In this section we load all the files of the subjects and modalities, and
% compute an SNR table. We also plot their normalized magnitudes comparing
% our 4 modalities.
if setup.compute_ComputeSNR == 1
    config.store_data           = 1;            % Set to 1 to store the data in a .mat file. Use 0 otherwise (useful for testing purposes).
    config.channel_modality     = 'temporal';   % These are the channels where we can find the cleanest representation of the SNR.
    
    [snr_results] = proc_ComputeSNR(trigger_list, subject_list, config);
end


%% Compute the PSD results for the Alpha waves.
% Similar to the previous analysis, although in this case we compute the
% average magnitude of the PSD previously extracted for the Alpha waves.
% The considered bandwith is f = [8.5, 12] Hz.
if setup.compute_ComputeAlpha == 1
    config.store_data           = 1;            % Set to 1 to store the data in a .mat file. Use 0 otherwise (useful for testing purposes).
    config.channel_modality     = 'occipital';  
    
    [alpha_results] = proc_ComputeAlpha(trigger_list, subject_list, config);
end

%% Compute Repeated Measures ANOVA (1-way and 2-way) for the SNR table.
% Credit of the RM-ANOVA functions: 
% * anova_rm    : Arash Salarian, mailto://arash.salarian@ieee.org. https://uk.mathworks.com/matlabcentral/fileexchange/22088-repeated-measures-anova
% * rm_anova2.m : Aaron Schurger. https://uk.mathworks.com/matlabcentral/fileexchange/6874-two-way-repeated-measures-anova
if setup.compute_AnovaSNR == 1
    % We load the data that will be used to compute the statistics of the
    % SNR.
    selected_subjects = [1:12,14]; % Subjects 1 to 15 discarding number 14.
    load(fullfile('..','Results', config.out_folder,'Numerical_results', sprintf('SNR')), 'snr_results');
    
    % We compute the Repeated Measures 1-way ANOVA for the different 
    % combinatios of modalities (RAND-REG) vs (LONG-SHORT).
    clear data_stats_1way
    data_stats_1way.anova_snr_short_RANDvsREG = anova_rm( snr_results([1,3],selected_subjects)','off');    
    data_stats_1way.anova_snr_long_RANDvsREG = anova_rm(snr_results([2,4],selected_subjects)','off');      
    data_stats_1way.anova_snr_RAND_LONGvsSHORT = anova_rm( snr_results([1,2],selected_subjects)','off');
    data_stats_1way.anova_snr_REG_LONGvsSHORT = anova_rm(snr_results([3,4],selected_subjects)','off');
        
    config.store_data   = 1;
    config.data_name    = 'SNR';    
    if config.store_data == 1
        mkdir(fullfile('..','Results',config.out_folder,'Numerical_results'));
        save(fullfile('..','Results',config.out_folder,'Numerical_results', sprintf('RM1WayANOVA_%s', config.data_name)), 'data_stats_1way');
        close all
        
    end
    
    % We compute then the Repeated Measures 2-way anova for all the data.
    [SNR_stats] = proc_RMTwoWayAnova(snr_results(:,selected_subjects), config);

end

%% Compute Repeated Measures ANOVA (1-way and 2-way) for the ALPHA WAVES table.
% Credit of the RM-ANOVA functions: 
% * anova_rm    : Arash Salarian, mailto://arash.salarian@ieee.org. https://uk.mathworks.com/matlabcentral/fileexchange/22088-repeated-measures-anova
% * rm_anova2.m : Aaron Schurger. https://uk.mathworks.com/matlabcentral/fileexchange/6874-two-way-repeated-measures-anova
if setup.compute_AnovaAlpha == 1
    % We load the data that will be used to compute the statistics of the
    % alpha waves.
    selected_subjects = [1:12,14];
    load(fullfile('..','Results', config.out_folder,'Numerical_results', sprintf('Alpha_magnitude')), 'alpha_results');
    
    % We compute the Repeated Measures 1-way ANOVA for the different 
    % combinatios of modalities (RAND-REG) vs (LONG-SHORT).
    clear data_stats_1way
    data_stats_1way.anova_alpha_short_RANDvsREG = anova_rm(alpha_results([1,3],selected_subjects)','off');
    data_stats_1way.anova_alpha_long_RANDvsREG = anova_rm(alpha_results([2,4],selected_subjects)','off');
    data_stats_1way.anova_alpha_RAND_LONGvsSHORT = anova_rm(alpha_results([1,2],selected_subjects)','off');
    data_stats_1way.anova_alpha_REG_LONGvsSHORT = anova_rm(alpha_results([3,4],selected_subjects)','off');
    
    config.store_data   = 1;
    config.data_name    = 'Alpha';            
    if config.store_data == 1
        mkdir(fullfile('..','Results',config.out_folder,'Numerical_results'));
        save(fullfile('..','Results',config.out_folder,'Numerical_results', sprintf('RM1WayANOVA_%s', config.data_name)), 'data_stats_1way');
        close all
        
    end
    
    % We compute then the Repeated Measures 2-way anova for all the data.
    [Alpha_stats] = proc_RMTwoWayAnova(alpha_results(:,selected_subjects), config);
    
end