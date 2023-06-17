% File used to get the PSD from the DSS transformed files. If they are
% already computed, it plots the averaged results for the selected
% participants.
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
% Last update: 18/June/2018
clc
trigger_list            =   [10, 20];               % [5 10 15 20], 5=RAND SHORT, 10=RAND LONG, 15=REG SHORT, 20=REG LONG. Put them in couples (5,15) and (10,20). Indicate the trigger list.
subject_list            =   [2];                   % [2:13, 15]. Indicate the subjects to analyze.
config.hpfreq           =   2;                     
config.lpfreq           =   30; 
config.baseline         =   'silence';              % Options: 'tone', 'silence' (original prestimuli data), 'activity' (2 to 2.5 seconds)
config.n_components     =   3;                      % DSS components to keep.
config.out_folder       =   sprintf('Trigger_analysis_PRE_HP%d_LP%d', config.hpfreq, config.lpfreq);  % Output data folder. High passed data (2Hz). Detrended.
% config.channels_path    =   fullfile('..','Results','Trigger_analysis_Module','Channels');
config.channels_path    =   fullfile('..','Results','Channels_DSS');
config.fs               =   600;                      % Sampling frequency

% fieldtrip_path      =   fullfile('..','fieldtrip'); % Fieldtrip path.  
fieldtrip_path      =   'C:\Users\Chait Lab\Documents\MATLAB\Toolboxes\fieldtrip-20170831';
addpath(genpath(fieldtrip_path)); 

% We define which parts of the code should be (re)computed.
setup.compute_PSD               =   1;

%% PSD estimation for and individual subject and condition (trigger).
% In this section we make use of the preprocessed (and stored) information
% and the precomputed channels for a set of subjects.
if setup.compute_PSD == 1
    config.store_data           = 1;            % Set to 1 to store the data in a .mat file. Use 0 otherwise (useful for testing purposes).
    config.channel_modality     = 'temporal';   % 'temporal', 'occipital'

    [psd_data] = proc_DSS_PSD(trigger_list, subject_list, config);

end