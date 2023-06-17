% Script used to perform the tone analysis of the LONG sequence.
% Consequently, the triggers used should be the couple (10,20).
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
% Last update: 08/August/2018
clc
trigger_list            =   [10, 20];               % [5 10 15 20], 5=RAND SHORT, 10=RAND LONG, 15=REG SHORT, 20=REG LONG. Put them in couples (5,15) and (10,20). Indicate the trigger list.
subject_list            =   [2:13,15];                   % [2:13, 15]. Indicate the subjects to analyze.
config.hpfreq           =   0;         % Standard: 2            
config.lpfreq           =   30;        % Standard: 30
config.baseline         =   'silence';              % Options: 'tone', 'silence' (original prestimuli data), 'activity' (2 to 2.5 seconds)
config.temporal_frame   =   [8,14]; % We determine the temporal frame in seconds that we want to process.
config.n_components     =   3;                      % DSS components to keep.


config.out_folder       =   sprintf('Trigger_analysis_PRE_HP%d_LP%d', config.hpfreq, config.lpfreq);  % Output data folder. High passed data (2Hz). Detrended.
% config.channels_path    =   fullfile('..','Results','Trigger_analysis_Module','Channels');
config.channels_path    =   fullfile('..','Results','Channels_DSS');
config.fs               =   600;                      % Sampling frequency

% fieldtrip_path      =   fullfile('..','fieldtrip'); % Fieldtrip path.  
fieldtrip_path      =   'C:\Users\Chait Lab\Documents\MATLAB\Toolboxes\fieldtrip-20170831';
addpath(genpath(fieldtrip_path)); 


setup.compute_ToneButterfly     =   1; % Must be performed for HP=2.



%%
if setup.compute_ToneButterfly == 1

    config.reject_visual    = 1; % Set to 1 to reject trials visually. Use 0 otherwise.
    config.store_data       = 0; % Set to 1 to store the data in a .mat file. Use 0 otherwise (useful for testing purposes).
    config.load_channels    = 1;
    config.DSS              = 1; % Set to 1 to use DSS data. Use 0 for data that is not transformed.

    
%     config.temporal_frame   = [0, 2.5    ]; % We determine the temporal frame in seconds that we want to process.

    source_ToneSTORAGE(trigger_list, subject_list, config);
%     source_SequenceSTORAGE(trigger_list, subject_list, config);
    
end







