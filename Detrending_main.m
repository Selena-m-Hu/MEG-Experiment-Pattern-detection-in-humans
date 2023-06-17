%% Detrending scripts

clc
clear all

subject_list            = [2 3 4 5 6 7 8 9 10 11 12 13 15 16 17 18 19 20 21 22 23 24];                   % [2:13, 15]. Indicate the subjects to analyze.
oldsubjects = 1;
% subject_list            =   [7]; %trying just one subject for now
% subject_list = [16:24];
config.hpfreq           =   2;                     
config.lpfreq           =   30; 
% config.baseline         =   'activity';           % Options: 'tone', 'silence' (original prestimuli data), 'activity' (2 to 2.5 seconds)
config.n_components     =   3;                      % DSS components to keep. We usually keep 3 components
config.single           =   0;                      % Set to 1 to compute a DSS matrix for each one of the conditions.
config.out_folder       =   sprintf('Trigger_analysis_PRE_HP%d_LP%d', config.hpfreq, config.lpfreq);  % Output data folder. High passed data (2Hz). Detrended.
% config.channels_path    =   fullfile('D:\Results','Trigger_analysis_PRE_HP0_LP30','Channels_DSS');  % channels selected based on the data after proprocessing 
% config.channels_path    =   fullfile('..','Results','Channels_DSS');
% %channels selected after DSS
config.fs               =   600;                      % Sampling frequency for MEG is 600, for EEG is 2048
% config.load_channels    =   0;                      % Sampling frequency
%   trigger_list            =   [10, 20];                 % old triggers[5 10 15 20], 5=RAND SHORT, 10=RAND LONG, 15=REG SHORT, 20=REG LONG. Put them in couples (5,15) and (10,20). Indicate the trigger list.
%  trigger_list            =   [5 15]; 
 trigger_list            =   [10, 20];                  % we use the old triggers in our following analysis for simplicity
%  trigger_list2           =   [40, 80];                % new triggers[40 60 80 100] 40 = RAND SHORT, 60 RAND LONG, 80 = REG SHORT, 100 = REG LONG
%  trigger_list2           =   [60, 100];  
addpath('D:\fieldtrip-20220707'); 
addpath('D:\NoiseTools');

% We define which parts of the code should be (re)computed.
setup.preprocessingDetrend      =   0;   % Detrending and smoothing
setup.preprocessingMEM          =   1;   % Filtering and more (epoching etc)
% setup.channel_selection         =   0; % Set to 1 to select automatically the optimal channels.
setup.compute_DSScomputation    =   0;
% setup.compute_ToneButterfly     =   0; % Must be performed for HP=2.
% setup.compute_PSD               =   0;

% setup.to_SPM                    =   0;
% nicolas_procedure               =   0; % 1: Used to compute DSS maximizing differences between conditions. 0: Computes DSS after grouping the trials of both conditions.

%% Detrending and smoothing
if setup.preprocessingDetrend== 1
    % We have considered the standard configuration:
    % * Allows to reject trials visually.
    % * It stores the output structure that contains the merged info of all the
    %   blacks of the subject.
     config.store_data       = 1; % Set to 1 to store the data in a .mat file. Use 0 otherwise (useful for testing purposes).

    
    for subject = subject_list
        [data_subject] = pre_Detrend(subject, config);
    end
    
end
