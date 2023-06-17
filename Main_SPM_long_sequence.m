% Script used to perform the SPM conversion of the LONG sequence.
% Consequently, the triggers used should be the couple (10,20). 
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
trigger_list            =   [10, 20];               % [10, 20], 5=RAND SHORT, 10=RAND LONG, 15=REG SHORT, 20=REG LONG. Put them in couples (5,15) and (10,20). Indicate the trigger list.
% subject_list            =   [3:13,15];                   % [2:13, 15]. Indicate the subjects to analyze.
subject_list            =   [18:24]; 
config.hpfreq           =   0;   % Standard: 2                   
config.lpfreq           =   30;  % Standard: 30
% config.baseline         =   'silence';              % Options: 'tone', 'silence' (original prestimuli data), 'activity' (2 to 2.5 seconds)
config.temporal_frame   =   [8, 14]; % We determine the temporal frame in seconds that we want to process. Standard: [-0,.5]
config.n_components     =   3;                      % DSS components to keep.
config.out_folder       =   sprintf('Trigger_analysis_PRE_HP%d_LP%d', config.hpfreq, config.lpfreq);  % Output data folder. High passed data (2Hz). Detrended.
% config.channels_path    =   fullfile('..','Results','Trigger_analysis_Module','Channels');
config.channels_path    =   fullfile('..','Results','Channels_DSS');
config.fs               =   600;                      % Sampling frequency

% fieldtrip_path      =   fullfile('..','fieldtrip'); % Fieldtrip path.  
fieldtrip_path      =   'C:\Users\Chait Lab\Documents\MATLAB\Toolboxes\fieldtrip-20170831';
addpath(genpath(fieldtrip_path)); 

% setup.compute_ToneButterfly     =   1; % Must be performed for HP=2.

% %%
if trigger_list(1) == 10

%   config.store_data       = 0; % Set to 1 to store the data in a .mat file. Use 0 otherwise (useful for testing purposes).
    config.DSS              = 1; % Set to 1 to use DSS data. Use 0 for data that is not transformed.

    source_SequenceSTORAGE(trigger_list, subject_list, config);
    
end


%Load raw data
% out_folder = config.out_folder;
% T = 10;
% S=3;
% 
% load(fullfile('..','Results_Antonio_S1_S15',out_folder, 'Preprocessed_data_AllChannels',...
% sprintf('Long_timelock-TRIG_%d-SUBJ_%d.mat',T, S)));
% if T == 5 | T == 15
%     cfg.toilim = [-0.2, 4];   % fast trial
% else
%     cfg.toilim = [-0.2, 16];  % slow trial
% end
% data_subject = ft_redefinetrial(cfg,data_subject);
% 
% short_data_subject = data_subject; 
% 
% save(fullfile('..','Results',out_folder, 'Preprocessed_data_AllChannels','short_timelock',...
% sprintf('Short_timelock-TRIG_%d-SUBJ_%d.mat',T, S)),'short_data_subject'); 
%                




