 
function visualRejection(trigger, subject, config)
%% Visual rejection
 
 % This script is for visually inspecting and rejecting bad trials of DSSed
 % data or raw data
 % for subject 16 to 24
 % Written by Mingyue Hu, sep, 2022

% Inputs required: 
%   path
%   Trigger
%   Subject number
%   raw data 
%   DSS projected data
%---------------------------------------------------------------------------
    dssLoad = 1;       %load DSSed data
    out_folder = config.out_folder;
    n_components = config.n_components; 
    
    %% Trigger data averaging
    % In order to select the proper channels, we average the information from all
    % the triggers.     
       
       %Load raw data
       load(fullfile('..','Results',out_folder, 'Preprocessed_data_AllChannels',...
       sprintf('Long_timelock-TRIG_%d-SUBJ_%d.mat',trigger, subject)));
        if trigger == 5 | trigger == 15
            cfg.toilim = [-0.2, 4];   % fast trial
        else
            cfg.toilim = [-0.2, 16];  % slow trial
        end
        data_subject = ft_redefinetrial(cfg,data_subject);
     
      
       %Load DSS data
       load(fullfile('..','Results',out_folder,'DSS_components','New_Transformed',...
       sprintf('Xdss-TRIG_%d-SUBJ_%d-SINGLE-COMP_%d.mat',trigger, subject, n_components)));    
       cfg = [];
       
        
       short_data_subject = data_subject; 
        
        % We save the shortened raw data structure, so we can directly load
        % it for further analysis
        
        if exist(fullfile('..','Results',out_folder,'Preprocessed_data_AllChannels','short_timelock')) == 0
        mkdir(fullfile('..','Results',out_folder,'Preprocessed_data_AllChannels','short_timelock'));
        else 
        end 
        save(fullfile('..','Results',out_folder, 'Preprocessed_data_AllChannels','short_timelock',...
        sprintf('Short_timelock-TRIG_%d-SUBJ_%d.mat',trigger, subject)),'short_data_subject'); 
               
        
        if dssLoad
        %% replace the data of preprocessing data with the DSSed data
        data_subject.trial = x_dss2;
        else 
        end
    
        %% Visually reject bad trials
        cfg.channel = 'all';
        dss_data_subject = ft_rejectvisual(cfg,data_subject);
        
        %% Save the data structure adapted for FT requirement
        if exist(fullfile('..','Results',out_folder,'DSS_components','normalizedDSS_cfg')) == 0
        mkdir(fullfile('..','Results',out_folder,'DSS_components','normalizedDSS_cfg'));
        else 
        end 
        save(fullfile('..','Results',out_folder,'DSS_components','normalizedDSS_cfg',...
        sprintf('Xdss-TRIG_%d-SUBJ_%d-SINGLE-Clean-COMP_%d.mat',trigger, subject, n_components)),'dss_data_subject');  
        
        
 
  