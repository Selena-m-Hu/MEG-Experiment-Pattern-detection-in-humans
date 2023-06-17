%% Temporal plotting for single subject 
% For newly aquired subjects
% By Mingyue Hu, Aug, 2022

out_folder = config.out_folder;
n_components = config.n_components; 

load_DSS_chann = 1;  % load selected channels based on DSSed signal
plotDSS = 1; %1:plot DSSed data; 0: plot the raw data(preprocessed)
labelBL = ['BL_SilentInterval'];
for subject_ind = 1:length(subject_list)
    
%% load the dss data computed across conditions (single = 1)
  for trigger_ind = 1:length(trigger_list)

load(fullfile('..','Results',out_folder,'Preprocessed_data_AllChannels',...
sprintf('Long_timelock-TRIG_%d-SUBJ_%d', trigger_list(trigger_ind) ,subject_list(subject_ind))));

load(fullfile('..','Results',out_folder,'DSS_components','New_Transformed',...
sprintf('Xdss-TRIG_%d-SUBJ_%d-SINGLE-COMP_%d.mat',trigger_list(trigger_ind), subject_list(subject_ind), n_components)));

%% specify the time interval that we are interested about
cfg = [];
if trigger_list(trigger_ind) == 5 | trigger_list(trigger_ind) == 15
    cfg.toilim = [-0.2, 4];  % fast trial
else
    cfg.toilim = [-0.2, 16];  % slow trial
end
data_subject = ft_redefinetrial(cfg,data_subject);


%% replace the data of preprocessing data with the DSSed data
dssData = data_subject;
dssData.trial = x_dss2;
dssData.time = data_subject.time;

%% Averaging the data across trials 
cfg = [];
cfg.trials = 'all';
cfg.covariance         = 'no';
cfg.covariancewindow   = 'all';
cfg.keeptrials         = 'no';
cfg.removemean         = 'no';
cfg.vartrllength       = 0;
timelock_dss = ft_timelockanalysis(cfg,dssData);
timelock_raw = ft_timelockanalysis(cfg,data_subject);

all_RAW_data{trigger_ind} = timelock_raw;   %preprocessed data
all_Dss_data{trigger_ind} = timelock_dss;   %put the dssed timelock data into cell matrix

% if tirgger_ind == 2 
% mkdir(fullfile('..','Results',out_folder,'DSS_timelock'));
% save(fullfile('..','Results',out_folder,'DSS_timelock',...
% sprintf('Xdss--SUBJ_%d-SINGLE-timelock.mat',trigger_list(trigger_ind), subject_list(subject_ind))), 'all_Dss_data');
% end 
  end 

% Load selected channels
if load_DSS_chann == 1  % load dssed data based channels
 load(fullfile('..','Results',out_folder,'Channels_DSS',...
     sprintf('Channels-TRIG_%d_%d-SUBJ_%d',trigger_list(1), trigger_list(2),subject_list(subject_ind))),'channels', 'channels_num');
  if plotDSS  
    mkdir(fullfile('..','Results',out_folder,'DSS_components','RMSplot','DSS_seleChann'))
    f1 = figure(subject_ind);
    plot(all_Dss_data{1}.time,rms(all_Dss_data{1}.avg(channels_num,:)),'Color',[0.8500, 0.3250, 0.0980],'Linewidth', 1.5);
    legend('RAND')
    hold on
    plot(all_Dss_data{2}.time,rms(all_Dss_data{2}.avg(channels_num,:)),'Color', [0, 0.4470, 0.7410], 'Linewidth', 1.5);
    legend('REG')
    title(['Subject' num2str(subject_list(subject_ind))])
    mkdir(fullfile('..','Results',out_folder,'DSS_components','RMSplot','DSS_seleChann',labelBL))
    plotFolder = fullfile('..','Results',out_folder,'DSS_components','RMSplot','DSS_seleChann',labelBL...
    ,sprintf('Power-TRIG_%d_%d-SUBJ_%d',trigger_list(1),trigger_list(2),subject_list(subject_ind)));
    savefig(f1,plotFolder)
  else % we do not apply the dssed data based channels on raw data, so end here
  end 
else % load raw data based channels
 load(fullfile('..','Results',out_folder,'Channels',...
     sprintf('Channels-TRIG_%d_%d-SUBJ_%d',trigger_list(1), trigger_list(2),subject_list(subject_ind))),'channels', 'channels_num');
   if plotDSS
    mkdir(fullfile('..','Results',out_folder,'DSS_components','RMSplot','RAWdata_seleChann',labelBL))   
    f1 = figure(subject_ind);
    plot(all_Dss_data{1}.time,rms(all_Dss_data{1}.avg(channels_num,:)),'Color',[0.8500, 0.3250, 0.0980],'Linewidth', 1.5);
    legend('RAND')
    hold on
    plot(all_Dss_data{2}.time,rms(all_Dss_data{2}.avg(channels_num,:)),'Color', [0, 0.4470, 0.7410], 'Linewidth', 1.5);
    legend('REG')
    title(['Subject' num2str(subject_list(subject_ind))])
    plotFolder = fullfile('..','Results',out_folder,'DSS_components','RMSplot','RAWdata_seleChann',labelBL...
    ,sprintf('Power-TRIG_%d_%d-SUBJ_%d',trigger_list(1),trigger_list(2),subject_list(subject_ind)));
    savefig(f1,plotFolder)
   else     
    mkdir(fullfile('..','Results',out_folder,'Preprocessed_data_AllChannels','RMSplot',labelBL))
    f1 = figure(subject_ind);
    plot(all_RAW_data{1}.time,rms(all_RAW_data{1}.avg(channels_num,:)),'Color',[0.8500, 0.3250, 0.0980],'Linewidth', 1.5);
    legend('RAND')
    hold on
    plot(all_RAW_data{2}.time,rms(all_RAW_data{2}.avg(channels_num,:)),'Color', [0, 0.4470, 0.7410], 'Linewidth', 1.5);
    legend('REG')
    title(['Subject' num2str(subject_list(subject_ind))])
    plotFolder = fullfile('..','Results',out_folder,'Preprocessed_data_AllChannels','RMSplot',labelBL...
    ,sprintf('Power-TRIG_%d_%d-SUBJ_%d',trigger_list(1),trigger_list(2),subject_list(subject_ind)));
    savefig(f1,plotFolder)   
   end 
     
end    
end