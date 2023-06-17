
%% Tone analysis results plot
% for condition of slow sequence, (trigger 10 and trigger 20)
% This script is adapted by Mingyue Hu, Sep, 2022

 clear all;
 clc;
 
 addpath('D:\Antonio\fieldtrip'); 
 %parameter
 subject_list = [2 3 4 5 6 7 8 9 10 11 12 13 15 16 17 18 19 20 21 22 23 24]; 
 out_folder = 'Trigger_analysis_PRE_HP2_LP30'; 
 n_components = 3; 
%output_appendix = '_silence'; %baseline methods: 'tone', 'silence', 'activity' 
 T_init = 8;  
 T_end = 14;
 hpfreq= 2;
 lpfreq = 30;
 dss = 0; % 1:load dss data; 0:load raw data
 firstCycle = 0; 

%% if we need to concatenate data matrix from different subjects
%  out_folder = 'Trigger_analysis_PRE_HP2_LP30'; 
%  n_components = 3; 
%  output_appendix = ''; %baseline methods: '_tone', '_silence', '_activity'
%  cat_Tone_Data_Matrix(out_folder,n_components,output_appendix); 


%% We load the data from all subjects 
% those data matrix are the output of cat_Tone_Data_Matrix function

output_appendix = '_tone'; %baseline methods: 'tone', 'silence', 'activity'
baseline = output_appendix;

if strcmp(output_appendix,'_tone')
  bl_tone = '-30'; 
else
  bl_tone = '-';
end 

% load 40 channels data matrix
if dss
    if firstCycle
    load(fullfile('..','Results',out_folder,'tone_analysis',...
    sprintf('Tone_40Channels_allSub-FirstCycle-COMP_%d_BL_%s_BLwindow_%s', n_components, output_appendix, bl_tone)), 'stable_average_shape');

    % load all channels data matrix
    load(fullfile('..','Results',out_folder,'tone_analysis',...
    sprintf('Tone_275Channels_allSub-FirstCycle-COMP_%d_BL_%s_BLwindow_%s', n_components, output_appendix, bl_tone)), 'stable_average_shape275');

    else    
    load(fullfile('..','Results',out_folder,'tone_analysis',...
    sprintf('Tone_40Channels_allSub-COMP_%d_BL_%s_BLwindow_%s', n_components, output_appendix, bl_tone)), 'stable_average_shape');

    % load all channels data matrix
%     load(fullfile('..','Results',out_folder,'tone_analysis',...
%     sprintf('Tone_275Channels_allSub-COMP_%d_BL_%s_BLwindow_%s', n_components, output_appendix, bl_tone)), 'stable_average_shape275');
    end
    
else
    
   if firstCycle
    load(fullfile('..','Results',out_folder,'tone_analysis',...
    sprintf('Tone_40Channels_allSub-FirstCycle-RAWdata_BL_%s_BLwindow_%s', output_appendix, bl_tone)), 'stable_average_shape');

   else 
    load(fullfile('..','Results',out_folder,'tone_analysis',...
    sprintf('Tone_40Channels_allSub-RAWdata_BL_%s_BLwindow_%s', output_appendix, bl_tone)), 'stable_average_shape');

   end
   
end
    
%we need channel information from timelock for plotting topography, here we only load one file from the last subject
load(fullfile('..','Results',out_folder,'Preprocessed_data_AllChannels','short_timelock',...
sprintf('Short_timelock-TRIG_%d-SUBJ_%d.mat',20,24)),'short_data_subject')                                     
timelock = ft_timelockanalysis([],short_data_subject);  

all_subject_channel_40 = stable_average_shape;  
% allChann_allSub_Tone_Data = stable_average_shape275;

%% Subject curve-topography graphs

    for subject_ind = 1:length(subject_list)
        
        subject_data.time = (1:150)/150*250-50;
%       subject_data.label = timelock.label;
        subject_data.dimord = 'chan_time';       
        %We load and plot the temporal data.
        channels_num = 40;
        mean_subjects_data(:, subject_ind,:) = squeeze(rms(all_subject_channel_40(:,:, subject_ind, :),1));
    end
    
    
       
    %% Plot parameters
    timeseries = 0;
    topography = 1; 
    subject_list_short = 1:length(subject_list);

    %% Global results (average of all the subjects)
    
    if timeseries
    Diff=mean_subjects_data(:,subject_list_short,1) - mean_subjects_data(:,subject_list_short,2); 

    % We use bootstrap and compute if there is a significant difference
    % between conditions.
    
    perc = 0.05;
    perc2 = 0.01;
    dataB=bootstrap(Diff'); 
    s=findSigDiff(dataB, perc);
    s2=findSigDiff(dataB, perc2);
    k1 = 1;
    k2 = 3;
    
    condi = squeeze(mean(mean_subjects_data(:,subject_list_short,:),2));
    REG = condi(:,2);
    RAND = condi(:,1);
    REGall = mean_subjects_data(:,:,2);
    RANall = mean_subjects_data(:,:,1);
    REGstd = std(REGall')/sqrt(size(REGall,2));
    RANstd = std(RANall')/sqrt(size(RANall,2));
    time = 30:150;
    figure;
    shadedErrorBar(subject_data.time(time),RAND(time),RANstd(time),'lineProps','k');
    hold on 
    shadedErrorBar(subject_data.time(time),REG(time),REGstd(time),'lineProps','r');
    xlim([-0.5,200])
    
    hold on
    plot(subject_data.time(time), k1*abs(s(time)),'Linewidth', 12);
    hold on
    plot(subject_data.time(time), k2*abs(s2(time)),'Linewidth', 12);
    xlabel('Time (ms)')
    ylabel('Magnitude (fT)')
    title(sprintf('All subjects. Baseline: %s', baseline));
    end 
    
   %% Topography 
    if topography
    subject_data.label = timelock.label;
    subject_data.avg = mean(mean(all_subject_channel_40(:,:,subject_list_short,:),4),3);
    
    for ind_T = 1:2
        switch ind_T
            case 1
                T_plot = [70,71]; %in ms 
            case 2
                T_plot = [124,125];
%             case 3
%                 T_plot = [50,100];
%             case 4
%                 T_plot = [100,150];
%             case 5
%                 T_plot = [150,200];
        end
%         subplot(1,5,ind_T)
      figure(ind_T)
        cfg = [];
        cfg.parameter = 'avg';
        cfg.layout='CTF275.lay';
        cfg.xlim=[T_plot(1), T_plot(2)]';
        cfg.marker = 'none';
        cfg.interactive = 'yes';
        cfg.colorbar = 'no';
        ft_topoplotER(cfg, subject_data); title(sprintf('%dms -- %dms', T_plot(1), T_plot(2)));

    end
    end 
   
%%
%     if store_output == 1
        mkdir(fullfile('..','Results',out_folder,'ToneTopography_NewSubjects'))
        savefig(fullfile('..','Results',out_folder,'ToneTopography_NewSubjects',sprintf('GLOBAL_%s_Topography_%.2f-T%.2f_HP%d_LP%d.fig', baseline, T_init, T_end,hpfreq,lpfreq)));
%     end

 
%% Look at the data subject by subject
 mean_subjects_data(:, subject_list_short,:) = squeeze(rms(all_subject_channel_40(:,:, subject_list_short, :),1));   
 for subject_ind = 1:length(subject_list)
 
        subject_data.time = (1:150)/150*250-50;
%       subject_data.label = timelock.label;
        subject_data.dimord = 'chan_time';       
        %We load and plot the temporal data.
        figure(subject_ind);
%         subplot(1,9,subject_ind)
     
        Diff=mean_subjects_data(:,subject_ind,1) - mean_subjects_data(:,subject_ind,2); 

        % We use bootstrap and compute if there is a significant difference
        % between conditions.
        perc = 0.05;
        perc2 = 0.01;
        dataB=bootstrap(Diff'); 
        s=findSigDiff(dataB, perc);
        s2=findSigDiff(dataB, perc2);
        k1 = 1;
        k2 = 2;

        condi = mean_subjects_data(:,subject_ind,:)*100;
        plot(subject_data.time, condi(:,1), 'color', 'k', 'Linewidth',3);   % RAND
        hold on
        plot(subject_data.time, condi(:,2), 'color', 'r', 'Linewidth',3);   % REG

        hold on
        plot(subject_data.time, k1*abs(s),'Linewidth', 3);
        hold on
        plot(subject_data.time, k2*abs(s2),'Linewidth', 3);
        xlabel('Time (ms)')
        ylabel('RMS magnitude (fT)')
        title(sprintf('subjects_%d. Baseline: %s',subject_list(subject_ind), baseline));
        legend('RAND','REG','p=0.05','p=0.01')
 
 end
 
    