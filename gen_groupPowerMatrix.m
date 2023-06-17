%% This script is to generate the power matrix of group from the ft data
% compute the mean response within the time interval of interest
% Edited by Mingyue Hu, 14/04/2023

clear all; clc;
% addpath('D:\fieldtrip-20220707'); 
DSS = 0;
trigger_list = [5, 15]; 
fs = 600; 

%% compute the group data matrix
computeMat = 1; % 1: compute the data matrix

if computeMat
  subject_list = [2 3 4 5 6 7 8 9 10 11 12 13 15 16 17 18 19 20 21 22 23 24];

    for subject_ind = 1:length(subject_list)
       
        % load channels
        load(fullfile('D:\MEGGAP\Channels_DSS',...
        sprintf('Channels-SUBJ_%d',subject_list(subject_ind))),'channels', 'channels_num');
            
        for trigger_ind = 1:length(trigger_list)
         
         if DSS
         input_folder = 'D:\Results\Trigger_analysis_PRE_HP0_LP2\DSS\DSStransformed\visualRejection';
         n_components = 3; 
    
         load(fullfile(input_folder,sprintf('DSSdata_subject-TRIG_%d-SUBJ_%d-COMP_%d.mat',...
         trigger_list(trigger_ind), subject_list(subject_ind), n_components)),'dss_data_subject');
         
         cfg = [];
         % We compute the average and other statistics from the trials using
         % ft_timelockanalysis.
         timelock_dss = ft_timelockanalysis(cfg, dss_data_subject);
         timelock_dss.avg = timelock_dss.avg(:,1:9720); 
%          timelock_dss.avg = timelock_all(:,:,trigger_ind,subject_ind);
         powerMatrix(:,trigger_ind,subject_ind) = rms(timelock_dss.avg(channels_num,:),1)*1e15;  %generate and save rms power for each subject in each condition
         
         clear timelock_dss
         clear dss_data_subject
         clear channels_num

         else 
         
         input_folder = 'D:\Results\Trigger_analysis_PRE_HP0_LP30\Preprocessed_data_AllChannels';

         load(fullfile(input_folder, sprintf('data_subject-TRIG_%d-SUBJ_%d.mat',trigger_list(trigger_ind),subject_list(subject_ind))),'timelock');
%          timelock.avg = timelock.avg(:,1:9720); 
         powerMatrix(:,trigger_ind,subject_ind) = rms(timelock.avg(channels_num,:),1)*1e15;  %generate and save rms power for each subject in each condition

         clear timelock
        
         end 

        end
     
          clear channels_num
    end

%% save the group data matrix
    if DSS
    save(fullfile(input_folder, sprintf('groupDSS_N22_power-TRIG_%d_%d.mat',trigger_list(1), trigger_list(2))), 'powerMatrix')
    else
    save(fullfile(input_folder, sprintf('groupRAW_N22_power-TRIG_%d_%d.mat',trigger_list(1), trigger_list(2))), 'powerMatrix')
    end 

end 
%% make quick plots
plotData = 0; 

if plotData
all_subjects_rms = powerMatrix;
time = [-200:1000/600:4000];
% time = [-200:1000/600:15999.9];
% time = [-200:1000/600:16000];
timelock.time = time/1000;

 for trigger_ind = 1:length(trigger_list)    
     switch trigger_list(trigger_ind)
        case 5
            out_5 = all_subjects_rms(:,1,:); %rms for selected channels for each subject
            time_5 = timelock.time;
            out_5 = squeeze(out_5);
        case 10
            out_10 = all_subjects_rms(:,1,:); %rms for selected channels for each subject
            time_10 = timelock.time;
            out_10 = squeeze(out_10);

        case 15
            out_15 = all_subjects_rms(:,2,:); %rms for selected channels for each subject
            time_15 = timelock.time;
            out_15 = squeeze(out_15);

        case 20
            out_20 = all_subjects_rms(:,2,:); %rms for selected channels for each subject
            time_20 = timelock.time;
            out_20 = squeeze(out_20);
     end
  end  
     
  % Plot data
    figure;
        shadedErrorBar(time_10(1,:),mean(out_10,2),(std(out_10')/sqrt(size(out_10,2))),'lineProps','k');
        hold on 
        shadedErrorBar(time_20(1,:),mean(out_20,2),(std(out_20')/sqrt(size(out_20,2))),'lineProps','r');

        xlabel('Time (s)');
        ylabel('Magnitude (fT)');  
        xlim([min(time_10(1,:)), max(time_10(1,:))])
%         legend('RAND','REG')
        title(sprintf('Timelock. LONG. GLOBAL, N=22'))
                
        hold on
        %Bootstrap analysis
          diffREGRAN= out_10 - out_20; 
          dataB=bootstrap(diffREGRAN'); 
          s=findSigDiff(dataB, 0.01);
          s1=findSigDiff(dataB, 0.05);
          
          s_out = proc_DiffPruning(s, trigger_list);
          s_out1 = proc_DiffPruning(s1, trigger_list);
         
          plot(timelock.time,40*abs(s),'color',[0.9290 0.6940 0.1250],'Linewidth', 12, 'color','k');
          plot(timelock.time,60*abs(s1),'color',[0.8500 0.3250 0.0980],'Linewidth', 12, 'color','r');

end 

%% Compute the mean response within certain time interval
    
    computeMean = 0; %1: compute the mean response
    
    if computeMean
       
        %load the group data matrix
        load(fullfile(input_folder, sprintf('groupDSS_N22_power-TRIG_%d_%d.mat',trigger_list(1), trigger_list(2))), 'powerMatrix')
    
        timeWindow = [3.27,7.69]; % in second
    
        powerMatrix = mean(powerMatrix(timeWindow(1)*fs:timeWindow(2)*fs,:,:),1);
        avgPower = squeeze(powerMatrix);
        avgPower = avgPower';
        
        xlswrite(fullfile(input_folder, sprintf('meanResponse_3.27_7.69sec_longSequence_RAN_REG.xlsx')),avgPower); 
    end 


%% quick plot

    % if short
    addpath('C:\Users\i7 System\Desktop\Online_studies\ONLINE_DATA\beeswarm');
    num = avgPower;
    % Beeswarm spread dots plot
    %how many groups of data do you have?
    y = num(:); %stack all columns of data in one column
    NSub= size(num);
    x = zeros(NSub(1),NSub(2));
    x = zeros(NSub(1),1);
    for i = 1:NSub(2)
        x(:,i) = i
    end
    x = x(:);  %stack the group numbers in on column 
    
    figure;
    %plot the violin
    violin(num,'edgecolor','',...
    'mc','k--',...
    'medc','')
    hold on;
    %plot the individuals
    beeswarm(x,y,'dot_size',1,'sort_style','up'); hold on
    ax = gca;
    ax.LineWidth = 1;
    xlabel('');   %input the x label
    xticks([1 2]); hold on
    ylim([-30 280]);
    yticks([0 100 200 300 400]); %for step reaction time
    
    tick = {'RAND' 'REG'};
    set(gca,'XTickLabel',{'RAND' 'REG'});
    set(gca,'XTickLabelRotation',25,'FontSize',12,'FontWeight','bold'); hold on  %Set x label with certain angle
    ylabel('Magnitude(fT)','FontSize',14)
    box off
    title('Slow');
  
    hold on;
    mean_all_TOI = avgPower;
    % compute the mean and error bar
    mRAND = mean(mean_all_TOI(:,1));
    sRAND = std(mean_all_TOI(:,1));
    mREG  = mean(mean_all_TOI(:,2));
    sREG  = std(mean_all_TOI(:,2));
    
    bmRAN = plot(1.3,[mRAND],'.','LineWidth',2,'MarkerSize',20,'color',[0.15,0.15,0.15]); hold on;
    bmREG = plot(2.3,[mREG],'.','LineWidth',2,'MarkerSize',20,'color',[0.15,0.15,0.15]); hold on;
    eRAN = errorbar(1.3,mRAND,sRAND,'.','LineWidth',2,'color',[0.15,0.15,0.15]);
    eREG = errorbar(2.3,mREG,sREG,'.','LineWidth',2,'color',[0.15,0.15,0.15]);

