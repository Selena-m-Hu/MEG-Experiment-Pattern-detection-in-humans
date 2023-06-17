    %% Temporal plot for new subjects (Subject16 to subject 24)
    % Scripts designed to plot the timelock information for all the
    % triggers (5, 10, 15, 20) for a specific subject (or group of 
    % subjects). Esentially, it plots two separate graphs for LONG and 
    % SHORT, considering the modalities REG and RAND for each one of the
    % configurations. NO TIMELOCK DATA IS COMPUTED IN THIS SCRIPT.
    % 
    % trigger_list: formed by the elements {5, 10, 15, 20]. 
    % * 5 : 3 second RAND sequences.
    % * 10: 15 second RAND sequences.
    % * 15: 3 second REG sequences.
    % * 20: 15 second REG sequences.
    % {Subject 9, the 8th subject is the outlier, whose performance (RMS) in long
    % REG is extremly high during the last half time of the trial.}
    % Last update: 07/June/2018
    % Adapted by Mingyue Hu, Aug, 2022
    
    clear all;
    clc
%     in_folder = 'Trigger_analysis_PRE_HP0_LP30';  % for loading/retriving channels information
    out_folder = 'Trigger_analysis_PRE_HP0_LP30';  % for loading time_lock data information
%     store_output = config.store_data;
    subject = [2 3 4 5 6 7 8 9 10 11 12 13 15 16 17 18 19 20 21 22 23 24]; 
%     trigger_list = [5,15];
    trigger_list = [5,15];
%     plotDSS = 1; %plot DSS data
    
    %% We create some matrices with the timelock information.
%    %load the timelock data from all new subjects
%     load(fullfile('..','Results',out_folder,'DSS timelock',...
%     sprintf('Channels-TRIG_%d-allSUBJ',trigger_list(1))), 'all_RAND_data');     
%     load(fullfile('..','Results',out_folder,'DSS timelock',...
%     sprintf('Channels-TRIG_%d-allSUBJ',trigger_list(2))), 'all_REG_data');
     load('D:\Results\Trigger_analysis_PRE_HP0_LP30\Preprocessed_data_AllChannels\data_subject-TRIG_5-SUBJ_2.mat');
     time = data_subject.time{1};
     NKEEP = 3;
     for subject_ind = 1:length(subject)    
         %load selected channels
      
       %load DSS channels
        load(fullfile('D:\MEGGAP\Channels_DSS',...
        sprintf('Channels-SUBJ_%d',subject(subject_ind)))...
        ,'channels', 'channels_num');

        load(fullfile('D:\Results',out_folder,'DSS_components','DSS_data',sprintf('Xdss-TRIG_%d-SUB_%d-COMP_%d_NBprocedure.mat',...
            trigger_list(1),subject(subject_ind), NKEEP)),'DSS_timelock1');
        DSS_timelock1=permute(DSS_timelock1,[2,1]);
          
        load(fullfile('D:\Results',out_folder,'DSS_components','DSS_data',sprintf('Xdss-TRIG_%d-SUB_%d-COMP_%d_NBprocedure.mat',...
            trigger_list(2),subject(subject_ind), NKEEP)),'DSS_timelock2');
        DSS_timelock2=permute(DSS_timelock2,[2,1]);
 
        out_5(subject_ind,:) = rms(DSS_timelock1(channels_num,:)); %rms for selected channels for each subject
        time_5(subject_ind,:) = time;
  
        out_15(subject_ind,:) = rms(DSS_timelock2(channels_num,:));
        time_15(subject_ind,:) = time;

        clear DSS_timelock1
        clear DSS_timelock2

                           
    end

    %% We plot data.
    
        % SHORT condition (3 seconds long).
        figure;
        shadedErrorBar(time_5(1,:),mean(out_5,1)*1e15,(std(out_5,1)/sqrt(size(out_5,1)))*1e15,'lineProps','r');
        hold on 
        shadedErrorBar(time_15(1,:),mean(out_15,1)*1e15,(std(out_15,1)/sqrt(size(out_15,1)))*1e15,'lineProps','k');
        
%         f1 = figure;
%         plot(time_5(1,:),mean(out_5,1)*1e15,'color','k','Linewidth', 1.5); hold on;
%         plot(time_15(1,:),mean(out_15,1)*1e15,'color','r','Linewidth', 1.5)
        xlabel('Time (s)');
        ylabel('Magnitude (fT)');
        xlim([min(time_5(1,:)), max(time_5(1,:))])
        legend('RAND','REG')
        if length(subject) == 1
            title(sprintf('Timelock. SHORT. SUBJ: %d',subject))
        else
            title(sprintf('Timelock, 0-30 Hz. SHORT. GLOBAL N=22'))
        end
        %% Bootstrap analysis
        
         switch trigger_list(1) 
           case 5
           diffREGRAN= out_5 - out_15; 
           case 10
           diffREGRAN= out_10 - out_20; 
         end
          hold on
          dataB=bootstrap(diffREGRAN); 
          s=findSigDiff(dataB, 0.01);
          s1=findSigDiff(dataB, 0.05);
          plot(timelock.time,3*abs(s),'color',[0.9290 0.6940 0.1250],'Linewidth', 3);
          plot(timelock.time,15*abs(s1),'color',[0.8500 0.3250 0.0980],'Linewidth', 3);

%         % LONG condition (15 seconds long).
%         subplot(212)
%         f1 = figure;
%         plot(time_10(1,:),mean(out_10,1)*1e15,'color','k','Linewidth', 1.5); hold on;
%         plot(time_20(1,:),mean(out_20,1)*1e15,'color','r','Linewidth', 1.5)
%         xlabel('Time (s)');
%         ylabel('Magnitude (fT)');  
%         xlim([min(time_10(1,:)), max(time_10(1,:))])
%         legend('RAND','REG')
%         if length(subject) == 1
%             title(sprintf('Timelock. LONG. SUBJ: %d',subject))
%         else
%             title(sprintf('Timelock. LONG. GLOBAL'))
%         end      
%     
%         plot(time_10(1,:),mean(out_10,1)*1e15,'color','k','Linewidth', 1.5); hold on;
%         plot(time_20(1,:),mean(out_20,1)*1e15,'color','r','Linewidth', 1.5)
%         xlabel('Time (s)');
%         ylabel('Magnitude (fT)');  
%         xlim([min(time_10(1,:)), max(time_10(1,:))])
%         legend('RAND','REG')
%         if length(subject) == 1
%             title(sprintf('Timelock. LONG. SUBJ: %d',subject))
%         else
%             title(sprintf('Timelock. LONG. GLOBAL'))
%         end      
%     end
%    mkdir(fullfile('..','Results',out_folder,'DSS timelock','allRMSplot'));    
%    plotFolder = fullfile('..','Results',out_folder,'DSS timelock','allRMSplot'...
%     ,sprintf('Power-TRIG_%d_%d-allSUBJ',trigger_list(1),trigger_list(2)));
%     savefig(f1,plotFolder)
    
   
   

    %% We can store the block information in a .fig file.
%     if store_output == 1
%         mkdir(fullfile('..','Results',out_folder,'Timelock_graphs'))
%         if length(subject) == 1
%             savefig(fullfile('..','Results',out_folder,'Timelock_graphs',sprintf('Timelock-SUBJ_%d.fig',subject)))
%         else
%             savefig(fullfile('..','Results',out_folder,'Timelock_graphs',sprintf('Timelock-GLOBAL.fig')))
%         end
%     end
%     

    
% end
for subject_ind = 1:length(subject)    
         %load selected channels
      
       %% Trigger 1
        load(fullfile('D:\Results',out_folder,'DSS_components','DSS_data',sprintf('Xdss-TRIG_%d-SUB_%d-COMP_%d_NBprocedure.mat',...
        trigger_list(1),subject(subject_ind), NKEEP)),'DSSdata_subject1','DSS_timelock1');
     
        DSS_timelock1 = DSS_timelock1(1:2520,:);  
        for i = 1:length(DSSdata_subject1(1,:))
        DSSdata_subject1{i} = DSSdata_subject1{i}(:,1:2520);
        end 
        save(fullfile('D:\Results',out_folder,'DSS_components','DSS_data',sprintf('Xdss-TRIG_%d-SUB_%d-COMP_%d_NBprocedure.mat',...
        trigger_list(1),subject(subject_ind), NKEEP)),'DSSdata_subject1','DSS_timelock1');
        clear DSSdata_subject1
        clear DSS_timelock1
     %% Trigger 2

        load(fullfile('D:\Results',out_folder,'DSS_components','DSS_data',sprintf('Xdss-TRIG_%d-SUB_%d-COMP_%d_NBprocedure.mat',...
        trigger_list(2),subject(subject_ind), NKEEP)),'DSSdata_subject2','DSS_timelock2');
        DSS_timelock2 = DSS_timelock2(1:2520,:);  

        for i = 1:length(DSSdata_subject2(1,:))
        DSSdata_subject2{i} = DSSdata_subject2{i}(:,1:2520);
        end 
        save(fullfile('D:\Results',out_folder,'DSS_components','DSS_data',sprintf('Xdss-TRIG_%d-SUB_%d-COMP_%d_NBprocedure.mat',...
        trigger_list(2),subject(subject_ind), NKEEP)),'DSSdata_subject2','DSS_timelock2');

        clear DSSdata_subject2
        clear DSS_timelock2
       
 end        
      
