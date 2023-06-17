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
    in_folder = 'Trigger_analysis_PRE_HP0_LP30';  % for loading/retriving channels information
    out_folder = 'Trigger_analysis_PRE_HP0_LP2';  % for loading time_lock data information
%     store_output = config.store_data;
    subject = [2 3 4 5 6 7 8 9 10 11 12 13 15 16 17 18 19 20 21 22 23 24]; 
%     trigger_list = [5,15];
    trigger_list = [5,15];
    plotDSS = 1; %plot DSS data
    
    %% We create some matrices with the timelock information.
%    %load the timelock data from all new subjects
%     load(fullfile('..','Results',out_folder,'DSS timelock',...
%     sprintf('Channels-TRIG_%d-allSUBJ',trigger_list(1))), 'all_RAND_data');     
%     load(fullfile('..','Results',out_folder,'DSS timelock',...
%     sprintf('Channels-TRIG_%d-allSUBJ',trigger_list(2))), 'all_REG_data');
     load('D:\Results\Trigger_analysis_PRE_HP0_LP30\Preprocessed_data_AllChannels\data_subject-TRIG_5-SUBJ_2.mat');
     time = 
   
     for subject_ind = 1:length(subject)    
         %load selected channels
           if plotDSS  
           %load DSS channels
            load(fullfile('D:\MEGGAP\Channels_DSS',...
            sprintf('Channels-TRIG_%d_%d-SUBJ_%d',trigger_list(1), trigger_list(2),subject(subject_ind)))...
            ,'channels', 'channels_num');

            load(fullfile('D:\Results',out_folder,'DSS_components','DSS_data',sprintf('Xdss-TRIG_%d-SUB_%d-COMP_%d_NBprocedure.mat',...
                trigger_list(1),subject_list(subject_ind), NKEEP)),'DSS_timelock1');
              
            load(fullfile('D:\Results',out_folder,'DSS_components','DSS_data',sprintf('Xdss-TRIG_%d-SUB_%d-COMP_%d_NBprocedure.mat',...
                trigger_list(2),subject_list(subject_ind), NKEEP)),'DSS_timelock2');
           
           else
           %load channels
            load(fullfile('..','Results',in_folder,'Channels',...
            sprintf('Channels-TRIG_%d_%d-SUBJ_%d',trigger_list(1), trigger_list(2),subject(subject_ind))),...
            'channels', 'channels_num');   
           end
           switch trigger_list(trigger_ind) 
                case 5
                    timelock = all_RAND_data{subject_ind};
                    out_5(subject_ind,:) = rms(timelock.avg(channels_num,:)); %rms for selected channels for each subject
                    time_5(subject_ind,:) = timelock.time;
                case 10
                    timelock = all_RAND_data{subject_ind};
                    out_10(subject_ind,:) = rms(timelock.avg(channels_num,:));
                    time_10(subject_ind,:) = timelock.time;
                case 15
                    timelock = all_REG_data{subject_ind};
                    out_15(subject_ind,:) = rms(timelock.avg(channels_num,:));
                    time_15(subject_ind,:) = timelock.time;
                case 20
                    timelock = all_REG_data{subject_ind};
                    out_20(subject_ind,:) = rms(timelock.avg(channels_num,:));
                    time_20(subject_ind,:) = timelock.time;
            end      
          
    end

    %% We plot data.
    
        % SHORT condition (3 seconds long).
        f1 = figure;
        plot(time_5(1,:),mean(out_5,1)*1e15,'color','k','Linewidth', 1.5); hold on;
        plot(time_15(1,:),mean(out_15,1)*1e15,'color','r','Linewidth', 1.5)
        xlabel('Time (s)');
        ylabel('Magnitude (fT)');
        xlim([min(time_5(1,:)), max(time_5(1,:))])
        legend('RAND','REG')
        if length(subject) == 1
            title(sprintf('Timelock. SHORT. SUBJ: %d',subject))
        else
            title(sprintf('Timelock. SHORT. GLOBAL'))
        end

        % LONG condition (15 seconds long).
        subplot(212)
        f1 = figure;
        plot(time_10(1,:),mean(out_10,1)*1e15,'color','k','Linewidth', 1.5); hold on;
        plot(time_20(1,:),mean(out_20,1)*1e15,'color','r','Linewidth', 1.5)
        xlabel('Time (s)');
        ylabel('Magnitude (fT)');  
        xlim([min(time_10(1,:)), max(time_10(1,:))])
        legend('RAND','REG')
        if length(subject) == 1
            title(sprintf('Timelock. LONG. SUBJ: %d',subject))
        else
            title(sprintf('Timelock. LONG. GLOBAL'))
        end      
    else
        plot(time_10(1,:),mean(out_10,1)*1e15,'color','k','Linewidth', 1.5); hold on;
        plot(time_20(1,:),mean(out_20,1)*1e15,'color','r','Linewidth', 1.5)
        xlabel('Time (s)');
        ylabel('Magnitude (fT)');  
        xlim([min(time_10(1,:)), max(time_10(1,:))])
        legend('RAND','REG')
        if length(subject) == 1
            title(sprintf('Timelock. LONG. SUBJ: %d',subject))
        else
            title(sprintf('Timelock. LONG. GLOBAL'))
        end      
    end
   mkdir(fullfile('..','Results',out_folder,'DSS timelock','allRMSplot'));    
   plotFolder = fullfile('..','Results',out_folder,'DSS timelock','allRMSplot'...
    ,sprintf('Power-TRIG_%d_%d-allSUBJ',trigger_list(1),trigger_list(2)));
    savefig(f1,plotFolder)
    
   
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
  plot(time,3*abs(s),'color',[0.9290 0.6940 0.1250],'Linewidth', 3);
  plot(time,15*abs(s1),'color',[0.8500 0.3250 0.0980],'Linewidth', 3);

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
