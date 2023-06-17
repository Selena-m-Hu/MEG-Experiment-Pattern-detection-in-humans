    %% Temporal plot for all subjects
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
    in_folder = 'Trigger_analysis_PRE_HP0_LP30';  %channels are saved in the 0-30Hz folder
    out_folder = 'Trigger_analysis_PRE_HP0_LP2';
%     store_output = config.store_data;
    subject = [2 3 4 5 6 7 8 9 10 11 12 13 15 16 17 18 19 20 21 22 23 24];
    % Look at subject 9, whose performance is the outlier in long trial REG
    % condition
%    trigger_list = [5,15];
    trigger_list = [10,20];
    plotDSS = 1; %plot DSS data
    
    %% We create some matrices with the timelock information.
   %load the timelock data from all new subjects
   
   %We load the old subjects file first
   load(fullfile('..','Results',out_folder,'DSS timelock',...
    sprintf('Channels-TRIG_%d-oldAllSUBJ',trigger_list(1))), 'all_RAND_data');     
    load(fullfile('..','Results',out_folder,'DSS timelock',...
    sprintf('Channels-TRIG_%d-oldAllSUBJ',trigger_list(2))), 'all_REG_data');
   
  % We use N_22_XX_Data variable to save all data
   N_22_RAND_Data = all_RAND_data;
   N_22_REG_Data = all_REG_data;
   
   
   %We load the new subjects data  
    load(fullfile('..','Results',out_folder,'DSS timelock',...
    sprintf('Channels-TRIG_%d-allSUBJ',trigger_list(1))), 'all_RAND_data');     
    load(fullfile('..','Results',out_folder,'DSS timelock',...
    sprintf('Channels-TRIG_%d-allSUBJ',trigger_list(2))), 'all_REG_data');
   
    %We combine the old subjects with the new subjects  
    N_22_RAND_Data = [N_22_RAND_Data all_RAND_data];
    N_22_REG_Data = [N_22_REG_Data all_REG_data];
    
    
   
   for trigger_ind = 1:length(trigger_list)
     for subject_ind = 1:length(subject)    
         %load selected channels
         if subject(subject_ind) < 16
           if plotDSS    
           %load DSS channels
            load(fullfile('..','Results_Antonio_S1_S15','Channels_DSS',...
            sprintf('Channels-SUBJ_%d',subject(subject_ind)))...
            ,'channels', 'channels_num');
           else
           %load channels
%             load(fullfile('..','Results_Antonio_S1_S15',in_folder,'Channels',...
%             sprintf('Channels-TRIG_%d_%d-SUBJ_%d',trigger_list(1), trigger_list(2),subject(subject_ind))),...
%             'channels', 'channels_num');   
           end
         else 
           if plotDSS  
           %load DSS channels
            load(fullfile('..','Results',in_folder,'Channels_DSS',...
            sprintf('Channels-TRIG_%d_%d-SUBJ_%d',trigger_list(1), trigger_list(2),subject(subject_ind)))...
            ,'channels', 'channels_num');
           else
           %load channels
            load(fullfile('..','Results',in_folder,'Channels',...
            sprintf('Channels-TRIG_%d_%d-SUBJ_%d',trigger_list(1), trigger_list(2),subject(subject_ind))),...
            'channels', 'channels_num');   
           end
         end             
             
           switch trigger_list(trigger_ind) 
                case 5
                    timelock = N_22_RAND_Data{subject_ind};
                    out_5(subject_ind,:) = rms(timelock.avg(channels_num,:)); %rms for selected channels for each subject
                    time_5(subject_ind,:) = timelock.time;
                case 10
                    timelock = N_22_RAND_Data{subject_ind};
                    out_10(subject_ind,:) = rms(timelock.avg(channels_num,:));
                    time_10(subject_ind,:) = timelock.time;
                case 15
                    timelock = N_22_REG_Data{subject_ind};
                    out_15(subject_ind,:) = rms(timelock.avg(channels_num,:));
                    time_15(subject_ind,:) = timelock.time;
                case 20
                    timelock = N_22_REG_Data{subject_ind};
                    out_20(subject_ind,:) = rms(timelock.avg(channels_num,:));
                    time_20(subject_ind,:) = timelock.time;
            end      
      end    
    end

    %% We plot data.
    
        %% SHORT condition (3 seconds long).        
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

        %% LONG condition (15 seconds long). 
        figure;
        shadedErrorBar(time_10(1,:),mean(out_10,1)*1e15,(std(out_10,1)/sqrt(size(out_10,1)))*1e15,'lineProps','r');
        hold on 
        shadedErrorBar(time_20(1,:),mean(out_20,1)*1e15,(std(out_20,1)/sqrt(size(out_20,1)))*1e15,'lineProps','k');
%         f1 = figure;
%         plot(time_10(1,:),mean(out_10,1)*1e15,'color','k','Linewidth', 1.5); hold on;
%         plot(time_20(1,:),mean(out_20,1)*1e15,'color','r','Linewidth', 1.5)
        xlabel('Time (s)');
        ylabel('Magnitude (fT)');  
        xlim([min(time_10(1,:)), max(time_10(1,:))])
        legend('RAND','REG')
        if length(subject) == 1
            title(sprintf('Timelock. LONG. SUBJ: %d',subject))
        else
            title(sprintf('Timelock. LONG. GLOBAL, N=21'))
        end      
       
        %% Plot performance of individuals 
        
        figure;
        
        plot(time_20(1,:),out_20*1e15)
        
   
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
  plot(timelock.time,0.3*abs(s),'color',[0.9290 0.6940 0.1250],'Linewidth', 3);
  plot(timelock.time,0.15*abs(s1),'color',[0.8500 0.3250 0.0980],'Linewidth', 3);
  
  
  %% Save figure
  mkdir(fullfile('..','Results',out_folder,'DSS timelock','allRMSplot'));    
   plotFolder = fullfile('..','Results',out_folder,'DSS timelock','allRMSplot'...
    ,sprintf('Power-TRIG_%d_%d-allSUBJ',trigger_list(1),trigger_list(2)));
    savefig(f1,plotFolder)

  
