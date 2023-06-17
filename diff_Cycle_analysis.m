%% This script is to analyse the response differences of cycle 2 to cycle 6 related to cycle 1
% Edited by Mingyue Hu, 27/02/2023

%    define the time interval of interest
clear all; clc;
TOI = [0 2.5;
       2.5 5;
       5 7.5;
       7.5 10;
       10 12.5;
       12.5 15]; %time interval of interest
%parameters for the files
output_appendix = '_toneOnset';
compute_data = 0;
plot_data = 1;

if compute_data

       for t_ind = 1:length(TOI)
          cycle           = TOI(t_ind,:);
          T_init          = cycle(1);
          T_end           = cycle(2);
        
          %load the data matrix of 40 channels
            load(fullfile('D:\Results','Trigger_analysis_PRE_HP2_LP30','Tone_Trials_BLrawdata','Cycle1-Cycle6',...
            sprintf('allsub_40chann_toneresponse_TRIG_10_20-Time_%s_%s-BL%s.mat',mat2str(T_init),mat2str(T_end),...
            output_appendix)),'stable_average_shape'); 

%             cycleData =  squeeze(rms(stable_average_shape,1));
%             save(fullfile('D:\Results','Trigger_analysis_PRE_HP2_LP30','Tone_Trials_BLrawdata','Cycle1-Cycle6',...
%             sprintf('allsub_rms_toneresponse_TRIG_10_20-Time_%s_%s-BL%s.mat',mat2str(T_init),mat2str(T_end),...
%             output_appendix)),'cycleData'); 
%        end 


            
            if t_ind == 1
                cycle1_RAN = squeeze(stable_average_shape(:,:,:,1));
                avgCycle1_RAN = squeeze(rms(cycle1_RAN,1));
                cycle1_REG = squeeze(stable_average_shape(:,:,:,2));
                avgCycle1_REG = squeeze(rms(cycle1_REG,1));
                clear stable_average_shape
            else
                cycle_RAN = squeeze(stable_average_shape(:,:,:,1));
                avgCycle_RAN = squeeze(rms(cycle_RAN,1));
                diffCycle_RAN = avgCycle_RAN-avgCycle1_RAN; %Calculate the differences between cycle 1 with cycle 2/3/4/5/6 in RAN
                cycle_REG = squeeze(stable_average_shape(:,:,:,2));
                avgCycle_REG = squeeze(rms(cycle_REG,1));  
                diffCycle_REG = avgCycle_REG-avgCycle1_REG; %Calculate the differences between cycle 1 with cycle 2/3/4/5/6 in REG

                clear stable_average_shape
                clear avgCycle_RAN
                clear avgCycle_REG
                clear cycle_RAN
                clear cycle_REG
                %save the data into one matrix
                diffCycle(t_ind-1,:,:,1) = diffCycle_RAN;
                diffCycle(t_ind-1,:,:,2) = diffCycle_REG;

                clear diffCycle_RAN
                clear diffCycle_REG
            end                   
       end 
     save(fullfile('D:\Results','Trigger_analysis_PRE_HP2_LP30','Tone_Trials_BLrawdata/',...
     sprintf('meanDiff_cycle1_with_restCycles_2-3-4-5-6_rms_TRIG_10_20.mat')),'diffCycle');
end 

%% we quickly plot the data 

load(fullfile('D:\Results','Trigger_analysis_PRE_HP2_LP30','Tone_Trials_BLrawdata/',...
     sprintf('meanDiff_cycle1_with_restCycles_2-3-4-5-6_rms_TRIG_10_20.mat')),'diffCycle');
    if plot_data
    
        TOI_diff = [2.5 5;
               5 7.5;
               7.5 10;
               10 12.5;
               12.5 15]; %time interval of interest
    
        time = (1:150)/150*250-50;
        time_length = 30:150;  
        control = zeros(1,121);
    
       for time_ind = 1:length(TOI_diff) %we plot from cycle 2
          cycle           = TOI_diff(time_ind,:);
          T_init          = cycle(1);
          T_end           = cycle(2);
      
        dataCycle = squeeze(diffCycle(time_ind,:,:,:));
        dataCycle_RAN = dataCycle(:,:,1);
        dataCycle_REG = dataCycle(:,:,2);
        
        REGDdiffavg = squeeze(mean(dataCycle_REG,2));
        RANdiffavg = squeeze(mean(dataCycle_RAN,2));
    
        REGdiffstd = std(dataCycle_REG,0,2)/sqrt(size(dataCycle_REG,2));
        RANdiffstd = std(dataCycle_RAN,0,2)/sqrt(size(dataCycle_RAN,2));
    
        figure;
    %     plot(time(time_length), dataCycle_RAN(time_length),'Linewidth', 5);
    %     hold on
    %     plot(time(time_length), dataCycle_REG(time_length),'Linewidth', 5);
    %     hold on 
    
        shadedErrorBar(time(time_length),RANdiffavg(time_length),RANdiffstd(time_length),'lineProps','c');
        hold on 
        shadedErrorBar(time(time_length),REGDdiffavg(time_length),REGdiffstd(time_length),'lineProps','m'); hold on;
        ylim([-5,5]);
        plot(time(time_length), control,'--');
    
        xlabel('Time (ms)')
        ylabel('Magnitude differences(fT)')
        title(sprintf('Average. Baselined. cycle-%s relative to cycle 1. Time interval %s-%s',mat2str(time_ind+1), mat2str(T_init), mat2str(T_end)));
        legend('RANdiff','REGdiff');
        box off
       end 

end 