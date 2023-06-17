%% Analyse and calculate the mean response of sustained response
% Edited by Mingyue Hu, 02/03/2023

clc;
clear all;

subject_list = [2 3 4 5 6 7 8 9 10 11 12 13 15 16 17 18 19 20 21 22 23 24]; 
% trigger_list = [5 15];
trigger_list = [10 20];
fs = 600; 

for trigger_ind = 1:length(trigger_list)

    for subject_ind = 1:length(subject_list)

     %load 40 selected channels
     load(fullfile('D:\MEGGAP\Channels_DSS',sprintf('Channels-SUBJ_%d',subject_list(subject_ind))),...
    'channels', 'channels_num');

     %load data
     % we analyse DSS data for long sequence; raw data for short sequence
     load(fullfile('D:\Results\0-2_slow_sequence_DSS',...
     sprintf('dataDSS_subject-TRIG_%d-SUBJ_%d.mat',trigger_list(trigger_ind), subject_list(subject_ind))),'dataDSS_subject');
     
     %compute the rms of selected 49 channels
     data_power = rms(dataDSS_subject(channels_num,:),1);
     all_power(:,subject_ind,trigger_ind) = data_power; 
    
    end 

end 

%% compute the mean response of time interval of interest
timeWindow = [3.27,7.69]; % in second

TOI_avgpower = mean(all_power(timeWindow(1)*fs:timeWindow(2)*fs,:,:),1);
avgPower = squeeze(TOI_avgpower);

save('meanResponse_8-15sec_longSequence.mat','avgPower'); 




