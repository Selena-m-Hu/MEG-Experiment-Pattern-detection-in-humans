%% Compute the mean power of certain time window for individuals
% 0-2 Hz, [10 20] long sequence, time window [8-15S];
% 0-30 Hz, [5 15] short sequence, time window [2-3S];

clc;
clear all;

subject_list = [2 3 4 5 6 7 8 9 10 11 12 13 15 16 17 18 19 20 21 22 23 24]; 
trigger_list_short = [5 15];
trigger_list_long = [10 20];
fs = 600; 

%% Load the data   
short = 1; % 1: load the short data; 2: load the long data
if short
load('D:\Results\0-30_fast_slow_rawdata\allsub_power_RAN.mat');
load('D:\Results\0-30_fast_slow_rawdata\allsub_power_REG.mat');

timeWindow = [2,3];  % in second
else
load('D:\Results\Trigger_analysis_PRE_HP0_LP2\DSS_components\normalizedDSS_cfg\ALLchann_allSubj_timelock-TRIG_10-20-COMP_3.mat')
timeWindow = [8,15]; % in second
end 
mean_all_TOI = []; 
for s = 1:length(subject_list)
    
%load channels
load(fullfile('D:\Antonio\Results\Trigger_analysis_PRE_HP0_LP30\Channels_DSS',sprintf('Channels-SUBJ_%d',subject_list(s))),...
    'channels', 'channels_num');

%% Compute the mean response of the timewindow 
data_subject = squeeze(timelock_all(:,:,:,s)); 
data_subject = squeeze(rms(data_subject(channels_num,:,:),1)); %compute the rms based on the selected channels

TOI_data_subject = data_subject(timeWindow(1)*fs:timeWindow(2)*fs,:);
mean_TOI_subject = mean(TOI_data_subject(:,:),1)*1e15; 

%% Allocate mean data from all subject in one matrix

mean_all_TOI = [mean_all_TOI; mean_TOI_subject]; %column 1: RAND; column 2: REG

end 

%% Plot the data

% if short
addpath('C:\Users\i7 System\Desktop\Online_studies\ONLINE_DATA\beeswarm');
num = mean_all_TOI;
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
if short
title('Fast, 0-30 Hz');
else
title('Slow');
end 
hold on;
% compute the mean and error bar
mRAND = mean(mean_all_TOI(:,1));
sRAND = std(mean_all_TOI(:,1));
mREG  = mean(mean_all_TOI(:,2));
sREG  = std(mean_all_TOI(:,2));

bmRAN = plot(1.3,[mRAND],'.','LineWidth',2,'MarkerSize',20,'color',[0.15,0.15,0.15]); hold on;
bmREG = plot(2.3,[mREG],'.','LineWidth',2,'MarkerSize',20,'color',[0.15,0.15,0.15]); hold on;
eRAN = errorbar(1.3,mRAND,sRAND,'.','LineWidth',2,'color',[0.15,0.15,0.15]);
eREG = errorbar(2.3,mREG,sREG,'.','LineWidth',2,'color',[0.15,0.15,0.15]);



