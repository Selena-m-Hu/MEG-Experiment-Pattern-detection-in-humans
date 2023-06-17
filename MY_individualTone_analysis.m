%% Individual Tone response analysis

%------The purpose of this script-------
% find and epoch each individual tone
% compute the mean of the response of M100 time interval
%----Edited by Mingyue Hu, 22/02/2023

clear all;clc;
%% Parameter
hpfreq          = 2;
lpfreq          = 30;
baseline        = '_toneOnset'; %select the baseline scheme： 'toneOnset'； 'none'
fs              = 600; % Sampling frequency
pre_estim       = 0.2; % Prestimuli time (in second). In our case, 200ms.
window_size     = 0.250*fs; % 250 ms, the length of the stimuli (50ms signal + 200ms silence).
trigger_list    =[10 20];
subject_list = [2 3 4 5 6 7 8 9 10 11 12 13 15 16 17 18 19 20 21 22 23 24];
addpath('D:\fieldtrip-20220707'); 
compute = 1;
plot = 0;

if compute
 for trigger_ind = 1:length(trigger_list)
   
    for subject_ind = 1:length(subject_list)

    % load the data
        load(fullfile('D:','Results','Trigger_analysis_PRE_HP2_LP30','Preprocessed_data_AllChannels',sprintf('data_subject-TRIG_%d-SUBJ_%d.mat',...
        trigger_list(trigger_ind), subject_list(subject_ind))),'data_subject'); 
        timelock = ft_timelockanalysis([],data_subject);

    % load the selected channels
        load(fullfile('D:\MEGGAP\Channels_DSS',sprintf('Channels-SUBJ_%d',subject_list(subject_ind))),...
        'channels', 'channels_num');
         counter = 1;
                
            for t_ind = 1:window_size:length(timelock.time)
                % We fill the first one (which will not be used) with
                % zeros, since it will be incomplete. This way, the
                % structure of the data is easier and we consider a 50ms
                % to baseline and 200ms of signal.
                % (so the tone epoch trial structure is pre-stim(50 ms) + tone(50ms) + silence(150ms) =
                % 150 ms). 
                % Here we use the pre-tone interval activity to deliver the 'tone' activity baseline 
                
                if t_ind == 1                 
                    trigger_shape(:,:,counter) = zeros(length(channels_num),150); %first counter include no information
                else  %[-0.2,0] include no information,the tone onset corresponding to time point of 120, 90 is -0.5 sec(no signal), followed with 200 ms signal
                    trigger_shape(:,:,counter) = timelock.avg(channels_num,t_ind-.1*fs:t_ind+0.15*fs-1); %second counter(first tone), information should start from -0.05
                end
                counter = counter + 1;
            end
            
%            Multiply by a constant in order to get the units into
            % femtoTeslas.
            trigger_shape = trigger_shape*1e15;
            
            %% Compute the individual tone activity
            % The first tone information was carried by the second counter
            % We have 60 tones in total, the first tone start from the
            % second counter
            target_trigger = trigger_shape(:,:,2:61);
                      
            %% We baseline the data according to the criteria we want to use.
            % 'tone onset': each tone is baselined using the activity of the onset of the tone.                       
            bl_window = 30;  % define the time window you want to use for tone baseline
            bl_tone = '-30'; %useful for saving file name 
            baseline_data = target_trigger(:,bl_window,:); %generate the baseline matrix

            for bl_ind = 1:length(baseline_data(1,1,:)) %loop the baseline tone by tone 
            baselined_tone_allsub(:,:,bl_ind,subject_ind,trigger_ind) = target_trigger(:,:,bl_ind) - repmat(baseline_data(:,:,bl_ind),1,150);
            end 
           

%             clear trigger_shape
            clear target_trigger
            clear baseline_data
            clear data_subject
            clear timelock
            clear channels_num
    
    end  
  
 end                
            save(fullfile('D:\Results','Trigger_analysis_PRE_HP2_LP30','Tone_Trials_BLrawdata/',...
            sprintf('allsub_40chann_individualTone_TRIG_10_20-Baselined.mat')),'baselined_tone_allsub'); 
            clear baselined_tone_allsub  
end 

%% quick plot
TI = [90,200]; %M100 response time interval; in ms [90 ms to 200 ms]
fs = 600;
TIcorr = (TI(1)*fs/1000):(TI(2)*fs/1000);%convert to MEG sample unit

if plot

   load(fullfile('D:\Results','Trigger_analysis_PRE_HP2_LP30','Tone_Trials_BLrawdata/',...
   sprintf('allsub_40chann_individualTone_TRIG_10_20-Baselined.mat')),'baselined_tone_allsub'); 
   % baselined_tone_allsub:(d1:channels; d2:data points; d3:tone position;
   % d4:subject number; d5:condition type)
   
   %condition RND 
    RAN_data = baselined_tone_allsub(:,:,:,:,1);
    powerRAN = squeeze(rms(RAN_data,1));
    M100RAN = squeeze(mean(powerRAN(TIcorr,:,:),1));
   %condition REG
    REG_data = baselined_tone_allsub(:,:,:,:,2);
    powerREG = squeeze(rms(REG_data,1));
    M100REG = squeeze(mean(powerREG(TIcorr,:,:),1));

    m100_tone = [M100RAN' M100REG'];
    xlswrite('tone_M100_90-200ms_RAN-REG.xlsx', m100_tone); %save the data

    

    %stats/bootstrap
    perc = 0.05;
    perc2 = 0.01;
    dataB=bootstrap(Diff'); 
    s=findSigDiff(dataB, perc);
    s2=findSigDiff(dataB, perc2);
    k1 = 1;
    k2 = 3;
    
    %mean and standard error
    REG = mean(REGall,2);
    RAND = mean(RANall,2);
   
    REGstd = std(REGall')/sqrt(size(REGall,2));
    RANstd = std(RANall')/sqrt(size(RANall,2));
    
    %Making plot
    time = (1:150)/150*250-50;
    timeind = 30:150;    
    figure;
    shadedErrorBar(time(timeind),RAND(timeind),RANstd(timeind),'lineProps','k');
    hold on 
    shadedErrorBar(time(timeind),REG(timeind),REGstd(timeind),'lineProps','r');
    xlim([-0.5,200])
    
    hold on
    plot(time(timeind), k1*abs(s(timeind)),'Linewidth', 12);
    hold on
    plot(time(timeind), k2*abs(s2(timeind)),'Linewidth', 12);
    xlabel('Time (ms)')
    ylabel('Magnitude (fT)')
    title(sprintf('All subjects. Baseline: %s. Time interval %s-%s', 'baselined', mat2str(tonePosition(1)), mat2str(tonePosition(10))));


    
%    %save the data mat
    M100_TONE(:,:,1) = permute(M100RAN,[2,1]);
    M100_TONE(:,:,2) = permute(M100REG,[2,1]);
    save(fullfile('D:\Results','Trigger_analysis_PRE_HP2_LP30','Tone_Trials_BLrawdata/',...
    sprintf('M100_Tone_response_allSubjects_1-RAN_2-REG.mat')),'M100_TONE'); 

    %compute the average tone response in each cycle
    %save the data matrix
    ranCycle(:,1) = mean(M100_TONE(:,2:10,1),2);
    ranCycle(:,2) = mean(M100_TONE(:,11:20,1),2);
    ranCycle(:,3) = mean(M100_TONE(:,21:30,1),2);
    ranCycle(:,4) = mean(M100_TONE(:,31:40,1),2);
    ranCycle(:,5) = mean(M100_TONE(:,41:50,1),2);
    ranCycle(:,6) = mean(M100_TONE(:,51:60,1),2);

    regCycle(:,1) = mean(M100_TONE(:,2:10,2),2);
    regCycle(:,2) = mean(M100_TONE(:,11:20,2),2);
    regCycle(:,3) = mean(M100_TONE(:,21:30,2),2);
    regCycle(:,4) = mean(M100_TONE(:,31:40,2),2);
    regCycle(:,5) = mean(M100_TONE(:,41:50,2),2);
    regCycle(:,6) = mean(M100_TONE(:,51:60,2),2);
    
%% compute the diff response between rest cycles with cycle 1
   %compute RAND
    RANdiffCycle(:,1) = ranCycle(:,2)-ranCycle(:,1);
    RANdiffCycle(:,2) = ranCycle(:,3)-ranCycle(:,1);
    RANdiffCycle(:,3) = ranCycle(:,4)-ranCycle(:,1);
    RANdiffCycle(:,4) = ranCycle(:,5)-ranCycle(:,1); 
    RANdiffCycle(:,5) = ranCycle(:,6)-ranCycle(:,1);

   %compute REG
    REGdiffCycle(:,1) = regCycle(:,2)-regCycle(:,1);
    REGdiffCycle(:,2) = regCycle(:,3)-regCycle(:,1);
    REGdiffCycle(:,3) = regCycle(:,4)-regCycle(:,1);
    REGdiffCycle(:,4) = regCycle(:,5)-regCycle(:,1); 
    REGdiffCycle(:,5) = regCycle(:,6)-regCycle(:,1);


    %% compute the mean response differences (averaged across subjects)
%     diff_power = mean(M100_TONE(:,:,1),1)-mean(M100_TONE(:,:,2),1);
%     M100_diff = mean(diff_power(TIcorr,:),1);

    M100_diff = M100_TONE(:,1)-M100_TONE(:,2);
%     stdM100 = std(M100_diff',0);
%     meanM100 = mean(M100_diff,2);
%     errorbar(meanM100,stdM100,'x'); hold on;
    %% Run linear regression and stats
    b = 1:60;
    M100_diff = M100_diff(14:60);
    mdl=fitlm(b,M100_diff);
    time =(1:150)/150*250-50;

    figure;
    plot(mdl);
    box off;
    xlabel('Tone (order)')
        ylabel('Power differences (RAN-REG)');
        title(sprintf('M100 mean response difference(baselined), all subjects(mean)'));
        xticks([14:60]); hold on;
        b=zeros(1,61);
        plot(0:60,b);  %plot the line where indicate no differences
        anova(mdl,'summary')
   

    %% Run two-way anova
    a = M100_TONE(:,14:60,1); %select the tones after 13 tones(1.3 cycle)
    a(23:44,:) = M100_TONE(:,14:60,2); % build the anova input data matrix
    nmbsub = 22; % Number of subjects from each condition, i.e., number of replications
    [~,~,stats] = anova2(a,nmbsub);
    c = multcompare(stats);
    tbl = array2table(c,"VariableNames", ...
    ["Group A","Group B","Lower Limit","A-B","Upper Limit","P-value"]);
    
end 



%% Run one sample t test for relative mean response from cycle 2 to cycle 6
% include the multiple correction
xlsread('diffPower_c2c3c4c5c6_RAN_REG_M100_100-200ms.xlsx');

[h,p,ci,stats] = ttest(data,0)






