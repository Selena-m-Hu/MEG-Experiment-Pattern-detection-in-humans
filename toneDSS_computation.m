function toneDSS_computation(trigger_list, subject_list,time_frame,dsscomp,cut_epoch)
%% Tone DSS components computation 
% --edited by Mingyue Hu, 8/11/2022
% ----------parameter--------------
    T_init          = time_frame(1);
    T_end           = time_frame(2);
    hpfreq          = 2;
    lpfreq          = 30;
    out_folder      = sprintf('Trigger_analysis_PRE_HP%d_LP%d',hpfreq,lpfreq);
    fs              = 600; % Sampling frequency
    window_size     = 0.250*fs; % 250 ms, the length of the stimuli (50ms signal + 200ms silence).
    computeEpoch    = cut_epoch;
    cycle           = '_TOI';  %input the time interval of interest ('_First';'_TOI')
    baseline        = '_none'; %select the baseline scheme： 'toneOnset'； 'none'
    
 for subject_ind = 1:length(subject_list)
    
    if computeEpoch

      for trigger_ind = 1:length(trigger_list)
          
           %% Load the raw data 
            load(fullfile('..','Results',out_folder,'Preprocessed_data_AllChannels','short_timelock',...
            sprintf('Short_timelock-TRIG_%d-SUBJ_%d.mat',trigger_list(trigger_ind),subject_list(subject_ind))),'short_data_subject')
                                             
%           timelock = ft_timelockanalysis([],short_data_subject);             

            %% we cut the sequence into tone epoch trial by trial first
         for ti = 1:length(short_data_subject.trial)
              
            % Now, cut the sequence into individual epoch
            counter = 1;
            timelock_trial = short_data_subject.trial(ti);
            timelock_trial = cell2mat(timelock_trial);
            for t_ind = 1:window_size:length(cell2mat(short_data_subject.time(ti)))
                % We fill the first one (which will not be used) with
                % zeros, since it will be incomplete. This way, the
                % structure of the data is easier and we consider a 50ms
                % to baseline and 200ms of signal. (so the trial structure is pre-stim(50 ms) + tone(50ms) + silence(150ms) =
                % 150 ms). 
                %-----------------------------------------------------------
                if t_ind == 1             
                    trigger_shape275(:,:,counter) = zeros(length(short_data_subject.label),150);                  
                else  
                    trigger_shape275(:,:,counter) = timelock_trial(:,t_ind-.1*fs:t_ind+0.15*fs-1);                  
                end
                counter = counter + 1;
            end
            %% select the epochs from time interval of interest
                      
            % We get the value of the initial time (0 seconds).
            t0 = 0.25/(window_size/fs);
            
            % We get the temporal index of the beginning and the end of our
            % data from the global matrix.
            t_stable = max(t0,1) + ((T_init/(window_size/fs)):(T_end/(window_size/fs)));
            % Calculate the average epoch across all epochs of this time
            % interval 
            tone_average  = mean(trigger_shape275(:,:,t_stable),3);
            
            % Allocate the averaged data into ft data structure
            tone_subject.trial{ti} = tone_average;
                            
         end    
            % Transform the epoch data structure for DSS analysis                
%         tone_subject.trial = nt_mat2trial(trigger_shape275);         
           x_orig{trigger_ind} = cat(3,tone_subject.trial{:});          
           tone_subject.time = (1:150)/150*250-50;                           
    
      end
       %% save the epoched raw data of each subject        
        save(fullfile('..','Results',out_folder,'ToneDSS','DSS_components',...
        sprintf('toneData-TRIG_%d_%d_Time_%d_%d-SUBJ_%d.mat',trigger_list(1),trigger_list(2),round(T_init),round(T_end),...
        subject_list(subject_ind))),'x_orig'); 4

        clear trigger_shape275
        clear tone_subject
         
     else
        load(fullfile('D:\Results','Trigger_analysis_PRE_HP2_LP30','Tone_Trials_rawdata/',...
        sprintf('toneTrials-Cycle%s_TRIG_%d_%d_Time_%d_%d-SUBJ_%d_BL%s.mat',cycle,trigger_list(1),trigger_list(2),round(T_init),round(T_end),...
        subject_list(subject_ind),baseline)),'data_merged'); 
        data1 = data_merged.trial(:,data_merged.trialinfo==10);%put the RAND trial in the first cell of x_orig
        x_orig{1} = cat(3,data1{:}); 
        data2 = data_merged.trial(:,data_merged.trialinfo==20);%REG in second 
        x_orig{2} = cat(3,data2{:}); 
        tone_subject.time = (1:150)/150*250-50;
     end 
   
    if dsscomp
       %% Compute DSS components         
          %We reshape the data to fit the structure of DSS functions(NoiseTools)
                t = tone_subject.time;
                
                % In case of single, we compute it twice.
                x = cat(3,x_orig{:});
                x = permute(x,[2,1,3]);
                
                % c0: baseline covariance
                % c1: biased covariance
                c0 = nt_cov(x);
                c1 = nt_cov(mean(x,3));
                % DSS
                [todss,pwr0,pwr1] = nt_dss0(c0,c1);
                % todss: matrix to convert data to normalized DSS components
                % pwr0: power per component (baseline)
                % pwr1: power per component (biased)

                z=nt_mmat(x,todss); %z is the output data(components) of DSS analysis, 
                clear x
                % We plot the power of the components and the top10 components.
                f1 = figure(1); clf; plot(pwr1./pwr0,'.-'); ylabel('score'); xlabel('component');
                f2 = figure(2); clf;
                for iComp=1:10
                    subplot(2,5,iComp);
                    nt_bsplot(z(:,iComp,:),[],[],t);
                    title(iComp);
                    xlim([t(1) t(end)]); xlabel('s');
                end
                
                % We store the DSSed data into a Fieldtrip variable.
                z_timelock.avg = z;
                z_timelock.time = t;
                z_timelock.label = short_data_subject.label;
                z_timelock.samples_cond1 = size(x_orig{1},3);
                z_timelock.samples_cond2 = size(x_orig{2},3);
                
                
                mkdir(fullfile('D:\Results','Trigger_analysis_PRE_HP2_LP30','Tone_Trials_DSS','Plots'));

                savefig(f1, fullfile('D:\Results','Trigger_analysis_PRE_HP2_LP30','Tone_Trials_DSS','Plots',...
                sprintf('Power-TRIG_%d_%d_Time_%d_%d-SUBJ_%d', trigger_list(1),trigger_list(2),round(T_init),round(T_end),...
                subject_list(subject_ind))));
            
                savefig(f2, fullfile('D:\Results','Trigger_analysis_PRE_HP2_LP30','Tone_Trials_DSS','Plots',...
                sprintf('Components-TRIG_%d_%d_Time_%d_%d-SUBJ_%d', trigger_list(1),trigger_list(2),round(T_init),round(T_end),...
                subject_list(subject_ind))));
            
                % We store the transformed data so we can use it later.
                mkdir(fullfile('D:\Results','Trigger_analysis_PRE_HP2_LP30','Tone_Trials_DSS','DSS_components'));
                save(fullfile('D:\Results','Trigger_analysis_PRE_HP2_LP30','Tone_Trials_DSS','DSS_components',...
                sprintf('toneDSS-TRIG_%d_%d_Time_%d_%d-SUBJ_%d-COMP_%d.mat',trigger_list(1), trigger_list(2),round(T_init),round(T_end),...
                subject_list(subject_ind), 274)),'z_timelock','todss');     
            
                 
                 clear z_timelock
                 clear todss
                
    end 
                 
                 clear x_orig
                 clear tone_subject
  end
 
 
end 