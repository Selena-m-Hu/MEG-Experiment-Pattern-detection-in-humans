function post_TriggerAnalysis(trigger_list, subject, config)
    % Function designed to plot the timelock information for all the
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
    %
    % subject: we can choose a single subject or a list of them.
    % * 2-15, the index of the subject.
    %
    % config: allows for certain configurations.
    %   .out_folder: name of the folder where data will be stored.
    %   .store_data: set to 1 to store in a .mat file the whole data of the
    %       subject once that it has been processed.
    %
    % Visitor: 
    % Antonio Rodriguez Hidalgo 
    % Dept. of Signal Theory and Communications
    % Universidad Carlos III de Madrid
    % arodh91@gmail.com
    %
    % Principal Investigator:
    % Maria Chait 
    % Ear Institute
    % University College London
    % m.chait@ucl.ac.uk
    %
    % Last update: 28/June/2018
    

    out_folder = config.out_folder;
    store_output = config.store_data;
    
    fs          = 600;
    pre_estim   = 0.5;
    window_size = .250*fs; % 250 ms.
                
    trigger_length = 0.25; % 250ms.
    event_length = .05;
    spacing = window_size*.05/trigger_length;
    T_init = config.temporal_frame(1);
    T_end = config.temporal_frame(2);
    hpfreq = config.hpfreq;
%     verbose = 1;
    
    %% We create some matrices with the timelock information.
    for trigger_ind = 1:length(trigger_list)
        for subject_ind = 1:length(subject)
            % We load the timelock information that we computed in a
            % previous stage.
            load(fullfile('..','Results',out_folder,'Timelock',sprintf('Timelock-TRIG_%d-SUBJ_%d.mat',trigger_list(trigger_ind),subject(subject_ind))),'timelock')
            counter = 1;

            for t_ind = 1:window_size:length(timelock.time)
                % We fill the first one (which will not be used) with
                % zeros, since it will be incomplete. This way, the
                % structure of the data is easier and we consider a 50ms
                % to baseline and 200ms of signal.
                if t_ind == 1
                    trigger_shape(:,:,counter) = zeros(1,150);
                else
                    classical_baseline = 1;
                    if classical_baseline == 1
                        trigger_shape(:,:,counter) = rms(timelock.avg(:,t_ind-.05*fs:t_ind+0.2*fs-1));
                        % We baseline the signal.
                        trigger_shape(:,:,counter) = trigger_shape(:,:,counter) - mean( trigger_shape(:,1:window_size*.05/trigger_length,counter));

%                     else 
                        % NEW BASELINE (nope nope)
                        trigger_shape2(:,:,counter) = rms(timelock.avg(:,t_ind-.05*fs:t_ind+0.2*fs-1) - repmat(mean(timelock.avg(:,t_ind-.05*fs:(t_ind-1)),2),1,150));
                    end
                end
                counter = counter + 1;
            end
            trigger_shape = trigger_shape*1e15;
            trigger_shape2 = trigger_shape2*1e15;


            %% Baselined data.
            % This is useful to get the AVERAGE shape of the signal, and
            % event to extract some features from this feature.
            t0 = pre_estim/(window_size/fs);
            
            % Stable data (6-14 seconds after stimuli)
            t_stable = t0 + ((T_init/(window_size/fs)):(T_end/(window_size/fs)));
            % Unstable data (1-4 seconds after stimuli)
%             t_unstable = t0 + ((1/(window_size/fs)):(4/(window_size/fs)));
            
            % We don't need to baseline again since we did it in a previous
            % stage.
            stable_shape(:,:,subject_ind, trigger_ind)      = squeeze(trigger_shape(:,:,t_stable))';
            stable_average_shape(:,:,subject_ind, trigger_ind)      = mean(trigger_shape(:,:,t_stable),3);
            stable_std_shape(:,:,subject_ind, trigger_ind)      = std(bootstrp(1000, @mean, squeeze(trigger_shape(:,:,t_stable))'));
%             unstable_average_shape(:,:,subject_ind, trigger_ind)    = mean(trigger_shape(:,:,t_unstable),3);
               
            %% SLOPES
            % This feature should be useful to get how steep is the slope
            % of each one of the modalities for all the peaks. The way to
            % compute it is using the average magnitude per modality.
%             slopes(:,subject_ind, trigger_ind) = mean(permute(diff(rms(trigger_shape(:,:,t_stable),1),1,2),[2,3,1]),2);
            slopes_stable(:,subject_ind, trigger_ind) = diff(stable_average_shape(:,:,subject_ind, trigger_ind),1)';
%             slopes_unstable(:,subject_ind, trigger_ind) = diff(unstable_average_shape(:,:,subject_ind, trigger_ind),1)';
            
            %% Gaussian Mixture results (NON CONSIDERED NOW)
            % In this case, we generate several mixtures of two Gaussians
            % with different deviations. Our goal is to correlate these
            % kernels with the trigger data in order to determine whether
            % their shape is similar to the mixture ones. Both Gaussians
            % are spaced 50ms, according to the nature of our experiment.
            % We also compute the kurtosis of each one of the triggers.
            
%             spacing = window_size*[.05:.01:.12]./trigger_length; % Spacing from 30ms to 90ms.
%             prestim_delay = window_size*event_length/trigger_length; % 50ms (30 samples).
%             gauss_std = [4]; % x/150*.25
%             init_delay = prestim_delay; % We suppose the first gaussian appears 50ms post-stimuli. This parameter will affect our estimation of the delay.
%             x = 1:window_size;
%             ind_gauss = 1;
%             for ind_spacing = 1:length(spacing)
%                 for ind_std = 1:length(gauss_std)
%                     gauss1 = gaussmf(x,[gauss_std(ind_std) prestim_delay+init_delay]);
%                     gauss2 = gaussmf(x,[gauss_std(ind_std) prestim_delay+init_delay+spacing(ind_spacing)]);
%                     kernel(:,ind_gauss) = mat2gray((gauss1+gauss2));
%                     kernel_config(ind_gauss).spacing = spacing(ind_spacing);
%                     kernel_config(ind_gauss).std = gauss_std(ind_std);
%                     
%                     ind_gauss = ind_gauss + 1;
%                     
%                 end
%             end
%             
%             % We check the alignment of each one of the triggers of the
%             % subject.
%             clear val
%             for ind_kernel = 1:length(gauss_std)*length(spacing)
%                 for ind_trigg = 1:size(trigger_shape,3)
%                     aux = trigger_shape(:,:,ind_trigg);
%                     % Since we only care about the shape of the signal and
%                     % not the magnitude, we normalize the data.
%                     aux_zeros = aux; 
%                     [max_val(ind_kernel,ind_trigg),val(ind_kernel,ind_trigg)] = max(xcorr(mat2gray(aux_zeros), kernel(:,ind_kernel)));
%                     kurtosis_val(ind_trigg, subject_ind, trigger_ind) = kurtosis(aux);
%                 end
%             end
%             delay_vec (:,:, subject_ind, trigger_ind) = val;
%             max_vec (:,:, subject_ind, trigger_ind) = max_val;
          
            clear trigger_shape
  
        end    
    end
    aux = squeeze(stable_average_shape);
    
    mkdir(fullfile('..','Results','Global_trigger'))
    save(fullfile('..','Results','Global_trigger',sprintf('GLOBAL_T%.2f-T%.2f_HP%d.mat', T_init, T_end, hpfreq)), 'aux', 'stable_shape');
      
    
    keyboard
    %% Gaussian Mixture results (NON CONSIDERED NOW)
%     ind_ker = 1;
%     % Cross-correlation DELAY
%     post_SubjectTemporalStats(delay_vec(ind_ker,t_stable,:,:), subject)
%     subplot(161)
%     ylabel('Average delay with GM')
%     subplot(1,6,2:6)
%     xlabel('Subjects')
%     title('Delay')
%     
%     % MAXIMUM cross-correlation values.
%     post_SubjectTemporalStats(max_vec(ind_ker,t_stable,:,:), subject)
%     subplot(161)
%     ylabel('Average Cross-correlation')
%     subplot(1,6,2:6)
%     xlabel('Subjects')
%     title('Max Cross-correlation with GM')
    

    
    

    %% Analysis of signal peaks considering 2-GM with different spacing.
    % Here we intend to measure the spacing of each one of the subjects
    % responses using a Gaussian Mixture with two elements. By doing so, we
    % intend to verify if a subject shows a different spacing depending on
    % the modality.
    pre_stim_time = 30; % 50ms - 30 samples.
    for subject_ind = 1:length(subject) 
%         plot(ks_stable(:,subject_ind,1)); hold on;
%         plot(ks_stable(:,subject_ind,2))
        for trigger_ind = 1:length(trigger_list)
            % We start processing the data from the stable temporal area.
            % We get both the magnitudes (normalized to [0,1]) and their
            % slopes. We are removing the initial 50ms window, since we
            % used them to baseline the data previously. Consequently, we
            % have a temporal length of 200ms (120 samples).
            magnitudes = mat2gray(squeeze(stable_average_shape(:,pre_stim_time+1:end,subject_ind,trigger_ind)));
            slope = slopes_stable(pre_stim_time+1:end,subject_ind,trigger_ind);
            max_indexes = find(diff(slope > 0) == -1)+1;
            min_indexes = [1; find(diff(slope > 0) == 1)+1; length(slope)];
            
            % We get the position of each one of the maxima from each
            % averaged trigger signal (for subject and condition).
            max_activity = zeros(size(magnitudes));
            for ind_min = 1:length(min_indexes)-1
                % We make use chunks of the normalized signal (from 0 to 1,
                % chunk) and the one with the original magnitude values 
                % (chunk_orig). We get chunks going from minimum to
                % minimum, including the initial and final indexes.
                chunk = magnitudes(min_indexes(ind_min):min_indexes(ind_min+1));
                chunk_orig = squeeze(stable_average_shape(:,pre_stim_time+1:end,subject_ind,trigger_ind));
                chunk_orig = chunk_orig(min_indexes(ind_min):min_indexes(ind_min+1));
                
                % We determine the maxima region considering the following
                % criteria: a maximum value exists from the extreme point
                % and also its surroundings, considering that they belong
                % to the maximum point when their magnitude is (as a limit
                % value) 90 of the maximum's one .
                % We process and normalize each chunk using the maximum
                % value contained, and determine the 90% region area.
                for ind_max = 1:length(max_indexes)
                    norm_chunk = chunk/magnitudes(max_indexes(ind_max));
                    if sum(norm_chunk > 1) == 0 & isempty(find(norm_chunk == 1)) == 0
                        % We extract a vector that will contain the
                        % position of each one of the maxima areas along
                        % the temporal values (stored in a binary vector).
                        % This will be useful to extract the onset and
                        % offset times of each maximum.
                        max_activity(find(norm_chunk > .9)+min_indexes(ind_min)-1) = 1;
                        
                        % We also store the max values (non-normalized in 
                        % order to keep their sign). This gives us a list
                        % with the average magnitude of the 90% area.
%                         mean_mag(ind_max) = mean(chunk(norm_chunk > .9));
                        mean_mag(ind_max) = mean(chunk_orig(norm_chunk > .9));

                    end
                end
            end
            
            % Once we have the maximum data and their onset-offset times,
            % we store them in a structure.
            try
                % Original latency of the detected maxima.
                peak(subject_ind, trigger_ind).detected_centers = max_indexes' + pre_stim_time;
                
                % Onset latency considering the 90% area.
                peak(subject_ind, trigger_ind).t_onset = find(diff(max_activity) == 1) + pre_stim_time;
                % We consider that an initial peak might happen with
                % latency 0ms. For that case:
                if max_activity(1) == 1
                    peak(subject_ind, trigger_ind).t_onset  = [1, peak(subject_ind, trigger_ind).t_onset];
                end
                
                % Offset latency considering the 90% area.
                peak(subject_ind, trigger_ind).t_offset = find(diff(max_activity) == -1) + pre_stim_time;
                
                % Magnitude considering the 90% area.
                peak(subject_ind, trigger_ind).ave_magnitude = mean_mag;
                
                % Normalized magnitude for the 90% area, using the latency
                % onset-offset difference to normalize.
                peak(subject_ind, trigger_ind).norm_magnitude = mean_mag./(peak(subject_ind, trigger_ind).t_offset - peak(subject_ind, trigger_ind).t_onset);
                
                % We also store the original signal used to get the maxima
                % and their estimated (binary) positions.
                peak(subject_ind, trigger_ind).signal = squeeze(stable_average_shape(:,pre_stim_time+1:end,subject_ind,trigger_ind));
                peak(subject_ind, trigger_ind).peaks_bin = max_activity;
                
                % We get an estimation of the central latency value of the
                % peaks after the 90% processing.
                peak(subject_ind, trigger_ind).t_central = (peak(subject_ind, trigger_ind).t_offset + peak(subject_ind, trigger_ind).t_onset)/2;
                
            catch
                fprintf('Something went wrong while computing the data from the peaks!\n');
                keyboard
            end
            clear mean_mag

        end
    end
    
    
    %% Peak plotting
    if store_output == 1
        for subject_ind = [1:length(subject)]
            for trigger_ind = 1:length(trigger_list)
                signal_aux = peak(subject_ind,trigger_ind).signal;
                plot(   (1:length(signal_aux))*50/30, signal_aux, 'Linewidth',3); hold on;
                % Bootstrap std value from the mean.
                std_val = (squeeze(stable_std_shape(:,pre_stim_time+1:end,subject_ind,trigger_ind)));
                fill([1:length(std_val), fliplr(1:length(std_val))]*50/30,...
                    [signal_aux + std_val, fliplr(signal_aux - std_val)],'b', 'FaceAlpha',.2, 'LineStyle','none' )
                
                
                bin_aux = peak(subject_ind,trigger_ind).peaks_bin;
                fill([1:length(bin_aux), fliplr(1:length(bin_aux))]*50/30,...
                        [bin_aux*max(signal_aux), zeros(size(bin_aux))],'r','FaceAlpha',.3)
%                 plot((1:length(peak(subject_ind,trigger_ind).signal))*50/30,...
%                         peak(subject_ind,trigger_ind).peaks_bin)
                legend('RMS magnitude', '', 'Selected peaks')
                xlabel('Time (ms)')
                ylabel('Baselined RMS magnitude (fT)')
                title(sprintf('Subj: %d. Trig: %d\n', subject(subject_ind), trigger_list(trigger_ind)))
                mkdir(fullfile('..','Results',out_folder,'Peak_selection'))
                savefig(fullfile('..','Results',out_folder,'Peak_selection',sprintf('Peaks_T%d-T%d_SUBJ-%d_TRIG-%d.fig', T_init, T_end, subject(subject_ind), trigger_list(trigger_ind))))
                
                clf
            end
        end
    end
    
    
    %% Peak analysis
    % Once that we have the peak descriptors, we analyze them.
    for subject_ind = [1:length(subject)]
        for trigger_ind = 1:length(trigger_list)
            
            % We sort the peaks using the normalized magnitude, since it is
            % able to detect properly the magnitudes of maxima peaks
            % against the magnitude of other artifact peaks. Probably other
            % peak components should work too, although with different
            % peak_idx values.
            [peak_val, peak_idx] = sort(peak(subject_ind,trigger_ind).norm_magnitude,'descend');

            % We usually keep two peaks for each one of the subjects and
            % conditions. However, there are certain exceptions for case 8
            % (subject #9) and case 10 (subject #11), whose maxima indexes
            % were not properly sorted in the previous line of code. The
            % procedure to select those peaks was performed manually
            % (visually).
            
            
            if T_init == 6 & T_end == 14
                switch subject_ind
                    case 8
                        if trigger_ind == 2
                            aux_ind = [1,4];
                        else
                            aux_ind = 1:2;
                        end

                    case 10
                        if trigger_ind == 1
                            aux_ind = [1,3];
                        else
                            aux_ind = 1;                        
                        end

                    otherwise
                        aux_ind = [1:2];

                end
            elseif T_init == 1 & T_end == 4
                
                switch subject_ind
                    case 1
                        aux_ind = [1,3];
                    case 5  
                        aux_ind = [1,3];
                    case 10
                        if trigger_ind == 1
                            aux_ind = [1,3];
                        else
                            aux_ind = 1:2;                        
                        end
                    otherwise
                        aux_ind = [1:2];

                end
            end
            % Once we have the indexes of the peaks for each subject and
            % condition, we extract their information.
            % First, we average the absolute NORMALIZED magnitude of the 
            % selected peaks for a subject and condition.
            out_peaks_NormalizedMagnitude(subject_ind, trigger_ind) = mean(abs(peak_val(aux_ind)));
            
            % We repeat the procedure for the NON-NORMALIZED magnitudes.
            out_peaks_AverageMagnitude(subject_ind, trigger_ind) = mean(abs(peak(subject_ind,trigger_ind).ave_magnitude(peak_idx(aux_ind))));
            
            % We store the number of peaks selected for each subject and
            % condition. We should expect to have at least two peaks for
            % every combo. If a subject does not comply, we remove it from
            % the analysis since the average value might be unfair.
            num_peaks(subject_ind, trigger_ind) = length(peak_val(aux_ind));
            
            % We select the central latency obtained using the 90% area
            % criteria, converted to milliseconds (50ms - 30 samples).
            t_central{subject_ind, trigger_ind} = sort(peak(subject_ind, trigger_ind).t_central(peak_idx(aux_ind)))*50/30;
            
            % We determine the delay existing between the two peaks that
            % are analyzed.
            if num_peaks(subject_ind, trigger_ind) == 2
                % Using the 90% criteria.
                peak_delay90(subject_ind, trigger_ind) = diff(t_central{subject_ind, trigger_ind});
                % Using the exact position of the maximum value.
                peak_delayMAX(subject_ind, trigger_ind) = diff(peak(subject_ind, trigger_ind).detected_centers(peak_idx(aux_ind)));
            end
            clear aux_ind
                
        end
    end

    % Some initial statistics, removing Subject #11 (index 11), since it
    % only had one peak for condition 2 (REG).
    selected_subjects = [1:9,11:13];
    [~, table_delay90] = anova_rm(peak_delay90(selected_subjects,:), 'off');
    [~, table_delayMAX] = anova_rm(abs(peak_delayMAX(selected_subjects,:)), 'off');
    [~, table_NormalizedMagnitude] = anova_rm(out_peaks_NormalizedMagnitude(selected_subjects,:), 'off');
    [~, table_AverageMagnitude] = anova_rm(out_peaks_AverageMagnitude(selected_subjects,:), 'off');
    
    % We store the results.
    if store_output == 1
        results.peak_delay90 = peak_delay90;
        results.peak_delayMAX = abs(peak_delayMAX);
        results.out_peaks_NormalizedMagnitude = out_peaks_NormalizedMagnitude;
        results.out_peaks_AverageMagnitude = out_peaks_AverageMagnitude;

        tables.table_AverageMagnitude = table_AverageMagnitude;
        tables.table_NormalizedMagnitude = table_NormalizedMagnitude;
        tables.table_delayMAX = table_delayMAX;
        tables.table_delay90 = table_delay90;
        mkdir(fullfile('..','Results',out_folder,'Peak_analysis'))
        save(fullfile('..','Results',out_folder,'Peak_analysis',sprintf('Peak-Results_T%d-T%d.mat', T_init, T_end)), 'tables', 'results', 'peak')

    end
    
    
    %% We compute the AUC for the ||.|| of data and the peak curves for each subject an condition.
    if store_output == 1
        for subject_ind = 1:length(subject)
         
            % Modality 10 (RAND)
            signal_aux = stable_average_shape(1,:,subject_ind,1);
            plot((1:150)*50/30-50,signal_aux , 'Linewidth',3); hold on; 
            % Bootstrap std value from the mean.
            std_val = (squeeze(stable_std_shape(1,:,subject_ind,1)));
            fill([1:length(std_val), fliplr(1:length(std_val))]*50/30-50,...
                [signal_aux + std_val, fliplr(signal_aux - std_val)],'b', 'FaceAlpha',.2, 'LineStyle','none' )
            hold on
            
            % Modality 20 (REG)
            signal_aux = stable_average_shape(1,:,subject_ind,2);
            plot((1:150)*50/30-50, signal_aux, 'Linewidth',3);
            % Bootstrap std value from the mean.
            std_val = (squeeze(stable_std_shape(1,:,subject_ind,2)));
            fill([1:length(std_val), fliplr(1:length(std_val))]*50/30-50,...
                [signal_aux + std_val, fliplr(signal_aux - std_val)],'r', 'FaceAlpha',.2, 'LineStyle','none' )
            hold on
            
       
            
            
            
            title(sprintf('Event average. Subject: %d', subject(subject_ind)))
            auc(subject_ind, 1) = trapz(abs(stable_average_shape(1,:,subject_ind,1)));
            auc(subject_ind, 2) = trapz(abs(stable_average_shape(1,:,subject_ind,2)));
            legend(sprintf('RAND - AUC:%.2f',auc(subject_ind, 1)),'',...
                   sprintf('REG - AUC:%.2f',auc(subject_ind, 2)),'');
            xlabel('Time (ms)');
            ylabel('Baselined RMS magnitude (fT)')
            mkdir(fullfile('..','Results',out_folder,'Peak_analysis'))
            savefig(fullfile('..','Results',out_folder,'Peak_analysis',sprintf('Fig-AUC_T%d-T%d_SUBJ-%d.fig', T_init, T_end, subject(subject_ind))))
      
            clf
        end
        save(fullfile('..','Results',out_folder,'Peak_analysis',sprintf('AUC_T%d-T%d_SUBJ-%d.mat', T_init, T_end, subject(subject_ind))), 'auc')
    end
    

    %% We extract chunks from the signal.
    if store_output == 1
%         figure('units','normalized','outerposition',[0 0 1 1])   
        % RAND modality
        signal_aux  = mean(stable_average_shape(:,:,:,1),3);
        plot((1:window_size)/fs-.05, signal_aux, 'Linewidth',3); hold on;
        % Bootstrap std value from the mean.
        std_val = mean(stable_std_shape(:,:,:,1),3);        
        fill([1:length(std_val), fliplr(1:length(std_val))]/fs-0.050,...
            [signal_aux + std_val, fliplr(signal_aux - std_val)],'b', 'FaceAlpha',.2, 'LineStyle','none' );
        hold on
        
        % REG modality
        signal_aux = mean(stable_average_shape(:,:,:,2),3);
        plot((1:window_size)/fs-.05, signal_aux, 'Linewidth',3);
        % Bootstrap std value from the mean.
        std_val = mean(stable_std_shape(:,:,:,2),3);        
        fill([1:length(std_val), fliplr(1:length(std_val))]/fs-0.050,...
            [signal_aux + std_val, fliplr(signal_aux - std_val)],'r', 'FaceAlpha',.2, 'LineStyle','none' );
        hold on
        
        
        xlabel('Time (s)')
        ylabel('Baselined RMS magnitude (fT)')
        legend('RAND','','REG','')
        title(sprintf('GLOBAL data from %d to %d seconds', T_init, T_end))

        mkdir(fullfile('..','Results',out_folder,'Peak_analysis'))
        savefig(fullfile('..','Results',out_folder,'Peak_analysis',sprintf('GLOBAL_T%d-T%d.fig', T_init, T_end)))
                
    end

  


    
    %% KURTOSIS
%     % In this section we analyze and plot the kurtosis results for each one
%     % of the available triggers.
%     % Kurtosis values
%     post_SubjectTemporalStats(reshape(kurtosis_val(t_stable,:,:),1, 33,13,2), subject)
%     subplot(161)
%     ylabel('Average Kurtosis')
%     subplot(1,6,2:6)
%     xlabel('Subjects')
%     title('Kurtosis')
    
  
%     


    
end
