function post_TriggerAnalysisCHANNEL(trigger_list, subject, config)
    % Function useful to study the response observed after each one of the
    % tones.
    % It reads the timelock data obtained for each subject and
    % trigger and divides it into tones of 250ms. Then, it baselines them
    % according to three different criteria, and store the results in a
    % matrix.
    %
    % trigger_list: 
    % * [5, 15]     : SHORT sequences.
    % * [10, 20]    : long sequences.
    %
    % subject:
    % * 2-15, the index of the subject.
    %
    % config: allows for certain configurations.
    %   .out_folder: name of the folder where data will be stored.
    %   .store_data: set to 1 to store in a .mat file the whole data of the
    %       subject once that it has been processed.
    %   .temporal_frame: temporal region where the tones of interest are
    %       located. It gets divided into T_init and T_end, which indicate
    %       the beginning and the end of the data of interest. Is a vector
    %       with [T_init, T_end] in seconds.
    %   .hpfreq: in this case, useful to store the results in the proper
    %       folder.
    %   .lpfreq: in this case, useful to store the results in the proper 
    %       folder.
    %   .baseline: baseline modality that we want to use. We have three
    %       different ones:
    %       * 'tone': baselined independently for each tone using its
    %           initial 50ms (from the total of 250ms per tone).
    %       * 'silence': baseline data is computed using the prestimuli 
    %           time. In our case, 500ms.
    %       * 'activity': baseline data is acquired from 2 to 2.5 seconds
    %           poststimuli, just 500ms before the second repetition of the
    %           sequence occurs (in case of REG).
    %
    % OUTPUT FOLDER (one of the following):
    %   * ../Results/***/Global_trigger,
    %   * ../Results/***/Global_trigger_silence
    %   * ../Results/***/Global_trigger_activity
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
    % Last update: 16/July/2018
    

    out_folder      = config.out_folder;
    store_output    = config.store_data;
    T_init          = config.temporal_frame(1);
    T_end           = config.temporal_frame(2);
    hpfreq          = config.hpfreq;
    lpfreq          = config.lpfreq;
    baseline        = config.baseline;
    fs              = 600; % Sampling frequency
    pre_estim       = 0.5; % Prestimuli time. In our case, 500ms.
    window_size     = .250*fs; % 250 ms, the length of the stimuli (50ms signal + 200ms silence).
    channels_path   = config.channels_path;            
%     trigger_length = 0.25; % 250ms.
   

    
    %% We create some matrices with the timelock information.
    for trigger_ind = 1:length(trigger_list)
        for subject_ind = 1:length(subject)
            % We load the timelock information that we computed in a
            % previous stage.
%             load(fullfile('..','Results',out_folder,'Timelock',sprintf('Timelock-TRIG_%d-SUBJ_%d.mat',trigger_list(trigger_ind),subject(subject_ind))),'timelock')
%             channels_num = 1:40;
            load(fullfile('..','Results',out_folder,'Preprocessed_data_AllChannels',sprintf('Long_timelock-TRIG_%d-SUBJ_%d.mat',trigger_list(trigger_ind),subject(subject_ind))),'timelock')
            % We load the channels that we obtained during the PSD analysis.
            load(fullfile(channels_path, sprintf('Channels-SUBJ_%d', subject(subject_ind))));
%             
            counter = 1;

            for t_ind = 1:window_size:length(timelock.time)
                % We fill the first one (which will not be used) with
                % zeros, since it will be incomplete. This way, the
                % structure of the data is easier and we consider a 50ms
                % to baseline and 200ms of signal.
                if t_ind == 1
                    trigger_shape(:,:,counter) = zeros(40,150);
                else
                    trigger_shape(:,:,counter) = timelock.avg(channels_num,t_ind-.05*fs:t_ind+0.2*fs-1);


                end
                counter = counter + 1;
            end
            
            % Multiply by a constant in order to get the units into
            % femtoTeslas.
            trigger_shape = trigger_shape*1e15;

            % We get the value of the initial time (0 seconds).
            t0 = pre_estim/(window_size/fs);
            
            % We get the temporal index of the beginning and the end of our
            % data from the global matrix.
            t_stable = t0 + ((T_init/(window_size/fs)):(T_end/(window_size/fs)));
            stable_average_shape(:,:,subject_ind, trigger_ind)      = mean(trigger_shape(:,:,t_stable),3);
            
            
            
            % We baseline the data according to the criteria we want to
            % use.
            % 'tone'    :  each tone is baselined using its initial 50ms.
            % 'silence' :  each tone is baselined using the 500ms 
            %              prestimuli data.  
            % 'activity':  each tone is baselined using the data from 2 to
            %              2.5 seconds poststimuli, just before the tone
            %              sequence starts repeating.
            switch baseline
                case 'tone'
                    baseline_data = mean(stable_average_shape(:,1:30,subject_ind, trigger_ind) ,2);
                    output_appendix = '';
                case 'silence'
                    baseline_data = mean(timelock.avg(channels_num, 1:fs*.5),2)*1e15;
                    output_appendix = '_silence';
                case 'activity'
                    baseline_data = mean(timelock.avg(channels_num, fs*(.5+2.0):fs*(.5+2.5)),2)*1e15;
                    output_appendix = '_activity';
            end
                    
            stable_shape(:,:,:,subject_ind, trigger_ind) = trigger_shape - repmat(mean(stable_average_shape(:,1:30,subject_ind, trigger_ind) ,2), 1, 150, size(trigger_shape,3));
            stable_average_shape(:,:,subject_ind, trigger_ind)  =  stable_average_shape(:,:,subject_ind, trigger_ind) - repmat(baseline_data, 1, 150);
            
            % NON-CHUNKED PATH (security check)
            timelock_global(:,:,subject_ind, trigger_ind) = timelock.avg(channels_num,:)*1e15 -repmat(baseline_data,1,size(timelock.avg(channels_num,:),2));
            
            counter = 1;
            for t_ind = 1:window_size:length(timelock.time)
                % We fill the first one (which will not be used) with
                % zeros, since it will be incomplete. This way, the
                % structure of the data is easier and we consider a 50ms
                % to baseline and 200ms of signal.
                if t_ind == 1
                    trigger_global(:,:,counter) = zeros(40,150);
                else
                    trigger_global(:,:,counter) = timelock_global(:,t_ind-.05*fs:t_ind+0.2*fs-1, subject_ind, trigger_ind);


                end
                counter = counter + 1;
            end
            verbose = 0;
            if verbose == 1
                plot(rms(mean(trigger_global(:,:,t_stable),3)),'Linewidth',4); hold on
                plot(rms(stable_average_shape), 'Linewidth',2)
                legend('Chunked then baselined', 'Baselined then chunked')
%                 plot(rms(mean(trigger_global(:,:,t_stable),3))-rms(stable_average_shape))

            end
            clear trigger_shape
  
        end    
    end
    
    % This is the final data that we will store.
    aux = squeeze(rms(stable_average_shape,1));
    tone_aux = squeeze(stable_shape);
    out = mean(rms(timelock_global(:,:,:,1),1),3);
    if verbose == 1
        Diff =  squeeze(rms(timelock_global(:,:,:,1),1) - rms(timelock_global(:,:,:,2),1) );

        dataB=bootstrap(Diff'); 

        perc = .05;
        s=findSigDiff(dataB, perc);

        perc2 = .001;
        s2=findSigDiff(dataB, perc2);


        %% Removal of noise from the differences.
        % ind_search = 1:3*600;
        ind_search = (2+0.5)*600:(2.5+0.5)*600;
        s_aux = abs(s);
        offset1 = find(diff(s_aux(ind_search) == 1) == -1);
        onset1 = find(diff(s_aux(ind_search) == 1) == 1);
        diff1 = (offset1 - onset1);

        offset2 = find(diff(s_aux == 1) == -1);
        onset2 = find(diff(s_aux == 1) == 1);
        diff2 = offset2 - onset2;

        output = ones(size(s))*NaN;
        pos = find(diff2 > max(diff1));
        for ind = 1:length(pos)
            output(onset2(pos(ind)):offset2(pos(ind))) = 1;
        end

        s2_aux = abs(s2);
        offset1 = find(diff(s2_aux(ind_search) == 1) == -1);
        onset1 = find(diff(s2_aux(ind_search) == 1) == 1);
        diff1 = (offset1 - onset1);

        offset2 = find(diff(s2_aux == 1) == -1);
        onset2 = find(diff(s2_aux == 1) == 1);
        diff2 = offset2 - onset2;

        output2 = ones(size(s2_aux))*NaN;
        if isempty(diff1) == 1
            pos = 1:length(onset2);
        else
            pos = find(diff2 > max(diff1));
        end
        for ind = 1:length(pos)
            output2(onset2(pos(ind)):offset2(pos(ind))) = 1;
        end

        output = abs(output).*sign(s);
        output2 = abs(output2).*sign(s2);
        
        output_pos = output;
        output_pos(output_pos <0) = NaN;        
        output_neg = output;
        output_neg(output_neg >0) = NaN;
        
        output_pos2 = output2;
        output_pos2(output_pos2 <0) = NaN;        
        output_neg2 = output2;
        output_neg2(output_neg2 >0) = NaN;
        %%
        
        t_index = (2.+0.5)*600+1:(15+0.5)*600;
        plot(timelock.time(t_index), squeeze(mean(rms(timelock_global(:,(t_index),:,:),1),3)))
        hold on;
        plot(timelock.time(t_index), 1*(s(t_index)), 'Linewidth', 3);
        hold on;
        plot(timelock.time(t_index), 4*(s2(t_index)), 'Linewidth', 3);
        
        hold on;
        plot(timelock.time(t_index), 2*output(t_index), 'Linewidth', 3);
        hold on;
        plot(timelock.time(t_index), 5*output2(t_index), 'Linewidth', 3);
        
        
        xlabel('Time (s)')
        ylabel('RMS magnitude (fT)')
        legend('RAND', 'REG', 'p=0.05','p=0.001')
        grid on
        xlim([2, 15])
        
        %%
        figure
        subplot(211)
        s_pos = s;
        s_pos(s_pos <0) = NaN;        
        s_neg = s;
        s_neg(s_neg >0) = NaN;
        
        s_pos2 = s2;
        s_pos2(s_pos2 <0) = NaN;        
        s_neg2 = s2;
        s_neg2(s_neg2 >0) = NaN;
        
        t_index = (2.+0.5)*600+1:(15+0.5)*600;
        plot(timelock.time(t_index), squeeze(mean(rms(timelock_global(:,(t_index),:,:),1),3)))
        hold on;
        plot(timelock.time(t_index), 1*(s_pos(t_index)), 'Linewidth', 3);
        hold on;
        plot(timelock.time(t_index), 4*(s_pos2(t_index)), 'Linewidth', 3);
        
        hold on;
        plot(timelock.time(t_index), 2*output_pos(t_index), 'Linewidth', 3);
        hold on;
        plot(timelock.time(t_index), 5*output_pos2(t_index), 'Linewidth', 3);
        
        
        xlabel('Time (s)')
        ylabel('RMS magnitude (fT)')
        legend('RAND', 'REG', 'p=0.05','p=0.001')
        title('Positive. RAND > REG')
        grid on
        xlim([2, 15])
        
        subplot(212)   
        t_index = (2.+0.5)*600+1:(15+0.5)*600;
        plot(timelock.time(t_index), squeeze(mean(rms(timelock_global(:,(t_index),:,:),1),3)))
        hold on;
        plot(timelock.time(t_index), 1*(s_neg(t_index)), 'Linewidth', 3);
        hold on;
        plot(timelock.time(t_index), 4*(s_neg2(t_index)), 'Linewidth', 3);
        
        hold on;
        plot(timelock.time(t_index), 2*output_neg(t_index), 'Linewidth', 3);
        hold on;
        plot(timelock.time(t_index), 5*output_neg2(t_index), 'Linewidth', 3);
        
        
        xlabel('Time (s)')
        ylabel('RMS magnitude (fT)')
        legend('RAND', 'REG', 'p=0.05','p=0.001')
        grid on
        xlim([2, 15])
        title('Negative. RAND < REG')
        
        savefig(fullfile('..','Results','ToneGraphs',sprintf('GLOBAL%s-HP_%d-LP_%d.fig', output_appendix, hpfreq, lpfreq)));
        
        
    end
    
    %%
    if store_output == 1
        mkdir(fullfile('..','Results',sprintf('Global_trigger%s', output_appendix)))
        save(fullfile('..','Results',sprintf('Global_trigger%s', output_appendix),sprintf('GLOBAL_T%.2f-T%.2f_HP%d_LP%d.mat', T_init, T_end, hpfreq, lpfreq)), 'aux', 'tone_aux');
      
    end


    
end
