function pre_TestGraph(trigger_list, subject_list, config)
    % Function used to apply ft_timelockanalysis to the data of a subject
    % for a specific trigger. It is useful to reduce computational time,
    % since we precompute all the information for later plots. In order to
    % do that, it uses the pre-processed data. It is computed only over the
    % channels of interest in order to reduce the computational time.
    %
    % trigger: 
    % * 5 : 3 second RAND sequences.
    % * 10: 15 second RAND sequences.
    % * 15: 3 second REG sequences.
    % * 20: 15 second REG sequences.
    %
    % subject:
    % * 2-15, the index of the subject.
    %
    % config: allows for certain configurations.
    %   .out_folder: name of the folder where data will be stored.
    %   .store_data: set to 1 to store in a .mat file the whole data of the
    %       subject once that it has been processed.
    %   .channels_path: path indicating where to find the channel
    %       information file.
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
    % Last update: 18/June/2018
    

    out_folder = config.out_folder;
    store_output = config.store_data;
    channels_path = config.channels_path;
    
    compute =1;
    n_components = 274;
    %%   
    if compute == 1 
        for subject_ind = 1:length(subject_list)
            % We read the data from both triggers and use it to compute the
            % DSS transformation matrix.
            for trigger_ind = 1:length(trigger_list)
    %         for subject_ind = 1:length(subject_list)
    %             load(fullfile('..','Results',out_folder,'Preprocessed_data', sprintf('Long_timelock-TRIG_%d-SUBJ_%d', trigger_list(trigger_ind) ,subject_list(subject_ind))));

                load(fullfile(channels_path, sprintf('Channels-SUBJ_%d', subject_list(subject_ind))));
                load(fullfile('..','Results',out_folder,'Preprocessed_data_AllChannels', sprintf('Long_timelock-TRIG_%d-SUBJ_%d', trigger_list(trigger_ind) ,subject_list(subject_ind))));
                
%                 cfg = [];
%                 if trigger_list(trigger_ind) == 5 | trigger_list(trigger_ind) == 15
%                     cfg.toilim = [-0.2, 4];
%                 else
%                     cfg.toilim = [-0.2,16];
%                 end
%                 data_subject = ft_redefinetrial(cfg,data_subject );
                % We baseline using the average value of the pre-stimuli
                % information.
                cfg = [];
                cfg.demean = 'yes'; % Necessary to baseline.
                cfg.baselinewindow = [-0.2 0];% in seconds
                data_subject_BS = ft_preprocessing(cfg, data_subject);
                timelock274 = ft_timelockanalysis([], data_subject_BS);
                
                cfg.channel = channels;
                timelock40 = ft_timelockanalysis(cfg, data_subject_BS);
                
                
                
                out274(:,:, subject_ind, trigger_ind) = timelock274.avg;
                out40(:,:, subject_ind, trigger_ind) = timelock40.avg;
                fprintf(sprintf('SUBJ_%d-TRIG_%d', subject_ind, trigger_ind))
%                 x_orig{trigger_ind} = cat(3,data_subject_BS.trial{:});
            end
            
          

        end
% 
%         dss_comp = dss_comp*1e15;
        % Plot showing the average cumulative sum of power in DSS. We keep 3
        % components, which represents the 80% of the power.
       
    else
        load(fullfile('..','Results',out_folder,'DSS_components',sprintf('DSS_TRIG-%d-%d-COMP_%d.mat',trigger_list(1), trigger_list(2), n_components)))
    end
    %% Results
    perc = 0.05;
    perc2 = 0.01;
%     dss_comp = dss_comp;
%     load(fullfile(channels_path, sprintf('Channels-SUBJ_%d', subject_list(subject_ind))));
    for n_signals =1:2
        switch n_signals
            case 1

                components = n_components;

                for subject_ind = 1:length(subject_list)
                    load(fullfile(channels_path, sprintf('Channels-SUBJ_%d', subject_list(subject_ind))));
                    out_signal(:,subject_ind, :) = squeeze(rms(dss_comp(:,:,subject_ind,:),1)*1e15);
                end
                k1 = 5;
                k2 = 10;
            case 2
                % Original signal
                for subject_ind = 1:length(subject_list)
                    load(fullfile(channels_path, sprintf('Channels-SUBJ_%d', subject_list(subject_ind))));
                    out_signal(:,subject_ind, :) = squeeze(rms(orig_comp(channels_num,:,subject_ind,:),1)*1e15);
%                      out_signal2(:,subject_ind, :) = squeeze(rms(orig_comp(:,:,subject_ind,:),1)*1e15);
                end
                k1 = 5;
                k2 = 10;
        end

        % We use bootstrap and compute if there is a significant difference
        % between conditions.  
        Diff = out_signal(:,:,1) - out_signal(:,:,2);

        dataB=bootstrap(Diff'); 
        s=findSigDiff(dataB, perc);
        s2=findSigDiff(dataB, perc2);
        subplot(2,1,n_signals)
        plot(time_data, squeeze(mean(out_signal,2)), 'Linewidth',3);
        hold on;
        plot(time_data, k1*s, 'Linewidth',3)
        hold on
        plot(time_data, k2*s2,'Linewidth',3)
        legend('RAND','REG','p=0.05','p=0.01')    
        if trigger_list(1) == 5
            xlim([-0.5, 4])
        else
            xlim([-0.5, 16])
        end
            
        xlabel('Time (s)')
        ylabel('RMS magnitude (fT)')
        grid on
        switch n_signals
            case 1
                title(sprintf('DSS signal. Components: %d', components))
            case 2
                title('Original signal')
        end
        clear out_signal
    end
    
%%    

    if store_output == 1
        mkdir(fullfile('..','Results',out_folder,'DSS_components'))
        savefig(fullfile('..','Results',out_folder,'DSS_components',sprintf('2DSS-Comp%d-TRIG%d-%d',components, trigger_list(1), trigger_list(2))))
                
%         save(fullfile('..','Results',out_folder,'DSS_components',sprintf('2DSS-Comp%d-TRIG%d-%d',components, trigger_list(1), trigger_list(2))), 'out_signal')
    end
    

    %% Re-baselined data
    fs = 600;
    clear out_comp
    for subject_ind = 1:length(subject_list)
        % We read the data from both triggers and use it to compute the
        % DSS transformation matrix.
        for trigger_ind = 1:length(trigger_list)
%             baseline_data = zeros(components,1);
            baseline_data = mean(dss_comp(1:components, fs*(0.5+2.0):fs*(0.5+2.5), subject_ind, trigger_ind),2);
%             baseline_data(1) = 1;
            out_comp(:,:,subject_ind, trigger_ind) = dss_comp(1:components, :, subject_ind, trigger_ind) - repmat(baseline_data, 1, size(dss_comp,2));
%             min_constant(subject_ind, trigger_ind) = abs(min(reshape(out_comp(:,:,subject_ind, trigger_ind),prod(size(out_comp(:,:,subject_ind, trigger_ind))),1)));
%             out_comp(:,:,subject_ind, trigger_ind) = out_comp(:,:,subject_ind, trigger_ind) +min_constant(subject_ind, trigger_ind);
            clear baseline_data

        end
    end
    output_appendix = '_activity';
%        plot(time_indexes,squeeze(mean(rms(out_comp(:,:,:,:),1),3)), 'Linewidth',1)
 
%     plot((1:11100)/fs-0.5,squeeze(mean(orig_comp(:,:,:,1),3))')
       close all
    % Bootstrapping    
    Diff =  squeeze(rms(out_comp(:,:,:,1),1)-rms(out_comp(:,:,:,2),1));
    dataB=bootstrap(Diff'); 
    perc = .05;
    s=findSigDiff(dataB, perc);
    perc2 = .01;
    s2=findSigDiff(dataB, perc2);
    % Plotting
    time_indexes = (1:size(out_comp,2))/fs-0.5;
    plot(time_indexes,squeeze(mean(rms(out_comp(:,:,:,:),1),3)), 'Linewidth',1)
    hold on
    plot(time_indexes, 2*s, 'Linewidth',3);
    hold on
    plot(time_indexes, 4*s2, 'Linewidth',3);
    xlabel('Time (s)')
    ylabel('RMS magnitude (fT)')
    xlim([2, 15])
    legend('RAND','REG','p=0.05','p=0.01')    
    if store_output == 1
        mkdir(fullfile('..','Results',out_folder,'DSS_components'))
        savefig(fullfile('..','Results',out_folder,'DSS_components',sprintf('22DSS-BS%s-Comp%d-TRIG%d-%d', output_appendix,components, trigger_list(1), trigger_list(2))))
                
    end   
        
end
