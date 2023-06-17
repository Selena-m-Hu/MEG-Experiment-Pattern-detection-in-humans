function Main_statsDSS()
    % This function plots and computes the statistics from the DSS
    % (projected into channel space) data.
    %
    % trigger_list: couples (5,15) or (10, 20).
    % * 5 : 3 second RAND sequences.
    % * 10: 15 second RAND sequences.
    % * 15: 3 second REG sequences.
    % * 20: 15 second REG sequences.
    %
    % subject_list:
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
    % Last update: 06/August/2018
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%% EDIT THE FOLLOWING CODE %%%%%%%%%%%%%%%%%%%%%
    hpfreq           =  2;
    lpfreq           =  30;
    n_components     =  3;   % DSS components to use (they are precomputed in pre_DSSprojection.m).
    SHORT            =  0;   % Set to 1 to work with SHORT sequences. 0 to LONG sequences.
    store_output     =  0;   % Set to 1 to store the resultant graphs.
    single           =  0;   % Set to 1 to compute a DSS matrix for each one of the conditions.
    % ____________________________________________________________________________
    % SUBJECT           | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 | 10 | 11 | 12 | 13 | 15 |
    % subject_selection | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 |  9 | 10 | 11 | 12 | 13 |
    subject_selection = [1:13]; 
    perc_vector = [0.05, 0.01];   % P-values vector.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%% NO MORE PARAMETERS BELOW %%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %% Results
    out_folder       =  sprintf('Trigger_analysis_PRE_HP%d_LP%d', hpfreq, lpfreq);  % Output data folder. High passed data (2Hz). Detrended.
    channels_path    =  fullfile('..','Results','Channels_DSS');
    fs               =  600; % Sampling frequency
    
    % Depending on the condition, we pick proper trigger values.
    if SHORT == 1  % SHORT sequences
        trigger_list = [5,15];
    else           % LONG sequences
        trigger_list = [10, 20];
    end

    if isdir(fullfile('..','Results',out_folder,'DSS_components','Transformed')) == 1 | ...
            (single == 0 & exist(fullfile('..','Results',out_folder,'DSS_components','Transformed',sprintf('Xdss-TRIG_%d-%d-COMP_%d.mat',trigger_list(1), trigger_list(2), n_components))) == 1)...
            (single == 1 & exist(fullfile('..','Results',out_folder,'DSS_components','Transformed',sprintf('Xdss-TRIG_%d-%d-SINGLE-COMP_%d.mat',trigger_list(1), trigger_list(2), n_components))) == 1)
        if single == 0
         load(fullfile('..','Results',out_folder,'DSS_components','Transformed',sprintf('Xdss-TRIG_%d-%d-COMP_%d.mat',trigger_list(1), trigger_list(2), n_components)),'dss_comp','x_comp');
         
        else
            load(fullfile('..','Results',out_folder,'DSS_components','Transformed',sprintf('Xdss-TRIG_%d-%d-SINGLE-COMP_%d.mat',trigger_list(1), trigger_list(2), n_components)),'dss_comp','x_comp');
        end
    end
    subject_list = [2:13,15]; % DO NOT EDIT THIS VECTOR.
    %%
    
%     subject_selection = subject_list;
    figure
    time_data = (1:size(dss_comp,2))/fs-.2;

    for n_signals =1%:2
        for subject_ind = 1:length(subject_selection)
%             fprintf('SUBJ: %d. TRIG: %d', subject_list(subject_ind), trigger_list(trigger_ind))
            load(fullfile(channels_path, sprintf('Channels-SUBJ_%d', subject_list(subject_selection(subject_ind)))));
%             channels_num = 1:274;
            switch n_signals
                case 1
                    components = n_components;         
                    out_signal(:,subject_ind, :) = squeeze(rms(dss_comp(channels_num,:,subject_selection(subject_ind),:),1)*1e15);
                    k1 = 1;
                    k2 = .2;
                                               
                case 2
                    % Original signal
                    out_signal(:,subject_ind, :) = squeeze(rms(x_comp(channels_num,:,subject_selection(subject_ind),:),1)*1e15);
                    k1 = 5;
                    k2 = 10;
                                        
            end
                
        end
    
        if length(subject_selection) > 1
            % We use bootstrap and compute if there is a significant difference
            % between conditions.  
            Diff = out_signal(:,:,1) - out_signal(:,:,2);

            % We compute the statistics from the data.
            s = post_computeStats(Diff, perc_vector, trigger_list)  ;
        else
            s = ones(size(out_signal,1),2)*NaN;
        end
        std_vals = squeeze(std(out_signal,[],2))/sqrt(size(out_signal,2));
        colors = [0    0.4470    0.7410; ...
                  0.8500    0.3250    0.0980; ...
                  0.9290    0.6940    0.1250; ...
                  0.4940    0.1840    0.5560];

        
        for trigger_ind = 1:length(trigger_list)
           % Security check: the Stdev should be symmetrical around a ones
            % vector.
            % signal = ones(size(squeeze(mean(out_signal(:,:,trigger_ind),2))));
            signal = squeeze(mean(out_signal(:,:,trigger_ind),2));
            h = plot(time_data,signal,'Color',colors(trigger_ind,:), 'Linewidth',3); hold on;


        end
        
       

        if trigger_list(1) == 10
            hold on;
            h = plot(time_data, k1*s(:,1), 'Linewidth',3, 'Color',colors(3,:));
            hold on
            h =plot(time_data, k2*s(:,2),'Linewidth',3, 'Color',colors(4,:));
            
            if length(subject_selection) > 1
                legend('RAND','REG',sprintf('p=%.2f',perc_vector(1)),sprintf('p=%.2f', perc_vector(2)))     
            else
                legend('RAND','REG')    
            end
        else
            hold on
            h =plot(time_data, k2*s(:,2),'Linewidth',3, 'Color',colors(4,:));
            
            if length(subject_selection) > 1
                legend('RAND','REG',sprintf('p=%.2f',perc_vector(2)))  
            else
                legend('RAND','REG')  
            end
        end

        trigger_ind = 1;
        signal = squeeze(mean(out_signal(:,:,trigger_ind),2));            
        fill([time_data, fliplr(time_data)], [signal+std_vals(:,trigger_ind); flipud(signal-std_vals(:,trigger_ind))]',...
                 colors(trigger_ind,:), 'edgecolor','none','facealpha',0.2)
        set(gca,'children',flipud(get(gca,'children'))); % Send shaded areas in background
        trigger_ind = 2;
        signal = squeeze(mean(out_signal(:,:,trigger_ind),2));            
        fill([time_data, fliplr(time_data)], [signal+std_vals(:,trigger_ind); flipud(signal-std_vals(:,trigger_ind))]',...
                 colors(trigger_ind,:), 'edgecolor','none','facealpha',0.2)
        set(gca,'children',flipud(get(gca,'children'))); % Send shaded areas in background
        
        if trigger_list(1) == 5
            xlim([-0.2, 4])
        else
            xlim([-0.2, 16])
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
        mkdir(fullfile('..','Results',out_folder,'DSS_components','Signal_curves'))
        savefig(fullfile('..','Results',out_folder,'DSS_components','Signal_curves',sprintf('DSS-Comp%d-TRIG%d-%d_ALLCHANNELS',components, trigger_list(1), trigger_list(2))))
                
%         save(fullfile('..','Results',out_folder,'DSS_components',sprintf('2DSS-Comp%d-TRIG%d-%d',components, trigger_list(1), trigger_list(2))), 'out_signal')
    end
    
    

    