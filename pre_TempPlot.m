function pre_TempPlot(trigger_list, subject, config)
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
    % OUTPUT FOLDER:
    % ../Results/***/Timelock_graphs
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
    % Last update: 07/June/2018
    
    out_folder = config.out_folder;
    store_output = config.store_data;
    
    
    %% We create some matrices with the timelock information.
    for trigger_ind = 1:length(trigger_list)
        for subject_ind = 1:length(subject)
            % We load the timelock information that we computed in a
            % previous stage.
            load(fullfile('..','Results',out_folder,'Timelock',sprintf('Timelock-TRIG_%d-SUBJ_%d.mat',trigger_list(trigger_ind),subject(subject_ind))),'timelock')

            switch trigger_list(trigger_ind) 
                case 5
                    out_5(subject_ind,:) = rms(timelock.avg);
                    time_5(subject_ind,:) = timelock.time;
                case 10
                    out_10(subject_ind,:) = rms(timelock.avg);
                    time_10(subject_ind,:) = timelock.time;
                case 15
                    out_15(subject_ind,:) = rms(timelock.avg);
                    time_15(subject_ind,:) = timelock.time;
                case 20
                    out_20(subject_ind,:) = rms(timelock.avg);
                    time_20(subject_ind,:) = timelock.time;
            end      
        end    
    end

    %% We plot data.
    if exist('out_5') == 1
        figure('units','normalized','outerposition',[0 0 1 1])   
        % SHORT condition (3 seconds long).
        subplot(211)
        plot(time_5(1,:),mean(out_5,1)*1e15); hold on;
        plot(time_15(1,:),mean(out_15,1)*1e15)
        xlabel('Time (s)');
        ylabel('Magnitude (fT)');
        xlim([min(time_5(1,:)), max(time_5(1,:))])
        legend('RAND','REG')
        if length(subject) == 1
            title(sprintf('Timelock. SHORT. SUBJ: %d',subject))
        else
            title(sprintf('Timelock. SHORT. GLOBAL'))
        end

        % LONG condition (15 seconds long).
        subplot(212)
        plot(time_10(1,:),mean(out_10,1)*1e15); hold on;
        plot(time_20(1,:),mean(out_20,1)*1e15)
        xlabel('Time (s)');
        ylabel('Magnitude (fT)');  
        xlim([min(time_10(1,:)), max(time_10(1,:))])
        legend('RAND','REG')
        if length(subject) == 1
            title(sprintf('Timelock. LONG. SUBJ: %d',subject))
        else
            title(sprintf('Timelock. LONG. GLOBAL'))
        end      
    else
        plot(time_10(1,:),mean(out_10,1)*1e15); hold on;
        plot(time_20(1,:),mean(out_20,1)*1e15)
        xlabel('Time (s)');
        ylabel('Magnitude (fT)');  
        xlim([min(time_10(1,:)), max(time_10(1,:))])
        legend('RAND','REG')
        if length(subject) == 1
            title(sprintf('Timelock. LONG. SUBJ: %d',subject))
        else
            title(sprintf('Timelock. LONG. GLOBAL'))
        end      
    end
        

    
    
    %% We can store the block information in a .fig file.
    if store_output == 1
        mkdir(fullfile('..','Results',out_folder,'Timelock_graphs'))
        if length(subject) == 1
            savefig(fullfile('..','Results',out_folder,'Timelock_graphs',sprintf('Timelock-SUBJ_%d.fig',subject)))
        else
            savefig(fullfile('..','Results',out_folder,'Timelock_graphs',sprintf('Timelock-GLOBAL.fig')))
        end
    end
    

    
end
