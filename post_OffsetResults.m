function post_OffsetResults(trigger_list, subject, config)
    % This function reads the epoched (and preprocessed) data and changes
    % its temporal reference. Instead of 0 seconds, we set it to be placed
    % in the offset temporal area (15 seconds).
    %    
    % trigger_list: 
    % * [5,15] : for SHORT sequences.
    % * [10, 20] : for LONG sequences.
    %
    % subject: we can choose a single subject or a list of them.
    % * 2-15, the index of the subject.
    %
    % config: allows for certain configurations.
    %   .out_folder: name of the folder where data will be stored.
    %   .store_data: set to 1 to store in a .mat file the whole data of the
    %       subject once that it has been processed.
    %   .hpfreq: represents the high pass frequency. In this case is used
    %       to store the data.
    %   .lpfreq: represents the low pass frequency. In this case is used
    %       to store the data.
    %
    % OUTPUT FOLDER:
    % ../Results/ToneGraphs
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
    hpfreq          = config.hpfreq;
    lpfreq          = config.lpfreq;
    
   

    %% We create some matrices with the timelock information.
    for trigger_ind = 1:length(trigger_list)
        for subject_ind = 1:length(subject)
            % We load the timelock information that we computed in a
            % previous stage.
            load(fullfile('..','Results',out_folder,'Offset',sprintf('Timelock-TRIG_%d-SUBJ_%d.mat',trigger_list(trigger_ind),subject(subject_ind))),'timelock')
            data(:, :, subject_ind, trigger_ind) = timelock.avg;
        end    
    end
    data = data*1e15;
    aux = squeeze(rms(data,1));
   
    Diff =  aux(:,:,1) - aux(:,:,2);
    dataB=bootstrap(Diff'); 

    perc = .05;
    s=findSigDiff(dataB, perc);
    perc2 = .01;
    s2=findSigDiff(dataB, perc2);

    t_indexes = 301:901; % We plot from 14.5 to 15.5.
    figure('units','normalized','outerposition',[0 0 1 1])
    plot(timelock.time(t_indexes), squeeze(mean(aux(t_indexes,:,:),2)), 'Linewidth',3);
    hold on
    plot(timelock.time(t_indexes),abs(s(t_indexes)), 'Linewidth',3)
    hold on
    plot(timelock.time(t_indexes),2*abs(s2(t_indexes)), 'Linewidth',3)
    legend('RAND','REG','p=0.05','p=0.01')
    title(sprintf('Offset analysis. HP: %d Hz. LP: %d Hz.', hpfreq, lpfreq));
    grid on
    if store_output == 1
        mkdir(fullfile('..','Results','ToneGraphs'))
        savefig(fullfile('..','Results','ToneGraphs',sprintf('OFFSET-HP_%d-LP_%d.fig', hpfreq, lpfreq)))
        print(fullfile('..','Results','ToneGraphs',sprintf('OFFSET-HP_%d-LP_%d.emf', hpfreq, lpfreq)), '-dmeta','-opengl')
%         keyboard

    end
    
end
