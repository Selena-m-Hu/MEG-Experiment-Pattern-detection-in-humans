function [] = PAR_Trigger_analysis(trigger, store)

        dataset = 2; % 2 - Long dataset ; 1 - Short dataset;
%         filenum = [2,5,13]; % Subject identifier. Currently: 2, 15.
        filenum = [11,12,14,15];

%     store = 1;
    out_folder = 'Trigger_analysis_Post16_Pre500';
    erroneous_block = [];
    if (abs(store)== 1) | (abs(store) == 3)
        cfg = [];
        cfg.feedback = 'no';
        cfg.channel = 'MEG';


        error_counter = 1;
    %     trigger = [5]; % Trigger signal to read: (5,15) for short events - (10,20) for long events.
        subjects = [];
        blocks = [];
        for ind_subjects = 1:length(filenum)
           switch filenum(ind_subjects)
                case 2
                    n_blocks = 6;
                case 9 
                    n_blocks = 8;
               case  10
                    n_blocks = 6; 
                otherwise
                    n_blocks = 7;
            end
            [subjects_aux, blocks_aux ] = ndgrid( filenum(ind_subjects), 2:n_blocks);
            subjects = [subjects; subjects_aux(:)];
            blocks = [blocks; blocks_aux(:)];
        end
                 
        for trigger_ind = 1:length(trigger)   
            data_subject = [];
            
            if store == -1
                store_counter = 0;
                data_subject = [];
            end

%             parfor ind_subjects = 1:length(subjects)
%                 data_block(ind_subjects) = Trigger_reader(dataset, subjects(ind_subjects), blocks(ind_subjects), trigger(trigger_ind));
%             end 
%             
%             for ind_subjects = 1:length(subjects)
%                 if iscell(data_block(ind_subjects).trial) == 0
%                     error(error_counter) = ind_subjects;
%                     error_counter = error_counter +1;
%                 else                    
%                     if isempty(data_subject) 
%                         data_subject = data_block(ind_subjects);
%                     else
% 
%                         data_subject = ft_appenddata(cfg, data_subject, data_block(ind_subjects));
%                     end
%                 end
%             end
            counter = 1;
            data_out = [];
            for ind_subjects = 1:length(subjects)
                if ind_subjects == 1
                    sub_num = subjects(ind_subjects);
                end
                 
                if sub_num ~= subjects(ind_subjects)
                    aux = ft_rejectvisual(cfg,data_subject);
                    aux_sum = aux.trial{1};
                    for ind_aux = 2:length(aux.trial)
                        aux_sum = aux_sum + aux.trial{ind_aux};
                    end
                    len_sum = length(aux.trial);
                    mkdir(fullfile('..','Results',out_folder));
                    save(fullfile('..','Results',out_folder,sprintf('Long_timelock-DATA_%d-TRIG_%d-SUBJ_%d.mat',dataset,trigger,sub_num)),'aux_sum','len_sum','aux')
%                     save(fullfile('..','Results',out_folder,sprintf('Long_timelock-DATA_%d-TRIG_%d-SUBJ_%d.mat',dataset,trigger,sub_num)),'aux_sum','len_sum')
                     
                    clear aux len_sum aux_sum
                    data_subject = [];

                    sub_num = subjects(ind_subjects);
                end
                
                data_block = Trigger_reader(dataset, subjects(ind_subjects), blocks(ind_subjects), trigger(trigger_ind));
%                 counter = counter +1;
                if iscell(data_block.trial) == 0
                    error(error_counter) = ind_subjects;
                    error_counter = error_counter +1;
                else                    
                    if isempty(data_subject) 
                        data_subject = data_block;
                    else
                        data_subject = ft_appenddata(cfg, data_subject, data_block);
                    end
                end
                
                if  ind_subjects == length(subjects) % For the last iteration...
                                      
                    aux = ft_rejectvisual(cfg,data_subject);
                    aux_sum = aux.trial{1};
                    for ind_aux = 2:length(aux.trial)
                        aux_sum = aux_sum + aux.trial{ind_aux};
                    end
                    len_sum = length(aux.trial);
%                     save(fullfile('..','Results',out_folder,sprintf('Long_timelock-DATA_%d-TRIG_%d-SUBJ_%d.mat',dataset,trigger,sub_num)),'aux_sum','len_sum')
                    save(fullfile('..','Results',out_folder,sprintf('Long_timelock-DATA_%d-TRIG_%d-SUBJ_%d.mat',dataset,trigger,sub_num)),'aux_sum','len_sum','aux')
                    %{
                        We store a file with the structure, that will be used
                        during the visual rejection of the channels (not the
                        trials, since that rejection has already been
                        performed).
                    %}                    
                    save(fullfile('..','Results',out_folder,sprintf('BAIT-DATA_%d-TRIG_%d.mat',dataset,trigger)),'aux')
                     
                    clear aux len_sum aux_sum
                    data_subject = [];

                    sub_num = subjects(ind_subjects);
                end
                
                clear data_block 
            end 
            clear data_subject
         


        end
    
        % Extraemos algunos estadísticos de los trials. Por ejemplo, el promedio de
        % la magnitud.
    

%         config.dataset = dataset;
%         config.filenum = filenum;
%         config.trigger = trigger;
%         config.errors = erroneous_block;
        mkdir(fullfile('..','Results',out_folder));
% %         save(fullfile('..','Results','Trigger_analysis',sprintf('Long_timelock-DATA:%d-TRIG:%d',dataset,trigger)),'timelock','config')
%         save(fullfile('..','Results','Trigger_analysis',sprintf('Long_timelock-DATA_%d-TRIG_%d.mat',dataset,trigger)),'timelock','config')
        
                    
    elseif store == 4
        
        for ind = 1:length(filenum)
           load(fullfile('..','Results',out_folder,sprintf('Long_timelock-DATA_%d-TRIG_%d-SUBJ_%d.mat',dataset,trigger,filenum(ind))))
           if ind == 1
              trial_mean = aux_sum;
              len_mean = len_sum;
           else
               trial_mean = trial_mean + aux_sum;
               len_mean = len_mean + len_sum;
           end

        end
        trial_mean = trial_mean/len_mean;

        load(fullfile('..','Results',out_folder,sprintf('BAIT-DATA_%d-TRIG_%d.mat', dataset,trigger)))
        for ind = 1:length(aux.trial)
            aux.trial{ind} = trial_mean;
        end

        % We reject data visually, and store the resutls for this trigger.
        cfg=[];
        cfg.method = 'summary';
        cfg.channel = {'MEG'};
        data_out = ft_rejectvisual(cfg,aux);


        timelock = ft_timelockanalysis(cfg, data_out);
        clear data_subject data_block data_subject

                
        config.dataset = dataset;
        config.filenum = filenum;
        config.trigger = trigger;
%         config.errors = erroneous_block;
        mkdir(fullfile('..','Results',out_folder));
%         save(fullfile('..','Results','Trigger_analysis',sprintf('Long_timelock-DATA:%d-TRIG:%d',dataset,trigger)),'timelock','config')
        save(fullfile('..','Results',out_folder,sprintf('Long_timelock-DATA_%d-TRIG_%d.mat',dataset,trigger)),'timelock','config')
        
                    
                
                
                
                
        
    else
        dataset = 2;
        figure('units','normalized','outerposition',[0 0 1 1])
        for trigger_ind = 1:length(trigger)
            load(fullfile('..','Results',out_folder,sprintf('Long_timelock-DATA_%d-TRIG_%d',dataset,trigger(trigger_ind))))
            plot(timelock.time, rms(timelock.avg));
            hold on;
            legend_string{trigger_ind} = sprintf('Trig: %d',trigger(trigger_ind));
        end
        legend(legend_string, 'Location','northwest')
        ylabel('RMS magnitude');
        xlabel('Time (s)');
        title('RMS magnitude across ALL channels')
        if (trigger(trigger_ind) == 5) | (trigger(trigger_ind) == 15)
            if length(trigger) > 1
                savefig(fullfile('..','Results',out_folder,'Short_analysis.fig'))
            else
                savefig(fullfile('..','Results',out_folder,sprintf('Short_analysis_TRIG-%d.fig',trigger(trigger_ind))))
            end
        else
            if length(trigger) > 1
                savefig(fullfile('..','Results',out_folder,'Long_analysis.fig'))
            else
                savefig(fullfile('..','Results',out_folder,sprintf('Long_analysis_TRIG-%d.fig',trigger(trigger_ind))))
            end
        end
        
        close all
        %% Manual selection of channels
        figure('units','normalized','outerposition',[0 0 1 1])
        for trigger_ind = 1:length(trigger)
            load(fullfile('..','Results',out_folder,sprintf('Long_timelock-DATA_%d-TRIG_%d',dataset,trigger(trigger_ind))))
            
            % Channel selection
    %         find(strcmp(timelock.label,'MLP34'))
%             switch trigger(trigger_ind)     
%                 case 20
%                     chns_L_down = [123,117,124,130,131,116,122,115,69,72];
%                     chns_R_up = [232,233,238,239,245,246,231,188,181,187];
%                     chns_L_up = [96,97,90,89,85,84,7,12,14,57];
%                     chns_R_down = [200,199,204,203,198,196,256,249,242,263];
%                 case 10
%                     chns_L_down = [117 124 123 116 109 110 130 122 69 72];
%                     chns_R_up = [232 239 233 238 246 188 187 181 138 229];
%                     chns_L_up = [112, 105,98,104,111,119,113,106,99,33];
%                     chns_R_down = [242 249 256 262 255 263 200 196 191 204];
%                 case 5
%                     chns_L_down = [106 105 102 104 110 111 109 101 69 66];
%                     chns_R_up = [209 207 206 203 212 202 205 160 166 167];
%                     chns_L_up = [7 12 16 87 82 86 85 81 77 74];
%                     chns_R_down = [216 221 215 220 180 176 183 179 173 170];
%                     
%                 case 15
%                     chns_L_down = [111 110 68 115 116 71 106 109 67 114];
%                     chns_R_up = [222 223 227 219 218 232 228 215 233 224];
%                     chns_L_up = [7 12 16 6 11 14 57 84 88 89];
%                     chns_R_down = [185 181 226 231 184 177 189 180 176 237];
%             end


            % Ordenamos las amplitudes, para lo cual escogemos un rango de instantes de
            timelock.fsample = 600;
            M100dat=timelock.avg'; % Obsérvese que aquí se transpone la media, por lo que tenemos (tiempo x canal).
            t0 = 0.2*timelock.fsample;
            % tiempo, y todos los canales.
             
            amps=mean(M100dat((t0+0.09*timelock.fsample):(t0+0.11*timelock.fsample), :),1);
            [ampsSorted,idx]= sort(amps,2,'descend');

            chnsSorted = timelock.label(idx);

            %selecting channels:

            chns_selectedLpos=[];
            chns_selectedRpos=[];
            chns_selectedLneg=[];
            chns_selectedRneg=[];
            leftChansCountPos=0;
            rightChansCountPos=0;
            leftChansCountNeg=0;
            rightChansCountNeg=0;

            for count=1:length(ampsSorted)
                strPos=chnsSorted(count);
                strNeg=chnsSorted(end-count+1);
                % Lookup in Left Hemisphere
                if  ~isempty(strfind(strPos{1},'ML'))
                    if(leftChansCountPos<10)
                        leftChansCountPos=leftChansCountPos+1;
                        chns_selectedLpos=[chns_selectedLpos strPos];
                    end
                end
                if ~isempty(strfind(strNeg{1},'MLT'))
                    if(leftChansCountNeg<10)
                        leftChansCountNeg=leftChansCountNeg+1;
                        chns_selectedLneg=[chns_selectedLneg strNeg];
                    end

                end

                % Lookup in Right Hemisphere
                if  ~isempty(strfind(strPos{1},'MRT'))
                    if(rightChansCountPos<10)
                        rightChansCountPos=rightChansCountPos+1;
                        chns_selectedRpos=[chns_selectedRpos strPos];
                    end
                end
                if ~isempty(strfind(strNeg{1},'MR'))
                    if(rightChansCountNeg<10)
                        rightChansCountNeg=rightChansCountNeg+1;
                        chns_selectedRneg=[chns_selectedRneg strNeg];
                    end
                end
            end

            chns_selectedL = [chns_selectedLpos chns_selectedLneg];
            chns_selectedR = [chns_selectedRpos chns_selectedRneg];

            chnsL_num=[];
            for count1=1:length(timelock.label)
                for count2=1:length(chns_selectedL)
                    if (strcmp(timelock.label{count1},chns_selectedL{count2}) ~= 0)
                        chnsL_num=[chnsL_num count1];
                    end
                end
            end

            chnsR_num=[];
            for count1=1:length(timelock.label)
                for count2=1:length(chns_selectedR)
                    if (strcmp(timelock.label{count1},chns_selectedR{count2}) ~= 0)
                        chnsR_num=[chnsR_num count1];
                    end
                end
            end

            chns_selectedL = timelock.label(chnsL_num); 
            chns_selectedR = timelock.label(chnsR_num);

%             chns_selectedL = unique([chns_L_up chns_L_down]);
%             chns_selectedR = unique([chns_R_up chns_R_down]);

            channel = ft_channelselection(unique([chns_selectedL, chns_selectedR]), timelock.label);

            verbose = 0
            if verbose == 1
                figure(3);
                cfg = [];
                cfg.parameter = 'avg';
                cfg.layout='CTF275.lay';
                cfg.xlim=[0.04, 0.08]';
                cfg.marker = 'labels';
                cfg.interactive = 'yes';
                cfg.colorbar = 'yes';

                cfg.markerfontsize = 8;
                cfg.highlight='on';
                cfg.highlightchannel=channel;
                cfg.highlightfontsize=20;
                subplot(121)
                ft_topoplotER(cfg, timelock); title ('M100 response');
                subplot(122)
                t_chunk = timelock.time;
                t_chunk = t_chunk(t_chunk < .5);
                plot(t_chunk,timelock.avg([chnsL_num, chnsR_num],t_chunk < .5))
            end
            
            out_signal(trigger_ind,:) = rms(timelock.avg(unique([chns_selectedL, chns_selectedR]),:));

            plot(timelock.time, out_signal(trigger_ind,:));
            hold on;
            legend_string{trigger_ind} = sprintf('Trig: %d',trigger(trigger_ind));
            legend(legend_string, 'Location','northwest')
            ylabel('RMS magnitude');
            xlabel('Time (s)');
            title('RMS magnitude across 40 channels')             
            
        end
        
        if (trigger(trigger_ind) == 5) | (trigger(trigger_ind) == 15)
            if length(trigger) > 1
                savefig(fullfile('..','Results',out_folder,'Short_analysis_TOP-channels.fig'))
            else
                savefig(fullfile('..','Results',out_folder,sprintf('Short_analysis_TOP-channels_TRIG-%d.fig',trigger(trigger_ind))))
            end
        else
            if length(trigger) > 1
                savefig(fullfile('..','Results',out_folder,'Long_analysis_TOP-channels.fig'))
            else
                savefig(fullfile('..','Results',out_folder,sprintf('Long_analysis_TOP-channels_TRIG-%d.fig',trigger(trigger_ind))))
            end
        end
        close all
         %% PSD estimation
        contador_subplot = 1;
        figure('units','normalized','outerposition',[0 0 1 1])        
        for trigger_ind = 1:length(trigger)
            fs = 600;
            window = fs;% 1 second window.
            noverlap = fs/2;
            nfft = 1024;
            % Welch.
            [pxx_welch,f_welch] = pwelch(out_signal(trigger_ind,:),window,noverlap,nfft, fs);

            % Periodogram
            nfft = 2^nextpow2(length(out_signal(trigger_ind,:)));
            pxx = abs(fft(out_signal(trigger_ind,:),nfft)).^2/length(out_signal(trigger_ind,:))/fs;
            f_per = (1:nfft)/nfft*fs;
            
            % Fieildtrip
            cfg = [];
            cfg.method = 'mtmfft';
            cfg.output = 'pow';
            cfg.taper = 'hanning'
%             cfg.tapsmofrq = 5
% cfg.keeptrials = 'yes'
            cfg.channel = channel;
            cfg.foilim = [.1, 30]
             [freq] = ft_freqanalysis(cfg, aux)
             semilogy(freq.freq ,rms(freq.powspctrm))
%             plot(timelock.time, rms(timelock.avg));
            hold on;
            
            
            
            subplot(2, max(1,length(trigger)), contador_subplot)
            loglog(f_per,pxx)
            grid on
            hold on
            loglog(f_welch, pxx_welch,'Linewidth',3)
            grid on
            ylabel('Magnitude (log)')
            title(sprintf('PSD estimation. TRIG: %d',trigger(trigger_ind)))
            legend('Periodogram','Welch method')
            contador_subplot = contador_subplot + length(trigger);
            
            subplot(2, max(1,length(trigger)), contador_subplot)
            loglog(f_per((f_per > 1) & (f_per < 30)) ,pxx((f_per > 1) & (f_per < 30)))
            grid on
            hold on
            loglog(f_welch((f_welch > 1) & (f_welch < 30)), pxx_welch((f_welch > 1) & (f_welch < 30)),'Linewidth',3)
            grid on
            ylabel('Magnitude (log)')
            xlabel('Frequency (Hz)')
            title(sprintf('Cutoff freq: 30 Hz. TRIG: %d',trigger(trigger_ind)))
            legend('Periodogram','Welch method')
            contador_subplot = contador_subplot - length(trigger)/2;
        
        end
        if (trigger(trigger_ind) == 5) | (trigger(trigger_ind) == 15)
            if length(trigger) > 1
                savefig(fullfile('..','Results',out_folder,'Short_PSD_TOP-channels.fig'))
            else
                savefig(fullfile('..','Results',out_folder,sprintf('Short_PSD_TOP-channels_TRIG-%d.fig',trigger(trigger_ind))))
            end
        else
            if length(trigger) > 1
                savefig(fullfile('..','Results',out_folder,'Long_PSD_TOP-channels.fig'))
            else
                savefig(fullfile('..','Results',out_folder,sprintf('Long_PSD_TOP-channels_TRIG-%d.fig',trigger(trigger_ind))))
            end
        end
    close all
        
    end
    
    
end



function [data_block] = Trigger_reader(dataset, filenum, block, trigger)
% clear all;
% close all;
addpath(genpath(fullfile('..','fieldtrip')));
if nargin == 2
    block = [];
elseif nargin == 1
    data_folder = dataset.data_folder;
    username = dataset.username;
    dataset = -1;
    filenum = 1;
end
sampling = 600;
switch dataset 
    case 1
        data_folder = fullfile('..','data');
        switch filenum
            case 1
                username = 'ac040981_MChait2_20120417_01.ds';
            case 2  
                username = 'cp071189_MChait2_20120416_01.ds'; %% FILE NOT WORKING
            case 3
                username = 'lt290589_MChait2_20120417_01.ds';
            case 4
                username = 'ns070382_MChait2_20120423_01.ds';
            case 5
                username = 'rs120685_MChait2_20120417_01.ds';
            otherwise
                if filenum < 0
                    username = fullfile('Loc',...
                                        sprintf('Subj%d', abs(filenum)),...
                                        sprintf('Subj%d_loc', abs(filenum)));
                end
        end
    case 2
        data_folder = fullfile('..','data_longshort', sprintf('Subj%d/',filenum));
        filelist_bad = dir([data_folder,'*.ds']);
        contador = 1;
        for ind = 1:length(filelist_bad)
            if filelist_bad(ind).name(1)~= '.'
                filelist{contador} = filelist_bad(ind).name;
                contador = contador+1;
            end
        end
        fprintf('** Number of blocks: %d **\n', length(filelist));
        if block <= length(filelist)
            fprintf('** Reading block %d **\n', block);
            username = filelist{block};
        else
            fprintf('** BLOCK NUMBER OUT OF BOUNDS **\n');
            keyboard
        end
    
end

cfg = [];
cfg.feedback = 'no';
cfg.channel = 'MEG'; % Le decimos que queremos capturar todos los canales de MEG.
% filenum = 4;
switch dataset
    case 1
        signature = 'Analysis_Wave4'; % 'params'
    case 2
%         signature = 'Dataset2-Long_first';
        signature = 'Dataset2-Long_first-DETREND';
    case -1
        signature = 'null';
end

if filenum > 0
    cfg.dataset = fullfile(data_folder,username);  % Fichero a leer.
else
    load(fullfile('..','data',username));
end

% We create the results directory if it doesn't exist.
if isdir(fullfile('..','Results')) == 0
    mkdir(fullfile('..','Results'))    
end

if isdir(fullfile('..','Results', signature)) == 0
    mkdir(fullfile('..','Results', signature))    
end

if isdir(fullfile('..','Results', signature, username)) == 0
    mkdir(fullfile('..','Results', signature, username))    
end

% cfg.dataset = fullfile(dataset,username);
% data = ft_preprocessing(cfg);
% plot((1:length(data.trial{1}))/600,data.trial{1}(250:255,:)')



hdr   = ft_read_header(cfg.dataset);
event = ft_read_event(cfg.dataset);
tr_channel = 'UPPT001';

counter = 1;
for i = 1:length(event)
    if strcmp(event(i).type, tr_channel)% & event(i).value <= 40
        sample(counter) = event(i).sample;
        value(counter) = event(i).value;
        counter = counter+1;
    end
        
end
clear hdr event

%% Trigger signals
verbose = 0;
if verbose == 1
    subplot(211)
    sample_time = sample/sampling;
    stem(sample_time, value)
    xlabel('Time (s)');
    ylabel('Trigger val');
    mean((sample(1:end-1)-sample(2:end))/600)
    title(sprintf('Time diff.(s) Mean: %.2f, Std: %.2f', mean(diff(sample_time)), std(diff(sample_time))));
    diff(sample_time)
    subplot(212)
    channel_val =5;
    sample_time = sample(value == channel_val)/sampling;
    stem( sample_time, value(value == channel_val)) 
    xlabel('Time (s)');
    ylabel('Trigger val')
    title(sprintf('Time diff.(s) Mean: %.2f, Std: %.2f', mean(diff(sample_time)), std(diff(sample_time))));
    % fprintf('Mean delay: %.2f\n', mean(diff(sample_time)))
    diff(sample_time)
end
%% We read data according to each trigger value.
% We begin with number 5.
% data = ft_preprocessing(cfg);

% cfg.dataset = fullfile('am061187_MChait2_20120510_01.ds');
cfg.trialdef.eventtype  = 'UPPT001'; % Aquí indicamos el itpo de evento que nos ayudará a separar en epochs.
cfg.trialdef.eventvalue = trigger;

if (trigger == 5) | (trigger == 15)
    cfg.trialdef.prestim    = 0.5; % Espacio de preestímulo (antes del estímulo). Puede usarse para extraer la media.
    cfg.trialdef.poststim   = 3.5%.5; % Espacio de señal que nos interesa, después del evento.
else
    cfg.trialdef.prestim    = 0.5; % Espacio de preestímulo (antes del estímulo). Puede usarse para extraer la media.
    cfg.trialdef.poststim = 16;
end
try
    cfg = ft_definetrial(cfg);
    % cfg = [];
    cfg.demean      = 'yes';
    % cfg.detrend = 'yes';
    cfg.baselinewindow = [-cfg.trialdef.prestim 0];% in seconds
    data = ft_preprocessing(cfg);
    cfg=[];
    cfg.lpfilter ='yes';
    cfg.lpfreq = 30;
    data_block=ft_preprocessing(cfg, data);
catch
%     data_block = struct('hdr',{},'fsample',{},'sampleinfo',{},'trialinfo',{},'grad',{},'trial',-1,...
%                         'time',{},'label',{},'cfg',{});   
%                     data_block.trial = -1:;
    data_block.hdr = -1;
    data_block.fsample = -1;
    data_block.sampleinfo = -1;
    data_block.trialinfo = -1;
    data_block.grad = -1;
    data_block.trial = -1;    
    data_block.time = -1;
    data_block.label = -1;
    data_block.cfg = -1;    
    
    
    
    
    
    
    
end
end
% We append data to the final storage structure.
% ft_appenddata(cfg, data_pre, data)





