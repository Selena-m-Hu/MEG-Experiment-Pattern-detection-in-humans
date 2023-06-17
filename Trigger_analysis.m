function [] = Trigger_analysis(trigger, store)

%     store = 1;
    erroneous_block = [];
    if (abs(store)== 1) | (abs(store) == 3)
        cfg = [];
        cfg.feedback = 'no';
        cfg.channel = 'MEG';

        dataset = 2; % 2 - Long dataset ; 1 - Short dataset;
%         filenum = [2,3,4,5,6,9,15]; % Subject identifier. Currently: 2, 15.
        filenum = [2:15]; % Subject identifier. Currently: 2, 15.


        error_counter = 1;
    %     trigger = [5]; % Trigger signal to read: (5,15) for short events - (10,20) for long events.
        for trigger_ind = 1:length(trigger)   
            data_subject = [];
            
            if store == -1
                store_counter = 0;
                data_subject = [];
            end

            for ind_subjects = 1:length(filenum)
                switch filenum(ind_subjects)
                    case 2
                        n_blocks = 6;
                    case 9 
                        n_blocks = 8;
                    otherwise
                        n_blocks = 7;
                end
                if abs(store) == 1
                    for block = 2:n_blocks
                        try
                            data_block = Trigger_reader(dataset, filenum(ind_subjects), block, trigger(trigger_ind));
                            if store == 1
                                if isempty(data_subject) 
                                    data_subject = data_block;
                                else
                                    data_subject = ft_appenddata(cfg, data_subject, data_block);
                                end
                            else
                                if isempty(data_subject) 
                                    aux = data_block.trial;
                                    data_subject = zeros(size(aux{1}));
                                end

                                for ind_store = 1:length(aux)
                                    data_subject = data_subject + aux{ind_store};

                                end
                                store_counter = store_counter + length(aux);

                            end
    %                         clear data_block
                        catch
                            erroneous_block(error_counter).trigger = trigger(trigger_ind);
                            erroneous_block(error_counter).subject = filenum(ind_subjects);
                            erroneous_block(error_counter).block = block;
                            error_counter = error_counter + 1;
        %                     keyboard
                        end
                    end
               
                    
                elseif store == 3
                    parfor block = 2:n_blocks
%                         try
                            data_block(block) = Trigger_reader(dataset, filenum(ind_subjects), block, trigger(trigger_ind));
%                             if store == 1
%                                 if isempty(data_subject) 
%                                     data_subject = data_block;
%                                 else
%                                     data_subject = ft_appenddata(cfg, data_subject, data_block);
%                                 end
%                                 
%                             end
%                         end
                    end
                    
                    
                    for block = 2:n_blocks
                        if isempty(data_subject) 
                            data_subject = data_block(block);
                        else
                            data_subject = ft_appenddata(cfg, data_subject, data_block(block));
                        end
                                
                    end
                    clear data_block
                end
            end

            
            if (store == 1) | (store == 3)
                % We reject data visually, and store the resutls for this trigger.
                cfg=[];
                cfg.method = 'summary';
                cfg.channel = {'MEG'};
                data_out(trigger_ind) = ft_rejectvisual(cfg,data_subject);


                timelock(trigger_ind) = ft_timelockanalysis(cfg, data_out(trigger_ind));
                clear data_subject data_block data_subject
                
            else
%                 keyboard
                cfg=[];
                cfg.method = 'summary';
                cfg.channel = {'MEG'};
                timelock(trigger_ind) = ft_timelockanalysis(cfg, data_block);
               
                timelock(trigger_ind).avg = data_subject/store_counter;
                
                
                
                
                
                
                
            end
        end
        % Extraemos algunos estadísticos de los trials. Por ejemplo, el promedio de
        % la magnitud.
    

        config.dataset = dataset;
        config.filenum = filenum;
        config.trigger = trigger;
        config.errors = erroneous_block;
        mkdir(fullfile('..','Results','Trigger_analysis'));
%         save(fullfile('..','Results','Trigger_analysis',sprintf('Long_timelock-DATA:%d-TRIG:%d',dataset,trigger)),'timelock','config')
        save(fullfile('..','Results','Trigger_analysis',sprintf('Long_timelock-DATA_%d-TRIG_%d.mat',dataset,trigger)),'timelock','config')
        
        
        
    else
        dataset = 2;
        figure('units','normalized','outerposition',[0 0 1 1])
        for trigger_ind = 1:length(trigger)
            load(fullfile('..','Results','Trigger_analysis',sprintf('Long_timelock-DATA:%d-TRIG:%d',dataset,trigger(trigger_ind))))
            plot(timelock.time, rms(timelock.avg));
            hold on;
            legend_string{trigger_ind} = sprintf('Trig: %d',trigger(trigger_ind));
        end
        legend(legend_string, 'Location','northwest')
        ylabel('RMS magnitude');
        xlabel('Time (s)');
        if (trigger(trigger_ind) == 5) | (trigger(trigger_ind) == 15)
            savefig(fullfile('..','Results','Trigger_analysis','Short_analysis.fig'))
        else
            savefig(fullfile('..','Results','Trigger_analysis','Long_analysis.fig'))
        end
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
cfg.trialdef.prestim    = 0.2; % Espacio de preestímulo (antes del estímulo). Puede usarse para extraer la media.
if (trigger == 5) | (trigger == 15)
    cfg.trialdef.poststim   = 3%.5; % Espacio de señal que nos interesa, después del evento.
else
    cfg.trialdef.poststim = 15;
end
cfg = ft_definetrial(cfg);
% cfg = [];
cfg.demean      = 'yes';
% cfg.detrend = 'yes';
cfg.baselinewindow = [-.2 0];% in seconds
data = ft_preprocessing(cfg);
cfg=[];
cfg.lpfilter ='yes';
cfg.lpfreq = 30;
data_block=ft_preprocessing(cfg, data);
end
% We append data to the final storage structure.
% ft_appenddata(cfg, data_pre, data)





