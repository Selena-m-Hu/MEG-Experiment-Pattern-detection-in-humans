%this script analyzes the M100 response in the localizer script and finds
%the 20 strongest channels in each hemisphere.

clear all;close all;
addpath(fullfile('..','fieldtrip'));
data_folder = fullfile('..','data');
cfg = [];
cfg.channel = 'MEG'; % Le decimos que queremos capturar todos los canales de MEG.
filenum = 5;

signature = 'Pre_analysis_Whole4';

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
% cfg.dataset = fullfile('am061187_MChait2_20120510_01.ds');

% 
%%
if filenum > 0
    cfg.lpfilter ='yes';
    cfg.lpfreq = 30;
    data = ft_preprocessing(cfg);
%     data.trial{1} = ones(size(data.trial{1}));
    %% Aquí podemos preprocesar la señal a lo bruto.
    pre_mode = 4;
    % We load the preselected channels.
    load(fullfile('..','Results', 'params', username,'Channel_selection.mat'));
    
    data_denoised = data;
    data_out = zeros(size(data.trial{1}));%zeros(length(params.chnsL_num)+length(params.chnsR_num),size(data.trial{1},2));
    selection = abs([params.chnsL_num, params.chnsR_num]);
    sign_selection = [params.sign_L, params.sign_R];
    switch pre_mode
        case 1
            % Wavelets.
            wname = 'coif4'; % db4, coif4, sym8, fk14, bior3.7, rbio1.5 (desfasa)
            lev = 4;
            % We denoise using wavelets.
            % TRABAJAR EXCLUSIVAMENTE CON LOS CANALES QUE HEMOS
            % PRESELECCIONADO EN EL ANÁLISIS ANTERIOR. ES IMPORTANTE
            % HACERLO DE ESTA FORMA, YA QUE EN LUGAR DE 180-240 CANALES
            % PASAMOS A PROCESAR EXCLUSIVAMENTE 40. TENIENDO EN CUENTA QUE
            % DSS TARDA BASTANTE TIEMPO, E  STO SERÁ BENEFICIOSO.
            
            contador = 1;
            for ind =abs([params.chnsL_num, params.chnsR_num])%1:size(data.trial{1},1)% 
%                 [data_out(ind, :),~,~] = wden(data.trial{1}(ind,:),'minimaxi','s','mln',lev,wname);
                [data_out(ind, :),~,~] = wdenoise(data.trial{1}(ind,:),lev,'Wavelet',wname);
%                 contador = contador+1;
            end
            data_out = {data_out};
            
%             ind = 10;             
%             plot(data.trial{1}(ind,12000:14000)); hold on; 
%             plot(data_out(ind,12000:14000));
        case 2 % EMD
            
            data_denoised.trial{1} = data_denoised.trial{1}(selection,:);
            [data_emd] = EMD_denoise(data_denoised);
            
            data_out = data.trial;
            data_out{1}(selection,:) = data_emd.trial{1};
            clear data_denoised data_emd
        case 3 % DSS - Dudoso, ya que su coste computacional es elevado.
        case 4 
            % Adaptative filtering - Aún falta por aclarar qué modalidad 
            % del mismo se va a utilizar. Por una parte, LMS tarda bastante
            % tiempo en converger pero no presenta problemas con la
            % varianza temporal. En cambio, RLS es muy eficiente pero
            % requiere un mayor coste computacional y es inestable frente a
            % cambios temporales.
             
            n_elem = floor(length(selection)*.25);
            n_rep = 50;
            
            out_sig = [];
            for ind = 1:n_rep
                rand_weight = rand(1,n_elem); rand_weight = rand_weight/sum(rand_weight);
                rand_ind = randperm(length(selection), n_elem);
                if ind == 1
                    out_sig = ((rand_weight.*sign_selection(rand_ind))*data.trial{1}(selection(rand_ind),:))/n_rep;
                else
                    out_sig = out_sig + ((rand_weight.*sign_selection(rand_ind))*data.trial{1}(selection(rand_ind),:))/n_rep;
                end
            end
            
            clear x y
            for ind = 1:length(selection)
                  sample_size = 100;
                % With the Length of the filter we can control the
                % smoothness of the output signal. The bigger the length,
                % the smoother the output.
                lmsfilt2 = dsp.LMSFilter('Length',ceil(2*log2(sample_size)),'Method','LMS', ...
                    'StepSize',0.25);        
                
             
                for ind_iter = 1:sample_size:length(data.trial{1})
                    x(ind,ind_iter:ind_iter+sample_size-1) = (data.trial{1}(selection(ind),ind_iter:ind_iter+sample_size-1)*sign_selection(ind))';
                    [y(ind,ind_iter:ind_iter+sample_size-1),e(ind_iter:ind_iter+sample_size-1),w] = lmsfilt2(x(ind,ind_iter:ind_iter+sample_size-1)',...
                                       out_sig(ind_iter:ind_iter+sample_size-1)');  
                end
                
                
                
            end
            
            verbose = 0
            if verbose == 1
                ind_plot = 6000;
                offset = 1000;
                subplot(211)
                plot(mean(y(:,ind_plot:ind_plot+offset),1));
                subplot(212)
                plot(mean(x(:,ind_plot:ind_plot+offset),1));
            end         
            data_out = {y};
            data_denoised.cfg.channel = data_denoised.cfg.channel(selection);
            data_denoised.label = data_denoised.label(selection);
            
    end
    
    data_denoised.trial = data_out;
    
    %% Aquí troceamos la señal, una vez procesada.
    
%     cfg = []
    cfg.trialdef.eventtype  = 'UPPT001'; % Aquí indicamos el itpo de evento que nos ayudará a separar en epochs.
    cfg.trialdef.eventvalue = 10;

    cfg.trialdef.prestim    = 0.2; % Espacio de preestímulo (antes del estímulo). Puede usarse para extraer la media.
    cfg.trialdef.poststim   = .5; % Espacio de señal que nos interesa, después del evento.
    cfg = ft_definetrial(cfg);
    
    cfg.demean      = 'yes';
    
    cfg.baselinewindow = [-.2 0];% in seconds
    data = ft_redefinetrial(cfg,data);
    data_denoised = ft_redefinetrial(cfg,data_denoised);
    
    % We filtered before, so this is not necessary.
%     cfg=[];
%     cfg.hpfilter ='yes';
%     cfg.hpfreq = 0.5;
%     cfg.lpfilter ='yes';
%     cfg.lpfreq = 30;
%     data=ft_preprocessing(cfg, data);

    % save data
    cfg=[];
    cfg.method = 'summary';
    cfg.channel = 'MEG';
    % data = ft_rejectvisual(cfg,data);

end

% A CONTINUACIÓN, REPRESENTAMOS GRÁFICAMENTE LOS RESULTADOS FILTRADOS. HAY
% QUE PROMEDIAR LOS TRIALS QUE SE HAN EXTRAÍDO, POR LO QUE ESTA SECCIÓN DE
% GRÁFICOS DEBERÍA SER COMÚN PARA EMD, DSS Y WAVELET. AQUÍ PODEMOS
% INTRODUCIR TAMBIÉN EL FILTRADO ADAPTATIVO.

%% Left-right channel processing
% Aquí se extraen los canales, considerando que pueden tener magnitudes
% opuestas: esto es, el mismo gráfico pero con el signo cambiado a negativo
% o positivo. En cualquier caso, el resultado promedio debería ser el
% adecuado.
cfg = [];
cfg.channel = 'MEG';

% Extraemos algunos estadísticos de los trials. Por ejemplo, el promedio de
% la magnitud.
timelock = ft_timelockanalysis(cfg, data);
timelock.fsample = 600; 

timelock_denoised = ft_timelockanalysis(cfg, data_denoised);
timelock_denoised.fsample = 600; 



output_original = timelock.avg([abs(params.chnsR_num), abs(params.chnsL_num)],:).*[params.sign_R params.sign_L]';
if pre_mode == 4
    output_denoised = timelock_denoised.avg;
else
    output_denoised = timelock_denoised.avg([abs(params.chnsR_num), abs(params.chnsL_num)],:).*[params.sign_R params.sign_L]';  
end

time = timelock.time;
plot(time, mean(output_original),'Linewidth',3);
hold on; 
plot(time, mean(output_denoised), 'Linewidth',3);
% legend({'Original',sprintf('Wavelet (%s - lev %d)',wname,lev)});
% title(sprintf('Wavelet (%s - lev %d). Subj: %d',wname, lev, filenum));

figure(2)
aux_ori  = mean(output_original); aux_ori = (aux_ori - min(aux_ori))/(max(aux_ori)-min(aux_ori));
aux_den  = mean(output_denoised); aux_den = (aux_den - min(aux_den))/(max(aux_den)-min(aux_den));
time = timelock.time;
plot( aux_ori,'Linewidth',3); hold on;
plot( aux_den, 'Linewidth',3); hold on;
plot( aux_den(2*log2(sample_size):end),'Linewidth',3)
legend('Original','Denoised','Denoised shift-corrected')
% legend({'Original',sprintf('Wavelet (%s - lev %d)',wname,lev)});
