    %this script analyzes the M100 response in the localizer script and finds
%the 20 strongest channels in each hemisphere.
function Partial(dataset, filenum, block)
% clear all;
% close all;
addpath(fullfile('..','fieldtrip'));
if nargin < 3
    block = [];
end

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
cfg.channel = 'MEG'; % Le decimos que queremos capturar todos los canales de MEG.
% filenum = 4;
switch dataset
    case 1
        signature = 'Analysis_Wave4'; % 'params'
    case 2
        signature = 'Dataset2-Long_first';
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

data = ft_preprocessing(cfg);
plot((1:length(data.trial{1}))/600,data.trial{1}(250:255,:)')



hdr   = ft_read_header(cfg.dataset);
event = ft_read_event(cfg.dataset);

%%
if filenum > 0
    
%     cfg.trialdef.eventtype  = '?'; % Aquí indicamos el itpo de evento que nos ayudará a separar en epochs.
    cfg.trialdef.eventtype  = 'UPPT001'; % Aquí indicamos el itpo de evento que nos ayudará a separar en epochs.
    cfg.trialdef.eventvalue = 10;

    cfg.trialdef.prestim    = 1; % Espacio de preestímulo (antes del estímulo). Puede usarse para extraer la media.
    cfg.trialdef.poststim   = 15; % Espacio de señal que nos interesa, después del evento.

    cfg = ft_definetrial(cfg);
    cfg.demean      = 'yes';

    cfg.baselinewindow = [-1 0];% in seconds
    data = ft_preprocessing(cfg);
    cfg=[];
%     cfg.hpfilter ='yes';
%     cfg.hpfreq = 0.5;
    cfg.lpfilter ='yes';
    cfg.lpfreq = 30;

    data=ft_preprocessing(cfg, data);

    % save data
    cfg=[];
    cfg.method = 'summary';
    cfg.channel = 'MEG';
    % data = ft_rejectvisual(cfg,data);

end


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

%finding the 20 strongest channels in each hemisphere

M100dat=timelock.avg'; % Obsérvese que aquí se transpone la media, por lo que tenemos (tiempo x canal).

t0 = 0.2*timelock.fsample;
% Ordenamos las amplitudes, para lo cual escogemos un rango de instantes de
% tiempo, y todos los canales.
amps=mean(M100dat((t0+0.09*timelock.fsample):(t0+0.11*timelock.fsample), :),1);
[ampsSorted,idx]=sort(amps,2,'descend');

chnsSorted = timelock.label(idx);

%selecting channels:ft_timelockanalysis

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
% chns_selectedL = [chns_selectedLpos];
chns_selectedR = [chns_selectedRpos chns_selectedRneg];
% chns_selectedR = [chns_selectedRpos];

chnsL_num=[];
signo = +1;

for count1=1:length(timelock.label)
    for count2=1:length(chns_selectedL)
        if (strcmp(timelock.label{count1},chns_selectedL{count2}) ~= 0)
            if (count2 > 10) & (length(chns_selectedL) == 20)
                signo = -1;
            end
            chnsL_num=[chnsL_num signo*count1];
            signo = +1;
%             chnsL_num=[chnsL_num count1];
        end
    end
end

chnsR_num=[];
signo = +1;
for count1=1:length(timelock.label)
    for count2=1:length(chns_selectedR)
        if (strcmp(timelock.label{count1},chns_selectedR{count2}) ~= 0)
            if (count2 > 10) & (length(chns_selectedR) == 20)
                signo = -1;
            end
            chnsR_num=[chnsR_num signo*count1];
            signo = +1;
%             chnsR_num=[chnsR_num count1];
        end
    end
end

chns_selectedL = timelock.label(abs(chnsL_num)); 
chns_selectedR = timelock.label(abs(chnsR_num));
sign_L = sign(chnsL_num);
sign_R = sign(chnsR_num);
% We store the channels and its signs.
params.chns_selectedL = chns_selectedL;
params.chns_selectedR = chns_selectedR;
params.sign_L = sign_L;
params.sign_R = sign_R;
params.chnsL_num = chnsL_num;
params.chnsR_num = chnsR_num;
% mkdir(fullfile('..','Results', 'params', username));
% save (fullfile('..','Results', 'params', username,'Channel_selection.mat'), 'params');
selection = [abs(chnsR_num), abs(chnsL_num)];
selection_sign = [sign_R, sign_L];

selected_dataL = timelock.avg(abs(chnsL_num),:);
selected_dataR = timelock.avg(abs(chnsR_num),:);

RMS_L = rms(selected_dataL,1);
RMS_R = rms(selected_dataR,1);
RMS_results = rms([selected_dataL; selected_dataR],1);
figure;
plot(timelock.time,    RMS_results)
% keyboard
%%
% denoiser(data, [-2], selection.*selection_sign)
% data = denoiser(data, [1,-1,2,-2,3,4,5], selection.*selection_sign,1);
% savefig(fullfile('..','Results', signature, username,'Denoiser_results.fig'))
% data = denoiser(data, [1,2,3,5], selection.*selection_sign,1);
% savefig(fullfile('..','Results', signature, username,'Denoiser_selection.fig'))
%%
% %% Rossy DSS
% output_original = timelock.avg([abs(chnsR_num), abs(chnsL_num)],:).*[sign_R sign_L]';
% addpath(genpath(fullfile('..','..','Antonio')));
% 
% % Rossy's code
% x = cat(3,data.trial{:}); % x is simply a matrix of the data from all trials
% % DSS wants it in time * channels * trials
% x = permute(x,[2 1 3]);
% c1 = zeros(size(x,2));
% c0 = nt_cov(x);
% c1 = nt_cov(mean(x,3));
% t = data.time{1};
% % baseline/demean
% x = nt_demean2(x,find(t<0)); % use nt_demean() instead if there's a strong trend
% % 2. prepare for DSS for evoked repeatability
% % calculate covariances
% c0 = zeros(size(x,2));
% c1 = zeros(size(x,2));
% c0 = nt_cov(x);
% c1 = nt_cov(mean(x,3));
% 
% % DSS
% [todss,pwr0,pwr1] = nt_dss0(c0,c1);
% 
% % 2. do DSS
% z = nt_mmat(x,todss); % z is components timeseries
% % f1 = figure(1); clf; plot(pwr1./pwr0,'.-'); ylabel('score'); xlabel('component');
% % Elegimos las componentes adecuadas, que sumen un porcentaje importante de
% % la energía.
% ind_comp = 5;%min(find((cumsum(pwr1./pwr0)/sum(pwr1./pwr0)) > .99));
% out_dss_rossy = mean(z(:,1:ind_comp,:),3)';
% 
% rmpath(genpath(fullfile('..','..','Antonio')));
% 
% %% Wavelet
% % output_original = timelock.avg([abs(chnsR_num), abs(chnsL_num)],:).*[sign_R sign_L]';
% % addpath(genpath(fullfile('wave')));
% % qmf = MakeONFilter('Coiflet',4);
% % lev = 1;
% % ind = 39;
% % y = output_original(ind,:);
% % % We introduce padding with the initial and last values of the vector.
% % padding = 2^ceil(log2(length(y)))-length(y);
% % init_padding = y(1)*ones(1,floor(padding/2)); end_padding = y(end)*ones(1,floor(padding/2)+mod(padding,2));
% % y_long = [init_padding, y, end_padding];
% % 
% % ydwt=FWT_PO(y_long,lev,qmf);
% % % yr_long=ebayesthresh_wavelet(ydwt, qmf);
% % yr_long=ebayesthresh_wavelet(ydwt,qmf,'independent',ceil(log2(length(y_long)))-2,'cauchy',1);
% % yr = yr_long(length(init_padding)+1:end-length(end_padding));
% % plot(y); hold on; plot((yr))
% 
% 
% %% Wavelet denoise
% 
% output_original = timelock.avg([abs(chnsR_num), abs(chnsL_num)],:).*[sign_R sign_L]';
% wname = 'coif4'; % db4, coif4, sym8, fk14, bior3.7, rbio1.5 (desfasa)
% lev = 4;
% % We denoise using wavelets.
% out_wavelet = zeros(size(output_original));
% for ind = 1:(length(chnsR_num)+length(chnsL_num))
% %     [out_wavelet(ind, :),~,~] = wden(output_original(ind,:),'sqtwolog','s','mln',lev,wname);
%     [out_wavelet(ind, :),~,~] = wdenoise(output_original(ind,:),lev,'Wavelet',wname);
%         
% end
% 
% % 
% % % We store the results.
% % savefig(fullfile('..','Results', signature, username,'An_GLOB.fig'))
% % close all
% 
% 
% %% DSS denoise
% data_out = dss(data);
% out_dss = data_out.avg(1,:);
% 
% 
% 
% %% ICA
% ica_on = 0;
% if ica_on == 1
%     ica_out = ica_denoise(data);
%     data_ica = data;
%     data_ica.trial = ica_out.trial;
%     
%     % Aquí se extraen los canales, considerando que pueden tener magnitudes
%     % opuestas: esto es, el mismo gráfico pero con el signo cambiado a negativo
%     % o positivo. En cualquier caso, el resultado promedio debería ser el
%     % adecuado.
%     cfg = [];
%     cfg.channel = 'MEG';
% 
%     % Extraemos algunos estadísticos de los trials. Por ejemplo, el promedio de
%     % la magnitud.
%     timelock_ica = ft_timelockanalysis(cfg, data_ica);
%     timelock_ica.fsample = 600;
%     selected_data_ICA = timelock_ica.avg([abs(chnsR_num), abs(chnsL_num)],:).*[sign_R sign_L]';
%     
%     time = timelock_ica.time;
%     
%     figure('units','normalized','outerposition',[0 0 1 1])
%     plot(time, mean(selected_data_ICA),'Linewidth',3);
% 
% 
% 
%     glob2 = mean(selected_data_ICA); glob_mean2 = (glob2-min(glob2))/(max(glob2)-min(glob2));
%     figure('units','normalized','outerposition',[0 0 1 1])
%     plot(time, glob_mean,'Linewidth',3);
%     hold on
%     plot(time, glob_mean2,'Linewidth',3);
%     hold on; 
%     plot(time, wave_mean, 'Linewidth',3);
%     % hold on;
%     % plot(time, dss_mean, 'Linewidth',3);
%     % hold on;
%     % plot(time, -dss_mean+1, 'Linewidth',3);
%     title(sprintf('Normalized curves. Subject: %d', filenum));
%     legend('Original (mean)','ICA+wavelet',sprintf('Wavelet (%s)',wname),'DSS (pos)','DSS (neg)');
%     savefig(fullfile('..','Results', signature, username,'An_Comparison_DSS-wave-ICA.fig'))
% end
% 
% 
% 
% %% EMD
% % data_emd = EMD_denoise(data);
% data_denoised = data;
% data_out = data.trial;
% for ind = 1:length(data_denoised.trial)
%     data_denoised.trial = {data.trial{ind}(selection,:)};
%     [data_emd] = EMD_denoise(data_denoised);
% %     clear data_denoised
% 
% 
%     data_out{ind} = data_emd.trial{1};
% end
% data_denoised.trial = data_out;
% % We tune the previous label and channel information to suit our current
% % configuration (reduced number of channels).
% data_denoised.cfg.channel = data_denoised.cfg.channel(selection);
% data_denoised.label = data_denoised.label(selection);
% 
% cfg = [];
% cfg.channel = 'MEG';
% 
% % Extraemos algunos estadísticos de los trials. Por ejemplo, el promedio de
% % la magnitud.
% timelock_emd = ft_timelockanalysis(cfg, data_denoised);
% clear data_out
% timelock_emd.fsample = 600;
% 
% 
% 
% %% Adaptative Filter (LMS)
% %    sample_size = 100;
% %                 % With the Length of the filter we can control the
% %                 % smoothness of the output signal. The bigger the length,
% %                 % the smoother the output.
% %                 lmsfilt2 = dsp.LMSFilter('Length',ceil(2*log2(sample_size)),'Method','LMS', ...
% %                     'StepSize',0.25);        
% %                 
% %              
% %                 for ind_iter = 1:sample_size:length(data.trial{1})
% %                     x(ind,ind_iter:ind_iter+sample_size-1) = (data.trial{1}(selection(ind),ind_iter:ind_iter+sample_size-1)*sign_selection(ind))';
% %                     [y(ind,ind_iter:ind_iter+sample_size-1),e(ind_iter:ind_iter+sample_size-1),w] = lmsfilt2(x(ind,ind_iter:ind_iter+sample_size-1)',...
% %                                        out_sig(ind_iter:ind_iter+sample_size-1)');  
% %                 end
% % keyboard
% 
% clear y error w x
% n_trial = length(data.trial); 
% % RLS: Hasta sample_size=29 se cumple que Y = floor((X-1)/2);
% % LMS: Hasta sample_size=12 se cumple que Y = floor((X-1)/2);  
% n_samples = 15;
% sample_size_rls = n_samples;%min(29,n_samples);   
% est_delay_rls = [];
% est_delay_rls = floor((sample_size_rls-1)/2);
% 
% % n_samples = 12;
% sample_size_lms = n_samples;%min(29,n_samples);   
% est_delay_lms = [];
% est_delay_lms = floor((sample_size_lms-1)/2);       
%                   
%     
%     
% lmsfilt2 = dsp.LMSFilter('Length',sample_size_lms,...%ceil(10*log2(sample_size)),...
%                          'Method','LMS', ...
%                          'StepSize',0.25);       
%                              
% RLSfilt2 = dsp.RLSFilter('Length', sample_size_rls);
% 
% clear x y_lms y_rls w_lms
% for ind_channel = 1:length(selection)     
%     contador = 1;
%     for ind_trial = 1:2:n_trial             
%     
%     
%         % Noisy signal
%         x(ind_channel,:,contador) =     selection_sign(ind_channel)*data.trial{ind_trial}(selection(ind_channel),:)';
%         % Expected signal
%         if (ind_trial+1) > n_trial
%             e = selection_sign(ind_channel)*data.trial{1}(selection(ind_channel),:)';
%         else
%             e = selection_sign(ind_channel)*data.trial{ind_trial+1}(selection(ind_channel),:)';
%         end
%         
%         [y_lms(ind_channel,:, contador), error_lms(ind_channel,:, ind_trial), w_lms(ind_channel,:, ind_trial)] = lmsfilt2(x(ind_channel,:,contador)', e);
%         [y_rls(ind_channel,:, contador), error_rls(ind_channel,:, ind_trial)] = RLSfilt2(x(ind_channel,:,contador)', e);
%         
%         contador = contador+1;
%     end
%     
% end
% % We correct the delay that appears after filtering.
% 
% % Instead of estimating the delay, we compute it using cross-correlation
% % from all the channels and trials.
% % Cross-correlation delay
% contador_corr = 1;
% for ind_channel = 1:length(selection)     
%     for ind_trial = 1:contador-1      
%         correl = xcorr(mean(y_lms(ind_channel,:,:),3), mean(x(ind_channel,:,:),3));
%         [~,max_correl_lms(contador_corr) ]= max(correl);
%         
%         
%         correl = xcorr(mean(y_rls(ind_channel,:,:),3), mean(x(ind_channel,:,:),3));
%         [~,max_correl_rls(contador_corr) ]= max(correl);
%         
%         contador_corr = contador_corr + 1;
%     end
% end
% if isempty(est_delay_lms) == 1
%     est_delay_lms = round(mean(max_correl_lms - size(y_lms,2)));
% end
% if isempty(est_delay_rls) == 1
%     est_delay_rls = round(mean(max_correl_rls - size(y_rls,2)));
% end
% fprintf('LMS delay (avg: %.2f, std: %.2f)\n', mean(est_delay_lms), std(est_delay_lms));
% 
% fprintf('RLS delay (avg: %.2f, std: %.2f)\n', mean(est_delay_rls), std(est_delay_rls));
% 
% y_corrected_lms = zeros(size(y_lms));
% y_corrected_lms(:, 1:end-est_delay_lms,:) = y_lms(:,est_delay_lms+1:end,:);
% clear('aux_af_lms'); aux_af_lms  =  mean(mean(y_corrected_lms, 3),1); aux_af_lms = (aux_af_lms - min(aux_af_lms))/(max(aux_af_lms)-min(aux_af_lms));
% 
% y_corrected_rls = zeros(size(y_rls));
% y_corrected_rls(:, 1:end-est_delay_rls,:) = y_rls(:,est_delay_rls+1:end,:);
% clear('aux_af_rls'); aux_af_rls  =  mean(mean(y_corrected_rls, 3),1); aux_af_rls = (aux_af_rls - min(aux_af_rls))/(max(aux_af_rls)-min(aux_af_rls));
% % 
% % % Adaptative filters
% % time = timelock.time;
% % aux_orig = mean(output_original); aux_orig = (aux_orig - min(aux_orig))/(max(aux_orig)-min(aux_orig));
% % time = timelock.time;
% % subplot(211)
% % plot(time, (aux_orig),'Linewidth',3);
% % hold on
% % plot(time, (aux_af_lms),'Linewidth',3);
% % legend('Original','LMS Filter (shift-corrected)');
% % title(sprintf('Adaptative filtering (LMS). Subj: %d',filenum));
% % 
% % subplot(212)
% % plot(time, (aux_orig),'Linewidth',3);
% % hold on
% % plot(time, (aux_af_rls),'Linewidth',3);
% % legend('Original (normalized)','RLS Filter (shift-corrected)');
% % title(sprintf('Adaptative filtering (RLS). Subj: %d',filenum));
% % fprintf('LMS %.2f - RLS %.2f\n',norm(aux_orig - aux_af_lms),norm(aux_orig - aux_af_rls))
% % % savefig(fullfile('..','Results', signature, username,'An_Comparison_Global.fig'))
% % % close all
%  
% %% Graphs
% % Wavelet vs original
% time = timelock.time;
% figure('units','normalized','outerposition',[0 0 1 1])
% plot(time, mean(output_original),'Linewidth',3);
% hold on; 
% plot(time, mean(out_wavelet), 'Linewidth',3);
% legend({'Original',sprintf('Wavelet (%s)',wname)});
% 
% % DSS vs origina
% figure('units','normalized','outerposition',[0 0 1 1])
% % sign_change = 1;%sign(trapz(mean(output_original)) /trapz(out_dss));
% aux_orig = mean(output_original); aux_orig = (aux_orig - min(aux_orig))/(max(aux_orig)-min(aux_orig));
% aux_out  = out_dss; aux_out = (aux_out - min(aux_out))/(max(aux_out)-min(aux_out));
% [sign_change] = sign_criteria(aux_orig, aux_out);
% if sign_change < 0
%     aux_out = (aux_out*sign_change) + 1;
% end
% 
% subplot(211);
% plot(time, mean(output_original), 'Linewidth',3);
% title('Original');
% subplot(212)
% plot(time,out_dss*sign_change,'r', 'Linewidth',3)
% title('DSS (Denoising Source Separation)');
% 
% % EMD vs original
% selected_data_EMD = timelock_emd.avg(:,:).*[sign_R sign_L]';
% figure('units','normalized','outerposition',[0 0 1 1])
% plot(time, mean(output_original),'Linewidth',3);
% hold on
% plot(time, mean(selected_data_EMD),'Linewidth',3);
% % title(sprintf('Normalized curves. Subject: %d', filenum));
% legend('Original (mean)','EMD');
% % savefig(fullfile('..','Results', signature, username,'An_Comparison_wave-EMD.fig'))
% 
% % Original vs wavelet vs EMD vs DSS
% close all
% figure('units','normalized','outerposition',[0 0 1 1]);
% subplot(311);
% plot(time, mean(output_original),'Linewidth',3);
% hold on; 
% plot(time, mean(out_wavelet), 'Linewidth',3);
% legend({'Original',sprintf('Wavelet (%s - lev %d)',wname,lev)});
% title(sprintf('Wavelet (%s - lev %d). Subj: %d',wname, lev, filenum));
% 
% subplot(312);
% 
% plot(time, aux_orig, 'Linewidth',3); 
% hold on;
% plot(time,aux_out, 'Linewidth',3)
% title(sprintf('DSS (Deno    ising Source Separation). Subj: %d',filenum));
% legend('Original (normalized)','DSS (normalized)');
% % clear aux_out aux_orig
% 
% subplot(313)
% plot(time, mean(output_original),'Linewidth',3);
% hold on
% plot(time, mean(selected_data_EMD),'Linewidth',3);
% legend('Original','EMD');
% title(sprintf('EMD (Empirical Mode Decomposition). Subj: %d',filenum));
% savefig(fullfile('..','Results', signature, username,'An_Comparison_Global.fig'))
% close all
% 
% %% Adaptative filters
% aux_orig = mean(output_original); aux_orig = (aux_orig - min(aux_orig))/(max(aux_orig)-min(aux_orig));
% time = timelock.time;
% subplot(211)
% plot(time, (aux_orig),'Linewidth',3);
% hold on
% sign_lms = sign_criteria(aux_orig, aux_af_lms);
% if sign_lms == -1
%     aux_af_lms = (-1*aux_af_lms)+1;
% end
% plot(time, (aux_af_lms),'Linewidth',3);
% legend('Original','LMS Filter (shift-corrected)');
% title(sprintf('Adaptative filtering (LMS). Subj: %d',filenum));
% 
% subplot(212)
% plot(time, (aux_orig),'Linewidth',3);
% hold on
% sign_rls = sign_criteria(aux_orig, aux_af_rls);
% if sign_rls == -1
%     aux_af_rls = (-1*aux_af_rls)+1;
% end
% plot(time, (aux_af_rls),'Linewidth',3);
% legend('Original (normalized)','RLS Filter (shift-corrected)');
% title(sprintf('Adaptative filtering (RLS). Subj: %d',filenum));
% savefig(fullfile('..','Results', signature, username,'An_Adaptative_filter.fig'))
% close all
% 
% 
% %% DSS - Rossy
% 
% clear('aux_orig'); aux_orig = mean(output_original); aux_orig = (aux_orig - min(aux_orig))/(max(aux_orig)-min(aux_orig));
% clear('aux_out_rossy'); aux_out_rossy  = mean(out_dss_rossy); aux_out_rossy = (aux_out_rossy - min(aux_out_rossy))/(max(aux_out_rossy)-min(aux_out_rossy));
% [sign_change] = sign_criteria(aux_orig, aux_out_rossy);
% if sign_change < 0
%     aux_out_rossy = (aux_out_rossy*sign_change) + 1;
% end
% 
% clear('aux_out'); aux_out  = out_dss; aux_out = (aux_out - min(aux_out))/(max(aux_out)-min(aux_out));
% [sign_change] = sign_criteria(aux_orig, aux_out);
% if sign_change < 0
%     aux_out = (aux_out*sign_change) + 1;
% end
% 
% % figure('units','normalized','outerposition',[0 0 1 1])
% % plot(time, aux_orig, 'Linewidth',3); 
% % clear('aux_orig'); aux_orig = mean(output_original); aux_orig = (aux_orig - min(aux_orig))/(max(aux_orig)-min(aux_orig));
% % clear('aux_out'); aux_out_rossy  = mean(out_dss_rossy); aux_out_rossy = (aux_out_rossy - min(aux_out_rossy))/(max(aux_out_rossy)-min(aux_out_rossy));
% % clear('aux_out'); aux_out  = out_dss; aux_out = (aux_out - min(aux_out))/(max(aux_out)-min(aux_out));
% figure('units','normalized','outerposition',[0 0 1 1])
% plot(time, aux_orig, 'Linewidth',3); 
% hold on;
% plot(time,aux_out, 'Linewidth',3)
% hold on;
% plot(time,aux_out_rossy, 'Linewidth',3)
% title(sprintf('DSS (Denoising Source Separation). Subj: %d',filenum));
% legend('Original (normalized)','DSS (normalized)','DSS - Rossy (normalized)');
% savefig(fullfile('..','Results', signature, username,'An_Comparison_DSS.fig'))
% % clear aux_out aux_orig
% % close all
% % 
% % hold on;
% % plot(time,aux_out, 'Linewidth',3)
% % hold on;
% % plot(time,aux_out_rossy, 'Linewidth',3)
% % title(sprintf('DSS (Denoising Source Separation). Subj: %d',filenum));
% % legend('Original (normalized)','DSS (normalized)','DSS - Rossy (normalized)');
% % savefig(fullfile('..','Results', signature, username,'An_Comparison_DSS.fig'))
% % clear aux_out aux_orig
% % close all
% 
% fprintf('*********** FINISHED ************');
% close all
%     