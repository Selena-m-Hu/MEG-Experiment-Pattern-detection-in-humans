function [data] = denoiser(data, denoise_selection, channel_selection, verbose)
% Function designed to denoise the channels contained in a Fieltrip data
% structure. It allows to choose among different denoising techniques,
% which can be processed in a single call. It also allows to choose certain
% channels of the dataset if prior knowledge is available.
%
% - data : Data structure provided by Fieldtrip.
% - denoise_selection : Integer vector indicating the chosen denoising
%   techniques.
% - channel_selection : Channels that were selected previously for the
%   analysis. If none introduced, the whole channel set is processed 
%   (slow AF). Each channel can have a positive or negative sign, which
%   indicates if the magnitude of the channel should be changed (useful for
%   the Adaptative Filter).
% 
% %% Antonio Rodríguez Hidalgo - 1 May 2018.

%%

    if isempty(denoise_selection) % No denoise selection. All of them will be computed.
        % TBD
    end
    if isempty(channel_selection) % No channel selection
        % TBD
        channel_selection = 1:size(data.trial{1},1);
    end
    % First, we check the number of trials. If it's one, data is not
    % chunked.
    n_trials = length(data.trial);
    x = cat(3,data.trial{:}); % All the data in a single matrix.
    selection_sign = sign(channel_selection);
    selection_num = abs(channel_selection);
    for ind_denoise = denoise_selection
        switch ind_denoise
            case 1 % Wavelet analysis - wdenoise(Matlab2018)                
                % It produces several output channels.
                % In this case, we average throught trials to get a more
                % stable signal.
                x_mean = mean(x,3);
                wname = 'coif4'; % db4, coif4, sym8, fk14, bior3.7, rbio1.5 (desfasa)
                lev = 4;
                
%                 out_wavelet = zeros(size(x_mean));
                for ind = 1:length(channel_selection)
                %     [out_wavelet(ind, :),~,~] = wden(output_original(ind,:),'sqtwolog','s','mln',lev,wname);
                    [out_wavelet((ind), :),~,~] = wdenoise(x_mean(selection_num(ind),:),lev,'Wavelet',wname);

                end
                data.denoised.wavelet_fine = selection_sign'.*out_wavelet;
                clear out_wavelet x_mean lev wname
                
                
            case -1 % Wavelet analysis - wden(Matlab2017 and older with wavelet toolbox)                
                % In this case, we average throught trials to get a more
                % stable signal.
                % It produces several output channels.
                x_mean = mean(x,3);
                wname = 'coif4'; % db4, coif4, sym8, fk14, bior3.7, rbio1.5 (desfasa)
                lev = 4;
               
%                 out_wavelet = zeros(size(x_mean));
                for ind = 1:length(channel_selection)
                    [out_wavelet(ind, :),~,~] = wden(x_mean(selection_num(ind),:),'sqtwolog','s','mln',lev,wname);
%                     [out_wavelet(channel_selection(ind), :),~,~] = wdenoise(x_mean(channel_selection(ind),:),lev,'Wavelet',wname);

                end
                data.denoised.wavelet_rough = selection_sign'.*out_wavelet;
                clear out_wavelet x_mean lev wname
                
                
            case 2 % Denoising Source Separation (DSS) - Cheveigné et al. 2014                       
                % It produces several output channels.
                data_out = dss(data);
                data.denoised.dss = data_out.avg(1,:);
                
                clear data_out
                
                
            case -2 % Denoising Source Separation (DSS) - NoiseTools Toolbox. Code provided by Rossy.
                addpath(genpath(fullfile('..','..','Antonio')));

                % Rossy's code
               
                % DSS wants it in time * channels * trials
                x_prime = permute(x,[2 1 3]);
%                 c1 = zeros(size(x_prime,2));
%                 c0 = nt_cov(x_prime);
%                 c1 = nt_cov(mean(x_prime,3));
                t = data.time{1};
                % baseline/demean
                x_prime = nt_demean2(x_prime,find(t<0)); % use nt_demean() instead if there's a strong trend
                % 2. prepare for DSS for evoked repeatability
                % calculate covariances
                c0 = zeros(size(x_prime,2));
                c1 = zeros(size(x_prime,2));
                c0 = nt_cov(x_prime);
                c1 = nt_cov(mean(x_prime,3));

                % DSS
                [todss,pwr0,pwr1] = nt_dss0(c0,c1);

                % 2. do DSS
                z = nt_mmat(x_prime,todss); % z is components timeseries
                % f1 = figure(1); clf; plot(pwr1./pwr0,'.-'); ylabel('score'); xlabel('component');
                % Elegimos las componentes adecuadas, que sumen un porcentaje importante de
                % la energía.
                ind_comp = 5;%min(find((cumsum(pwr1./pwr0)/sum(pwr1./pwr0)) > .99));
                data.denoised.dss_noisetools = mean(z(:,1:ind_comp,:),3)';

                rmpath(genpath(fullfile('..','..','Antonio')));
                clear ind_comp z todss pwr0 pwr1 c0 c1 x_prime t
                
                
            case 3 % Empirical Mode Decomposition (EMD).                
                data_denoised = data;
                data_out = data.trial;
                for ind = 1:n_trials
                    data_denoised.trial = {data.trial{ind}(selection_num,:)};
                    [data_emd] = EMD_denoise(data_denoised);
                    if ind == 1
                        data_out = data_emd.trial{1};
                    else
                        data_out = data_out + data_emd.trial{1};
                    end
                end
                data.denoised.emd = (selection_sign'.*data_out)/n_trials;
                clear data_denoised data_out data_emd
                
            case 4 % Adaptative filter (LMS).                 
                % LMS: Hasta sample_size=12 se cumple que Y = floor((X-1)/2);  
                n_samples = 30;
                sample_size_lms = n_samples;  
                est_delay_lms = floor((sample_size_lms-1)/2);                  
%                 x = (x-min(x(:)))/(max(x(:))-min(x(:)));
                
                % Filter object
                lmsfilt2 = dsp.LMSFilter('Length',sample_size_lms,...%ceil(10*log2(sample_size)),...
                         'Method','LMS', ...
                         'StepSize',0.25);           
                % Filtering stage. We choose a channel and filter
                % each consecutive trial. The same filter is used and
                % updated for each successive channel.
                for ind_channel = 1:length(selection_num)     
                    contador = 1;
                    for ind_trial = 1:2:n_trials                 

                        % Expected signal
                        if (ind_trial+1) > n_trials
                            e = selection_sign(ind_channel)*x(selection_num(ind_channel),:,1)';
                        else
                            e = selection_sign(ind_channel)*x(selection_num(ind_channel),:,ind_trial+1)';
                        end

                        [y_lms(ind_channel,:, contador), error_lms(ind_channel,:, ind_trial), w_lms(ind_channel,:, ind_trial)] = lmsfilt2(selection_sign(ind_channel)*x(selection_num(ind_channel),:,contador)',e);
                        contador = contador+1;
                    end

                end

                y_corrected_lms = zeros(size(y_lms));   
                y_corrected_lms(:, 1:end-est_delay_lms,:) = y_lms(:,est_delay_lms+1:end,:);
                
                % We store the resultant channels. They are not related
                % with the channel_selection parameter of the input.
                data.denoised.af_lms = mean(y_corrected_lms,3);                
                clear n_samples sample_size_lms est_delay_lms lmsfilt2 y_lms error_lms w_lms e y_corrected_lms

            case 5 % Adaptative filter (RLS).                 
                % RLS: Hasta sample_size=29 se cumple que Y = floor((X-1)/2);
                n_samples = 15;
                sample_size_rls = n_samples;
                est_delay_rls = floor((sample_size_rls-1)/2);   
                             
                RLSfilt2 = dsp.RLSFilter('Length', sample_size_rls);
                for ind_channel = 1:length(selection_num)     
                    contador = 1;
                    for ind_trial = 1:2:n_trials            
                        if (ind_trial+1) > n_trials
                            e = selection_sign(ind_channel)*x(selection_num(ind_channel),:,1)';
                        else
                            e = selection_sign(ind_channel)*x(selection_num(ind_channel),:,ind_trial+1)';
                        end
        
                        [y_rls(ind_channel,:, contador), error_rls(ind_channel,:, ind_trial)] = RLSfilt2(selection_sign(ind_channel)*x(selection_num(ind_channel),:,contador)', e);
        
                        contador = contador+1;
                    end

                end
                y_corrected_rls = zeros(size(y_rls));
                y_corrected_rls(:, 1:end-est_delay_rls,:) = y_rls(:,est_delay_rls+1:end,:);
                
                % We store the resultant channels. They are not related
                % with the channel_selection parameter of the input.
                data.denoised.af_rls = mean(y_corrected_rls,3);                
                clear y_corrected_rls est_delay_rls y_rls contador error_rls  e RLSfilt2 sample_size_rls
        end
    end
%     
%     
%                 % We tune the previous label and channel information to suit our current
%                 % configuration (reduced number of channels).
%                 data_denoised.cfg.channel = data_denoised.cfg.channel(selection);
%                 data_denoised.label = data_denoised.label(selection);
    data.denoised.x = x;
%     keyboard      
    if verbose == 1 % Use verbose to show graphical representations of the signals.
        %%
        x_mean = mean(x,3);
        x_selec = mean(selection_sign'.*x_mean(selection_num,:));
        x_norm = normalization_fun(x_selec);
        
        n_cols = 2;
        n_rows = ceil(length(denoise_selection)/n_cols);
        campos = fields(data.denoised);
        [campos_out, normalize] = reformat_fields(campos);
        
        datos = struct2cell(data.denoised);
        
        figure('units','normalized','outerposition',[0 0 1 1])
        for ind = 1:length(denoise_selection)
            subplot(n_rows, n_cols, ind)
            if normalize{ind} == 1
                % In this case, we need to normalize both the input data
                % and the output results to the range [0,1]. Additionally,
                % we might need to correct the sign of the output data.
                plot(x_norm, 'Linewidth',3);
                hold on;                
                norm_output = normalization_fun(mean(datos{ind},1));
                sign_output = sign_criteria(x_norm, norm_output);
                if sign_output == -1
                    norm_output = (norm_output*-1) + 1;
                end
                plot(norm_output, 'Linewidth',2)
                legend('Original','Denoised');
                
            else
                plot(x_selec, 'Linewidth',3);                
                hold on;
                plot(mean(datos{ind},1), 'Linewidth',2)
                legend('Original','Denoised');
                
            end
            title(campos_out{ind})            
        end
        
    end

end

function [campos_out, normalize] = reformat_fields(campos)
    for ind = 1:length(campos)
        switch campos{ind}
            case 'wavelet_fine'
                campos_out{ind} = 'Wavelet (Emp. Bayes)';
                normalize{ind} = 0;
            case 'wavelet_rough'
                campos_out{ind} = 'Wavelet (Univ.)';
                normalize{ind} = 0;
            case 'dss'
                campos_out{ind} = 'DSS (Normalized signals)';
                normalize{ind} = 1;
            case 'dss_noisetools'
                campos_out{ind} = 'DSS (noisetools, Normalized signals) ';
                normalize{ind} = 1;
            case 'emd'
                campos_out{ind} = 'EMD';
                normalize{ind} = 0;
            case 'af_lms'
                campos_out{ind} = 'LMS Filter (Normalized signals)';
                normalize{ind} = 1;
            case 'af_rls'
                campos_out{ind} = 'RLS Filter (Normalized signals)';
                normalize{ind} = 1;
            case 'x'
                campos_out{ind} = 'Input';
                normalize{ind} = 0;
        end
    end
    
end

function x_norm = normalization_fun(x)
    x_norm = (x-min(x))/(max(x)-min(x));
end