function [data_emd] = EMD_denoise(data)
    data_emd = data;
    data_emd.trial = {};

    % EMD decomposition is applied to each channel independently, where the
    % first set of IMF components is removed in order to substract the noisy
    % components (considering that the first components are high frequency and
    % so is the noise). The higher the order of the IMF, the lower is supposed 
    % to be its frequency.
    [n_channels, n_time] = size(data.trial{1});
    n_trials = length(data.trial);

    for ind_trial = 1:n_trials
        clear out aux imf residual
        aux = data.trial{ind_trial}; 
        for ind_channel = 1:n_channels
            [imf,residual] = emd(aux(ind_channel,:),'Display',0);
            out(ind_channel,:) = sum(imf(:,2:end),2)+residual;
            
        end

        data_emd.trial{ind_trial} = out;


        verbose = 0;
        if verbose == 1
            figure('units','normalized','outerposition',[0 0 1 1])
            n_subplots = 3+size(imf,2);
            subplot(n_subplots,1,1);
            plot(aux(ind_channel,:)); 
            title('Signal');
            for ind_subplot = 1:size(imf,2)
                subplot(n_subplots,1,ind_subplot+1);
                plot(imf(:,ind_subplot));
                title(sprintf('IMF %d',ind_subplot));
            end
            subplot(n_subplots,1,2+size(imf,2));
            plot(residual);
            title('Residual')        
            subplot(n_subplots,1,n_subplots);
            plot(aux(ind_channel,:),'Linewidth',2);
            hold on;
            plot(sum(imf(:,2:end),2)+residual, 'Linewidth',2);
            legend('Original','Denoised');
            title('Denoised');

        end
    end
end