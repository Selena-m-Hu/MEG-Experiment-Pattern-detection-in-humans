function output = ica_denoise(data)

    addpath(fullfile('..','FastICA_25'));
    
    for ind_trial = 1:length(data.trial)
        norm_data = data.trial{ind_trial}/max(data.trial{ind_trial}(:));
        [icasig, A, W] = fastica (norm_data); % A: mixing matrix  -  W: Reconstruction matrix.
%         [icasig, A, W] = fastica (norm_data,'maxNumIterations',50); % A: mixing matrix  -  W: Reconstruction matrix.


        wname = 'coif4'; % db4, coif4, sym8, fk14, bior3.7, rbio1.5 (desfasa)
        lev = 4;
        % We denoise using wavelets.
        aux_ICA = zeros(size(icasig));
        for ind = 1:size(icasig,1)
            [aux_ICA(ind, :),~,~] = wden(icasig(ind,:),'sqtwolog','s','mln',lev,wname);
        end
        % plot(mean(A*aux_ICA))
        output.trial{ind_trial} = A*aux_ICA*max(data.trial{ind_trial}(:));
    end
    % imagesc(output - norm_data);colorbar
    % W*norm_data-icasig

end