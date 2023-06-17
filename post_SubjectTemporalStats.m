function post_SubjectTemporalStats(data, subject_list)
%     subject_list = 1:length(subject);
%     data = delay_vec(ind_ker,t_stable,subject_list,:);
    data_mean       = squeeze(mean(data, 2));
    data_std   = squeeze(std(data, [], 2));

    % We get the distribution of the delays and maximum values from all the
    % subjects in order to estimate their distributions.
    data = squeeze(data);
    [n_trig, n_subj, n_mod] =size(data);
    for ind_modality = 1:n_mod
        data_long(:, ind_modality) = reshape(data(:,:,ind_modality), n_trig*n_subj, 1);
        if mean(data_long(:, ind_modality)) < 10
            points{ind_modality} = min(data_long(:, ind_modality))-10:1/4/mean(data_long(:, ind_modality)): max(data_long(:, ind_modality))+10;
        else
            points{ind_modality} = min(data_long(:, ind_modality)) - 10: max(data_long(:, ind_modality)) + 10;
        end
        pdf_estimation{ind_modality} = ksdensity(data_long(:, ind_modality), points{ind_modality});
        
    end
   
    % We extract the minimum and maximum values from data.
    min_y = max(min(data_mean(:) - data_std(:)) - 10, 0); 
    max_y = max(data_mean(:) + data_std(:))*1.2;
    
    %% Plots
    transp_level = .1;    
    for ind_modality = 1:n_mod
        subplot(1,6,1);     
        switch ind_modality
            case 1
                color = [0 0 1];
            case 2
                color = [1 0 0];
        end
        fill([-pdf_estimation{ind_modality}, zeros(size(pdf_estimation{ind_modality}))], ...
                                    [points{ind_modality},fliplr(points{ind_modality})],...
                                    color,'FaceAlpha',transp_level);
        ylim([min_y, max_y])
        set(gca,'xticklabel',{[]}) 
        hold on;
        
        subplot(1,6,2:6); 
        plot(data_mean(:, ind_modality),'Color',color, 'Linewidth',3);
        grid on
        hold on
        inBetween = [(data_mean(:, ind_modality) + data_std(:, ind_modality))', fliplr((data_mean(:, ind_modality) - data_std(:, ind_modality))')];
        fill([1:length(data_mean(:, ind_modality)), fliplr(1:length(data_mean(:, ind_modality)))]', inBetween,color, 'FaceAlpha',transp_level, 'LineStyle','none');
        ylim([min_y, max_y])
    end
%     title(sprintf('DELAYS. Kernel: %d',ind_ker))
%     xlabel('Subjects')
%     ylim([min_y, max_y])

end