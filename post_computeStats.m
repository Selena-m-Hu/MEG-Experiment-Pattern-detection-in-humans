function s_out = post_computeStats(Diff, perc_vector, trigger_list)
    % This function computes bootstrapping from two different conditions
    % using the percentiles defined in 'perc_vector'. The variable
    % 'trigger_list' is necessary to prune the outcoming differences.
    % 
    % Diff: variable that keeps the difference between conditions.
    % perc_vector: vector of p-values. For example: [0.05, 0.01].
    % trigger_list: couples [5,15] or [10, 20].
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
    % Last update: 06/August/2018



    for ind_perc = 1:length(perc_vector)
        dataB=bootstrap(Diff'); 
        s=findSigDiff(dataB, perc_vector(ind_perc));
        % S pruning
        s_out(:, ind_perc) = s; %proc_DiffPruning(s, trigger_list);
    end
end