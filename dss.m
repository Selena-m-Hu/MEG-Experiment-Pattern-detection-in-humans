function data_out = dss(data, selection)
    % we concatenate all the trials and select a part of them. This subset will
    % be used to compute the rotation and bias matrixes.
    

    n_trials = length(data.trial);
    [n_channels, t_trial] = size(data.trial{1});
    if nargin == 1
        selection = 1:n_channels;
    end
    rng(1492);
    ind_trials = randperm(n_trials,floor(n_trials*.2)); % 20% of the trials.
    X = [];
    L = [];
    for ind = ind_trials
        X = [X data.trial{ind}(selection,:)];
        L = [L diag(ones(1,t_trial))];
    end
    media = mean(X,2);
    X = X - media;

    [P,~,latent] = pca(X');
    % We select a percentaje of the components.
    n_features_th = .9;
    n_features = find(diff(cumsum(latent)/sum(latent) >= n_features_th) == 1);
    
    % We rotate data using the PCA matrix. We focus on some some chosen
    % particular components, the ones with higher variance.
    % We set a minimum number of three components.
    P = P(:,1:max(n_features,1));
    XP= X'*P;

    % We normalize data using a diagonal eigenvalues matrix.
    N = diag(1./sqrt(latent(1:n_features)));
    XPN = XP*N;
    % % Security check.
    % plot(std(XPN)) % It should be close to 1.

    % We compute the bias function.
    LXPN = L*XPN;

    % We obtain the PCA matrix from the normalized data, which gives us the final components.
    [Q, ~, latent_Q] = pca(LXPN);
    % We select a percentaje of the components.
    n_features_th = .9;
    n_features = find(diff(cumsum(latent_Q)/sum(latent_Q) >= n_features_th) == 1);
    Q = Q(:,1:max(1,n_features));
    % n_features = n_channels;
    % % Security check
    % plot(cumsum(latent_Q)/sum(latent_Q))

    % Output matrix.
    W = P*N*Q;

    %% Security check
    % Y = (L*X'*W)';
    % 
    % 
    % chunk_Y = Y(:,1:420);
    % chunk_X = X(:,1:420);
    % 
    % chunk_Y = chunk_Y/ max(chun_  k_Y(:));
    % chunk_X = chunk_X/max(chunk_X(:));
    % 
    % plot(mean(chunk_X,1));
    % hold on;
    % plot(mean(chunk_Y,1));
    % legend('X','Y (transformed)');


    %% Data transformation
    data_out.avg = zeros(n_features,t_trial);
    for ind = 1:n_trials
        X = data.trial{ind}(selection,:) - media;
        data_out.trial{ind} = (diag(ones(1,t_trial))*X'*W)';
        data_out.avg = data_out.avg + data_out.trial{ind};
    end
    data_out.avg = data_out.avg/n_trials;
    data_out.t_trial = t_trial;
    data_out.n_trials = n_trials;

end
