 function pre_DSScomputation(trigger_list, subject_list, config)
    % Here we compute and project the data into the DSS space keeping all
    % the components. The projection into the channel space is performed in
    % another script (pre_DSSprojection.m).
    %
    % trigger_list: couples (5,15) or (10, 20).
    % * 5 : 3 second RAND sequences.
    % * 10: 15 second RAND sequences.
    % * 15: 3 second REG sequences.
    % * 20: 15 second REG sequences.
    %
    % subject_list:
    % * 2-15, the index of the subject collected by Antonio.
    % * 16-24, the index of the subjects collected by MYH and RB
    %
    % config: allows for certain configurations.
    %   .out_folder: name of the folder where data will be stored.
    %   .store_data: set to 1 to store in a .mat file the whole data of the
    %       subject once that it has been processed.
    %   .channels_path: path indicating where to find the channel
    %       information file.
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
    % Last update: Mingyue Hu, July/2022
    

    out_folder = config.out_folder;
    compute = 1; % If it set to 1, we compute ALWAYS the DSS projection (slow).
    store_data = config.store_data;
    single = config.single;
       
    %% First we check if the transformed data exists. If it does, we do
    % nothing ( we save time).
    if  compute == 1
        
        for subject_ind = 1:length(subject_list)
            % We read the data from both triggers and use it to compute the
            % DSS transformation matrix.
            for trigger_ind = 1:length(trigger_list)

                load(fullfile('D:\Results',out_folder,'Preprocessed_data_AllChannels',...
                sprintf('data_subject-TRIG_%d-SUBJ_%d.mat',trigger_list(trigger_ind), subject_list(subject_ind))),'data_subject', 'timelock');
                

                % We baseline using the average value of the pre-stimuli
                % information. This is useful if we want to change the
                % baseline. 
                cfg = [];
                cfg.demean = 'yes'; % Necessary to baseline.
                cfg.baselinewindow = [-0.2 0];% in seconds
                data_subject = ft_preprocessing(cfg, data_subject);

                % We store the data for both conditions in a cell matrix.
                x_orig{trigger_ind} = cat(3,data_subject.trial{:}); %catenate the trial matrix in third dimension
         
            end

            
                % We reshape the data to fit the structure of DSS (NoiseTools)
                t = data_subject.time{1};
                
                % In case of single, we compute it twice.
                x = cat(3,x_orig{:});
                x = permute(x,[2,1,3]);
                
                % c0: baseline covariance
                % c1: biased covariance
                c0 = nt_cov(x);
                c1 = nt_cov(mean(x,3));
                
                % DSS
                [todss,pwr0,pwr1] = nt_dss0(c0,c1); %input the c1 biase
                % todss: matrix to convert data to normalized DSS components
                % pwr0: power per component (baseline)
                % pwr1: power per component (biased)

                z=nt_mmat(x,todss); %z is the output data(components) of DSS analysis,
%                 z1=nt_mmat(xx1,todss);
%                 z2=nt_mmat(xx2,todss);
                clear x
                % We plot the power of the components and the top10 components.
                f1 = figure(1); clf; plot(pwr1./pwr0,'.-'); ylabel('score'); xlabel('component');
                f2 = figure(2); clf;
                for iComp=1:10
                    subplot(2,5,iComp);
                    nt_bsplot(z(:,iComp,:),[],[],t);
                    title(iComp);
                    xlim([t(1) t(end)]); xlabel('s');
                end
                
                % We store the DSSed data into a Fieldtrip variable.
                z_timelock.avg = z;
                z_timelock.time = t;
                z_timelock.label = timelock.label;
                if single == 0
                    z_timelock.samples_cond1 = size(x_orig{1},3);
                    z_timelock.samples_cond2 = size(x_orig{2},3);
                else
                    z_timelock.samples_cond1 = size(x_orig{1},3);
                    z_timelock.samples_cond2 = -1;
                end
                % We store the output data.
                if store_data == 1
                    mkdir(fullfile('..','Results',out_folder,'DSS_components','Plots'))
                    if single == 0
                        savefig(f1, fullfile('..','Results',out_folder,'DSS_components','Plots',sprintf('Power-TRIG_%d_%d-SUBJ_%d', trigger_list(1),trigger_list(2), subject_list(subject_ind))));
                        savefig(f2, fullfile('..','Results',out_folder,'DSS_components','Plots',sprintf('Components-TRIG_%d_%d-SUBJ_%d', trigger_list(1),trigger_list(2), subject_list(subject_ind))))
                        % We store the transformed data so we can use it later.
                        mkdir(fullfile('..','Results',out_folder,'DSS_components'))
                        save(fullfile('..','Results',out_folder,'DSS_components',sprintf('DataDSS-TRIG_%d_%d-SUBJ_%d-COMP_%d.mat',trigger_list(1), trigger_list(2), subject_list(subject_ind), 274)),'z_timelock','todss', '-v7.3');
                    else                        
                        savefig(f1, fullfile('..','Results',out_folder,'DSS_components','Plots',sprintf('Power-TRIG_%d-SUBJ_%d', trigger_list(1), subject_list(subject_ind))));
                        savefig(f2, fullfile('..','Results',out_folder,'DSS_components','Plots',sprintf('Components-TRIG_%d-SUBJ_%d', trigger_list(1), subject_list(subject_ind))))
                        % We store the transformed data so we can use it later.
                        mkdir(fullfile('..','Results',out_folder,'DSS_components'))
                        save(fullfile('..','Results',out_folder,'DSS_components',sprintf('DataDSS-TRIG_%d-SUBJ_%d-COMP_%d.mat',trigger_list(1), subject_list(subject_ind), 274)),'z_timelock','todss', '-v7.3');
                    end
                end
                 
                clear z_timelock z
                clear x_orig
                clear data_subject
                
                     

        end

    end

        
end
