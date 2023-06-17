function pre_DSSprojection_Testing(trigger_list, subject_list, config)
    % This functions reads the DSSed data from pre_DSScomputation.m and
    % projects it back into channel space, keeping only the number of DSS
    % components that we want (usually 3).
    %
    % trigger: couples (5,15) or (10, 20).
    % * 5 : 3 second RAND sequences.
    % * 10: 15 second RAND sequences.
    % * 15: 3 second REG sequences.
    % * 20: 15 second REG sequences.
    %
    % subject:
    % * 2-15, the index of the subject collected by Antonio.
    % * 16-24, the new subjects aquired by RB, MYH
    %
    % config: allows for certain configurations.
    %   .out_folder: name of the folder where data will be stored.
    %   .store_data: set to 1 to store in a .mat file the whole data of the
    %       subject once that it has been processed.
    %   .n_components: number of components from DSS to keep.
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
    % Last update: Mingyue Hu, Sep/2022
    
    out_folder = config.out_folder;
    store_output = config.store_data;
    single = config.single;
    rejectVisual = config.reject_visual;   
    n_components = config.n_components;
    
 

        for subject_ind = 1:length(subject_list)
            % We read the data from both triggers and use it to compute the
            % DSS transformation matrix.
            if single == 0
                load(fullfile('..','Results',out_folder,'DSS_components',...
                    sprintf('DataDSS-TRIG_%d_%d-SUBJ_%d-COMP_%d.mat',trigger_list(1), trigger_list(2), subject_list(subject_ind), 274)));
                t_cond{1} = 1:z_timelock.samples_cond1;
                t_cond{2} = z_timelock.samples_cond1+(1:z_timelock.samples_cond2);       
            end

                           
            for trigger_ind = 1:length(trigger_list)
                if single == 1   %if you analysed the DSS across conditions
                    load(fullfile('..','Results',out_folder,'DSS_debugging',...
                        sprintf('DataDSS-TRIG_%d-SUBJ_%d-COMP_%d.mat',trigger_list(trigger_ind), subject_list(subject_ind), 274)));
                    t_cond{1} = 1:z_timelock.samples_cond1;
                    t_cond{2} =  1:z_timelock.samples_cond1;   
                end
                load(fullfile('..','Results',out_folder,'Preprocessed_data_AllChannels',...
                    sprintf('Long_timelock-TRIG_%d-SUBJ_%d', trigger_list(trigger_ind) ,subject_list(subject_ind))));

                % Inherited code: we limit the data that we use for DSS and
                % we baseline again.
                cfg = [];
                if trigger_list(trigger_ind) == 5 | trigger_list(trigger_ind) == 15
                    cfg.toilim = [-0.2, 4];
                else
                    cfg.toilim = [-0.2,16];
                end
                data_subject = ft_redefinetrial(cfg,data_subject);

                %% We baseline using the average value of the pre-stimuli
                % information.
                cfg = [];
                cfg.demean = 'yes'; % Necessary to baseline.
                cfg.baselinewindow = [-0.2 0];% in seconds
                data_subject = ft_preprocessing(cfg, data_subject);
          
                %% We reshape the data for the selected trigger into a cell.
%                 x_orig{trigger_ind} = cat(3,data_subject_BS.trial{:});
                x_orig{trigger_ind} = nt_trial2mat(data_subject.trial); 
                % We check if the folder where the projections
                % actually exist.
                dss_folder_flag = isdir(fullfile('..','Results',out_folder,'DSS_components')) == 0;


                %% Here we back-project the data into the original 274
                % channels, using the number of components that we want to
                % keep.
                
                 
                 mkdir(fullfile('..','Results',out_folder,'DSS_debugging'));

%                 c = nt_xcov(z_timelock.avg(:,:,t_cond{trigger_ind}),permute(x_orig{trigger_ind},[2,1,3]))/...
%                         size(z_timelock.avg,1)*size(z_timelock.avg,3); 
%                         c is cross-covariance between z(raw data) and
%                         x(DSS components)  %Antonio's script
                 c = nt_xcov(z_timelock.avg(:,:,t_cond{trigger_ind}),x_orig{trigger_ind}); % c is cross-covariance between z(raw data) and x(DSS components)

                    x_dss2 = nt_mmat(z_timelock.avg(:,1:n_components,t_cond{trigger_ind}),c(1:n_components,:)); % project from component to sensor space, only using the KEEP components
                    x_dss2 = nt_mat2trial(x_dss2);  %convert the dss projection into the configuration that fieldtrip expects
%                   x_dss = permute(x_dss(:,:,:), [2,1,3]); %Antonio's
%                   script
                                  
                 %% We reject outlier trials                    
%                     if rejectVisual
%                      subject = subject_list(subject_ind);
%                      trigger = trigger_list(trigger_ind);
%                     visualRejection(trigger, subject, config, data_subject, x_dss2)
%                     end 

                 %% We save the output matrix into folder
                    if store_output == 1
                        
                        if single == 0
                            save(fullfile('..','Results',out_folder,'DSS_debugging',...
                                sprintf('Trans-TRIG_%d-SUBJ_%d-COMP_%d.mat',trigger_list(trigger_ind), subject_list(subject_ind), n_components)),'x_dss2');
                        else                          
                            save(fullfile('..','Results',out_folder,'DSS_debugging',...
                                sprintf('Trans_dss-TRIG_%d-SUBJ_%d-SINGLE-COMP_%d.mat',trigger_list(trigger_ind), subject_list(subject_ind), n_components)),'x_dss2');
                        end
                    end
                

                
                % We also store the projected data into a variable that we
                % will be used in posterior analysis and plots.
                % Subject 16 was saved seperately as it only has 272
                % channels
                
                x_dss = cat(3,x_dss2{:}); %convert the dss data into 3 dimension and save it into one cell
                
                % save the variable as the way Antoniot did ('x_dss'), for adapting
                % other scripts he programmed 
                if single == 0
                save(fullfile('..','Results',out_folder,'DSS_debugging',...
                    sprintf('dss-TRIG_%d-SUBJ_%d-COMP_%d.mat',trigger_list(trigger_ind), subject_list(subject_ind), n_components)),'x_dss');
                else
                save(fullfile('..','Results',out_folder,'DSS_debugging',...
                    sprintf('dss-TRIG_%d-SUBJ_%d-SINGLE-COMP_%d.mat',trigger_list(trigger_ind), subject_list(subject_ind), n_components)),'x_dss');
                end
                
                
                %% Put data from all subjects into one mat file, for further analysis purpose
%                 dss_comp(:,:, subject_ind, trigger_ind) = mean(x_dss,3);  %Dss components              
%                 x_comp(:, :, subject_ind, trigger_ind) = mean(x_orig{trigger_ind},3); % mean of all trials from raw data after preprocessing

                clear x_dss c z
            end
            clear z_timelock


        end

        
    
end
      

        

