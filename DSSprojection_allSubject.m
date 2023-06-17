function DSSprojection_allSubject(trigger_list, subject_list, config)
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
        rejectVisual = config.reject_visual;   
        n_components = config.n_components;
   

for subject_ind = 1:length(subject_list)
            
        %load DSS components
            load(fullfile('..','Results',out_folder,'DSS_components',sprintf('DataDSS-TRIG_%d_%d-SUBJ_%d-COMP_%d.mat',trigger_list(1), trigger_list(2)...
             ,subject_list(subject_ind), 274)),'z_timelock');

                    t_cond{1} = 1:z_timelock.samples_cond1;
                    t_cond{2} = z_timelock.samples_cond1+(1:z_timelock.samples_cond2);       
           

                                         
       for trigger_ind = 1:length(trigger_list)
                
                load(fullfile('D:\Results',out_folder,'Preprocessed_data_AllChannels',...
                sprintf('data_subject-TRIG_%d-SUBJ_%d.mat',trigger_list(trigger_ind), subject_list(subject_ind))),'data_subject');
                      
          
                %% We reshape the data for the selected trigger into a cell.
%                 x_orig{trigger_ind} = cat(3,data_subject_BS.trial{:});
                x_orig{trigger_ind} = nt_trial2mat(data_subject.trial); 


                %% Here we back-project the data into the original 274
                % channels, using the number of components that we want to
                % keep.
                
                % **For newly aquired data by MYH,RB, there are only 273
                % channels**
                
                 c = nt_xcov(z_timelock.avg(:,:,t_cond{trigger_ind}),x_orig{trigger_ind}); % c is cross-covariance between z(raw data) and x(DSS components)

                    x_dss2 = nt_mmat(z_timelock.avg(:,1:n_components,t_cond{trigger_ind}),c(1:n_components,:)); % project from component to sensor space, only using the KEEP components
                    x_dss2 = nt_mat2trial(x_dss2);  %convert the dss projection into the configuration that fieldtrip can read
%                   x_dss = permute(x_dss(:,:,:), [2,1,3]); %Antonio's
%                   script
                    data_subject.trial = x_dss2;              
                 %% We reject outlier trials                    
                    if rejectVisual
                    cfg.channel = 'all';
                    dss_data_subject = ft_rejectvisual(cfg,data_subject);                                     
                    else
                    dss_data_subject  = data_subject;
                    end 

                   x_dss = cat(3,dss_data_subject.trial{:}); %convert the dss data into 3 dimension and save it into one cell
                 %% We save the output matrix into folder
                    if store_output == 1
                        
                      mkdir(fullfile('D:\Results',out_folder,'DSS','DSStransformed','visualRejection')) % data structure is saved as 3D 
                       
                      save(fullfile('D:\Results',out_folder,'DSS','DSStransformed','visualRejection',...
                      sprintf('DSSdata_subject-TRIG_%d-SUBJ_%d-COMP_%d.mat',trigger_list(trigger_ind), subject_list(subject_ind), n_components)),'dss_data_subject','x_dss','-v7.3');
                       
 

                   clear x_dss c z
                   clear x_dss2
                   clear x_orig
                    end

            clear z_timelock
            clear dss_data_subject
            clear data_subject


        end

        
end
