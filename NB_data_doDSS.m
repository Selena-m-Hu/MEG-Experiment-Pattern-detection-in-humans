% DSS Script, does in order:
% 1) suppress slow drift components (recognizable because they maximize
% variance between trials).
% 2) remove some outlier trials.
% 3) suppress 50 Hz and harmonics.
% 4) cut data to focus on transition part.
% 5) baseline correct and remove some more outlier trials.
% 6) apply DSS to find maximally repeatable components.
% 7) remove a few more outlier trials.
% 8) project back to sensor space.

clear all
subject_list = [13];

for kk=subject_list
    
    clearvars -except subject_list kk
    cd(['.\Data Preparation\subject' num2str(kk)])
    
    %load CTF274_FIL_lay.mat
    %cfg.layout=lay;
    cfg.emarkersize=1;
    fname='data_REG10_'; suffix='.mat';
    
    % h=ft_read_header(fname);
    
    idxMEG=33:306;
    idxREF=4:32;
    
    conds={'ran','ranreg','reg','regran'};
    
    % load and perform initial cleaning
    for iCond=1:4
        
        fn=[fname,conds{iCond},'_merged',suffix];
        disp(fn);
        disp('load data...')
        
        load(['.\Start\' fn]);
        x=nt_trial2mat(data_merged.trial);
        data_merged.trial=[];
        datas{iCond}=data_merged;
        DSR=1; % downsample ratio
        x=nt_dsample(x,DSR);
        sr=data_merged.fsample/DSR;
        disp('done')
        nt_whoss;
        
        t=data_merged.time{1}(1:DSR:end);
        t=t(1:size(x,1));
        
        meg=x(:,idxMEG,:);
        ref=x(:,idxREF,:);
        miscdata{iCond}=x(:,[1:3,307],:);
        refdata{iCond} = ref;
        x=meg; clear meg;
        
        [~,~,ntrials]=size(x);
        trial_list=1:ntrials;
        
        % zap bad channels
        prof=mean(nt_unfold(x.^2));
        BADFACTOR=20;
        badIdx=find(prof>BADFACTOR*mean(prof));
        x(:,badIdx,:)=0;
        ref(:,badIdx,:)=0;
        if ~isempty(badIdx);
            disp(['WARNING: zapping channels ',num2str(badIdx)]);
        end
        
        % merge meg and ref to derive estimates of noise components
        xx=nt_pca(nt_normcol(cat(2,x,ref)));
        THRESH=0.1;
        idx= mean(nt_unfold(xx.^2))>THRESH;
        xx=xx(:,idx,:); % retain only components shared between channels
        
        suppressed=[];
        
        % suppress components with large variance between trials (slow drifts)
        REMOVE=1:10;
        c0=nt_cov(nt_demean2(xx)); c1=nt_cov(xx);
        [todss,pwr0,pwr1]=nt_dss0(c0,c1);
        figure(10); clf; plot(pwr1./pwr0,'.-');
        hold on; plot(pwr1(REMOVE)./pwr0(REMOVE),'.r');
        title('slow trends'); drawnow
        z=nt_mmat(xx,todss);
        x=nt_tsr(x,z(:,REMOVE,:));
        xx=nt_tsr(xx,z(:,REMOVE,:));
        suppressed=cat(2,suppressed,z(:,REMOVE,:));
        
        % discard outliers
        disp([num2str(size(x,3)), ' trials']);
        THRESH=2;
        figure(11); clf; subplot 121; nt_find_outlier_trials(nt_demean2(x)); title('before');
        idx=nt_find_outlier_trials2(nt_demean2(x),THRESH);
        x=x(:,:,idx);
        trial_list=trial_list(idx);
        x=nt_demean(x);
        xx=xx(:,:,idx);
        xx=nt_demean(xx);
        
        suppressed=suppressed(:,:,idx);
        figure(11); subplot 122; nt_find_outlier_trials(nt_demean2(x)); title('after'); drawnow
        disp(['discard outlier trials, remains ', num2str(size(x,3))]);
        
        % suppress components with 50Hz and harmonics
        REMOVE=1:3;
        [c0,c1]=nt_bias_filter(xx,[50,100,150]/sr,512);
        [todss,pwr0,pwr1]=nt_dss0(c0,c1);
        figure(10); clf; plot(pwr1./pwr0,'.-');
        hold on; plot(pwr1(REMOVE)./pwr0(REMOVE),'.r');
        title('50 Hz & harmonics'); drawnow
        z=nt_mmat(xx,todss);
        x=nt_tsr(x,z(:,REMOVE,:));
        suppressed=cat(2,suppressed,z(:,REMOVE,:));
        
        % baseline correct
        x=nt_demean2(x,find(t<0)); % use nt_demean() instead if there's a strong trend
        
        % regress out linear trend
        %x=nt_tsr_nodemean(x,repmat((t'),[1,1,size(x,3)]));
        
        % discard some more outlier trials
        THRESH=3;
        figure(11); clf; subplot 121; nt_find_outlier_trials((x)); title('before')
        idx=nt_find_outlier_trials2(x,THRESH);
        x=x(:,:,idx);
        x=nt_demean(x);
        trial_list=trial_list(idx);
        suppressed=suppressed(:,:,idx);
        figure(11); subplot 122; nt_find_outlier_trials((x)); title('after'); drawnow
        disp(['discard outlier trials, remains ', num2str(size(x,3))]);
        
        % baseline correct
        x=nt_demean(x,find(t<0));
        
        data{iCond}=x;
        supp{iCond}=suppressed;
        trial_lists{iCond}=trial_list;
        nt_whoss;
    end
    
    % calculate covariances
    c0=zeros(numel(idxMEG));
    c1=zeros(numel(idxMEG));
    for iCond=1:4
        c0=c0+nt_cov(data{iCond});
        c1=c1+nt_cov(mean(data{iCond},3));
    end
    
    % DSS
    [todss,pwr0,pwr1]=nt_dss0(c0,c1);
    figure(10); clf; plot(pwr1./pwr0,'.-');
    hold on; plot(pwr1(1:6)./pwr0(1:6),'.r');
    title('repeatability')
    
    
    % apply DSS matrix
    c=zeros(size(todss'));
    for iCond=1:4
        trial_list=trial_lists{iCond};
        x=data{iCond};
        suppressed=supp{iCond};
        z=nt_mmat(x,todss);
        
        % discard a few more outlier trials
        THRESH=2;
        figure(11); clf; subplot 121; nt_find_outlier_trials((x)); title('before');
        idx=nt_find_outlier_trials2(x,THRESH);
        x=x(:,:,idx);
        z=z(:,:,idx);
        trial_list=trial_list(idx);
        suppressed=suppressed(:,:,idx);
        figure(11); subplot 122; nt_find_outlier_trials((x)); title('after'); drawnow
        disp(['discard outlier trials, remains ', num2str(size(x,3))]);
        
        data{iCond}=x;
        comp{iCond}=z;
        supp{iCond}=suppressed;
        trial_lists{iCond}=trial_list;
        
        c=c+nt_xcov(z,x);
    end
    
    % plot components 1:6
    figure(1); clf
    nt_banner([fname]);
    for iComp=1:6;
        subplot (3,2,iComp);
        zz=[];
        for iCond=1:4
            zz=[zz,mean(comp{iCond}(:,iComp,:),3)];
        end
        plot(t,zz);
        title(num2str(iComp));
    end
    
    % plot topographies
    figure(2); clf
    nt_banner([fname]);
    for iComp=1:6;
        subplot (3,2,iComp);
        %  topoplot(cfg,c(iComp,:)');
        title(num2str(iComp));
        drawnow
    end
    
    % project back to sensor space
    %R=input('how many components to keep?');
    KEEP=1:2;
    for iCond=1:4
        xxx=nt_mmat(comp{iCond}(:,KEEP,:),c(KEEP,:));
        misc=miscdata{iCond};
        misc=misc(:,:,trial_lists{iCond}); %remove outliers from misc data a
        ref=refdata{iCond};
        ref=ref(:,:,trial_lists{iCond}); %remove outliers from misc data a
        data_merged=datas{iCond};
        
        % put back into ft structure
        data_merged.trial=nt_mat2trial([misc(:,1:3,:) ref xxx misc(:,4,:)]);  % we may need to tweak other fields of data
        data_merged.time = data_merged.time(trial_lists{iCond});
        data_merged.label = data_merged.label;
        
        fn=[fname,conds{iCond},'_merged' suffix];
        
        if ~exist(['.\DSS'],'dir')
            mkdir(['.\DSS']);
        end
        save(['.\DSS\' fn], 'data_merged');
    end
    
    
    
    cd ../..
    
end