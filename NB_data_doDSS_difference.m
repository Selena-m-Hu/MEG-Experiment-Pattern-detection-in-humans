clear all
close all

subject_list = [1 2 3 4 5 6 7 8 9 10 11 12 13];

dss.z1 = [];
dss.z2 = [];
for kk=subject_list
    
    clearvars -except subject_list kk dss
    path = ['.\Data Preparation\subject' num2str(kk) '\Start\'];
    
    idxMEG=33:306;
    idxREF=4:32;
    
    load([path 'data_REG10_ranreg_merged']);
    x1=nt_trial2mat(data_merged.trial);
    
    x1=x1(:,idxMEG,:);
    
    load([path 'data_REG10_ran_merged']);
    x2=nt_trial2mat(data_merged.trial);
    
    x2=x2(:,idxMEG,:);
    
    %nt_whoss;
    sr=data_merged.fsample;
    %clear data_merged;
    %nt_whoss;
    
    x1_0=x1; x2_0=x2;
    
    x1=nt_demean(x1, 1:300); % remove mean over pre-stimulus
    x2=nt_demean(x2, 1:300);
    
    % remove 50Hz & harmonics
    [c0,c1]=nt_bias_fft(x1,[50,100,150]/sr,512);
    [todss,pwr0,pwr1]=nt_dss0(c0,c1); figure(1); clf; plot(pwr1./pwr0,'.-');
    z=nt_mmat(x1,todss);
    NDISCARD=8;
    x1=nt_tsr(x1,z(:,1:NDISCARD,:));
    [c0,c1]=nt_bias_fft(x2,[50,100,150]/sr,512);
    [todss,pwr0,pwr1]=nt_dss0(c0,c1); figure(1); clf; plot(pwr1./pwr0,'.-');
    z=nt_mmat(x2,todss);
    NDISCARD=8;
    x2=nt_tsr(x2,z(:,1:NDISCARD,:));
    
    %x1=nt_demean(x1, 1:300); % remove mean over pre-stimulus
    %x2=nt_demean(x2, 1:300);
    
    % DSS to maximize repeatability within each condition (and reduce
    % dimensionality)
    c0=nt_cov(x1)+nt_cov(x2);
    c1=nt_cov(mean(x1,3))+nt_cov(mean(x2,3));
    [todss,pwr0,pwr1]=nt_dss0(c0,c1); figure(2); clf; plot(pwr1./pwr0,'.-');
    NKEEP=12;
    xx1=nt_mmat(x1,todss(:,1:NKEEP));
    xx2=nt_mmat(x2,todss(:,1:NKEEP));
    figure(3); clf;
    nt_banner('DSS to maximize reproducibility');
    for c=1:4;
        subplot (4,2,(c-1)*2+1);
        nt_bsplot(xx1(:,c,:));
        title(['Transition, DSS', num2str(c)]);
        subplot (4,2,c*2);
        nt_bsplot(xx2(:,c,:));
        title(['Control, DSS', num2str(c)]);
    end
    
    %xx1=nt_demean(xx1, 1:300); % remove mean over pre-stimulus
    %xx2=nt_demean(xx2, 1:300);
    
    % DSS to maximize differences between RAN and RANREG
    c0=nt_cov(xx1)+nt_cov(xx2);
    c1=nt_cov(mean(xx1,3)-mean(xx2,3));
    [todss,pwr0,pwr1]=nt_dss0(c0,c1); figure(4); clf; plot(pwr1./pwr0,'.-');
    z1=nt_mmat(xx1,todss);
    z2=nt_mmat(xx2,todss);
    figure(5); clf;
    nt_banner('DSS to maximize the Trans-Ctrl difference');
    for c=1:4;
        subplot (4,3,(c-1)*3+1);
        nt_bsplot(z1(:,c,:));
        title(['Transition, DSS', num2str(c)]);
        subplot (4,3,(c-1)*3+2);
        nt_bsplot(z2(:,c,:));
        title(['Control, DSS', num2str(c)]);
        subplot (4,3,(c-1)*3+3);
        plot([mean(z1(:,c,:),3),mean(z2(:,c,:),3)])
        title(['both', num2str(c)]);
    end
    
    cfg.layout='CTF274_FIL_lay.mat';
    figure(6); clf
    C=nt_xcov(z1(:,1,:),x1) + nt_xcov(z2(:,1,:),x2);
    fTopoplot(cfg,C');
    
    zz=cat(3,z1,z2);
    x=cat(3,x1,x2);
    NKEEP=2;
    C=nt_regcov(nt_xcov(x,zz(:,1:NKEEP,:)),nt_cov(zz(:,1:NKEEP,:)));
    
    % project back
    %av1=nt_mmat(z1(:,1:NKEEP,:),C);
    %av2=nt_mmat(z2(:,1:NKEEP,:),C);
    %data1=data_merged;
    %data1.trial=nt_mat2trial(av1);
    %data2=data_merged;
    %data2.trial=nt_mat2trial(av2);
    
    dss.z1{kk} = z1(:,1:NKEEP,:);
    dss.z2{kk} = z2(:,1:NKEEP,:);
end

save _DSSdiff dss

return;

%%

z1 = [dss.z1];
z2 = [dss.z2];
z1avg = []; z2avg = [];

for i=1:numel(dss.z1) % compute mean across trials
    z1{i}=nt_demean(z1{i}, 1:300); % remove mean over pre-stimulus
    z2{i}=nt_demean(z2{i}, 1:300);
    z1avg{i} = mean(z1{i},3);
    z2avg{i} = mean(z2{i},3);
end

dim = ndims(z1avg{1});
z1avg= cat(dim+1,z1avg{:});

dim = ndims(z2avg{1});
z2avg= cat(dim+1,z2avg{:});

nt_banner('average DSS');
for c=1:2;
    subplot (4,3,(c-1)*3+1);
    nt_bsplot(z1avg(:,c,:),2,'zerobased',[],1,1);
    title(['Transition, DSS', num2str(c)]);
    subplot (4,3,(c-1)*3+2);
    nt_bsplot(z2avg(:,c,:),2,'zerobased',[],1,1);
    title(['Control, DSS', num2str(c)]);
    subplot (4,3,(c-1)*3+3);
    plot([mean(z1avg(:,c,:),3),mean(z2avg(:,c,:),3)])
    title(['both', num2str(c)]);
end