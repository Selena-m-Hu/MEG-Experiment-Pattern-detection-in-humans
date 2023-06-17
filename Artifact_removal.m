%this script analyzes the M100 response in the localizer script and finds
%the 20 strongest channels in each hemisphere.

clear all;close all;
addpath(fullfile('..','fieldtrip'));
data_folder = fullfile('..','data_longshort/Subj15/');
cfg = [];
cfg.channel = 'MEG'; % Le decimos que queremos capturar todos los canales de MEG.
username = 'nh120784_Theofilos_20150205_05.ds';
% cfg.dataset = fullfile(data_folder,'ac040981_MChait2_20120417_01.ds');  % Fichero a leer.
cfg.dataset = fullfile(data_folder,username);  % Fichero a leer.

cfg.lpfilter ='yes';
cfg.lpfreq = 30;
data = ft_preprocessing(cfg);
% cfg = [];
% cfg.layout='CTF275.lay';
% cfg.viewmode = 'vertical';
% cfg.continuous = 'yes';
% cfg.blocksize = 30;
% cfg.channels = [1:5];
% cfg = ft_databrowser(cfg,data);



% cfg.dataset = fullfile('am061187_MChait2_20120510_01.ds');
cfg.trialdef.eventtype  = 'UPPT001'; % Aquí indicamos el itpo de evento que nos ayudará a separar en epochs.
cfg.trialdef.eventvalue = 10;
cfg.trialdef.prestim    = 0.2; % Espacio de preestímulo (antes del estímulo). Puede usarse para extraer la media.
cfg.trialdef.poststim   = 15%.5; % Espacio de señal que nos interesa, después del evento.
cfg = ft_definetrial(cfg);

cfg.demean      = 'yes';
cfg.baselinewindow = [-.2 0];% in seconds
data = ft_preprocessing(cfg);
cfg=[];
% cfg.hpfilter ='yes';
% cfg.hpfreq = 0.5;
cfg.lpfilter ='yes';
cfg.lpfreq = 30;

data=ft_preprocessing(cfg, data);
% save data
cfg=[];
cfg.method = 'summary';
cfg.channel = {'MEG'};
data = ft_rejectvisual(cfg,data);
% 
% 
% cfg = [];
% cfg.method = 'fastica';
% cfg.fastica.maxNumIterations = 100;
% ic_data = ft_componentanalysis(cfg,data);
% 
% cfg = [];
% cfg.layout='CTF275.lay';
% cfg.viewmode = 'component';
% cfg.continuous = 'yes';
% cfg.blocksize = 30;
% cfg.channels = [1:5];
% ft_databrowser(cfg,ic_data);

%%
% Aquí se extraen los canales, considerando que pueden tener magnitudes
% opuestas: esto es, el mismo gráfico pero con el signo cambiado a negativo
% o positivo. En cualquier caso, el resultado promedio debería ser el
% adecuado.
cfg = [];
cfg.channel = 'MEG';

% Extraemos algunos estadísticos de los trials. Por ejemplo, el promedio de
% la magnitud.
timelock = ft_timelockanalysis(cfg, data);
timelock.fsample = 600;

%finding the 20 strongest channels in each hemisphere

M100dat=timelock.avg'; % Obsérvese que aquí se transpone la media, por lo que tenemos (tiempo x canal).

t0 = 0.2*timelock.fsample;

% Ordenamos las amplitudes, para lo cual escogemos un rango de instantes de
% tiempo, y todos los canales.
amps=mean(M100dat((t0+0.09*timelock.fsample):(t0+0.11*timelock.fsample), :),1);
[ampsSorted,idx]=sort(amps,2,'descend');

chnsSorted = timelock.label(idx);

%selecting channels:

chns_selectedLpos=[];
chns_selectedRpos=[];
chns_selectedLneg=[];
chns_selectedRneg=[];
leftChansCountPos=0;
rightChansCountPos=0;
leftChansCountNeg=0;
rightChansCountNeg=0;

for count=1:length(ampsSorted)
    strPos=chnsSorted(count);
    strNeg=chnsSorted(end-count+1);
    % Lookup in Left Hemisphere
    if  ~isempty(strfind(strPos{1},'ML'))
        if(leftChansCountPos<10)
            leftChansCountPos=leftChansCountPos+1;
            chns_selectedLpos=[chns_selectedLpos strPos];
        end
    end
    if ~isempty(strfind(strNeg{1},'MLT'))
        if(leftChansCountNeg<10)
            leftChansCountNeg=leftChansCountNeg+1;
            chns_selectedLneg=[chns_selectedLneg strNeg];
        end
        
    end
    
    % Lookup in Right Hemisphere
    if  ~isempty(strfind(strPos{1},'MRT'))
        if(rightChansCountPos<10)
            rightChansCountPos=rightChansCountPos+1;
            chns_selectedRpos=[chns_selectedRpos strPos];
        end
    end
    if ~isempty(strfind(strNeg{1},'MR'))
        if(rightChansCountNeg<10)
            rightChansCountNeg=rightChansCountNeg+1;
            chns_selectedRneg=[chns_selectedRneg strNeg];
        end
    end
end

chns_selectedL = [chns_selectedLpos chns_selectedLneg];
chns_selectedR = [chns_selectedRpos chns_selectedRneg];

chnsL_num=[];
for count1=1:length(timelock.label)
    for count2=1:length(chns_selectedL)
        if (strcmp(timelock.label{count1},chns_selectedL{count2}) ~= 0)
            chnsL_num=[chnsL_num count1];
        end
    end
end

chnsR_num=[];
for count1=1:length(timelock.label)
    for count2=1:length(chns_selectedR)
        if (strcmp(timelock.label{count1},chns_selectedR{count2}) ~= 0)
            chnsR_num=[chnsR_num count1];
        end
    end
end

chns_selectedL = timelock.label(chnsL_num); 
chns_selectedR = timelock.label(chnsR_num);

selected_dataL = timelock.avg(chnsL_num,:);
time = timelock.time;
figure(1);plot(time, selected_dataL)
title ('LH channels');
legend(chns_selectedL);

selected_dataR = timelock.avg(chnsR_num,:);
time = timelock.time;
figure(2);plot(time, selected_dataR)
title ('RH channels');
legend(chns_selectedR);


%%
plot(time,rms([selected_dataL; selected_dataR]))
% imagesc([abs(selected_dataL); abs(selected_dataR)])
aux.username = username;
aux.data_folder = data_folder;
[trig_sample,trig_value] = Trigger_analysis(aux);
stem(trig_sample/600, trig_value)


%%
cfg = [];
% cfg.parameter = 'avg';
cfg.layout='CTF275.lay'; % 274
cfg.xlim=[0.1 0.1]';
cfg.marker = 'labels';
cfg.interactive = 'yes';
cfg.markerfontsize = 2;


figure(3);
channel = ft_channelselection([chns_selectedL, chns_selectedR], timelock.label);
cfg = [];
% cfg.parameter = 'avg';
cfg.layout='CTF275.lay';
cfg.xlim=[0.1 0.1]';
cfg.marker = 'labels';
cfg.interactive = 'yes';
cfg.markerfontsize = 2;
 cfg.highlight='on';
 cfg.highlightchannel=channel;
 cfg.highlightfontsize=20;
ft_topoplotER(cfg, timelock); title ('M100 response');
%%
% save am061187SelectedChannels chns_selectedR chns_selectedL 
