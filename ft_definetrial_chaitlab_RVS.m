function [cfg] = ft_definetrial_chaitlab_RVS(cfg)
%Sort the trials in chronological order
% Definition of trials based on events at Ear Institute BIOSEMI system
% Trigger values are determined from the DURATION of the trigger pulses !
%
% FORMAT [trl, conditionlabels, S] = spm_eeg_definetrial(S)
% S                 - input structure (optional)
% (optional) fields of S:
%   S.D             - MEEG object or filename of M/EEG mat-file
%   S.timewin       - time window (in PST ms)
%   S.trialdef      - structure array for trial definition with fields (optional)
%       S.trialdef.conditionlabel - string label for the condition
%       S.trialdef.eventtype      - string (should be 'trigger_up')
%       S.trialdef.eventvalue     - trigger value for this condition
%       S.trialdef.trlshift       - shift the triggers by a fixed amount (sec)
%                                   (e.g. projector delay). One per
%                                   condition/trigegr type
%   cfg.trialdef.trlshiftpertrial
%   S.filename      - path of BDF file to read
% OUTPUT:
%   trl             - Nx3 matrix [start end offset]
%   conditionlabels - Nx1 cell array of strings, label for each trial
%   S               - modified configuration structure (for history)
%__________________________________________________________________________
% Nicolas Barascud
% Adapted code from Vladimir Litvak, Robert Oostenveld
% Edited by Rosy Southwell - May 2018:
% - increased tolerance for trig durations to +- 3 samples
% - allow a different trialshift per condition
% - allow a different trialshift per trial (cfg.trialdef.trlshiftpertrial)
% - ignore trigger values / conditions which arent passed in in the trialdef, don't include
% these in cfg.event

%% Parameters
%--------------------------------------------------------------------------
% data = varargin{1};
% read the header, required to determine the stimulus channels and trial specification
cfgh=[];
if isfield(cfg,'headerformat')
 hdr = ft_read_header(cfg.dataset,'headerformat',cfg.headerformat);   
else
    hdr = ft_read_header(cfg.dataset);
end



% find the STATUS channel and read the values from it
trigchani = find(strcmpi(hdr.label,'STATUS'));
trigdata = ft_read_data(cfg.dataset,'chanindx', trigchani);

pretrig  = cfg.trialdef.prestim;
posttrig = cfg.trialdef.poststim;
begsample = 1; 
%% Read trigger channel (code adapted from ft_read_event.m)
%--------------------------------------------------------------------------

% convert to 32-bit integer representation and only preserve the lowest 24 bits
trigger = bitand(int32(trigdata), 2^24-1);

byte1 = 2^8  - 1;
byte2 = 2^16 - 1 - byte1;
byte3 = 2^24 - 1 - byte1 - byte2;

% get the respective status and trigger bits
trigger = bitand(trigger, bitor(byte1, byte2)); % this is contained in the lower two bytes
trigger = bitand(trigger, byte1); % this is contained in the lower two bytes


% determine when the respective status bits go up or down
flank_trigger = diff([trigger]);

%% Create event structure
%--------------------------------------------------------------------------

event       = [];
pad         = 0;
trigshift   = 0;
% convert the trigger into an event with a value at a specific sample
for i=find(flank_trigger>0)
    event(end+1).type   = 'trigger_up';        % distinguish between up and down flank
    event(end  ).sample = i + begsample-1;      % assign the sample at which the trigger has gone down
    event(end  ).value  = double(trigger(i+trigshift));      % assign the trigger value just _after_ going up
end

% Sort events in chronological order
[tmp, ind] = sort([event.sample]);
event = event(ind);

% Find distance between consecutive events
distance = find(flank_trigger<0,length(find(flank_trigger>0))) - find(flank_trigger>0,length(find(flank_trigger<0)));

% Remove events that are too close to each other
idx = find(distance< 4);
event(idx) = [];
distance(idx) = [];

%% Safecheck. It's best to use MATLAB 2015a at this point.
if verLessThan('matlab','8.5.0')
    disp('Warning: No proper trigger safecheck was made.');
    disp('Matlab 2015a (or later) is required for optimal processing.');
    disp('Performing ''basic'' check instead (this is more likely to fail later...)');
    
    % Find unique trigger values in data (there should be as many as conditions)
    safecheck = unique(distance);
    
    if numel(safecheck) ~=numel(cfg.trialdef)
        disp([num2str(numel(safecheck)) ' unique trigger values were found in ' ]);
        disp(['trigger channel (' num2str(numel(cfg.trialdef)) ' required).']);
        disp('Attempting to continue anyway...');
    else
        disp(['A total ' num2str(numel(safecheck)) ' unique trigger values were found in ']);
        disp(['trigger channel (' num2str(numel(cfg.trialdef)) ' required). All OK !']);
    end
    
else % execute code for R2015a later
    % Find unique trigger values in data (there should be as many as conditions)
    tolerance = 2/max(distance); %% Tolerance decided
    safecheck = uniquetol(distance,tolerance);
    
    if numel(safecheck) ~=numel(cfg.trialdef.conditionlabel)
        disp([num2str(numel(safecheck)) ' unique trigger values were found in ']);
        disp(['trigger channel (' num2str(numel(cfg.trialdef.conditionlabel )) ' required). ']);
        %         return; %%Commented by Sijia: Retrive trials with only 2 out of 4
        %         triggers in EEG2 2016/02/09
    else
        disp(['A total ' num2str(numel(safecheck)) ' unique trigger values were found in']);
        disp(['trigger channel (' num2str(numel(cfg.trialdef.conditionlabel )) ' required). All OK !']);
    end
    
end

% Dirty fix for inaccurate trigger durations. Biosemi system codes event
% durations with ±1 time sample precision, so we manually correct these
% inaccuracies
tol = 3; % RVS changed from 2
for i=1:numel(cfg.trialdef.conditionlabel) % for all trial types
    % find trigger values ± 3 samples to be on the safe side
    idx = (distance >= (cfg.trialdef.eventvalue(i)-tol) ...
        & distance <= (cfg.trialdef.eventvalue(i)+tol)); %%%%% Tolerance
    distance(idx) = cfg.trialdef.eventvalue(i);
end

% % Remove trigger_down events as we don't need them anymore
% idx = strcmp({event.type},'trigger_down');
% event(idx) = [];

for i=1:numel(event) % re-populate value with DURATION instead of AMPLITUDE of trigger event
    event(i).value = distance(i);
end
% remove any events not requested in trialdef
ix_rem = find(~ismember([event(:).value],cfg.trialdef.eventvalue));
if ~isempty(ix_rem)
    event(ix_rem) = [];
    distance(ix_rem) = [];
    disp(['Ignoring ' num2str(length(ix_rem)) ' triggers not requested in trialdef; ' num2str(length(event)) ' remaining.']);
end
cfg.event = event;

%% Build trl matrix based on selected events
%--------------------------------------------------------------------------

for j=1:numel(cfg.trialdef.eventvalue)
    if ~isfield(cfg.trialdef,'trlshift')
        trlshift(j) = 0;
    elseif numel(cfg.trialdef.eventvalue) == numel(cfg.trialdef.trlshift)
        trlshift(j) = round(cfg.trialdef.trlshift(j) * hdr.Fs); % assume passed as s
    elseif numel(cfg.trialdef.trlshift) == 1
        trlshift(j) = cfg.trialdef.trlshift * hdr.Fs;
    else
        error('cfg.trialdef.trlshift is wrong dimensions!')
    end
end
trl = [];
conditionlabels = {};condis=[];
%% New
for i=1:numel(event)
    if ~strcmp(cfg.trialdef.eventtype,'trigger_up')
        disp('ERROR: S.trialdef.eventtype should be ''trigger_up''. Aborting!')
        return;
    end
    [icondition icondition]=find(cfg.trialdef.eventvalue==event(i).value);
    if isempty(icondition) % then we don't need this event as it is not in the list
        error('Wut. some unwanted trials are in the event structure...')
    end
    trloff = round(pretrig*hdr.Fs); % assume passed as s
    trlbeg = event(i).sample - trloff; % trloff is prestim; positive number
    trldur = round((pretrig+posttrig)*hdr.Fs);% assume passed as s
    trlend = trlbeg + trldur;
    
    % Added by Rik in case wish to shift triggers (e.g, due to a delay
    % between trigger and visual/auditory stimulus reaching subject).
    % (i) shift trigger by set amount per condition
    
    trlbeg = trlbeg + trlshift(icondition);
    trlend = trlend + trlshift(icondition);
    % (ii) shift trigger by set amount per event
    if isfield(cfg.trialdef,'trlshiftpertrial')
        if numel(event) ~= numel(cfg.trialdef.trlshiftpertrial)
            warning('***number of trigger events found is not equal to length of trlshiftpertrial, weird things may happen!')
        end
        trlbeg = trlbeg + round(cfg.trialdef.trlshiftpertrial(i)* hdr.Fs);
        trlend = trlend + round(cfg.trialdef.trlshiftpertrial(i)* hdr.Fs);
    end
    % Add the beginsample, endsample and offset of this trial to the list
    trl = [trl; trlbeg trlend -trloff distance(i)]; % !! Fieldtrip difference prestim positive value means trial begins before trigger
    conditionlabels{end+1} = cfg.trialdef.conditionlabel(icondition);
        condis(end+1) = icondition;

end

cfg.trl = trl;
cfg.conditionlabels = conditionlabels;
cfg.condi = condis;

%% plot
figure(99); clf;

samp = [event.sample];
val = [event.value];
cmap = parula(8);
for i = 1:length(cfg.trialdef.eventvalue)
    yyaxis left
plot(samp(val==cfg.trialdef.eventvalue(i)), val(val==cfg.trialdef.eventvalue(i)),'*','Color',cmap(i,:));
hold on
end
ylabel('Trigger pulse duration')

legend(cfg.trialdef.conditionlabel,'Location','best','AutoUpdate','off')
yyaxis right
plot(trigdata,'r')
hold on
title('Trigger duration detection')
xlabel('Time (samples)')
ylabel('raw trigger channel')

