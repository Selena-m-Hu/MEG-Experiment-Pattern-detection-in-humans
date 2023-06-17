function [cfg] = ft_definetrial_filMEG_RB(cfg, plt, nevents)

% Definition of trials based on DURATION-BASED triggers
% the triggers are sent as an audio output 
% so that we get very precise timing of the auditory stimulus :-)
% MEG At the FIL:
% the triggers are saved in an analog auxiliary channel
% ~This is how we do at the EI~
% ~H8rs gonna h8~
% 
% FORMAT [trl, conditionlabels, S] = spm_eeg_definetrial(S)
% cfg                - input structure (optional)
% fields of cfg:
% (cfg.header) - hdr structure from earlier preprocessing - can be used to
% get around Fieldtrip errors with reading CTF_res4
%   cfg.data            - MEEG object or filename of M/EEG mat-file
%   cfg.trialdef      - structure array for trial definition with fields (optional)
% cfg.trialdef.prestim // cfg.trialdef.poststim
% - these badboiiis determine how much time in seconds
% you want to cut before and after the trigger 
% 
%      cfg.trialdef.conditionlabel - string label for the condition
%       cfg.trialdef.eventtype      - string (should be 'trigger_up')
%       cfg.trialdef.eventvalue     - trigger value for this condition
%       cfg.trialdef.trlshift       - shift the triggers by a fixed amount (sec)
%                                   (e.g. projector delay). One per
%                                   condition/trigegr type
%   (cfg.trialdef.trlshiftpertrial - what it says on the can 
%                         - specify as a vector with one element per event in your dataset)
%   cfg.filename      - path of BDF file to read
% OUTPUT:
%   trl             - Nx3 matrix [start end offset]
%   conditionlabels - Nx1 cell array of strings, label for each trial
%   cfg               - modified configuration structure (for history)
%__________________________________________________________________________
% Adapted code from the original gangsters Vladimir Litvak, Robert Oostenveld
% Thusly edited by by NBara, then remixed by Rosy Southwell in Juuune 2018:
% - tolerance for trig durations to +- 3 samples
% - allow a different trialshift per condition
% - allow a different trialshift per trial (cfg.trialdef.trlshiftpertrial)
% - ignore trigger values / conditions which arent passed in the trialdef, don't include
% these in cfg.event

%% Parameters
%--------------------------------------------------------------------------
% data = varargin{1};
% read the header, required to determine the stimulus channels and trial specification
if ~isfield(cfg,'header')
    % sometimes there is problem with automatically treading MEG header
    % here - so it is an option to extract the header yourself and include
    % as a field of cfg
 hdr = ft_read_header(cfg.dataset,'headerformat',cfg.headerformat);   
else
    hdr=cfg.header;
end
begsample = 1;
endsample = hdr.nSamples*hdr.nTrials;
dataformat  = [];

% find the TRIGGER channel and read the values from it
trigchani = find(strcmpi(hdr.label,cfg.trigchan));
trigdata = ft_read_data(cfg.dataset,'chanindx', trigchani);


pretrig  = cfg.trialdef.prestim;
posttrig = cfg.trialdef.poststim;

%% Read trigger channel (specifically designed for audio triggers sent to analogue input of fil MEG system!!!)
thresh = 0.1; % binarise trigdata so that all values greater than this proportion are set to 1, otherwise 0
switch cfg.trialdef.eventtype
    case 'trigger_up'
        trigger = trigdata>thresh*max(trigdata);
    case 'trigger_down'
        trigger = trigdata<thresh*min(trigdata);
end

flank_trigger = 0.5*(diff([trigger]));
figure;
plot(trigdata, 'b'); hold on; plot(flank_trigger, 'r'); hold on; plot(0.5*(trigger), 'k');
nevents = find(flank_trigger==0.5);
title(['found events ' num2str(numel(nevents))]);
display(['found events ' num2str(numel(nevents))]);

%% Create event structure
%--------------------------------------------------------------------------
event       = [];
pad         = 0;
trigshift   = 0;
% convert the trigger into an event with a value at a specific sample
for i=find(flank_trigger>0)
    event(end+1).type   = cfg.trialdef.eventtype;        % distinguish between up and down flank
    event(end  ).sample = i + begsample-1;      % assign the sample at which the trigger has gone down
    event(end  ).value  = double(trigger(i+trigshift));      % assign the trigger value just _after_ going up
end

% Sort events in chronological order
[tmp, ind] = sort([event.sample]);
event = event(ind);

% Find distance between consecutive events
distance = find(flank_trigger<0) - find(flank_trigger>0);

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
tol = 4; % RVS changed from 2/ RB change from 4
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
    if ~(strcmp(cfg.trialdef.eventtype,'trigger_up') || strcmp(cfg.trialdef.eventtype,'trigger_down'))
        disp('ERROR: S.trialdef.eventtype should be ''trigger_up'' or ''trigger_down''. Aborting!')
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
if plt
    figure(99); clf;
    plot(trigdata,'k')
    hold on
    samp = [event.sample];
    val = [event.value];
    for i = 1:length(cfg.trialdef.eventvalue)
        plot(samp(val==cfg.trialdef.eventvalue(i)), val(val==cfg.trialdef.eventvalue(i)),'*');
        hold on
    end
    legend([{'raw trigger channel'},cfg.trialdef.conditionlabel],'Location','best')
    title('Trigger duration detection')
    xlabel('Time (samples)')
    ylabel('Trigger pulse duration')
end
