function [nevents] = ft_definetrial_events_MYH(cfg)

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
% nevents: the time point of the onset of each events based on the audio
% channel triggers
%__________________________________________________________________________
% Adapted code from the original gangsters Vladimir Litvak, Robert Oostenveld
% Thusly edited by by NBara, then remixed by Rosy Southwell in Juuune 2018:
% - tolerance for trig durations to +- 3 samples
% - allow a different trialshift per condition
% - allow a different trialshift per trial (cfg.trialdef.trlshiftpertrial)
% - ignore trigger values / conditions which arent passed in the trialdef, don't include
% these in cfg.event
%% --------------------------------------------------------------------------
% In this script, we just need to know the time point of the onset of the last event, based
% on the triggers, for the following detrending purpose
%-----------this script was adapted by Mingyue Hu, 07/07/2022

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

end
