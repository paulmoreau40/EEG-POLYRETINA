function [EEG_trial_data, EEG_baseline_data] = extract_segments_EEG_compute_spectrum(EEG, EEG_trial_data, EEG_baseline_data, participant_id, make_plot)
% Custom function to extract baselines of interest and return a structure with the segmented regions of interest and associated metadata
%
% Inputs:
% EEG           - EEG struct where Button Press doublons need to be removed
% segment_of_interst        - Defining what segment in time is of interest
%                             to us
% participant_id            - id of participant
% make_plot                 - whether or not to make spectoplot
%
% Outputs:
% EEG_trial_data            = Struct with EEG segments and their respective
%                             time trials and metadata

% TRIALS
begin_segment_trial = {'TrialStart'};
end_segment_trial = {'TrialEnd'};
current_line_trial = size({EEG_trial_data.metaInfo.TrialIndex},2);
spectrumTrials = [];

% BASELINES
begin_segment_baseline = {'BaselineStart'};
end_segment_baseline = {'BaselineEnd'};
current_line_baseline = size({EEG_baseline_data.metaInfo.TrialIndex},2);
spectrumBaseline = [];
events = EEG.event;

if ismember(participant_id,'P001') % dont know why but Coarse baseline present
    events(end) = [];
    events(257).latency = events(257-1).latency; % instead of 411461 to have more baseline than just 125
end

% if ismember(participant_id,'P008')
%     events(347).latency = 582817;
%     events(370).latency = 620000;
% end
% 
% if ismember(participant_id,'P007')
%     events(43).latency = 72833; % baseline -> 859 samples (instead of 76)
% end



% FIND LAST NON-ZERO VALUE OF BLOCK AND INDEX TO COMPUTE N TOT TRIALS
lastBlock = []; lastIndex = [];
for i = length(events):-1:1
    if events(i).BlockIndex ~= 0
        lastBlock = events(i).BlockIndex;
        break;
    end
end
for i = length(events):-1:1
    if events(i).TrialIndex ~= 0
        lastIndex = events(i).TrialIndex;
        break;
    end
end
n_tot_trials = lastBlock * lastIndex;

blockLength = 3;
i_trial = 1;
FoV = [];

for trial = 1:n_tot_trials
    % Resetting values to identify if missing tag:
    start_idx = -1;
    end_idx = -1;
    isBaseline = 0;
 
    bl = ceil(trial/blockLength);
    tr_inBl = trial - (bl-1)*blockLength;
    trial_idx = find([events.BlockIndex] == bl & [events.TrialIndex] == tr_inBl);
    
    if length(trial_idx) ~= 2
        error('This trial has ' + trial_idx + " events")
    elseif ismember(events(trial_idx(end)).TrialType, "Baseline")
        isBaseline = 1;
    else
        % save the angle value to update later the baselines angles
        FoV(bl) = str2double(regexp(events(trial_idx(end)).TrialType, '\d+', 'match'));
    end

    % Find start and end idx (Trial / Baseline Start/End) for this trial
    for i = 1:length(trial_idx)
        if ismember(events(trial_idx(i)).type, [begin_segment_trial, begin_segment_baseline])
            start_idx = trial_idx(i);
        elseif ismember(events(trial_idx(i)).type, [end_segment_trial, end_segment_baseline])
            end_idx = trial_idx(i);
        else
            error("Didn't find any start or end event for this trial")
        end
    end

    if any([start_idx, end_idx] == -1)
        error("The start or end trial has not been updated/found")
    end

    data = EEG.data(:,round(events(start_idx).latency):round(events(end_idx).latency));
    
    if size(data,2) < 500
        warning('data too small for the power spectrum computation ? (< 500)')
    end

    if make_plot
        close all; % to avoid spectopo bugs
        plot_option = 'on';
    else
        plot_option = 'off';
    end

    if isBaseline
        EEG_baseline_data.metaInfo(current_line_baseline + bl).participant_id = participant_id;
        EEG_baseline_data.metaInfo(current_line_baseline + bl).BlockIndex = bl;
        EEG_baseline_data.metaInfo(current_line_baseline + bl).TrialIndex = tr_inBl;

        [spectrumBaseline(:,:,bl), freqs, ~, ~, ~] =...
        spectopo(data, size(data,2), EEG.srate,...
        'freq', [10.0, 20.0], ... %freq_of_interest, ... %10,...                %'plotchan', [],...
        'chanlocs', EEG.chanlocs,...
        'freqfac', 2,...
        'winsize', 250, ... % fréquence d'échantillonnage
        'overlap', 125,... 
        'wintype','hamming',...
        'freqrange',[2 40], ...
        'plot', plot_option,...
        'plotmean', 'off',...
        'verbose','off');
    else
        EEG_trial_data.metaInfo(current_line_trial + i_trial).participant_id = participant_id;
        EEG_trial_data.metaInfo(current_line_trial + i_trial).BlockIndex = bl;
        EEG_trial_data.metaInfo(current_line_trial + i_trial).TrialIndex = tr_inBl;
        EEG_trial_data.metaInfo(current_line_trial + i_trial).FieldOfView = FoV(bl);

        [spectrumTrials(:,:,i_trial), freqs, ~, ~, ~] =...
        spectopo(data, size(data,2), EEG.srate,...
        'freq', [10.0, 20.0], ... %freq_of_interest, ... %10,...                %'plotchan', [],...
        'chanlocs', EEG.chanlocs,...
        'freqfac', 2,...
        'winsize', 250, ... % fréquence d'échantillonnage
        'overlap', 125,... 
        'wintype','hamming',...
        'freqrange',[2 40], ...
        'plot', plot_option,...
        'plotmean', 'off',...
        'verbose','off');

        i_trial = i_trial + 1;
    end

    disp(['Finished computing spectrum for trial/baseline ' num2str(trial) '/' num2str(n_tot_trials) ' for Participant ' participant_id])
    % if trial == 347
    %     trial;
    % end
end





len_fov = length(FoV);
if ismember(participant_id,'P009')
    len_fov = len_fov - 1;
end
    
for i=1:len_fov
    EEG_baseline_data.metaInfo(current_line_baseline + i).FieldOfView = FoV(i);
end




EEG.event = events;

% Saving frequencie and spectre for participant considered
EEG_trial_data.(participant_id).srate = EEG.srate;
EEG_trial_data.(participant_id).chanlocs = EEG.chanlocs;
EEG_trial_data.(participant_id).freqs = freqs;
EEG_trial_data.(participant_id).spectrum = 10.^(spectrumTrials./10);

EEG_baseline_data.(participant_id).srate = EEG.srate;
EEG_baseline_data.(participant_id).chanlocs = EEG.chanlocs;
EEG_baseline_data.(participant_id).freqs = freqs;
EEG_baseline_data.(participant_id).spectrum = 10.^(spectrumBaseline./10);

end


