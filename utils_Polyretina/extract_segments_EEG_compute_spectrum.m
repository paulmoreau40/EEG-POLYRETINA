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


if ismember(participant_id,'P001') % dont know why but Coarse baseline present
    EEG.event(end) = [];
    EEG.event(253).latency = 410800; % instead of 411461 to have more baseline than just 125
end

if ismember(participant_id,'P008') % dont know why but Coarse baseline present
    EEG.event(347).latency = 582817;
    EEG.event(370).latency = 620000;
end


n_tot_trials = EEG.event(end).BlockIndex * EEG.event(end).TrialIndex;
blockLength = 3;
n_blocks = 70;
i_trial = 1;

FoV = zeros(n_blocks,1);

for trial = 1:n_tot_trials
    % Resetting values to identify if missing tag:
    start_idx = -1;
    end_idx = -1;
    isBaseline = 0;
 
    bl = ceil(trial/blockLength);
    tr_inBl = trial - (bl-1)*blockLength;
    trial_idx = find([EEG.event.BlockIndex] == bl & [EEG.event.TrialIndex] == tr_inBl);
    
    if length(trial_idx) ~= 2
        error('This trial has ' + trial_idx + " events")
    elseif ismember(EEG.event(trial_idx(end)).TrialType, "Baseline")
        isBaseline = 1;
    else
        % CHANGE IN THE FUTURE: the angle is not written in TrialType for
        % baselines
        FoV(bl) = str2double(regexp(EEG.event(trial_idx(end)).TrialType, '\d+', 'match'));
    end

    % Find start and end idx (Trial / Baseline Start/End) for this trial
    for i = 1:length(trial_idx)
        if ismember(EEG.event(trial_idx(i)).type, [begin_segment_trial, begin_segment_baseline])
            start_idx = trial_idx(i);
        elseif ismember(EEG.event(trial_idx(i)).type, [end_segment_trial, end_segment_baseline])
            end_idx = trial_idx(i);
        else
            error("Didn't find any start or end event for this trial")
        end
    end

    if any([start_idx, end_idx] == -1)
        error("The start or end trial has not been updated/found")
    end

    data = EEG.data(:,round(EEG.event(start_idx).latency):round(EEG.event(end_idx).latency));
    
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
    if trial == 347
        trial;
    end
end

% COMPLETE FOV FOR BASELINES

if length(FoV) ~= 70
    error("Not 70 blocks")
end

for i=1:length(FoV)
    EEG_baseline_data.metaInfo(current_line_baseline + i).FieldOfView = FoV(i);
end

% FoV = num2cell(FoV);
% [EEG_baseline_data.metaInfo.FoV] = deal(FoV{:});


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


