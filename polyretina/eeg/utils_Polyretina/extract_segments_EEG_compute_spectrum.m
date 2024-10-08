function [EEG_trial_data, EEG_baseline_data, EEG_coarse_data] = extract_segments_EEG_compute_spectrum(EEG, EEG_trial_data, EEG_baseline_data, EEG_coarse_data, participant_id, make_plot)
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

current_line_trial = size({EEG_trial_data.metaInfo.TrialIndex},2);
current_line_baseline = size({EEG_baseline_data.metaInfo.TrialIndex},2);
current_line_coarse = size({EEG_coarse_data.metaInfo.idx},2);

spectrumTrials = []; spectrumBaseline = []; spectrumCoarse = [];
dataTrials = []; data20Trials = []; data45Trials = []; 
dataBaselines = []; data20Baselines = []; data45Baselines = [];
events = EEG.event;

if ismember(participant_id,'P001') % dont know why but Coarse baseline present
    events(end) = [];
    events(257).latency = events(257-1).latency; % instead of 411461 to have more baseline than just 120
end

if ismember(participant_id,'P008')
    events(351).latency = events(351-1).latency;
    events(374).latency = events(374).latency +700;
end

if ismember(participant_id,'P007')
    events(47).latency = 75546; % baseline -> 859 samples (instead of 76)
end



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
FoV = []; % save the angle value to update later the baselines angles

for trial = 1:n_tot_trials
    % Resetting values to identify if missing tag:
    start_idx = -1;
    end_idx = -1;
    isBaseline = 0;
 
    bl = ceil(trial/blockLength);
    tr_inBl = trial - (bl-1)*blockLength;
    trial_idx = find([events.BlockIndex] == bl & [events.TrialIndex] == tr_inBl);

    plot_option = 'off';
    if make_plot
        close all; % To avoid spectopo bugs
        plot_option = 'on';
    end
    
    if length(trial_idx) ~= 2
        error('This trial has ' + trial_idx + " events")
    elseif ismember(events(trial_idx(end)).TrialType, "Baseline")
        isBaseline = 1;
    else
        FoV(bl) = str2double(regexp(events(trial_idx(end)).TrialType, '\d+', 'match'));
    end

    % Find start and end idx (Trial / Baseline Start/End) for this trial
    for i = 1:length(trial_idx)
        if ismember(events(trial_idx(i)).type, [{'TrialStart'}, {'BaselineStart'}])
            start_idx = trial_idx(i);
        elseif ismember(events(trial_idx(i)).type, [{'TrialEnd'}, {'BaselineEnd'}])
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


    % EXTRACT THE LAST TWO SECONDS OF EACH TRIAL / BASELINE
    n_samples = 2 * EEG.srate;
    data_short = EEG.data(:,round(events(end_idx).latency) - n_samples:round(events(end_idx).latency));

    % COMPUTE POWER SPECTRUM
    [spectrum, freqs, ~, ~, ~] = spectopo(data, size(data, 2), EEG.srate, ...
        'freq', [10.0, 20.0], ...
        'chanlocs', EEG.chanlocs, ...
        'freqfac', 2, ...
        'winsize', 250, ...
        'overlap', 125, ...   
        'wintype', 'hamming', ...
        'freqrange', [2 40], ...
        'plot', plot_option, ...
        'plotmean', 'off', ...
        'verbose', 'off');
    
    if isBaseline
        EEG_baseline_data.metaInfo(current_line_baseline + bl).participant_id = participant_id;
        EEG_baseline_data.metaInfo(current_line_baseline + bl).BlockIndex = bl;
        EEG_baseline_data.metaInfo(current_line_baseline + bl).TrialIndex = tr_inBl;
    
        spectrumBaseline(:, :, bl) = spectrum;
        dataBaselines(:,:,bl) = data_short;

    else
        EEG_trial_data.metaInfo(current_line_trial + i_trial).participant_id = participant_id;
        EEG_trial_data.metaInfo(current_line_trial + i_trial).BlockIndex = bl;
        EEG_trial_data.metaInfo(current_line_trial + i_trial).TrialIndex = tr_inBl;
        EEG_trial_data.metaInfo(current_line_trial + i_trial).FieldOfView = FoV(bl);
    
        spectrumTrials(:, :, i_trial) = spectrum;
        dataTrials(:,:,i_trial) = data_short;
        i_trial = i_trial + 1;
    end

    disp(['Finished computing spectrum for trial/baseline ' num2str(trial) '/' num2str(n_tot_trials) ' for Participant ' participant_id])
end





% POWER SPECTRUM OF COARSE BASELINES
BaseCoarse_idx = find(strcmp({events.type},'BaselineCoarseStart'));

for i=1:length(BaseCoarse_idx)
    if ~ismember(events(BaseCoarse_idx(i)+1).type,'BaselineCoarseEnd')
        error("The BaselineCoarseStart is not followed by a BaselineCoarseEnd")
    end

    data = EEG.data(:,round(events(BaseCoarse_idx(i)).latency):round(events(BaseCoarse_idx(i)+1).latency));
    
    [spectrumCoarse(:, :, i), freqs, ~, ~, ~] = spectopo(data, size(data, 2), EEG.srate, ...
        'freq', [10.0, 20.0], ...
        'chanlocs', EEG.chanlocs, ...
        'freqfac', 2, ...
        'winsize', 250, ...
        'overlap', 125, ...   
        'wintype', 'hamming', ...
        'freqrange', [2 40], ...
        'plot', plot_option, ...
        'plotmean', 'off', ...
        'verbose', 'off');

    EEG_coarse_data.metaInfo(current_line_coarse + i).participant_id = participant_id;
    EEG_coarse_data.metaInfo(current_line_coarse + i).idx = i;

    if ismember(i, [1, 3])
        EEG_coarse_data.metaInfo(current_line_coarse + i).type = "black";
    else
        EEG_coarse_data.metaInfo(current_line_coarse + i).type = "edge";
    end

    disp(['Finished computing spectrum for coarse baseline ' num2str(i) '/' num2str(length(BaseCoarse_idx)) ' for Participant ' participant_id])
end
    
    




    




EEG.event = events;

% EEG TRIALS
EEG_trial_data.(participant_id).srate = EEG.srate;
EEG_trial_data.(participant_id).chanlocs = EEG.chanlocs;
EEG_trial_data.(participant_id).freqs = freqs;
EEG_trial_data.(participant_id).spectrum = 10.^(spectrumTrials./10);
EEG_trial_data.(participant_id).data2sec = dataTrials;

% EEG BASELINES
EEG_baseline_data.(participant_id).srate = EEG.srate;
EEG_baseline_data.(participant_id).chanlocs = EEG.chanlocs;
EEG_baseline_data.(participant_id).freqs = freqs;
EEG_baseline_data.(participant_id).spectrum = 10.^(spectrumBaseline./10);
EEG_baseline_data.(participant_id).data2sec = dataBaselines;


len_fov = length(FoV);
if ismember(participant_id,'P009')
    len_fov = len_fov - 1; 
end
for i=1:len_fov
    EEG_baseline_data.metaInfo(current_line_baseline + i).FieldOfView = FoV(i);
end

% EEG COARSE BASELINES
EEG_coarse_data.(participant_id).srate = EEG.srate;
EEG_coarse_data.(participant_id).chanlocs = EEG.chanlocs;
EEG_coarse_data.(participant_id).freqs = freqs;
EEG_coarse_data.(participant_id).spectrum = 10.^(spectrumCoarse./10);

end
