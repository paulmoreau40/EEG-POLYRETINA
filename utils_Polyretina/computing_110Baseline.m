function [EEG_baseline_data] = computing_110Baseline(EEG_baseline_data)
% Computing the baseline from the 110° Field of view trials

% Getting overall number of trials to compute:
total_num_trials = size({EEG_baseline_data.metaInfo.trial_id},2);

EEG_spectrum = [];

% Looping over every trial
for trial = 1:total_num_trials
    
    % Retrieving participant information
    participant_id = EEG_baseline_data.metaInfo(trial).participant_id; 
    current_trial = EEG_baseline_data.metaInfo(trial).trial_id; 
    field_of_view = EEG_baseline_data.metaInfo(trial).field_of_view;
    
    % Retrieving spectrum if the field of view is 110° and ignoring if
    % otherwise
    if field_of_view == 110
        EEG_spectrum = EEG_baseline_data.(['P' participant_id]).spectrum(:,:,current_trial);
    end    
end

% Taking average of the respective spectrum
mean_EEG_spectrum = mean(EEG_spectrum,3);
EEG_baseline_data.spectrum_std = std(EEG_spectrum, [], 3)/sqrt(size(EEG_spectrum, 3));

% Re-iterating over all of the participants and redefining their baseline
% as the average of the 110° spectrum
for trial = 1:total_num_trials
    
    % Retrieving participant information
    participant_id = EEG_baseline_data.metaInfo(trial).participant_id; 
    current_trial = EEG_baseline_data.metaInfo(trial).trial_id; 

    % Set to averaged baseline per FOV
    EEG_baseline_data.(['P' participant_id]).spectrum(:,:,current_trial) = mean_EEG_spectrum;
    
end


end

