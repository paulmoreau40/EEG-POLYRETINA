function [EEG_baseline_data] = computing_singleBaseline(EEG_baseline_data, EEG_coarse_data)

% Getting overall number of trials to compute:
total_num_trials = size({EEG_baseline_data.metaInfo.trial_id},2);

EEG_spectrum = [];

% Looping over every trial
for trial = 1:total_num_trials
    
    % Retrieving participant information
    participant_id = EEG_baseline_data.metaInfo(trial).participant_id; 
    current_trial = EEG_baseline_data.metaInfo(trial).trial_id; 
    
    % Save spectrum in overall matrix
    EEG_spectrum = EEG_baseline_data.(participant_id).spectrum(:,:,current_trial);
    
end

% Taking average of the spectrum
mean_EEG_spectrum = mean(EEG_spectrum,3);
EEG_baseline_data.spectrum_std = std(EEG_spectrum, [], 3)/sqrt(size(EEG_spectrum, 3));

% Re-iterating over all of the participants and redefining their baseline
for trial = 1:total_num_trials
    
    % Retrieving participant information
    participant_id = EEG_baseline_data.metaInfo(trial).participant_id; 
    current_trial = EEG_baseline_data.metaInfo(trial).trial_id; 

    % Set to averaged baseline
    EEG_baseline_data.(participant_id).spectrum(:,:,current_trial) = mean_EEG_spectrum;
    
end


end
