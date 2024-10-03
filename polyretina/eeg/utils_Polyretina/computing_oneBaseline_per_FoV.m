function [EEG_baseline_data] = computing_oneBaseline_per_FoV(EEG_baseline_data)

% Getting overall number of trials to compute:
total_num_trials = size({EEG_baseline_data.metaInfo.trial_id},2);

EEG_spectrum_20 = [];
EEG_spectrum_45 = [];
% EEG_spectrum_110 = [];

% Looping over every trial
for trial = 1:total_num_trials
    
    % Retrieving participant information
    participant_id = EEG_baseline_data.metaInfo(trial).participant_id; 
    current_trial = EEG_baseline_data.metaInfo(trial).trial_id; 
    field_of_view = EEG_baseline_data.metaInfo(trial).field_of_view;
    
    % According to field of view, save spectrum in correct matrix
    if field_of_view == 20
        EEG_spectrum_20 = EEG_baseline_data.(participant_id).spectrum(:,:,current_trial);
    elseif field_of_view == 45
        EEG_spectrum_45 = EEG_baseline_data.(participant_id).spectrum(:,:,current_trial);
    % elseif field_of_view == 110
    %     EEG_spectrum_110 = EEG_baseline_data.(participant_id).spectrum(:,:,current_trial);
    end    
end

% Taking average of the respective spectrum
mean_EEG_spectrum_20 = mean(EEG_spectrum_20,3);
mean_EEG_spectrum_45 = mean(EEG_spectrum_45,3);
% mean_EEG_spectrum_110 = mean(EEG_spectrum_110,3);

EEG_baseline_data.spectrum_std_20 = std(EEG_spectrum_20, [], 3)/sqrt(size(EEG_spectrum_20, 3));
EEG_baseline_data.spectrum_std_45 = std(EEG_spectrum_45, [], 3)/sqrt(size(EEG_spectrum_45, 3));
% EEG_baseline_data.spectrum_std_110 = std(EEG_spectrum_110, [], 3)/sqrt(size(EEG_spectrum_110, 3));

% Re-iterating over all of the participants and redefining their baseline
for trial = 1:total_num_trials
    
    % Retrieving participant information
    participant_id = EEG_baseline_data.metaInfo(trial).participant_id; 
    current_trial = EEG_baseline_data.metaInfo(trial).trial_id; 
    field_of_view = EEG_baseline_data.metaInfo(trial).field_of_view;

    % According to field of view, set to averaged baseline per FOV
    if field_of_view == 20
        EEG_baseline_data.(participant_id).spectrum(:,:,current_trial) = mean_EEG_spectrum_20;
    elseif field_of_view == 45
        EEG_baseline_data.(participant_id).spectrum(:,:,current_trial) = mean_EEG_spectrum_45;
    % elseif field_of_view == 110
    %     EEG_baseline_data.(participant_id).spectrum(:,:,current_trial) = mean_EEG_spectrum_110;
    end
    
end




end

