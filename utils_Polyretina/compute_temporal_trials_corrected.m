function [EEG_selected_spectrum_FoV] = compute_temporal_trials_corrected(EEG_trial_data, EEG_baseline_data, electrodes_of_interest, bool_corrected)




baseline_means = struct();
participant_ids = unique([EEG_baseline_data.metaInfo.participant_id]);

for i = 1:length(participant_ids)    
    participant_data = EEG_baseline_data.(participant_ids(i)).data2sec;
    mean_baseline = mean(participant_data, 2);  % Average across time points (dimension 2)
    baseline_means.(participant_ids(i)) = mean_baseline;
end







% Step 1: Precompute the baseline mean for each block

baseline_means = struct();

for participant = 1:num_participants
    participant_id = EEG_baseline_data.metaInfo(participant).participant_id;
    
    % Initialize a structure to store baseline means for each block of this participant
    baseline_means.(participant_id) = struct();
    
    % Loop through each block to compute the mean baseline for that block
    blocks = unique([EEG_baseline_data.metaInfo.BlockIndex]);
    
    for block = blocks
        % Get the baseline data for the current block
        baseline_trials = EEG_baseline_data.(participant_id).data(:,:, EEG_baseline_data.metaInfo.BlockIndex == block);
        
        % Compute the mean baseline across all electrodes and time points for this block
        mean_baseline = mean(baseline_trials, 3);  % Average across time points and electrodes
        baseline_means.(participant_id).(['block_' num2str(block)]) = mean_baseline;
    end
end

% Step 2: Use the precomputed baseline means to normalize the trials

EEG_corrected_trials = struct();

for participant = 1:num_participants
    participant_id = EEG_trial_data.metaInfo(participant).participant_id;
    
    % Loop through each trial for this participant
    for trial = 1:num_trials
        block_idx = EEG_trial_data.metaInfo(trial).BlockIndex;
        
        % Retrieve the precomputed baseline mean for this block
        mean_baseline = baseline_means.(participant_id).(['block_' num2str(block_idx)]);
        
        % Normalize the current trial using the mean baseline for this block
        trial_data = EEG_trial_data.(participant_id).data(:,:,trial);
        corrected_trial = trial_data - mean_baseline;
        
        % Store the corrected trial
        EEG_corrected_trials.(participant_id).data(:,:,trial) = corrected_trial;
    end
end













EEG_selected_spectrum_FoV = [];

total_num_trials = size({EEG_relative_spectrum.metaInfo.participant_id},2);

% Defining variables to know if we moved onto another participant
previous_participant = EEG_relative_spectrum.metaInfo(1).participant_id;
count_trial = 1;
count_trial_bis = 1;
count_meta = 1;

if bool_divide_by_FoV %  TRIALS AND BASELINE
    disp(['Keeping only trials of interst for FoV: ' num2str(wanted_FoV) '...'])
    % 2. Loop over all the trials to keep only with field of view of interest
    for trial = 1:total_num_trials

        % Seeing if we moved onto next participant
        current_participant = EEG_relative_spectrum.metaInfo(trial).participant_id;
        if ~strcmp(current_participant, previous_participant)
            count_trial = 1;
            count_trial_bis=1;
        end

        % Check if the field of view is the one wanted
        current_FoV = EEG_relative_spectrum.metaInfo(trial).FieldOfView;
        if current_FoV == wanted_FoV

            trial_id = EEG_relative_spectrum.metaInfo(trial).TrialIndex;

            % Saving spectrum and related information for that
            % participant/trial
            EEG_selected_spectrum_FoV.(current_participant).srate = EEG_relative_spectrum.(current_participant).srate;
            EEG_selected_spectrum_FoV.(current_participant).chanlocs = EEG_relative_spectrum.(current_participant).chanlocs;
            EEG_selected_spectrum_FoV.(current_participant).freqs = EEG_relative_spectrum.(current_participant).freqs;

            switch absolute_or_relative
                case 'absolute'
                    EEG_selected_spectrum_FoV.(current_participant).relative_spectrum(:,:,count_trial) = EEG_relative_spectrum.(current_participant).spectrum(:,:,count_trial_bis); % POL : before, trial_idx
                % case 'relative'
                %     EEG_selected_spectrum_FoV.(current_participant).relative_spectrum(:,:,count_trial) = EEG_relative_spectrum.(current_participant).relative_spectrum(:,:,count_trial_bis); % POL : before, trial_idx
                % otherwise
                %     error('Type of spectra computed incorrect: type either "absolute" or "relative"');
            end

            EEG_selected_spectrum_FoV.metaInfo(count_meta).participant_id = current_participant;
            EEG_selected_spectrum_FoV.metaInfo(count_meta).TrialIndex = trial_id;
            EEG_selected_spectrum_FoV.metaInfo(count_meta).FieldOfView = current_FoV;
            
            % Updating counting variables
            count_trial = count_trial + 1;
            count_meta = count_meta + 1;
        end
        previous_participant = current_participant;
        count_trial_bis = count_trial_bis + 1;
    end
    
else




















disp('Keeping only electrodes of brain region of interest...')
if ischar(electrodes_of_interest)
    if strcmp(electrodes_of_interest, 'all')
        disp('Keeping all of the electrodes...')
    else
        error('Please enter correct option')
    end
else
    % Looping over all participants
    for p = 1:length(participants)
        % Retrieving chanlocs information for that participant (its the same across all participants but still computed out of precaution...)
        chanlocs = EEG_selected_spectrum_FoV.(current_participant).chanlocs;

        % Retrieve the indices that correspond to electrodes that we don't want to keep
        OoI_electrode_indices = find(~ismember({chanlocs.labels}, electrodes_of_interest));

        % Removing all data from electrodes that are out of interst (OoI)
        EEG_selected_spectrum_FoV.(participants{p}).relative_spectrum(OoI_electrode_indices,:,:) = [];

    end
end