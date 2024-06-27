function [EEG_selected_spectrum_FoV] = extract_trials_according_to_brainregion_and_frequency(EEG_relative_spectrum, electrodes_of_interest, range_freqs_of_interest, bool_divide_by_FoV, wanted_FoV, absolute_or_relative)
% Custom function to extract segments of interest and return a structure with the segmented regions of interest and associated metadata
%
% Inputs:
% EEG_relative_spectrum       - Struct with relative spectra for all participants 
%                               and all conditions of FoV and all electrodes over all
%                               frequencies 
% electrodes_of_interest      - Defining electrodes of interest
% range_freqs_of_interest     - Defining frequencies of interest
% bool_divide_by_FoV          - Separate by field of view
% FoV                         - Field of View of Interest
% absolute_or_relative        - Whether the input structure has absolute spectrum or relative spectrum
%
%
% Outputs:
% EEG_selected_spectrum_FoV            = Struct with relative spectra computed for field of view of 
%                                        interest over the brain region of interest


% 1. Defining  structures needed and other need variables
EEG_selected_spectrum_FoV = [];
kept_participant_id = [];
kept_trial_id = [];
kept_field_of_view = [];
% kept_auditory_cue = [];
% kept_auditory_answer = [];

% Getting overall number of trials to compute:
total_num_trials = size({EEG_relative_spectrum.metaInfo.TrialIndex},2);
%total_num_trials = EEG_relative_spectrum.metaInfo(end).BlockIndex * 2;

% Defining variables to know if we moved onto another participant
previous_participant = EEG_relative_spectrum.metaInfo(1).participant_id;
count_trial = 1;
count_trial_bis = 1;

if bool_divide_by_FoV
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

            % Get trial_id
            trial_id = EEG_relative_spectrum.metaInfo(trial).TrialIndex;
            % auditory_cue = EEG_relative_spectrum.metaInfo(trial).auditory_cue;
            % auditory_answer = EEG_relative_spectrum.metaInfo(trial).auditory_answer;

            % Save metaInfo that is kept in new structure
            kept_participant_id = [kept_participant_id; current_participant];
            kept_trial_id = [kept_trial_id; trial_id];
            kept_field_of_view = [kept_field_of_view; current_FoV];
            % kept_auditory_cue = [kept_auditory_cue; auditory_cue];
            % kept_auditory_answer = [kept_auditory_answer; auditory_answer];

            % Saving spectrum and related information for that
            % participant/trial
            EEG_selected_spectrum_FoV.(current_participant).srate = EEG_relative_spectrum.(current_participant).srate;
            EEG_selected_spectrum_FoV.(current_participant).chanlocs = EEG_relative_spectrum.(current_participant).chanlocs;
            EEG_selected_spectrum_FoV.(current_participant).freqs = EEG_relative_spectrum.(current_participant).freqs;

            switch absolute_or_relative
                case 'absolute'
                    EEG_selected_spectrum_FoV.(current_participant).relative_spectrum(:,:,count_trial) = EEG_relative_spectrum.(current_participant).spectrum(:,:,count_trial_bis); % POL : before, trial_idx
                case 'relative'
                    EEG_selected_spectrum_FoV.(current_participant).relative_spectrum(:,:,count_trial) = EEG_relative_spectrum.(current_participant).relative_spectrum(:,:,count_trial_bis); % POL : before, trial_idx
                otherwise
                    error('Type of spectra computed incorrect: type either "absolute" or "relative"');
            end

            EEG_selected_spectrum_FoV.metaInfo(count_trial).participant_id = current_participant;
            EEG_selected_spectrum_FoV.metaInfo(count_trial).TrialIndex = trial_id;
            EEG_selected_spectrum_FoV.metaInfo(count_trial).FieldOfView = current_FoV;
            

            % Updating counting variables
            count_trial = count_trial + 1;
            previous_participant = current_participant;
        end
        count_trial_bis = count_trial_bis + 1;
    end
    
else
    disp('Consider all the Fields of View:...')
    % 2. Loop over all the trials
    for trial = 1:total_num_trials

        % Seeing if we moved onto next participant
        current_participant = EEG_relative_spectrum.metaInfo(trial).participant_id;
        if ~strcmp(current_participant, previous_participant)
            count_trial = 1;
        end

        % Retrieving current FoV
        current_FoV = EEG_relative_spectrum.metaInfo(trial).FieldOfView;

        % Get trial_id
        trial_id = EEG_relative_spectrum.metaInfo(trial).TrialIndex;
        % auditory_cue = EEG_relative_spectrum.metaInfo(trial).auditory_cue;
        % auditory_answer = EEG_relative_spectrum.metaInfo(trial).auditory_answer;

        % Save metaInfo that is kept in new structure
        kept_participant_id = [kept_participant_id; current_participant];
        kept_trial_id = [kept_trial_id; trial_id];
        kept_field_of_view = [kept_field_of_view; current_FoV];
        % kept_auditory_cue = [kept_auditory_cue; auditory_cue];
        % kept_auditory_answer = [kept_auditory_answer; auditory_answer];

        % Saving spectrum and related information for that
        % participant/trial
        EEG_selected_spectrum_FoV.(current_participant).srate = EEG_relative_spectrum.(current_participant).srate;
        EEG_selected_spectrum_FoV.(current_participant).chanlocs = EEG_relative_spectrum.(current_participant).chanlocs;
        EEG_selected_spectrum_FoV.(current_participant).freqs = EEG_relative_spectrum.(current_participant).freqs;

        switch absolute_or_relative
            case 'absolute'
                EEG_selected_spectrum_FoV.(current_participant).relative_spectrum(:,:,count_trial) = EEG_relative_spectrum.(current_participant).spectrum(:,:,trial); % POL : before, trial_idx
            case 'relative'
                EEG_selected_spectrum_FoV.(current_participant).relative_spectrum(:,:,count_trial) = EEG_relative_spectrum.(current_participant).relative_spectrum(:,:,trial); % POL : before, trial_idx
            otherwise
                error('Type of spectra computed incorrect: type either "absolute" or "relative"');
        end


        EEG_selected_spectrum_FoV.metaInfo(count_trial).participant_id = current_participant;
        EEG_selected_spectrum_FoV.metaInfo(count_trial).TrialIndex = trial_id;
        EEG_selected_spectrum_FoV.metaInfo(count_trial).FieldOfView = current_FoV;
        
        % Updating counting variables
        count_trial = count_trial + 1;
        previous_participant = current_participant;

    end
end

% Saving all of the metadata information in the structure
% EEG_selected_spectrum_FoV.metaInfo = struct('participant_id', kept_participant_id,...
%                 'TrialIndex', kept_trial_id, 'FieldOfView', kept_field_of_view); %, 'auditory_cue', kept_auditory_cue, 'auditory_answer', kept_auditory_answer

disp('... finished') 

% Get participants which are kept:
participants = unique({EEG_relative_spectrum.metaInfo(:).participant_id});

% 3. Keep only frequencies of interset:
if ischar(range_freqs_of_interest)
    if strcmp(range_freqs_of_interest, 'all')
        disp('Keeping all of the frequencies...')
    else
        error('Please enter correct range of frequency')
    end
elseif length(range_freqs_of_interest) ~= 2
    error('Range of frequencies given is incorrect: please enter [min_freq max_freq]')
else
    disp(['Keeping only frequencies in range [' num2str(range_freqs_of_interest(1)) ', ' num2str(range_freqs_of_interest(2)) '] ...'])
    % Loop over all of those participants:
    for p = 1:length(participants)

        % Retrieve indicies of frequencies which are out of the scope of
        % interest (meaning beyond filtering frequencies)
        OoI_frequencies_indices = find((EEG_relative_spectrum.(participants{p}).freqs < range_freqs_of_interest(1)) + (EEG_relative_spectrum.(participants{p}).freqs > range_freqs_of_interest(end)));

        EEG_selected_spectrum_FoV.(participants{p}).relative_spectrum(:,OoI_frequencies_indices,:) = [];
        EEG_selected_spectrum_FoV.(participants{p}).freqs(OoI_frequencies_indices) = [];
    end
    
end
disp('...finished')


% 4. Keep only brain regions of interest

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
disp('...finished')

            
end
