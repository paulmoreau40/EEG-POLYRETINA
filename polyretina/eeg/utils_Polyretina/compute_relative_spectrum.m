function [EEG_relative_spectrum] = compute_relative_spectrum(EEG_trial_data, EEG_baseline_data, make_plot)
% Custom function to compute relative spectra for every trial and
% corresponding baselines computed previously
%
% Inputs:
% EEG_trial_data           - Struct with trial data spectrum
% EEG_baseline_data        - Struct with corresponding baseline data
%                            spectrum
% make_plot                - whether or not to make spectoplot
%
% Outputs:
% EEG_relative_spectrum    = Struct with relative spectrum computed for
%                             every spectrum


% 1. Creating empty structure
EEG_relative_spectrum = [];
EEG_relative_spectrum.metaInfo = EEG_trial_data.metaInfo;

% Getting overall number of trials to compute:
total_num_trials = size({EEG_relative_spectrum.metaInfo.trial_id},2);

% 2. Looping over every trial computed

for trial = 1:total_num_trials
    
    % Get participant_id and trial_id in question
    participant_id = EEG_relative_spectrum.metaInfo(trial).participant_id;
    trial_id = EEG_relative_spectrum.metaInfo(trial).trial_id;
    
    % First Sanity check: verify that frequency vectors are the same for
    % trial and corresponding baseline data
    if EEG_trial_data.(['P' participant_id]).freqs ~= EEG_baseline_data.(['P' participant_id]).freqs
        error('Frequencies do not match up for EEG_trial_data and EEG_baseline_data... recompute with same spectopo parameters');
    else
        EEG_relative_spectrum.(['P' participant_id]).freqs = EEG_trial_data.(['P' participant_id]).freqs;
    end
        
    if length(EEG_trial_data.(['P' participant_id]).spectrum(:,:,trial_id)) ~= length(EEG_baseline_data.(['P' participant_id]).spectrum(:,:,trial_id))
        error('Size of of spectra not compatible between EEG trial and EEG baseline')
    else
        EEG_relative_spectrum.(['P' participant_id]).relative_spectrum(:,:,trial_id) = EEG_trial_data.(['P' participant_id]).spectrum(:,:,trial_id) ./ EEG_baseline_data.(['P' participant_id]).spectrum(:,:,trial_id);
        EEG_relative_spectrum.(['P' participant_id]).srate = EEG_trial_data.(['P' participant_id]).srate;
        EEG_relative_spectrum.(['P' participant_id]).chanlocs = EEG_trial_data.(['P' participant_id]).chanlocs;
    end
    
end

disp('Finished computing relative spectrum for all trials and all participants...')



end

