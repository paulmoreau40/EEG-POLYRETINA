function [spectrum_data_20,spectrum_data_45, spectrum_data_110] = format_data_for_multisubject_stats(EEG_data, freq_boundaries)
% Computing and formating the spectrums to give to permutation analysis
%
%  INPUT
% EEG_data                  - Input data containing the spectra for the
%                             participants
% freq_boundaries           - Boundaries of frequencies used for filtering
%                             and which are used as limitations for permutation anaysis on frequency
%                             domain
%  OUTPUT
% spectrum_data_20           - Data formated for permutation for FoV 20
% spectrum_data_45           - Data formated for permutation for FoV 450
% spectrum_data_110          - Data formated for permutation for FoV 110

% Reorganizing data into 3 structures and separating by field of view, as
% well as removing frequencies that are not of interst
[EEG_selected_spectrum_20] = extract_trials_according_to_brainregion_and_frequency(EEG_data, 'all', freq_boundaries, 1, 20, 'absolute');
[EEG_selected_spectrum_45] = extract_trials_according_to_brainregion_and_frequency(EEG_data, 'all', freq_boundaries, 1, 45, 'absolute');
[EEG_selected_spectrum_110] = extract_trials_according_to_brainregion_and_frequency(EEG_data, 'all', freq_boundaries, 1, 110, 'absolute');

% Now averaging over every trial per participant to get final structure of
% interest

spectrum_data_20 = [];
spectrum_data_45 = [];
spectrum_data_110 = [];

% Getting list of participants considered:
participant_list = unique({EEG_data.metaInfo(:).participant_id});

% Counting how many subjects to concatenate structures
subject_count = 1;

for participant = 1:length(participant_list)
    % Retrieving for each FoV the spectrums and averaging over all trials
    spectrum_data_20(:,:,subject_count) = mean(EEG_selected_spectrum_20.(['P' participant_list{participant}]).relative_spectrum, 3);
    spectrum_data_45(:,:,subject_count) = mean(EEG_selected_spectrum_45.(['P' participant_list{participant}]).relative_spectrum, 3);
    spectrum_data_110(:,:,subject_count) = mean(EEG_selected_spectrum_110.(['P' participant_list{participant}]).relative_spectrum, 3);
    
     subject_count = subject_count + 1;
end

end

