function [EEG_selected_spectrum_all_trials] = format_for_plotting_spectra(EEG_selected_spectrum_FoV_Brain_Area_freqs)
% Custom function to format the data and concatenate all trials over all
% participants
%
% Inputs:
% EEG_selected_spectrum_FoV_Brain_Area_freqs       - Structure containing relative spectra, where only the FoV, Brain region and
%                                                    frequencies of interest were kept
%
%
% Outputs:
% EEG_selected_spectrum_all_trials            = Struct with with all information of before, but all trials of every participant concatenated
% freqs_to_plot                               = frequencies of spectrum
%                                    

% Defining needed variables
EEG_selected_spectrum_all_trials = [];

% Retrieving information related to frequency
spectrum_concatenated_all_trials = [];

% Identifying how many participants there are:
participants = unique({EEG_selected_spectrum_FoV_Brain_Area_freqs.metaInfo(:).participant_id});

% Loop over every participant
for p = 1:length(participants)
    
    % For every participant: retrieve trials and place in final structure 
    spectrum_for_this_participant = EEG_selected_spectrum_FoV_Brain_Area_freqs.(participants{p}).relative_spectrum(:,:,:);
    
    % Concatenate it to final structure spectrum:
    spectrum_concatenated_all_trials = cat(3, spectrum_concatenated_all_trials, spectrum_for_this_participant);
    
end

EEG_selected_spectrum_all_trials.metaInfo = EEG_selected_spectrum_FoV_Brain_Area_freqs.metaInfo;
EEG_selected_spectrum_all_trials.srate = EEG_selected_spectrum_FoV_Brain_Area_freqs.(participants{p}).srate; % last participant 
EEG_selected_spectrum_all_trials.chanlocs = EEG_selected_spectrum_FoV_Brain_Area_freqs.(participants{p}).chanlocs;
EEG_selected_spectrum_all_trials.freqs = EEG_selected_spectrum_FoV_Brain_Area_freqs.(participants{p}).freqs;
EEG_selected_spectrum_all_trials.spectrum = spectrum_concatenated_all_trials;

end

