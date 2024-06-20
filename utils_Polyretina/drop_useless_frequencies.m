function [EEG_relative_spectrum] = drop_useless_frequencies(EEG_relative_spectrum, low_freq, high_freq, absolute_or_relative)
% Custom function to remove frequencies which are out of interest
%
% Inputs:
% EEG_relative_spectrum    - Struct with relative spectrum, where
%                            frequencies need to be removed
% low_freq                 - Cutoff frequency for high_pass filter
% high_freq                - Cutoff frequency for low_pass filter
% absolute_or_relative     - String to identify if the field is relative
%                            spectrum or absolute spectrum
%
% Outputs:
% EEG_relative_spectrum    = Struct with relative spectrum computed for
%                             every spectrum and only frequencies of
%                             interest

participants = unique({EEG_relative_spectrum.metaInfo(:).participant_id});

% Loop over every participant
for p = 1:length(participants)
    
    % Retrieve indicies of frequencies which are out of the scope of interest (meaning beyond filtering frequencies)
    OoI_frequencies_indices = find((EEG_relative_spectrum.(['P' participants{p}]).freqs < low_freq) + (EEG_relative_spectrum.(['P' participants{p}]).freqs > high_freq));
    
    if ~isempty(OoI_frequencies_indices)
        if strcmp(absolute_or_relative, 'absolute')
            EEG_relative_spectrum.(['P' participants{p}]).spectrum(:,OoI_frequencies_indices,:) = [];
            EEG_relative_spectrum.(['P' participants{p}]).freqs(OoI_frequencies_indices) = [];
        elseif strcmp(absolute_or_relative, 'relative')
            EEG_relative_spectrum.(['P' participants{p}]).relative_spectrum(:,OoI_frequencies_indices,:) = [];
            EEG_relative_spectrum.(['P' participants{p}]).freqs(OoI_frequencies_indices) = [];
        else
            error('Error in entering whether "absolute" or "relative"');
        end
    end
end



end

