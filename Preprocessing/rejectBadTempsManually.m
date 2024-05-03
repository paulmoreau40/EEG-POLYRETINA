function EEG_manual = rejectBadTempsManually(EEG, cfg)
subject = cfg.subjects(cfg.current_subject).id;

% 1. HP filter
lowcutoff = cfg.filterPreProc.low_cut_off;
highcutoff = cfg.filterPreProc.high_cut_off;
if ~isempty(lowcutoff)
    fprintf('Highpass Filtering (%.1f Hz)...\n', lowcutoff)
end
if ~isempty(highcutoff)
    fprintf('Lowpass Filtering (%.1f Hz)...\n', highcutoff)
end
[EEG_HP] = custom_filter(EEG, lowcutoff, highcutoff);

% Add path to prepPipeline subdirectories if not in the list
tmp = which('getPipelineDefaults');
if isempty(tmp)
    myPath = fileparts(which('prepPipeline'));
    addpath(genpath(myPath));
end

% 2. Remove Line Noise with PREP pipeline functions
disp('Removing Line Noise...')
[EEG_HP_noLN, lineNoiseOut] = removeLineNoise_custom(EEG_HP, cfg.lineNoiseRemoval_method, false);
clear EEG_HP
% If you want to save the filteredEEG_noLN struct with the LineNoiseRemoval information for later:
EEG_HP_noLN.etc.lineNoiseRemoval = lineNoiseOut;

% 3. Mark artefacts manually
EEG_manual = pop_eegplot(EEG_HP_noLN, 1, 0, 1);

EEG_manual.etc.APP.rejectedSamples = ~logical(clean_sample_mask);
end
