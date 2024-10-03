function EEG_asr = rejectBadTempsWithASR(EEG, cfg, plot)

subject = cfg.subjects(cfg.current_subject).id;

% 1. HP filter
lowcutoff = cfg.filterPreProc.low_cut_off;
highcutoff = cfg.filterPreProc.high_cut_off;
if ~isempty(lowcutoff)
    fprintf('Highpass Filtering (%.1f Hz) for ASR...\n', lowcutoff)
end
if ~isempty(highcutoff)
    fprintf('Lowpass Filtering (%.1f Hz) for ASR...\n', highcutoff)
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

% 3. Reject transient artifacts with ASR
[EEG_asr,~,~] = clean_artifacts(EEG_HP_noLN, 'Highpass', 'off', 'FlatlineCriterion', 'off',...
    'ChannelCriterion', 'off', 'LineNoiseCriterion', 'off', 'BurstCriterion', cfg.burst_crit, 'BurstRejection', 'on');

if plot
    vis_artifacts(EEG_asr,EEG_HP_noLN);
end

EEG_asr.etc.ASR.BurstCrit = cfg.burst_crit;
% Get portions of data which have been rejected:
EEG_asr.etc.ASR.rejectedSamples = ~EEG_asr.etc.clean_sample_mask;

EEG_asr.etc = rmfield(EEG_asr.etc, 'clean_sample_mask');
end
