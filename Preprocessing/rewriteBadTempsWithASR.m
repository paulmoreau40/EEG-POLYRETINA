function EEG_asr = rewriteBadTempsWithASR(EEG, cfg)

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

% 3. Correct transient artifacts with ASR
[EEG_asr,~,~] = clean_artifacts(EEG_HP_noLN, 'Highpass', 'off', 'FlatlineCriterion', 'off',...
    'ChannelCriterion', 'off', 'LineNoiseCriterion', 'off', 'BurstCriterion', cfg.burst_crit, 'BurstRejection', 'off');

% To enable later vizualization of ASR effects:
out_filepath = fullfile([cfg.study_folder cfg.preprocessing_folder lower(cfg.globalArchitecture)],...
    'ASR_corrected', 'ASR_output');
if ~exist(out_filepath, 'dir')
    mkdir(out_filepath)
end
pop_saveset(EEG_asr, 'filename', [subject '_' cfg.ASRout_filename], 'filepath', out_filepath);

EEG_asr.etc.ASR.BurstCrit = cfg.burst_crit;
% Get portions of data which have been rejected:
EEG_asr.etc.ASR.rejectedSamples = ~EEG_asr.etc.clean_sample_mask;
% Get portions of data which have changed:
EEG_asr.etc.ASR.modifiedSamples = sum(abs(EEG_asr.data-EEG_HP_noLN.data(:,EEG_asr.etc.clean_sample_mask)),1) > 1e-10;

EEG_asr.etc = rmfield(EEG_asr.etc, 'clean_sample_mask');
end
