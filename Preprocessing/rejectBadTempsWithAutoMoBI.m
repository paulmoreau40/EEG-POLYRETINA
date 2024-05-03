function EEG_autoMoBI = rejectBadTempsWithAutoMoBI(EEG, cfg)
subject = cfg.subjects(cfg.current_subject).id;
N = makeFolderFileNames(cfg, subject);

% 1. HP filter
lowcutoff = cfg.filterPreProc.low_cut_off;
highcutoff = cfg.filterPreProc.high_cut_off;
if ~isempty(lowcutoff)
    fprintf('Highpass Filtering (%.1f Hz) for autoMoBI...\n', lowcutoff)
end
if ~isempty(highcutoff)
    fprintf('Lowpass Filtering (%.1f Hz) for autoMoBI...\n', highcutoff)
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

%% ----------------3. BERLIN PIPELINE--------------------------
% The EEG set is saved at different steps in the process
%% Perform ICA on specific parts of the dataset to detect eye components
dataOfInterest = cfg.autoMoBI.DoI4AMICAeye;
OoI_segments_index = select_data_of_interest_EEGAFF2(EEG_HP_noLN, dataOfInterest);
EEG_sel_for_ICAeye = eeg_eegrej(EEG_HP_noLN, OoI_segments_index);

% save the dataset (to have the data locally, better for AMICA)
EEG_sel_for_ICAeye = pop_saveset(EEG_sel_for_ICAeye,...
    'filename', N.eyedetectionFile, 'filepath', N.searchFolder_2arch_rej);

[EEG_AMICAeye] = FindEyesWithICA(EEG_sel_for_ICAeye, cfg);
clear EEG_sel_for_ICAeye

% Transfer AMICA information and save the dataset
EEG_prep4auto_cleaning = EEG_HP_noLN;
EEG_prep4auto_cleaning.icawinv = EEG_AMICAeye.icawinv;
EEG_prep4auto_cleaning.icasphere = EEG_AMICAeye.icasphere;
EEG_prep4auto_cleaning.icaweights = EEG_AMICAeye.icaweights;
EEG_prep4auto_cleaning.icachansind = EEG_AMICAeye.icachansind;
EEG_prep4auto_cleaning.etc.spatial_filter = EEG_AMICAeye.etc.spatial_filter;
EEG_prep4auto_cleaning.etc.spatial_filter.preprocessing.filter = EEG_AMICAeye.etc.filter;
EEG_prep4auto_cleaning.etc.spatial_filter.preprocessing.lineNoiseRemoval = EEG_AMICAeye.etc.lineNoiseRemoval;
EEG_prep4auto_cleaning.etc.ic_classification = EEG_AMICAeye.etc.ic_classification;
EEG_prep4auto_cleaning = eeg_checkset(EEG_prep4auto_cleaning);
EEG_prep4auto_cleaning = pop_saveset(EEG_prep4auto_cleaning,...
    'filename', N.eyedetectionFile, 'filepath', N.searchFolder_2arch_rej);

%% Define eyes and remove
class_results = EEG_prep4auto_cleaning.etc.ic_classification.ICLabel.classifications;

%% Visual inspection
%{
eye_comps = DefineEyeComps(class_results, bemobil_config.ICdetection_thresholds, 'all', [], []);
freq_range = [1 100]; % in Hz
spec_opt = {'freqrange', freq_range}; % cell array of options which are passed to spectopo()
erp_opt = {}; % cell array of options which are passed to erpimage()
pop_viewprops(EEG_prep4auto_cleaning, 0, eye_comps', spec_opt , erp_opt, 1, 'ICLabel');
%}

%% Only the unique ones
eye_comps = DefineEyeComps(class_results, cfg.ICdetect_thresholds, 'unique',...
    fullfile(cfg.figures_folder, 'AutoBadSampsCleaning', 'EyeDetection', filesep),...
    sprintf('%s_Eye_ICs_distribution_%s_%s', subject, upper(cfg.ICAmethod), cell2str(dataOfInterest, '_')));

%% Subtract components:
EEG_prep4auto_cleaning = pop_subcomp(EEG_prep4auto_cleaning, eye_comps);
%{
% project out eyes and save
EEG_AMICA_no_eyes = pop_subcomp(EEG_AMICA_raw, eyes);
% save the dataset without eye components, overwritting preceeding one
pop_saveset(EEG_AMICA_no_eyes, 'filename', [bemobil_config.filename_prefix num2str(subject) '_'...
    bemobil_config.without_eyes], 'filepath', output_filepath);
%}

%% determine automatic time domain cleaning boundaries on the channel level
% Remove events to save RAM
EEG_prep4auto_cleaning = pop_editeventvals(EEG_prep4auto_cleaning,'delete', 1:length(EEG_prep4auto_cleaning.event));

disp('Determining continuous data cleaning boundaries...')
[~, automatic_cleaning] = autoClean_continuousEEG_custom(...
    EEG_prep4auto_cleaning, cfg.autoMoBI, subject);

EEG_HP_noLN.etc.spatial_filter_eyerejection = EEG_prep4auto_cleaning.etc.spatial_filter;
EEG_HP_noLN.etc.ic_classification_eyerejection = EEG_prep4auto_cleaning.etc.ic_classification;
EEG_HP_noLN.etc.ic_classification_eyerejection.removed_comps = eye_comps;
EEG_HP_noLN.etc.autoMoBI = automatic_cleaning;
clear EEG_prep4auto_cleaning

% LUKAS adding: Add buffers to bad epochs
disp('Adding buffers to data cleaning boundaries...')
buffer_length_sample = (cfg.autoMoBI.buffer_length*cfg.autoMoBI.wind_ms)*(EEG_HP_noLN.srate/1000);
EEG_HP_noLN = add_buffers_continousCleaning(EEG_HP_noLN, buffer_length_sample);

% Get cleaning indices
invalid_segments_index = EEG_HP_noLN.etc.autoMoBI.invalid_segments_final_start_stop_sample;
% Reject the bad segments in the dataset
EEG_autoMoBI = eeg_eegrej(EEG_HP_noLN, invalid_segments_index);
end