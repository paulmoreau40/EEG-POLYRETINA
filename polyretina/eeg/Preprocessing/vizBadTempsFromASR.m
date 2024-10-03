function vizBadTempsFromASR(cfg)
subject = cfg.subjects(cfg.current_subject).id;
searchFolder_2 = [cfg.study_folder cfg.preprocessing_folder];
searchFolder_2arch = [searchFolder_2 lower(cfg.globalArchitecture) filesep];
switch cfg.ASR_use
    case 'rewrite'
        searchFolder_2arch_asr = [searchFolder_2arch 'ASR_corrected' filesep];
    case 'reject'
        searchFolder_2arch_asr = [searchFolder_2arch 'ASR_rejected' filesep];
    otherwise
        error('Unknown ASR use')
end
nobadchansFile = [subject '_' cfg.BadChansRemoved_filename];
preICAFile = [subject '_' cfg.beforeICA_filename];

EEG = pop_loadset('filename', nobadchansFile, 'filepath', searchFolder_2arch);

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

EEG_asr = pop_loadset('filename', preICAFile, 'filepath', searchFolder_2arch_asr);
% Necessary for vis_artifacts to work correctly
EEG_asr.etc.clean_sample_mask = ~EEG_asr.etc.ASR.rejectedSamples;

vis_artifacts(EEG_asr, EEG_HP_noLN);

% rejIntervals = mask2intervals(EEG_asr.etc.ASR.rejectedSamples);
% winrej = rejIntervals;
% winrej(:,end+1:end+3) = repmat([1,0,0], size(rejIntervals,1),1);
% winrej(:,end+1:end+EEG_HP_noLN.nbchan) = 1;
% 
% % Transform to seconds
% rejIntervals = rejIntervals./EEG_asr.srate;
% % Add 1s buffers:
% rejIntervals(:,1) = rejIntervals(:,1)-1;
% rejIntervals(:,2) = rejIntervals(:,2)+1;
% for i = 1:size(rejIntervals,1)
%     vis_artifacts(EEG_asr, EEG_HP_noLN);
% end
% 
% remainingData = EEG_HP_noLN.data;
% remainingData(:,EEG_asr.etc.ASR.rejectedSamples) = 0;
% rejectedData = EEG_HP_noLN.data;
% rejectedData(:,~EEG_asr.etc.ASR.rejectedSamples) = 0;
% 
% eegplot(remainingData,...
%     'srate', EEG_HP_noLN.srate,...
%     'data2', rejectedData,...
%     'events', EEG_HP_noLN.event,...
%     'submean', 'off')

end
