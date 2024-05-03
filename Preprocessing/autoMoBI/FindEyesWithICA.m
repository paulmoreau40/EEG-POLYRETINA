function [EEG_ica_raw] = FindEyesWithICA(EEG, cfg)
% Automatic detection of eye components with AMICA and ICLabel.
subject = cfg.subjects(cfg.current_subject).id;
N = makeFolderFileNames(cfg, subject);

%% AMICA Loop
% to find eye components that should be projected out before
% automatic continuous data cleaning algorithm

% Compute rank
switch lower(cfg.globalArchitecture)
    case 'simple'
        N_removed_chans = sum(~EEG.etc.clean_channel_mask);
    case 'bemobil'
        N_removed_chans = length(EEG.etc.noisyChannelsDetection.noisyChannels.all);
end

disp('AMICA computation on raw data for projecting eye components out...');
EEG_ica_raw = bemobil_custom_signal_decomposition(EEG, cfg, EEG.nbchan - N_removed_chans);
% Rename the AMICA folder to AMICAeye
movefile(fullfile(N.searchFolder_2arch_rej,sprintf('%s_AMICA', subject)),...
    fullfile(N.searchFolder_2arch_rej,sprintf('%s_AMICAeye', subject)))
%ICA_num = numel(EEG_ica_raw.etc.icaweights_beforerms(:,1));
%pop_topoplot(EEG,0, [1:ICA_num]);

%% Determine (ICLabel) eye components and project out
% iclabel
disp('Determine (ICLabel) eye components and project out...');
EEG_ica_raw = pop_iclabel(EEG_ica_raw, cfg.iclabel);

% plot eyes and save
class_results = EEG_ica_raw.etc.ic_classification.ICLabel.classifications;
titles = cell(size(class_results,1),1);
for t = 1:numel(titles)
    titles{t,1} = ['IC#' num2str(t) ' - ' num2str(class_results(t,3)*100,2) '% eye'];
end

bemobil_plot_patterns(EEG_ica_raw.icawinv, EEG_ica_raw.chanlocs, 'titles', titles,...
    'weights', class_results(:,3), 'minweight', cfg.ICdetect_thresholds(3));
suptitle(sprintf('%s eye components selected by ICLabel (%s)', subject, upper(cfg.ICAmethod)));
saveCurrentFig(fullfile(cfg.figures_folder, 'AutoBadSampsCleaning', 'EyeDetection', filesep),...
    sprintf('%s_eye-comps',subject), {'fig'}, []);
end

