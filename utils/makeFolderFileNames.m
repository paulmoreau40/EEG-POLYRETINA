function names = makeFolderFileNames(cfg, subject)
names.searchFolder_1 = [cfg.study_folder cfg.raw_EEGLAB_data_folder subject filesep];
createFolder(names.searchFolder_1);
names.searchFolder_2 = [cfg.study_folder cfg.preprocessing_folder];
createFolder(names.searchFolder_2);
names.searchFolder_2arch = [names.searchFolder_2 lower(cfg.globalArchitecture) filesep];
createFolder(names.searchFolder_2arch);
names.searchFolder_3 = [cfg.study_folder cfg.single_subject_analysis_folder];
createFolder(names.searchFolder_3);
names.searchFolder_3arch = [names.searchFolder_3 lower(cfg.globalArchitecture) filesep];
createFolder(names.searchFolder_3arch);
names.searchFolder_4 = [cfg.study_folder cfg.multi_subject_analysis_folder];
createFolder(names.searchFolder_4);
names.searchFolder_4arch = [names.searchFolder_4 lower(cfg.globalArchitecture) filesep];
createFolder(names.searchFolder_4arch);
switch cfg.badSampsRejection
    case 'manual'
        names.searchFolder_2arch_rej = [names.searchFolder_2arch 'Manual' filesep];
        names.searchFolder_3arch_rej = [names.searchFolder_3arch 'Manual' filesep];
        names.searchFolder_4arch_rej = [names.searchFolder_4arch 'Manual' filesep];
    case 'app'
        names.searchFolder_2arch_rej = [names.searchFolder_2arch 'APP' filesep];
        names.searchFolder_3arch_rej = [names.searchFolder_3arch 'APP' filesep];
        names.searchFolder_4arch_rej = [names.searchFolder_4arch 'APP' filesep];
    case 'asr'
        switch cfg.ASR_use
            case 'rewrite'
                names.searchFolder_2arch_rej = [names.searchFolder_2arch 'ASR_corrected' filesep];
                names.searchFolder_3arch_rej = [names.searchFolder_3arch 'ASR_corrected' filesep];
                names.searchFolder_4arch_rej = [names.searchFolder_4arch 'ASR_corrected' filesep];
            case 'reject'
                names.searchFolder_2arch_rej = [names.searchFolder_2arch 'ASR_rejected' filesep];
                names.searchFolder_3arch_rej = [names.searchFolder_3arch 'ASR_rejected' filesep];
                names.searchFolder_4arch_rej = [names.searchFolder_4arch 'ASR_rejected' filesep];
            otherwise
                error('Unknown ASR use')
        end
            case 'autoMoBI'
        names.searchFolder_2arch_rej = [names.searchFolder_2arch 'autoMoBI' filesep];
        names.searchFolder_3arch_rej = [names.searchFolder_3arch 'autoMoBI' filesep];
        names.searchFolder_4arch_rej = [names.searchFolder_4arch 'autoMoBI' filesep];
    otherwise
        error('Unknown bad samples rejection')
end
createFolder(names.searchFolder_2arch_rej);
createFolder(names.searchFolder_3arch_rej);
createFolder(names.searchFolder_4arch_rej);

for c = 1:numel(cfg.cats2keep)
    if c == 1
        names.searchFolder_2arch_rej_ICcats = [names.searchFolder_2arch_rej cfg.cats2keep{c}];
        names.searchFolder_3arch_rej_ICcats = [names.searchFolder_3arch_rej cfg.cats2keep{c}];
        names.searchFolder_4arch_rej_ICcats = [names.searchFolder_4arch_rej cfg.cats2keep{c}];
    else
        names.searchFolder_2arch_rej_ICcats = [names.searchFolder_2arch_rej_ICcats '_' cfg.cats2keep{c}];
        names.searchFolder_3arch_rej_ICcats = [names.searchFolder_3arch_rej_ICcats '_' cfg.cats2keep{c}];
        names.searchFolder_4arch_rej_ICcats = [names.searchFolder_4arch_rej_ICcats '_' cfg.cats2keep{c}];
    end
end
names.searchFolder_2arch_rej_ICcats = [names.searchFolder_2arch_rej_ICcats filesep];
createFolder(names.searchFolder_2arch_rej_ICcats);
names.searchFolder_3arch_rej_ICcats = [names.searchFolder_3arch_rej_ICcats filesep];
createFolder(names.searchFolder_3arch_rej_ICcats);
names.searchFolder_4arch_rej_ICcats = [names.searchFolder_4arch_rej_ICcats filesep];
createFolder(names.searchFolder_4arch_rej_ICcats);

names.searchFolder_3arch_rej_ICcats_dParam = [sprintf('%sdepth%d',...
    names.searchFolder_3arch_rej_ICcats, cfg.recon.fwd.normalizeDepthParam*100) filesep];
createFolder(names.searchFolder_3arch_rej_ICcats_dParam);
names.searchFolder_4arch_rej_ICcats_dParam = [sprintf('%sdepth%d',...
    names.searchFolder_4arch_rej_ICcats, cfg.recon.fwd.normalizeDepthParam*100) filesep];
createFolder(names.searchFolder_4arch_rej_ICcats_dParam);

names.postimportFile = [subject '_EEG_' cfg.merged_filename];
names.preparedFile = [subject '_' cfg.prepared_filename];
names.nobadchansFile = [subject '_' cfg.BadChansRemoved_filename];
names.eyedetectionFile = [subject '_' cfg.EyeDetection_filename];
names.preICAFile = [subject '_' cfg.beforeICA_filename];
names.postICAFile = [subject '_' cfg.icaOutput_filename];
names.dipfitFile = [subject '_' cfg.dipolesFitted_filename];
names.IClabelledFile = [subject '_' cfg.icLabelled_filename];
names.postLabelingFile = [subject '_' cfg.icaSelect_filename];
names.finalFile = [subject '_' cfg.postICA_tempRej_filename];

% switch cfg.epochs.event
%     case 'Intersection'
%         prefix_epoch = sprintf('%s_inter',subject);        
%     otherwise
%         error('Unknown event type')
% end
% 
% switch lower(cfg.epochs.window)
%     case 'full'
%         names.epochedFile = sprintf('%s_%s_full_%s',prefix_epoch,lower(cfg.model),cfg.epoched_filename);
%     case 'fixed'
%         names.epochedFile = sprintf('%s_%s_bef%d-aft%d_%s',prefix_epoch,lower(cfg.model),...
%             abs(cfg.epochs.limits_wdw(1)), cfg.epochs.limits_wdw(2), cfg.epoched_filename);
%     otherwise
%         error('Unknown window type')
% end
% names.baselineEpochedFile = sprintf('%s_%s_%s',subject,lower(cfg.model),cfg.base_epoched_filename);

%% Source reconstruction:
names.bemFolder = fullfile(cfg.recon.path, subject, 'bem');
names.roiFolder = fullfile(cfg.recon.path, subject, 'roi');
   
names.transFile = sprintf('%s-trans.fif', subject);
names.bemFile = sprintf('%s_bem.fif', subject);
if ischar(cfg.recon.src_spacing)
    file = sprintf('%s-%s_src', subject, cfg.recon.src_spacing);
else
    file = sprintf('%s-%dmm_src', subject, cfg.recon.src_spacing);
end
names.srcFile = sprintf('%s.fif', file);

if strcmp(cfg.recon.downsamp_bem, 'None')
    file = sprintf('%s_bem', file);
else
    file = sprintf('%s-ico%d_bem', file, cfg.recon.downsamp_bem);
end
names.fwdFile = sprintf('%s_fwd.fif', file);

if isempty(cfg.recon.dist_limit_mm)
    names.roiSumFile = sprintf('%s-res%dmm-%s-nolim-roi_summary.mat',...
        file, cfg.recon.resolution, cfg.recon.roi_creation);
    names.covFile = sprintf('%s-res%dmm-%s-nolim-roi_cov.mat',...
        file, cfg.recon.resolution, cfg.recon.roi_creation);
    names.covCholFile = sprintf('%s-res%dmm-%s-nolim-roi_covchol.mat',...
        file, cfg.recon.resolution, cfg.recon.roi_creation);
else
    names.roiSumFile = sprintf('%s-res%dmm-%s-lim%dmm-roi_summary.mat',...
        file, cfg.recon.resolution, cfg.recon.roi_creation, cfg.recon.dist_limit_mm);
    names.covFile = sprintf('%s-res%dmm-%s-lim%dmm-roi_cov.mat',...
        file, cfg.recon.resolution, cfg.recon.roi_creation, cfg.recon.dist_limit_mm);
    names.covCholFile = sprintf('%s-res%dmm-%s-lim%dmm-roi_covchol.mat',...
        file, cfg.recon.resolution, cfg.recon.roi_creation, cfg.recon.dist_limit_mm);
end

    function createFolder(path)
        if ~exist(path, 'dir')
            mkdir(path);
        end
    end
end