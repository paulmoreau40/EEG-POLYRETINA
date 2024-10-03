configEEGAFF_Ainhoa;

subject_inds = 2:9;
for subject_ind = subject_inds
    if ~exist('ALLEEG','var')
        launchEEGLAB;
    end
    
    subject = study_config.subjects(subject_ind).id;
    disp(['Subject ' subject]);
    study_config.current_subject = subject_ind;
    
    %% Folders and Files names:
    N = makeFolderFileNames(study_config, subject);
    path_compAnalysis = 'C:\Users\Ainhoa\Documents\M2\MATLAB\codeM2\myCode\CompAnalysis';
    
    %% Initial preparation of the data (definitive changes)
    EEG_prepared = pop_loadset('filename', N.preparedFile,'filepath', N.searchFolder_2, 'loadmode', 'info');
    if subject_ind == 2
        % already done
    else
        EEG_prepared = eeg_checkset(EEG_prepared, 'makeur');
        EEG_prepared.etc.CompAnalysis = importdata(fullfile(path_compAnalysis, sprintf('CompAnalysis_%s.mat', subject)));
        EEG_prepared = pop_saveset(EEG_prepared, 'filename', N.preparedFile,'filepath', N.searchFolder_2);
    end
    
    %% Preprocessing
    EEG_interpAvRef = pop_loadset('filename', N.nobadchansFile, 'filepath', N.searchFolder_2arch, 'loadmode', 'info');
    EEG_interpAvRef.event = EEG_prepared.event;
    EEG_interpAvRef.urevent = EEG_prepared.urevent;
    EEG_interpAvRef.etc.CompAnalysis = importdata(fullfile(path_compAnalysis, sprintf('CompAnalysis_%s.mat', subject)));
    pop_saveset(EEG_interpAvRef, 'filename', N.nobadchansFile,'filepath', N.searchFolder_2arch);
    clear EEG_interpAvRef
    
    EEG_forICA = pop_loadset('filename', N.preICAFile,'filepath', N.searchFolder_2arch_rej, 'loadmode', 'info');
    %EEG_forICA.event = EEG_prepared.event;
    %EEG_forICA.urevent = EEG_prepared.urevent;
    EEG_forICA.etc.CompAnalysis = importdata(fullfile(path_compAnalysis, sprintf('CompAnalysis_%s.mat', subject)));
    if strcmp(study_config.badSampsRejection, 'autoMoBI') && length(EEG_forICA.etc.autoMoBI.rejectedSamples) ~= EEG_prepared.pnts
        buffered_bounds = EEG_forICA.etc.autoMoBI.invalid_segments_final_start_stop_sample;
        % boolean indicating rejected samples in the original dataset
        rejSamps = zeros(1,EEG_prepared.pnts);
        for i = 1:size(buffered_bounds,1)
            if buffered_bounds(i,1) < 1
                inter_rej = 1:buffered_bounds(i,2);
                rejSamps(1:buffered_bounds(i,2)) = ones(1,length(inter_rej));
            elseif buffered_bounds(i,2) > EEG_prepared.pnts
                inter_rej = buffered_bounds(i,1):EEG_prepared.pnts;
                rejSamps(buffered_bounds(i,1):EEG_prepared.pnts) = ones(1,length(inter_rej));
            else
                inter_rej = buffered_bounds(i,1):buffered_bounds(i,2);
                rejSamps(buffered_bounds(i,1):buffered_bounds(i,2)) = ones(1,length(inter_rej));
            end
        end
        EEG_forICA.etc.autoMoBI.rejectedSamples = rejSamps;
    end    
    pop_saveset(EEG_forICA, 'filename', N.preICAFile,'filepath', N.searchFolder_2arch_rej);
    clear EEG_forICA
    
    EEG_ica = pop_loadset('filename', N.postICAFile,'filepath', N.searchFolder_2arch_rej, 'loadmode', 'info');
    EEG_ica.event = EEG_prepared.event;
    EEG_ica.urevent = EEG_prepared.urevent;
    EEG_ica.etc.CompAnalysis = importdata(fullfile(path_compAnalysis, sprintf('CompAnalysis_%s.mat', subject)));
    if strcmp(study_config.badSampsRejection, 'autoMoBI') && length(EEG_ica.etc.autoMoBI.rejectedSamples) ~= EEG_prepared.pnts
        EEG_ica.etc.autoMoBI.rejectedSamples = rejSamps;
    end
    pop_saveset(EEG_ica, 'filename', N.postICAFile, 'filepath', N.searchFolder_2arch_rej);
    clear EEG_ica
    
    if study_config.doDipoleFitting
        EEG_dipfit = pop_loadset('filename', N.dipfitFile,'filepath', N.searchFolder_2arch_rej, 'loadmode', 'info');
        EEG_dipfit.event = EEG_prepared.event;
        EEG_dipfit.urevent = EEG_prepared.urevent;
        EEG_dipfit.etc.CompAnalysis = importdata(fullfile(path_compAnalysis, sprintf('CompAnalysis_%s.mat', subject)));
        if strcmp(study_config.badSampsRejection, 'autoMoBI') && length(EEG_dipfit.etc.autoMoBI.rejectedSamples) ~= EEG_prepared.pnts
            EEG_dipfit.etc.autoMoBI.rejectedSamples = rejSamps;
        end
        pop_saveset(EEG_dipfit, 'filename', N.dipfitFile,'filepath', N.searchFolder_2arch_rej);
    end
    clear EEG_dipfit
    
    EEG_labelled = pop_loadset('filename', N.IClabelledFile,'filepath', N.searchFolder_2arch_rej, 'loadmode', 'info');
    %EEG_labelled.event = EEG_prepared.event;
    %EEG_labelled.urevent = EEG_prepared.urevent;
    EEG_labelled.etc.CompAnalysis = importdata(fullfile(path_compAnalysis, sprintf('CompAnalysis_%s.mat', subject)));
    if strcmp(study_config.badSampsRejection, 'autoMoBI') && length(EEG_labelled.etc.autoMoBI.rejectedSamples) ~= EEG_prepared.pnts
        EEG_labelled.etc.autoMoBI.rejectedSamples = rejSamps;
    end
    pop_saveset(EEG_labelled, 'filename', N.IClabelledFile,'filepath', N.searchFolder_2arch_rej);
    clear EEG_labelled
    
    EEG_final = pop_loadset('filename', N.postLabelingFile,'filepath', N.searchFolder_2arch_rej, 'loadmode', 'info');
    EEG_final.event = EEG_prepared.event;
    EEG_final.urevent = EEG_prepared.urevent;
    EEG_final.etc.CompAnalysis = importdata(fullfile(path_compAnalysis, sprintf('CompAnalysis_%s.mat', subject)));
    if strcmp(study_config.badSampsRejection, 'autoMoBI') && length(EEG_final.etc.autoMoBI.rejectedSamples) ~= EEG_prepared.pnts
        EEG_final.etc.autoMoBI.rejectedSamples = rejSamps;
    end
    pop_saveset(EEG_final, 'filename', N.postLabelingFile,'filepath', N.searchFolder_2arch_rej);
    
    clear EEG_prepared EEG_final
end