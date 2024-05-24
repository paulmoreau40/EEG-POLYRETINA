%% processing loop
% addpath(genpath('C:\Users\Paul Moreau\Documents\EEG_project\EEG-POL')); %paul
% addpath('C:\Users\Paul Moreau\Documents\EEG_project\eeglab2021.0');
% addpath('C:\Users\Paul Moreau\Documents\EEG_project\ParforProgMon-master')

addpath(genpath('F:\PM_Polyretina\EEG_project\EEG-POL')); %paul
addpath('F:\PM_Polyretina\EEG_project\eeglab2021.0');

configEEGAFF2;

skipImport = false; % Boolean to avoid running Mobilab step
overwriteImport = false; % Boolean to force Mobilab step to happen anyway
skipPrep = false; % Boolean to avoid preparation
overwritePrep = false | overwriteImport; % Boolean to force the preparation to happen anyway
skipPreproc = false; % Boolean to avoid preprocessing
overwritePreproc = false | overwritePrep; % Boolean to force the preprocessing to happen anyway
overwriteBadTempOnly = false | overwritePrep; % Boolean to force the bad epochs search to happen anyway
skipICA = false; % Boolean to avoid the decomposition
overwriteICA = false | overwriteBadTempOnly | overwritePreproc; % Boolean to force the decomposition to happen anyway
skipDipoles = false; % Boolean to avoid the dipole fitting
overwriteDipoles = false | overwriteICA; % Boolean to force the dipole fitting to happen anyway
skipAutoLabeling = false; % Boolean to avoid the automatic labeling
overwriteAutoLabeling = false | overwriteICA | (overwriteDipoles & study_config.doDipoleFitting); % Boolean to force IClabel to run anyway
overwriteManualLabeling = false | overwriteAutoLabeling; % Boolean to force reviewing of IC labels manually
overwriteSecondTempRej = false | overwriteManualLabeling; % Boolean to force second temporal rejection
for subject_ind = subject_inds
    if ~exist('ALLEEG','var')
        launchEEGLAB;
    end
    
    % clear RAM
    %STUDY = []; CURRENTSTUDY = 0; ALLEEG = []; EEG=[]; CURRENTSET=[];
    
    % Overwrite subject for testing (COMMENT / DECOMMENT)
    subject_ind = 6;

    subject = study_config.subjects(subject_ind).id;
    disp(['Subject ' subject]);
    study_config.current_subject = subject_ind;
    
    %% Folders and Files names:
    N = makeFolderFileNames(study_config, subject);
    
    %% Importing XDF files
    if ~skipImport && (~exist([N.searchFolder_1 N.postimportFile],'file') || overwriteImport)
        % load xdf files and process them with mobilab, export to eeglab, split MoBI and merge all conditions for EEG
        [ALLEEG, EEG_merged, CURRENTSET] = xdf2set(ALLEEG, CURRENTSET, subject_ind, study_config, overwriteImport);
    elseif ~skipPrep && (~exist([N.searchFolder_2 N.preparedFile],'file') || overwritePrep)
        %continue
        EEG_merged = pop_loadset('filename', N.postimportFile,'filepath', N.searchFolder_1);
        [ALLEEG, EEG_merged, CURRENTSET] = eeg_store(ALLEEG, EEG_merged, CURRENTSET);
    end
    
    %continue
    %% Initial preparation of the data (definitive changes)
    if ~skipPrep && (~exist([N.searchFolder_2 N.preparedFile],'file') || overwritePrep)
        EEG_merged = changeUnit2MicroVolt(EEG_merged, study_config);
        EEG_merged = importCustomElectrodes(EEG_merged, study_config);
        
        % Fill NaNs (if comprised between valid samples:
        EEG_merged = fillNaNs(EEG_merged, study_config);
        
        % Check events and report for missing data
        % Quick check for tests of the EEG_ET subjects.
        NanInspection = isnan(EEG_merged.data);
        NanInspection = sum(NanInspection,1)>0;
        perc = 100*sum(NanInspection)/length(NanInspection);
        
        %EEG_merged = events_check(EEG_merged, study_config);
        EEG_merged = events_check_EEGPOL(EEG_merged, study_config); % PAUL DECOMMENT
        
        % Read distance Answers (now included in the report function)
        %Answers = ReadAnswers_EEGAFF(study_config);
        
        % Create Report for global stats
        %EEG_merged = Report_EEGAFF(EEG_merged, study_config); % PAUL DECOMMENT
        
        %Merge the DataBases saving merged DataBase into .csv
        %MergedDataBase = [MergedDataBase,DataBase];        
        %writetable(struct2table(MergedDataBase),[study_config.study_folder study_config.preprocessing_folder '/MergedDatabase.csv'],'Delimiter',';');
         
        %plot zones and start positions
        %check_zones_EEGAFF(DataBase, study_config)
        % plot answers
        %check_Answers_EEGAFF(DataBase, study_config)
        
        % Boxplot-groups_stat = some preliminary stats (for the included participants
          
        % Resampling
        EEG_merged = pop_resample(EEG_merged, study_config.resample_freq);
        EEG_merged = eeg_checkset(EEG_merged);
        
%         % Remove NaN regions
%         nanTimePoints = sum(isnan(EEG_merged.data),1)>0;
%         EEG_selected = pop_select(EEG_merged, 'nopoint', mask2intervals(nanTimePoints));
        
        % Select Complete Trials only
        trials_intervals = getIntervals(EEG_merged, 'fullTrial', study_config.trialBuffer, true);
        EEG_selected = pop_select(EEG_merged,'time',trials_intervals);
        
        % Check there is no more NaN
        if ~all(sum(isnan(EEG_selected.data),2)==0)
            error('Still NaNs in the selected data set')
        end
        EEG_selected = eeg_checkset(EEG_selected);
        
        EEG_prepared = pop_saveset(EEG_selected, 'filename', N.preparedFile,'filepath', N.searchFolder_2);
        [ALLEEG, EEG_prepared, CURRENTSET] = eeg_store(ALLEEG, EEG_prepared, CURRENTSET);
        % Prepare data:
        %[ALLEEG, EEG_prepared, CURRENTSET] = prepareData(ALLEEG, EEG_merged, CURRENTSET, subject, study_config);
    elseif ~skipPreproc && (~exist([N.searchFolder_2arch_rej N.preICAFile],'file') || overwritePreproc)
        EEG_prepared = pop_loadset('filename', N.preparedFile,'filepath', N.searchFolder_2);
        [ALLEEG, EEG_prepared, CURRENTSET] = eeg_store(ALLEEG, EEG_prepared, CURRENTSET);
    else % No need to load the data
     end
    % Save RAM
    clear EEG_merged
    
    %continue
    %% Preprocessing
    if ~skipPreproc && (~exist([N.searchFolder_2arch_rej N.preICAFile],'file') || overwritePreproc || overwriteBadTempOnly)
        if (overwritePreproc || ~exist([N.searchFolder_2arch N.nobadchansFile],'file'))
            doBadChans = true;
        else
            % This step was entered because of overwriteBadTempOnly
            doBadChans = false;
            EEG_prepared = pop_loadset('filename', N.nobadchansFile, 'filepath', N.searchFolder_2arch);
            [ALLEEG, EEG_prepared, CURRENTSET] = eeg_store(ALLEEG, EEG_prepared, CURRENTSET);
        end
        
        EEG_preproc = preprocess(EEG_prepared, study_config, doBadChans);
        
        fileID = fopen('summary.txt','a');
        fprintf(fileID, '%s:\n', subject);
        fprintf(fileID, 'Rank of the preprocessed EEG set: %d\n',...
            EEG_preproc.nbchan - length(EEG_preproc.etc.noisyChannelsDetection.noisyChannels.all));
        fprintf(fileID, '%d data points in the preprocessed EEG set.\n', EEG_preproc.pnts);
        switch study_config.badSampsRejection
            case 'app'
                fprintf(fileID, '%.1f%% of data rejected by preprocessing.\n',...
                    100*sum(EEG_preproc.etc.APP.rejectedSamples)/length(EEG_preproc.etc.APP.rejectedSamples));
            case 'asr'
                fprintf(fileID, '%.1f%% of data rejected by preprocessing.\n',...
                    100*sum(EEG_preproc.etc.ASR.rejectedSamples)/length(EEG_preproc.etc.ASR.rejectedSamples));
            case 'autoMoBI'
                fprintf(fileID, '%.1f%% of data rejected by preprocessing.\n',...
                    100*sum(EEG_preproc.etc.autoMoBI.rejectedSamples)/length(EEG_preproc.etc.autoMoBI.rejectedSamples));
        end
        fclose(fileID);
        
        EEG_forICA = pop_saveset(EEG_preproc, 'filename', N.preICAFile,'filepath', N.searchFolder_2arch_rej);
        [ALLEEG, EEG_forICA, CURRENTSET] = eeg_store(ALLEEG, EEG_forICA, CURRENTSET);
    elseif ~skipICA && (~exist([N.searchFolder_2arch_rej N.postICAFile],'file') || overwriteICA)
        %continue
        EEG_forICA = pop_loadset('filename', N.preICAFile,'filepath', N.searchFolder_2arch_rej);
        %EEG_forICA = pop_loadset('filename', N.preICAFile,'filepath', N.searchFolder_2arch, 'loadmode', 'info');
        %fprintf('%d data points in the preprocessed EEG set.\n', EEG_forICA.pnts);
        
        [ALLEEG, EEG_forICA, CURRENTSET] = eeg_store(ALLEEG, EEG_forICA, CURRENTSET);
    end
    
    % Save RAM
    clear EEG_preproc
    %continue
    
    if ~skipICA && (~exist([N.searchFolder_2arch_rej N.postICAFile],'file') || overwriteICA)
        switch lower(study_config.globalArchitecture)
            case 'simple'
                N_removed_chans = sum(~EEG_forICA.etc.clean_channel_mask);
            case 'bemobil'
                N_removed_chans = length(EEG_forICA.etc.noisyChannelsDetection.noisyChannels.all);
        end
        
        EEG_ica = bemobil_custom_signal_decomposition(EEG_forICA, study_config, EEG_forICA.nbchan - N_removed_chans);
        %ICA_num = numel(EEG_raw.etc.icaweights_beforerms(:,1));
        %pop_topoplot(EEG,0, [1:ICA_num]);
        
        % Transfer information on the noBadChannels dataset
        EEG_noBadCh = pop_loadset('filename', N.nobadchansFile, 'filepath', N.searchFolder_2arch);
        
        switch study_config.badSampsRejection
            case 'app'
                % Transfer APP information
                EEG_noBadCh.etc.APP = EEG_ica.etc.APP;
            case 'asr'
                % Transfer ASR information
                EEG_noBadCh.etc.ASR = EEG_ica.etc.ASR;
            case 'autoMoBI'
                                % Transfer autoMoBI information
                EEG_noBadCh.etc.autoMoBI = EEG_ica.etc.autoMoBI;
        end
        % Transfer AMICA information
        EEG_noBadCh.icawinv = EEG_ica.icawinv;
        EEG_noBadCh.icasphere = EEG_ica.icasphere;
        EEG_noBadCh.icaweights = EEG_ica.icaweights;
        EEG_noBadCh.icachansind = EEG_ica.icachansind;
        EEG_noBadCh.etc.spatial_filter = EEG_ica.etc.spatial_filter;
        EEG_noBadCh.etc.spatial_filter.preprocessing.filter = EEG_ica.etc.filter;
        EEG_noBadCh.etc.spatial_filter.preprocessing.lineNoiseRemoval = EEG_ica.etc.lineNoiseRemoval;
        
        clear EEG_ica
        EEG_ica = pop_saveset(EEG_noBadCh, 'filename', N.postICAFile, 'filepath', N.searchFolder_2arch_rej);
        [ALLEEG, EEG_ica, CURRENTSET] = eeg_store(ALLEEG, EEG_ica, 0);
    elseif (~study_config.doDipoleFitting && (~exist([N.searchFolder_2arch_rej_ICcats N.postLabelingFile],'file') || overwriteAutoLabeling || overwriteManualLabeling)) ||...
            (study_config.doDipoleFitting && ~skipDipoles && (~exist([N.searchFolder_2arch_rej N.dipfitFile],'file') || overwriteDipoles))
        %continue
        EEG_ica = pop_loadset('filename', N.postICAFile,'filepath', N.searchFolder_2arch_rej);
        [ALLEEG, EEG_ica, CURRENTSET] = eeg_store(ALLEEG, EEG_ica, 0);
    end
    
    clear EEG_forICA
    continue
    
    if study_config.doDipoleFitting
        if ~skipDipoles && (~exist([N.searchFolder_2arch_rej N.dipfitFile],'file') || overwriteDipoles)
            % HP filter
            lowcutoff = study_config.filterICLabel.low_cut_off;
            highcutoff = study_config.filterICLabel.high_cut_off;
            fprintf('Highpass Filtering (%.1f Hz)...\n', lowcutoff)
            [EEG_HP] = custom_filter(EEG_ica, lowcutoff, highcutoff);
            
            % Remove Line Noise
            disp('Removing Line Noise...')
            [EEG_HP, lineNoiseOut] = removeLineNoise_custom(EEG_HP, study_config.lineNoiseRemoval_method, false);
            % If you want to save the filteredEEG_noLN struct with the LineNoiseRemoval information for later:
            EEG_HP.etc.lineNoiseRemoval = lineNoiseOut;
            
            % Remove Bad Temps
            switch study_config.badSampsRejection
                case 'app'
                    EEG_HP = pop_select(EEG_HP, 'nopoint', mask2intervals(EEG_HP.etc.APP.rejectedSamples));
                case 'asr'
                    EEG_HP = pop_select(EEG_HP, 'nopoint', mask2intervals(EEG_HP.etc.ASR.rejectedSamples));
                case 'autoMoBI'
                    EEG_HP = pop_select(EEG_HP, 'nopoint', mask2intervals(EEG_HP.etc.autoMoBI.rejectedSamples));
            end
            
            %% Warping of locations and dipole fitting
            % renames the specified channels, warps the chanlocs on a standard head model and fits dipoles for
            % each IC below the threshold of residual variance
            
            % do the warp and dipfit
            disp('Dipole fitting...');
            EEG_dipfit = dipfit_EEGAFF(EEG_HP, study_config.residualVariance_threshold,...
                study_config.do_remove_outside_head, study_config.number_of_dipoles);
            
            EEG_ica.dipfit = EEG_dipfit.dipfit;
            EEG_ica.etc.preproc_dipfit.filter = EEG_HP.etc.filter;
            EEG_ica.etc.preproc_dipfit.lineNoiseRemoval = lineNoiseOut;
            clear EEG_HP EEG_dipfit
            
            %Save the ica dataset
            EEG_ica = pop_saveset(EEG_ica, 'filename', N.dipfitFile,'filepath', N.searchFolder_2arch_rej);
            [ALLEEG, EEG_ica, CURRENTSET] = eeg_store(ALLEEG, EEG_ica, CURRENTSET);
        elseif (~skipAutoLabeling && (~exist([N.searchFolder_2arch_rej N.IClabelledFile],'file') || overwriteAutoLabeling) ||...
                ~exist([N.searchFolder_2arch_rej_ICcats N.postLabelingFile],'file') || overwriteManualLabeling)
            %continue
            clear EEG_ica 
            EEG_ica = pop_loadset('filename', N.dipfitFile,'filepath', N.searchFolder_2arch_rej);
            [ALLEEG, EEG_ica, CURRENTSET] = eeg_store(ALLEEG, EEG_ica, CURRENTSET);
        end
    end
    
    %continue    
    if ~skipAutoLabeling && (~exist([N.searchFolder_2arch_rej N.IClabelledFile],'file') || overwriteAutoLabeling)
        % HP filter
        lowcutoff = study_config.filterICLabel.low_cut_off;
        highcutoff = study_config.filterICLabel.high_cut_off;
        fprintf('Highpass Filtering (%.1f Hz)...\n', lowcutoff)
        [EEG_HP] = custom_filter(EEG_ica, lowcutoff, highcutoff);
        
        % Remove Line Noise
        disp('Removing Line Noise...')
        [EEG_HP, lineNoiseOut] = removeLineNoise_custom(EEG_HP, study_config.lineNoiseRemoval_method, false);
        % If you want to save the filteredEEG_noLN struct with the LineNoiseRemoval information for later:
        EEG_HP.etc.lineNoiseRemoval = lineNoiseOut;
        
        % Remove Bad Temps
        switch study_config.badSampsRejection
            case 'app'
                EEG_HP = pop_select(EEG_HP, 'nopoint', mask2intervals(EEG_HP.etc.APP.rejectedSamples));
            case 'asr'
                EEG_HP = pop_select(EEG_HP, 'nopoint', mask2intervals(EEG_HP.etc.ASR.rejectedSamples));
            case 'autoMoBI'
                EEG_HP = pop_select(EEG_HP, 'nopoint', mask2intervals(EEG_HP.etc.autoMoBI.rejectedSamples));
        end
        
        %% Automatic IC categorization
        EEG_labelled = iclabel(EEG_HP, study_config.iclabel);
        clear EEG_HP
        EEG_labelled = IC_categorization(EEG_labelled, study_config.ICdetect_thresholds);
        pop_saveset(EEG_labelled, 'filename', N.IClabelledFile,'filepath', N.searchFolder_2arch_rej);
    elseif (~exist([N.searchFolder_2arch_rej_ICcats N.postLabelingFile],'file') || overwriteManualLabeling)
        %continue
        clear EEG_labelled
        EEG_labelled = pop_loadset('filename', N.IClabelledFile,'filepath', N.searchFolder_2arch_rej);
        [ALLEEG, EEG_labelled, CURRENTSET] = eeg_store(ALLEEG, EEG_labelled, CURRENTSET);
    end
    
    continue
    if (~exist([N.searchFolder_2arch_rej_ICcats N.postLabelingFile],'file') || overwriteManualLabeling)
        %% Manual ICs selection
        keptComponents_repaired % calls the script edited manually
        switch study_config.badSampsRejection
            case 'app'
                compsInspect = app_KC;
            case 'asr'
                compsInspect = asr_KC;
            case 'autoMoBI'
                compsInspect = autoMoBI_KC;
        end
        
        if sum(strcmp({compsInspect.ID},subject))==0
            % Subject never inspected
            pop_viewprops(EEG_labelled, 0, 1:size(EEG_labelled.icaact,1),...
                {'freqrange',[1 60]}, {}, 1, 'ICLabel'); % for component properties
            %Removing components--components to be removed are plotted
            kept_comp = input('Components to keep: '); % enter [1,3,4,8,10...]
        else
            % Check components for subject already inspected
            i = find(strcmp({compsInspect.ID},subject));
            
            kept_comp = [];
            for ct = 1:numel(study_config.cats2keep)
                kept_comp = union(kept_comp, compsInspect(i).(study_config.cats2keep{ct}));
            end
           
            %comp2review = union(union(compsInspect(i).Brain,compsInspect(i).BrainWithNoise),compsInspect(i).Doubts);
            comp2review = union(setdiff(compsInspect(i).NoiseWithBrain,compsInspect(i).Doubts),...
                setdiff(compsInspect(i).Doubts_Ainhoa,union(union(union(compsInspect(i).Brain,...
                compsInspect(i).BrainWithNoise),compsInspect(i).Checkerboard),compsInspect(i).Doubts)));
            %comp2review = compsInspect(i).Doubts;
            %comp2review = union(compsInspect(i).Brain,compsInspect(i).BrainWithNoise);
            userDecision1 = input(sprintf('Review subject %s? ',subject));
            if userDecision1
                for c = 1:length(comp2review)
                    pop_prop_extended(EEG_labelled, 0, comp2review(c), NaN, {'freqrange',[1 60]}, {}, 1, 'ICLabel');
                    userDecision2 = input(sprintf('IC%d To keep? ',comp2review(c)));
                    %close gcf;
                    if userDecision2
                        kept_comp = union(kept_comp, comp2review(c));
                    else
                        kept_comp = setdiff(kept_comp, comp2review(c));
                    end
                end
            end
        end
        close all
        
        fprintf('Selected %d Brain components\n',length(kept_comp));
        rem_comp = setdiff(1:size(EEG_labelled.icaact,1),kept_comp); % get the components to remove...
        EEG_compRej = pop_subcomp(EEG_ica, rem_comp, 0); % ...and remove them from EEG_ica.
        
        % Add metainfo:
        EEG_compRej.etc.preproc_labelling.filter = EEG_labelled.etc.filter;
        EEG_compRej.etc.preproc_labelling.lineNoiseRemoval = EEG_labelled.etc.lineNoiseRemoval;
        EEG_compRej.etc.ic_classification = EEG_labelled.etc.ic_classification;
        EEG_compRej.etc.ic_cleaning = struct('method', 'manual inspection',...
            'keptClasses', 1,'keptICs', kept_comp, 'thrownICs', rem_comp);
        clear EEG_labelled
        
        %Save the dataset
        pop_saveset(EEG_compRej, 'filename', N.postLabelingFile,'filepath', N.searchFolder_2arch_rej_ICcats);
    elseif study_config.do_second_tempRej && (~exist([N.searchFolder_2arch_rej_ICcats N.finalFile],'file') || overwriteSecondTempRej)
        clear EEG_labelled
        EEG_compRej = pop_loadset('filename', N.postLabelingFile,'filepath', N.searchFolder_2arch_rej_ICcats);
        [ALLEEG, EEG_compRej, CURRENTSET] = eeg_store(ALLEEG, EEG_compRej, CURRENTSET);
    end
    
    if study_config.do_second_tempRej && (~exist([N.searchFolder_2arch_rej_ICcats N.finalFile],'file') || overwriteSecondTempRej)
        EEG_copy = EEG_compRej;
        
        % Add a channel corresponding to already flagged bad samples
        EEG_copy.nbchan = EEG_compRej.nbchan + 1;
        switch study_config.badSampsRejection
            case 'app'
                EEG_copy.data = cat(1, EEG_compRej.data, EEG_compRej.etc.APP.rejectedSamples);
            case 'asr'
                EEG_copy.data = cat(1, EEG_compRej.data, EEG_compRej.etc.ASR.rejectedSamples);
        end
        EEG_copy.chanlocs(EEG_copy.nbchan) = struct('labels', 'rejSamps', 'type', 'MASK', 'ref', '', 'urchan', [],...
            'X', [], 'Y', [], 'Z', [], 'sph_theta', [], 'sph_phi', [], 'sph_radius', [], 'theta', [], 'radius', []);
        
        % Select data in Conditions only
        CondNames = unique({EEG_compRej.event.Condition},'stable');
        CondNames = CondNames(2:end);
        CondLats = zeros(length(CondNames),2);
        c1 = 1;
        c2 = 1;
        while c2 <= length(CondNames)
            ev = find(strcmp({EEG_compRej.event.Condition},CondNames{c2}),1);
            CondLats(c1,1) = floor(EEG_compRej.event(ev).latency)-EEG_compRej.srate;
            if c1>1 && (CondLats(c1,1) <= CondLats(c1-1,2))
                if c1 == length(CondNames)
                    CondLats(c1,1) = 0; 
                end
                c1 = c1-1;                
            end
            
            ev = find(strcmp({EEG_compRej.event.Condition},CondNames{c2}) &...
                contains({EEG_compRej.event.Phase},'Phase'),1,'last');
            CondLats(c1,2) = ceil(EEG_compRej.event(ev).latency)+EEG_compRej.srate;
            
            c1 = c1+1;
            c2 = c2+1;
        end
        CondLats(CondLats==0) = [];
        EEG_selected = pop_select(EEG_copy, 'point', CondLats);
        clear EEG_copy EEG_compRej
        
        % copy bad samples information and remove the corresponding channel
        switch study_config.badSampsRejection
            case 'app'
                EEG_selected.etc.APP.rejectedSamples = logical(EEG_selected.data(end,:));
            case 'asr'
                EEG_selected.etc.ASR.rejectedSamples = logical(EEG_selected.data(end,:));
        end
        EEG_selected.nbchan = EEG_selected.nbchan - 1;
        EEG_selected.data = EEG_selected.data(1:end-1,:);
        EEG_selected.chanlocs = EEG_selected.chanlocs(1:end-1);
        EEG_selected.icaact = [];
        EEG_selected = eeg_checkset(EEG_selected);
        
        EEG_tempRej = preprocess_4SC(EEG_selected, study_config, false);
        
        %pop_viewprops(EEG_tempRej, 0, 1:size(EEG_tempRej.icaact,1),...
        %    {'freqrange',[1 60]}, {}, 1); % for component properties
        
        EEG_final = EEG_selected;
        switch study_config.badSampsRejection
            case 'app'
                EEG_final.etc.APP2 = EEG_tempRej.etc.APP;
            case 'asr'
                EEG_final.etc.ASR2 = EEG_tempRej.etc.ASR;
        end       
        clear EEG_tempRej
        
        pop_saveset(EEG_final, 'filename', N.finalFile,'filepath', N.searchFolder_2arch_rej_ICcats);
    end
    
    if ~isempty(ALLEEG)
        ALLEEG = pop_delset(ALLEEG, 1:CURRENTSET);
        CURRENTSET=1;
    end
end