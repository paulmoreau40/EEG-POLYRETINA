%% Start %%
%clear all;
configEEGAFF_Ainhoa;
addpath(genpath('C:\Users\Ainhoa\Documents\M2\MATLAB\codeM2\mne-matlab'))
addpath(genpath('C:\Users\Ainhoa\Documents\M2\MATLAB\codeM2\regu'))

%% Parameters
% Configuration struct for ft_sourceanalysis
cfg = importConfig(study_config, 'ft_sourceanalysis');
redoSourceModelling = false;
redoGCV = false;

for subject_ind = subject_inds
    if ~exist('ALLEEG','var')
        launchEEGLAB;
    end
    
    %subject_ind = 7;
    subject = study_config.subjects(subject_ind).id;
    disp(['Subject ' subject]);
    study_config.current_subject = subject_ind;
    
    % Folders and Files names:
    N = makeFolderFileNames(study_config, subject);
    
    if ~isfile(fullfile(N.searchFolder_3arch_rej,...
            sprintf('%s_ERSPsources_ROIs.mat',N.epochedFile(1:end-12)))) || redoSourceModelling
        
        % Save some information to a summary file
        fileID = fopen('summary_epochs.txt','a');
        fprintf(fileID, '%s:\n', subject);
        
        %% Load EEG data
        EEG = pop_loadset('filename', N.epochedFile, 'filepath', N.searchFolder_3arch_rej);
        data = eeglab2fieldtrip(EEG, 'raw', 'none');
        data.trialinfo = struct2table(EEG.etc.epochInfo.Observations);        
        
        %%%%% MODIFY ACCORDING TO SELECTION RULE
        % Select only trials unflagged at preprocessing
        cfg_sel = struct();
        cfg_sel.trials = data.trialinfo.percCleanEEG>80;
        data = ft_selectdata(cfg_sel,data);
        fprintf(fileID, 'Number of observations kept: %d\n', sum(cfg_sel.trials));   
        %%%%%             
        
        %clear EEG;
        
        EEG_b = pop_loadset('filename', N.baselineEpochedFile, 'filepath', N.searchFolder_3arch_rej);
        data_b = eeglab2fieldtrip(EEG_b, 'raw', 'none');
        data_b.trialinfo = struct2table(EEG_b.etc.epochInfo.Baselines);        
        %%%%% MODIFY ACCORDING TO SELECTION RULE
        % Select only trials unflagged at preprocessing
        cfg_sel.trials = data_b.trialinfo.cleanEEG_base;
        data_b = ft_selectdata(cfg_sel,data_b);
        fprintf(fileID, 'Number of baselines kept: %d\n', sum(cfg_sel.trials));
        %%%%%        
        %clear EEG_b;
        
        fclose(fileID);
        
        labels = data.label;
        if isempty(study_config.subjects(subject_ind).badLocElectrodes)
            chans2remove = {};
        else
            chans2remove = study_config.subjects(subject_ind).badLocElectrodes;
        end
        
        %% Frequency analysis
        cfg_freq_b = importConfig(study_config, 'ft_freqanalysis-base');
        cfg_freq_b.channel = setdiff(labels, chans2remove, 'stable');
        cfg_freq_b.trials = 'all';
        data_freq_b = ft_freqanalysis(cfg_freq_b, data_b);
        %clear data_b
        % noisecov is not computed by ft_sourceanalysis when dealing with frequency data
        cfg.mne.noisecov = speye(numel(data_freq_b.label));
        
        cfg_freq = importConfig(study_config, 'ft_freqanalysis-main');
        cfg_freq.channel = setdiff(labels, chans2remove, 'stable');
        cfg_freq.trials = 'all';
        data_freq = ft_freqanalysis(cfg_freq, data);
        
        %% Load the source model
        %%%%%% Only required if the leadfield is not given as input: %%%%%
        % The volume conduction model of the head should be specified as
        %cfg.headmodel     = structure with volume conduction model, see FT_PREPARE_HEADMODEL
        % The EEG or MEG sensor positions can be present in the data or can be specified as
        %cfg.elec          = structure with electrode positions or filename, see FT_READ_SENS
        %%%%%%
        
        params = importConfig(study_config, 'load_sourcemodel');
        params.chans2remove = chans2remove;
        [cfg.sourcemodel, cov, cov_chol, dipsROI] = load_sourcemodel(N, labels, params);
        if params.cholesky
            error('Not ready: leadfield has already been modified by cov_chol so the optimization will not be correct')
            cfg.mne.sourcecov = speye(size(cov_chol));
        else
            cfg.mne.sourcecov = cov;
        end
        
        % Prepare for optimization of regularization parameter
        G = cell2mat(cfg.sourcemodel.leadfield);
        G_non0 = G(:,all(G,1));
        % PREWHITENING needs to happen before adjusting gain matrix with correlation matrix
        if strcmp(cfg.mne.prewhiten,'yes')
            error('To modify: no cov matrix in frequency domain')
            P = prewhiteningMatrix(data_freq.cov, 1e-0);
            G_non0 = P*G_non0;
        end
        
        % Always use the cov_chol matrix for the optimization
        % Scale cov so that tr(G*G') = tr(Id) in the end
        % (for comparison with the regularization term)
        if strcmp(cfg.mne.scalesourcecov, 'yes')
            scaled_cov = cov_chol.*(size(G_non0,1)/trace(G_non0*cov_chol*(G_non0*cov_chol)'));
            G_non0 = G_non0*scaled_cov;
        else
            G_non0 = G_non0*cov_chol;
        end
        [U,s,V] = csvd(G_non0);
        
        %% GCV for lambda
        if isnan(study_config.subjects(subject_ind).lambdaOptFreq) || redoGCV
            plotParams = struct('do', true, 'subject', subject,...
                'saveFolder', [study_config.figures_folder, 'SourceModelOptim', filesep]);
            if strcmp(cfg.mne.prewhiten,'yes')
                cfg.mne.lambda = doGCVoptimization(G_non0, data_freq, 'freq', P, plotParams);
            else
                cfg.mne.lambda = doGCVoptimization(G_non0, data_freq, 'freq', [], plotParams);
            end
        else
            cfg.mne.lambda = study_config.subjects(subject_ind).lambdaOptFreq;
        end
        
        %% Actual source reconstruction:
        % Compute activations trial by trial to avoid memory overthrow
        cfg.keeptrials = 'no'; % Even for 1 trial, cfg.keeptrials = 'yes' exceeds the maximum memory
        % (has to do with the difference between "freq = ft_checkdata(freq, 'cmbstyle', 'full');"
        % and "freq = ft_checkdata(freq, 'cmbstyle', 'fullfast');" in prepare_freq_matrices.)
        % This function is supposed to compute the csd from fourierspctrm so
        % both options are completely equivalent in terms of results (not combining info from multiple trials)
        
        ROIs = fieldnames(dipsROI);
        for tr = 1:size(data.trialinfo,1)
            fprintf('Trial %d...\n',tr);
            cfg_sel.trials = tr;
            % Compute ERSPs
            data_freq_tr = ft_selectdata(cfg_sel, data_freq);
            % Computation of source activity
            stc_freq_tr = ft_sourceanalysis(cfg, data_freq_tr);
            clear data_freq_tr % save RAM
            
            if tr == 1
                % create stc_freq_rois
                stc_freq_rois = stc_freq_tr;
                stc_freq_rois = rmfield(stc_freq_rois, {'inside', 'pos', 'tri','avg'});
                stc_freq_rois.trialinfo = data_freq.trialinfo;
                stc_freq_rois.method = 'rawtrial';
                % for each ROI:
                for r = 1:numel(ROIs)
                    stc_freq_rois.(sprintf('%s_mom',ROIs{r})) = cell2mat(stc_freq_tr.avg.mom(dipsROI.(ROIs{r})));
                    stc_freq_rois.(sprintf('%s_pow',ROIs{r})) = stc_freq_tr.avg.pow(dipsROI.(ROIs{r}),:,:,:);
                    stc_freq_rois.(sprintf('%s_noisecov',ROIs{r})) = cell2mat(stc_freq_tr.avg.noisecov(dipsROI.(ROIs{r})));
                end
            else
                % for each ROI:
                for r = 1:numel(ROIs)
                    stc_freq_rois.(sprintf('%s_mom',ROIs{r})) = cat(2,stc_freq_rois.(sprintf('%s_mom',ROIs{r})),...
                        cell2mat(stc_freq_tr.avg.mom(dipsROI.(ROIs{r}))));
                    stc_freq_rois.(sprintf('%s_pow',ROIs{r})) = cat(2,stc_freq_rois.(sprintf('%s_pow',ROIs{r})),...
                        stc_freq_tr.avg.pow(dipsROI.(ROIs{r}),:,:,:));
                    stc_freq_rois.(sprintf('%s_noisecov',ROIs{r})) = cat(2,stc_freq_rois.(sprintf('%s_noisecov',ROIs{r})),...
                        cell2mat(stc_freq_tr.avg.noisecov(dipsROI.(ROIs{r}))));
                end
            end
            clear stc_freq_tr % save RAM
        end
        
        %     % Operations
        %     for r = 1:numel(ROIs)
        %         % Average over time:
        %         singletrialbase = mean(stc_freq_rois.(sprintf('%s_pow',ROIs{r})),4);
        %         % Single trial normalization (for pow only)
        %         stc_freq_rois.(sprintf('%s_pow',ROIs{r})) = stc_freq_rois.(sprintf('%s_pow',ROIs{r}))./...
        %             repmat(singletrialbase, [1,1,1,length(stc_freq_rois.time)]);
        %         % dB transformation:
        %         stc_freq_rois.(sprintf('%s_pow',ROIs{r})) = 10*log10(stc_freq_rois.(sprintf('%s_pow',ROIs{r})));
        %     end
        
        % Save activations
        save(fullfile(N.searchFolder_3arch_rej,sprintf('%s_ERSPsources_ROIs',N.epochedFile(1:end-12))), 'stc_freq_rois', '-v7.3');
        clear stc_freq_rois
        
        % Baseline activity
        for tr = 1:size(data_b.trialinfo,1)
            fprintf('Trial %d...\n',tr);
            cfg_sel.trials = tr;
            % Compute ERSPs
            data_freq_tr = ft_selectdata(cfg_sel, data_freq_b);
            % Computation of source activity
            stc_freq_tr = ft_sourceanalysis(cfg, data_freq_tr);
            clear data_freq_tr % save RAM
            
            if tr == 1
                % create stc_freq_rois
                stc_freq_rois_b = stc_freq_tr;
                stc_freq_rois_b = rmfield(stc_freq_rois_b, {'inside', 'pos', 'tri','avg'});
                stc_freq_rois_b.trialinfo = data_freq_b.trialinfo;
                stc_freq_rois_b.method = 'rawtrial';
                % for each ROI:
                for r = 1:numel(ROIs)
                    stc_freq_rois_b.(sprintf('%s_mom',ROIs{r})) = cell2mat(stc_freq_tr.avg.mom(dipsROI.(ROIs{r})));
                    stc_freq_rois_b.(sprintf('%s_pow',ROIs{r})) = stc_freq_tr.avg.pow(dipsROI.(ROIs{r}),:,:,:);
                    stc_freq_rois_b.(sprintf('%s_noisecov',ROIs{r})) = cell2mat(stc_freq_tr.avg.noisecov(dipsROI.(ROIs{r})));
                end
            else
                % for each ROI:
                for r = 1:numel(ROIs)
                    stc_freq_rois_b.(sprintf('%s_mom',ROIs{r})) = cat(2,stc_freq_rois_b.(sprintf('%s_mom',ROIs{r})),...
                        cell2mat(stc_freq_tr.avg.mom(dipsROI.(ROIs{r}))));
                    stc_freq_rois_b.(sprintf('%s_pow',ROIs{r})) = cat(2,stc_freq_rois_b.(sprintf('%s_pow',ROIs{r})),...
                        stc_freq_tr.avg.pow(dipsROI.(ROIs{r}),:,:,:));
                    stc_freq_rois_b.(sprintf('%s_noisecov',ROIs{r})) = cat(2,stc_freq_rois_b.(sprintf('%s_noisecov',ROIs{r})),...
                        cell2mat(stc_freq_tr.avg.noisecov(dipsROI.(ROIs{r}))));
                end
            end
            clear stc_freq_tr % save RAM
        end
        
        %    % Operations
        %     for r = 1:numel(ROIs)
        %         % Average over time:
        %         singletrialbase = mean(stc_freq_rois_b.(sprintf('%s_pow',ROIs{r})),4);
        %         % Single trial normalization (for pow only)
        %         stc_freq_rois_b.(sprintf('%s_pow',ROIs{r})) = stc_freq_rois_b.(sprintf('%s_pow',ROIs{r}))./...
        %             repmat(singletrialbase, [1,1,1,length(stc_freq_rois_b.time)]);
        %         % dB transformation:
        %         stc_freq_rois_b.(sprintf('%s_pow',ROIs{r})) = 10*log10(stc_freq_rois_b.(sprintf('%s_pow',ROIs{r})));
        %     end
        
        % Save activations
        save(fullfile(N.searchFolder_3arch_rej,sprintf('%s_baselines_ERSPsources_ROIs',subject)), 'stc_freq_rois_b', '-v7.3');
        clear stc_freq_rois_b
        
        %     % Separate the dataset per condition to avoid memory overthrow
        %     conditions = unique(data.trialinfo.Condition);
        %     for cnd = 1:numel(conditions)
        %         trials_cnd = find(strcmp(data.trialinfo.Condition,conditions{cnd}));
        %         fprintf('Block %d; %s condition...\n', data.trialinfo.Block(trials_cnd(1)), conditions{cnd}(1:end-1));
        %
        %         % Compute activations trial by trial to avoid memory overthrow
        %         cfg.keeptrials = 'no'; % Even for 1 trial, cfg.keeptrials = 'yes'
        %         exceeds the maximum memory, see explanation above
        %         for tr = 1:length(trials_cnd)
        %             fprintf('Trial %d...\n',tr);
        %             cfg_sel.trials = trials_cnd(tr);
        %             % Compute ERSPs
        %             data_freq_tr = ft_selectdata(cfg_sel, data_freq);
        %             % Computation of source activity
        %             stc_freq_tr = ft_sourceanalysis(cfg, data_freq_tr);
        %             clear data_freq_tr % save RAM
        %             if tr == 1
        %                 % Create stc_freq_cnd
        %                 stc_freq_cnd = stc_freq_tr;
        %                 stc_freq_cnd.method = 'rawtrial';
        %                 stc_freq_cnd = rmfield(stc_freq_cnd, 'avg');
        %                 stc_freq_cnd.trial = stc_freq_tr.avg;
        %             else
        %                 % Copy relevant information to stc_freq_cnd
        %                 stc_freq_cnd.cumtapcnt = ones(size(stc_freq_cnd.cumtapcnt)+[0 1]);
        %                 for m = 1:numel(stc_freq_cnd.trial.mom)
        %                     if ~isempty(stc_freq_cnd.trial.mom{m})
        %                         stc_freq_cnd.trial.mom{m} = cat(2,stc_freq_cnd.trial.mom{m},stc_freq_tr.avg.mom{m});
        %                     end
        %                 end
        %                 stc_freq_cnd.trial.pow = cat(2,stc_freq_cnd.trial.pow,stc_freq_tr.avg.pow);
        %             end
        %             clear stc_freq_tr % save RAM
        %         end
        %         % Save the activations for the whole condition
        %         save(fullfile(N.searchFolder_3arch_rej,sprintf('%s_ERSPsources_B%d_%s',subject,...
        %             data.trialinfo.Block(trials_cnd(1)), conditions{cnd}(1:end-1))), 'stc_freq_cnd', '-v7.3');
        %         clear stc_freq_cnd % save RAM
        %     end
    end
end