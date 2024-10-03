%% Start %%
%clear all;
config4SC_Alex;
addpath(genpath('C:\Users\Ainhoa\Documents\M2\MATLAB\codeM2\mne-matlab'))

%% Parameters
% Configuration struct for ft_sourceanalysis
cfg = importConfig(study_config, 'ft_sourceanalysis');
usebase = false;
redoSourceModelling = false;
redoGCV = false;

for subject_ind = subject_inds(1:end)
    if ~exist('ALLEEG','var')
        launchEEGLAB;
    end
    
    %subject_ind = 7;
    subject = study_config.subjects(subject_ind).id;
    disp(['Subject ' subject]);
    study_config.current_subject = subject_ind;
    
    % Folders and Files names:
    N = makeFolderFileNames(study_config, subject);
    
    % Save some information to a summary file
    fileID = fopen('summary_epochs.txt','a');
    fprintf(fileID, '%s:\n', subject);
    
    % Load EEG data
    EEG = pop_loadset('filename', N.epochedFile, 'filepath', N.searchFolder_3arch_rej);
    data = eeglab2fieldtrip(EEG, 'raw', 'none');
    data.trialinfo = struct2table(EEG.etc.epochInfo.Intersections);
    %%%%% MODIFY ACCORDING TO SELECTION RULE
    % Select only trials unflagged at preprocessing
    cfg_sel = struct();
    cfg_sel.trials = data.trialinfo.percCleanEEG>80;
    data = ft_selectdata(cfg_sel,data);
    fprintf(fileID, 'Number of observations kept: %d\n', sum(cfg_sel.trials));
    %%%%%
    %clear EEG;
    
    if usebase
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
    end
    fclose(fileID);
    
    labels = data.label;
    if isempty(study_config.subjects(subject_ind).badLocElectrodes)
        chans2remove = {};
    else
        chans2remove = study_config.subjects(subject_ind).badLocElectrodes;
    end
    
    %% Time lock analysis
    cfg_time = importConfig(study_config, 'ft_timelockanalysis-main');
    cfg_time.channel = setdiff(labels, chans2remove, 'stable');
    cfg_time.trials = 'all';
    data_time = ft_timelockanalysis(cfg_time, data);
    data_time.cov = repmat(reshape(eye(numel(data_time.label)), size(data_time.cov(1,:,:))),...
        [size(data_time.trial,1),1,1]);
    
    if usebase
        cfg_time_b = importConfig(study_config, 'ft_timelockanalysis-base');
        cfg_time_b.channel = setdiff(labels, chans2remove, 'stable');
        cfg_time_b.trials = 'all';
        data_time_b = ft_timelockanalysis(cfg_time_b, data_b);
        %cfg.mne.noisecov = data_time_b.cov; % does not work
        % Replace data covariance matrix with baseline covariance matrix
        data_time.cov = repmat(reshape(data_time_b.cov, [1,size(data_time_b.cov)]),...
            [size(data_time.trial,1),1,1]);
    end
    
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
        P = prewhiteningMatrix(data_time.cov, 1e-0);
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
    if isnan(study_config.subjects(subject_ind).lambdaOptTime) || redoGCV
        plotParams = struct('do', true, 'subject', subject,...
            'saveFolder', [study_config.figures_folder, 'SourceModelOptim', filesep]);
        if strcmp(cfg.mne.prewhiten,'yes')
            cfg.mne.lambda = doGCVoptimization(G_non0, data_time, 'time', P, plotParams);
        else
            cfg.mne.lambda = doGCVoptimization(G_non0, data_time, 'time', [], plotParams);
        end
    else
        cfg.mne.lambda = study_config.subjects(subject_ind).lambdaOptTime;
    end
    
    %% Actual source reconstruction:
    cfg.rawtrial = 'yes';
    % Compute activations trial by trial to avoid memory overthrow
    ROIs = fieldnames(dipsROI);
    for tr = 1:size(data.trialinfo,1)
        fprintf('Trial %d...\n',tr);
        cfg_sel.trials = tr;
        % Compute ERSPs
        data_time_tr = ft_selectdata(cfg_sel, data_time);
        % Computation of source activity
        stc_time_tr = ft_sourceanalysis(cfg, data_time_tr);
        clear data_time_tr % save RAM
        
        if tr == 1
            % create stc_freq_rois
            stc_time_rois = stc_time_tr;
            stc_time_rois = rmfield(stc_time_rois, {'inside', 'pos', 'tri', 'trial'});
            stc_time_rois.trialinfo = data_time.trialinfo;
            % for each ROI:
            for r = 1:numel(ROIs)
                stc_time_rois.(sprintf('%s_mom',ROIs{r})) = reshape(cell2mat(stc_time_tr.trial.mom(dipsROI.(ROIs{r}))),...
                    [length(dipsROI.(ROIs{r})),1,length(stc_time_tr.time)]);
                stc_time_rois.(sprintf('%s_pow',ROIs{r})) = reshape(stc_time_tr.trial.pow(dipsROI.(ROIs{r}),:),...
                    [length(dipsROI.(ROIs{r})),1,length(stc_time_tr.time)]);
                stc_time_rois.(sprintf('%s_noisecov',ROIs{r})) = reshape(cell2mat(stc_time_tr.trial.noisecov(dipsROI.(ROIs{r}))),...
                    [length(dipsROI.(ROIs{r})),1]);
            end
        else
            % for each ROI:
            for r = 1:numel(ROIs)
                stc_time_rois.(sprintf('%s_mom',ROIs{r})) = cat(2,stc_time_rois.(sprintf('%s_mom',ROIs{r})),...
                    reshape(cell2mat(stc_time_tr.trial.mom(dipsROI.(ROIs{r}))),[length(dipsROI.(ROIs{r})),1,length(stc_time_tr.time)]));
                stc_time_rois.(sprintf('%s_pow',ROIs{r})) = cat(2,stc_time_rois.(sprintf('%s_pow',ROIs{r})),...
                    reshape(stc_time_tr.trial.pow(dipsROI.(ROIs{r}),:,:),[length(dipsROI.(ROIs{r})),1,length(stc_time_tr.time)]));
                stc_time_rois.(sprintf('%s_noisecov',ROIs{r})) = cat(2,stc_time_rois.(sprintf('%s_noisecov',ROIs{r})),...
                    reshape(cell2mat(stc_time_tr.trial.noisecov(dipsROI.(ROIs{r}))),[length(dipsROI.(ROIs{r})),1]));
            end
        end
        clear stc_time_tr % save RAM
    end
    
%     % Operations
%     for r = 1:numel(ROIs)
%         % Average over time:
%         singletrialbase = mean(stc_time_rois.(sprintf('%s_pow',ROIs{r})),3);
%         % Single trial normalization (for pow only)
%         stc_time_rois.(sprintf('%s_pow',ROIs{r})) = stc_time_rois.(sprintf('%s_pow',ROIs{r}))-...
%             repmat(singletrialbase, [1,1,length(stc_time_rois.time)]);
%     end
    
    % Save activations
    save(fullfile(N.searchFolder_3arch_rej,sprintf('%s_ERPsources_ROIs',N.epochedFile(1:end-12))), 'stc_time_rois', '-v7.3');
    clear stc_time_rois
    
    %load(fullfile(N.searchFolder_3arch_rej,sprintf('%s_ERPsources_ROIs',N.epochedFile(1:end-12))))
%     %% Plot
%     style = 'Condition-wise';
%     ylims = [-5,5];
%     n_samp_smooth = 25; % 250 samples = 1s
%     figure
%     subplot(2,3,1)
%     plotERProi(stc_time_rois,'PPAlh',style,ylims,n_samp_smooth);
%     subplot(2,3,2)
%     plotERProi(stc_time_rois,'RSClh',style,ylims,n_samp_smooth);
%     subplot(2,3,3)
%     plotERProi(stc_time_rois,'OPAlh',style,ylims,n_samp_smooth);
%     subplot(2,3,4)
%     plotERProi(stc_time_rois,'PPArh',style,ylims,n_samp_smooth);
%     subplot(2,3,5)
%     plotERProi(stc_time_rois,'RSCrh',style,ylims,n_samp_smooth);
%     subplot(2,3,6)
%     plotERProi(stc_time_rois,'OPArh',style,ylims,n_samp_smooth);
%     suptitle(sprintf('%s - %s - %dms moving window', subject, style,1000*n_samp_smooth/250))
    
    
    if usebase
        % Baseline activity
        for tr = 1:size(data_b.trialinfo,1)
            fprintf('Trial %d...\n',tr);
            cfg_sel.trials = tr;
            % Compute ERSPs
            data_time_tr = ft_selectdata(cfg_sel, data_time_b);
            % Computation of source activity
            stc_time_tr = ft_sourceanalysis(cfg, data_time_tr);
            clear data_time_tr % save RAM
            
            if tr == 1
                % create stc_freq_rois
                stc_time_rois_b = stc_time_tr;
                stc_time_rois_b = rmfield(stc_time_rois_b, {'inside', 'pos', 'tri', 'trial'});
                stc_time_rois_b.trialinfo = data_freq_b.trialinfo;
                % for each ROI:
                for r = 1:numel(ROIs)
                    stc_time_rois_b.(sprintf('%s_mom',ROIs{r})) = reshape(cell2mat(stc_time_tr.trial.mom(dipsROI.(ROIs{r}))),...
                        [length(dipsROI.(ROIs{r})),1,length(stc_time_tr.time)]);
                    stc_time_rois_b.(sprintf('%s_pow',ROIs{r})) = reshape(stc_time_tr.trial.pow(dipsROI.(ROIs{r}),:),...
                        [length(dipsROI.(ROIs{r})),1,length(stc_time_tr.time)]);
                    stc_time_rois_b.(sprintf('%s_noisecov',ROIs{r})) = reshape(cell2mat(stc_time_tr.trial.noisecov(dipsROI.(ROIs{r}))),...
                        [length(dipsROI.(ROIs{r})),1]);
                end
            else
                % for each ROI:
                for r = 1:numel(ROIs)
                    stc_time_rois_b.(sprintf('%s_mom',ROIs{r})) = cat(2,stc_time_rois_b.(sprintf('%s_mom',ROIs{r})),...
                        reshape(cell2mat(stc_time_tr.trial.mom(dipsROI.(ROIs{r}))),[length(dipsROI.(ROIs{r})),1,length(stc_time_tr.time)]));
                    stc_time_rois_b.(sprintf('%s_pow',ROIs{r})) = cat(2,stc_time_rois_b.(sprintf('%s_pow',ROIs{r})),...
                        reshape(stc_time_tr.trial.pow(dipsROI.(ROIs{r}),:,:),[length(dipsROI.(ROIs{r})),1,length(stc_time_tr.time)]));
                    stc_time_rois_b.(sprintf('%s_noisecov',ROIs{r})) = cat(2,stc_time_rois_b.(sprintf('%s_noisecov',ROIs{r})),...
                        reshape(cell2mat(stc_time_tr.trial.noisecov(dipsROI.(ROIs{r}))),[length(dipsROI.(ROIs{r})),1]));
                end
            end
            clear stc_time_tr % save RAM
        end
        
%         % Operations
%         for r = 1:numel(ROIs)
%             % Average over time:
%             singletrialbase = mean(stc_time_rois_b.(sprintf('%s_pow',ROIs{r})),3);
%             % Single trial normalization (for pow only)
%             stc_time_rois_b.(sprintf('%s_pow',ROIs{r})) = stc_time_rois_b.(sprintf('%s_pow',ROIs{r}))-...
%                 repmat(singletrialbase, [1,1,length(stc_time_rois_b.time)]);
%         end
        
        % Save activations
        save(fullfile(N.searchFolder_3arch_rej,sprintf('%s_baselines_ERPsources_ROIs',N.epochedFile(1:end-12))), 'stc_time_rois_b', '-v7.3');
        clear stc_time_rois_b
    end
    
    
    
    %     %% Plot on brain:
    %     %if iscell(stc.avg.mom)
    %     %    stc.avg.mom = cell2mat(stc.avg.mom);
    %     %end
    %     cfg_plot = struct();
    %     cfg_plot.funparameter = 'pow';
    %     %cfg_plot.maskparameter = 'mom';
    %     cfg_plot.method = 'surface';
    %     %cfg_plot.latency = 'all';
    %     cfg_plot.latency = [0,0.5];
    %     cfg_plot.avgovertime = 'yes';
    %     ft_sourceplot(cfg_plot, stc_avg_time);
end