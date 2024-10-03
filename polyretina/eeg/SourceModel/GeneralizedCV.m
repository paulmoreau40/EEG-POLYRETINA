%clear all;
config4SC_Alex;
addpath(genpath('C:\Users\Alexandre\Desktop\code\mne-matlab'))

%% Parameters
% Configuration struct for ft_sourceanalyis
cfg = struct();
cfg.method = 'mne';
cfg.keeptrials = 'no'; % Only useful to keep trials info in the output struct
cfg.keepleadfield = 'yes'; % Only useful to keep leadfield info in the output struct

% Operations on the forward model
% cfg.reducerank = 'no'; % Reduce rank of dipole orientations if the input leadfield is free.
% cfg.backproject = 'no'; % Back project the rank reduction of the leadfield
% % Depth normalization
% cfg.normalize = 'yes';
% cfg.normalizeparam = 0.8;
% cfg.weight = 1;

% Parameters for mne inverse modelling
cfg.mne.noiselambda = 0; % Regularization parameter for noise covariance - only used for pre-whitening.
%cfg.mne.lambda = 1/9; % regularization parameter
% lambda is estimated from snr if not specified.
cfg.mne.snr = 3;
cfg.mne.prewhiten = 'no';
cfg.mne.scalesourcecov = 'no';

for subject_ind = subject_inds
    if ~exist('ALLEEG','var')
        launchEEGLAB;
    end
    
    %subject_ind = 4;
    subject = study_config.subjects(subject_ind).id;
    disp(['Subject ' subject]);
    study_config.current_subject = subject_ind;
    
    %% Folders and Files names:
    N = makeFolderFileNames(study_config, subject);
    
    %% EEG data
    if isempty(study_config.subjects(subject_ind).badLocElectrodes)
        chans2remove = {};
    else
        chans2remove = study_config.subjects(subject_ind).badLocElectrodes;
    end
    
    EEG = pop_loadset('filename', N.epochedFile, 'filepath', N.searchFolder_3arch_rej);
    data = eeglab2fieldtrip(EEG, 'raw', 'none');
    
    %% EEG preprocessing
    cfg_preproc = struct();
    cfg_preproc.channel = setdiff(data.label, chans2remove, 'stable');
    cfg_preproc.trials = 'all';
    cfg_preproc.latency = 'all';
    cfg_preproc.covariance = 'yes';
    cfg_preproc.covariancewindow = 'all';
    cfg_preproc.keeptrials = 'yes';
    cfg_preproc.removemean = 'no'; % Remove mean for covariance computation
    data_preproc = ft_timelockanalysis(cfg_preproc, data);
    clear data
    
    EEG_b = pop_loadset('filename', N.baselineEpochedFile, 'filepath', N.searchFolder_3arch_rej);
    data_b = eeglab2fieldtrip(EEG_b, 'raw', 'none');
    
    %% Noise covariance matrix
    % Calculating the average covariance matrix from each baseline trial
    cfg_base = struct();
    cfg_base.channel = setdiff(data_b.label, chans2remove, 'stable');
    cfg_base.trials = 'all';
    cfg_base.latency = 'all';
    cfg_base.covariance = 'yes';
    cfg_base.covariancewindow = 'all';
    cfg_base.keeptrials = 'no';
    cfg_base.removemean = 'yes'; % Remove mean for covariance computation
    data_b_preproc = ft_timelockanalysis(cfg_base, data_b);
    %cfg.mne.noisecov = data_b_preproc.cov; % does not work
    % Replace data covariance matrix with baseline covariance matrix
    data_preproc.cov = repmat(reshape(data_b_preproc.cov, [1,size(data_b_preproc.cov)]),...
        [size(data_preproc.trial,1),1,1]);
    clear EEG_b data_b
    
    %% Params struct for loading the source model
    params = struct();
    params.chans2remove = chans2remove;
    params.fixed_ori = study_config.recon.fixed_ori;
    params.patch_space = study_config.recon.patch_space;
    params.cholesky = true;
    params.normalizeDepth = true;
    params.weightWithPatchSize = false;
    params.normalizeDepthParam = 0.75;
    [cfg.sourcemodel, cov, dipsROI] = load_sourcemodel(N, {EEG.chanlocs.labels}, params);
    
    G = cell2mat(cfg.sourcemodel.leadfield);
    G_non0 = G(:,all(G,1));
    prewhiten = false;
    scaling = false;
    % PREWHITENING needs to happen before adjusting gain matrix with
    % correlation matrix
    if prewhiten
        P = prewhiteningMatrix(data_b_preproc.cov, 1e-0);
        G_non0 = P*G_non0;
    end
    
    if params.cholesky
        % Scale cov so that tr(G*G') = tr(Id) in the end
        % (for comparison with the regularization term)
        if scaling
            scaled_cov = cov.*(size(G_non0,1)/trace(G_non0*cov*(G_non0*cov)'));
            G_non0 = G_non0*scaled_cov;
        else
            G_non0 = G_non0*cov;
        end
        cfg.mne.sourcecov = eye(size(cov));
    else
        % Scale cov so that tr(G*cov*G') = tr(Id) in the end
        % (for comparison with the regularization term)
        if scaling
            scaled_cov = cov.*(size(G_non0,1)/trace(G_non0*cov*G_non0'));
            cfg.mne.sourcecov = scaled_cov;
        else
            cfg.mne.sourcecov = cov;
        end
    end
    
    [U,s,V] = csvd(G_non0);
    % One sample target at a time
    lambdas_optim = nan(size(data_preproc.trial,1),size(data_preproc.trial,3));
    
    fprintf('Trial ')
    for tr = 2
        fprintf('%d..',tr)
        for t = 1:size(lambdas_optim,2)
            M = squeeze(data_preproc.trial(tr,:,t))' - mean(data_b_preproc.avg,2);
            if prewhiten
                [lambdas_optim(tr,t),~,~] = gcv(U,s,P*M,'Tikh');
            else
                [lambdas_optim(tr,t),~,~] = gcv(U,s,M,'Tikh');
            end
        end
    end
    fprintf('\n')
    
    figure;
    hold on;
    histogram(lambdas_optim(tr,:))
    xline(mean(lambdas_optim(tr,:)),'-.r',...
        'Label', sprintf('Mean = %.2f', mean(lambdas_optim(tr,:))),...
        'LabelOrientation', 'horizontal');
    xlabel('Optimal lambda');
    ylabel('Count');
    title({'Distribution of optimal lambdas','1 optimum per time point (1500 in total)'});
    
    figure;
    xlabel('Time (s)');
    plot(data_preproc.time, lambdas_optim(tr,:));
    ylabel('Lambda');
    title('Evolution of lambda optimization across trial time stamps');
    
%     M_full = squeeze(data_preproc.trial(1,:,:)) - repmat(mean(data_b_preproc.avg,2),1,size(data_preproc.trial,3));
%     if prewhiten
%         figure
%         [lambda_optim,~,~] = gcv(U,s,P*M_full,'Tikh');
%         figure
%         l_curve(U,s,P*M_full)
%     else
%         figure
%         [lambda_optim,~,~] = gcv(U,s,M_full,'Tikh');
%         figure
%         l_curve(U,s,M_full)
%     end

    %% Computation of source activity
    % Remove baseline from source activity
    data_preproc.trial = data_preproc.trial - repmat(mean(data_b_preproc.avg,2)',size(data_preproc.trial,1),1,size(data_preproc.trial,3));
    %stc_avg = ft_sourceanalysis(cfg, data_preproc);
    
    cfg.rawtrial = 'yes';
    stc_trials = ft_sourceanalysis(cfg, data_preproc);
    
    %cfg.rawtrial = 'no';
    %cfg.bootstrap = 'yes';
    %cfg.numbootstrap = 50;
    %stc_bootstrapped = ft_sourceanalysis(cfg, data_preproc);
    
    % If option 1 was taken multiply solution with corr_chol
    %     if study_config.recon.fixed_ori
    %         load(fullfile(N.bemFolder, N.corrCholFile), 'corr_chol');
    %         % Wrong below: multiply mom instead of pow
    %         stc.avg.pow = corr_chol * stc.avg.pow;
    %         clear corr_chol
    %     else
    %         error('Not implemented')
    %     end
    
    %% Plot
    %     n_samp_smooth = 100; % 1s = 250 samples
    %     figure
    %     subplot(1,3,1)
    %     hold on
    %     plot(stc_avg.time,movmean(mean(stc_avg.avg.pow(dipsROI.PPAlh,:),1),n_samp_smooth))
    %     plot(stc_avg.time,movmean(mean(stc_avg.avg.pow(dipsROI.PPArh,:),1),n_samp_smooth))
    %     legend(study_config.recon.hemispheres)
    %     title('PPA')
    %     subplot(1,3,2)
    %     hold on
    %     plot(stc_avg.time,movmean(mean(stc_avg.avg.pow(dipsROI.RSClh,:),1),n_samp_smooth))
    %     plot(stc_avg.time,movmean(mean(stc_avg.avg.pow(dipsROI.RSCrh,:),1),n_samp_smooth))
    %     legend(study_config.recon.hemispheres)
    %     title('RSC')
    %     subplot(1,3,3)
    %     hold on
    %     plot(stc_avg.time,movmean(mean(stc_avg.avg.pow(dipsROI.OPAlh,:),1),n_samp_smooth))
    %     plot(stc_avg.time,movmean(mean(stc_avg.avg.pow(dipsROI.OPArh,:),1),n_samp_smooth))
    %     legend(study_config.recon.hemispheres)
    %     title('OPA')
    
    %% Plot on brain:
    %if iscell(stc.avg.mom)
    %    stc.avg.mom = cell2mat(stc.avg.mom);
    %end
    %cfg_plot = struct();
    %cfg_plot.funparameter = 'pow';
    %cfg_plot.maskparameter = 'mom';
    %cfg_plot.method = 'surface';
    %cfg_plot.latency = 'all';
    %cfg_plot.latency = [0,0.5];
    %cfg_plot.avgovertime = 'yes';
    %ft_sourceplot(cfg_plot, stc_avg);
end