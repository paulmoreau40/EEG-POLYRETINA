%clear all;
configEEGAFF_Ainhoa;
addpath(genpath('C:\Users\Ainhoa\Documents\M2\MATLAB\codeM2\mne-matlab'))

%% Parameters
% Configuration struct for ft_sourceanalyis
cfg = struct();
cfg.method = 'mne';
cfg.keeptrials = 'no'; % Only useful to keep trials info in the output struct
cfg.keepleadfield = 'yes'; % Only useful to keep leadfield info in the output struct

% Operations on the forward model (only if the leadfield is to be computed on the fly)
%cfg.reducerank = 'no'; % Reduce rank of dipole orientations if the input leadfield is free.
%cfg.backproject = 'no'; % Back project the rank reduction of the leadfield
% Depth normalization
%cfg.normalize = 'no';
%cfg.normalizeparam = 0.5;
%cfg.weight = 1;

% Parameters for mne inverse modelling
cfg.mne.noiselambda = 0; % Regularization parameter for noise covariance - only used for pre-whitening.
cfg.mne.lambda = 0; % regularization parameter
% lambda is estimated from snr if not specified.
%cfg.mne.snr = 3;
cfg.mne.prewhiten = 'no';
cfg.mne.scalesourcecov = 'no';

for subject_ind = subject_inds
    if ~exist('ALLEEG','var')
        launchEEGLAB;
    end
    
    %subject_ind = 9;
    subject = study_config.subjects(subject_ind).id;
    disp(['Subject ' subject]);
    study_config.current_subject = subject_ind;
    
    %% Folders and Files names:
    N = makeFolderFileNames(study_config, subject);
    
    %% EEG data
    EEG_b = pop_loadset('filename', N.baselineEpochedFile, 'filepath', N.searchFolder_3arch_rej);
    data_b = eeglab2fieldtrip(EEG_b, 'raw', 'none');
    
    if isempty(study_config.subjects(subject_ind).badLocElectrodes)
        chans2remove = {};
    else
        chans2remove = study_config.subjects(subject_ind).badLocElectrodes;
    end
    
    %% Noise covariance matrix
    % Just use this to have the right structure for inverse model input
    cfg_time = struct();
    cfg_time.channel = setdiff(data_b.label, chans2remove, 'stable');
    cfg_time.trials = 'all';
    cfg_time.latency = 'all';
    cfg_time.covariance = 'yes';
    cfg_time.covariancewindow = 'all';
    cfg_time.keeptrials = 'no';
    cfg_time.removemean = 'no'; % Remove mean for covariance computation
    data_b_time = ft_timelockanalysis(cfg_time, data_b);
    data_b_time.cov = eye(size(data_b_time.cov));
    %clear EEG_b data_b
    
    %     cfg_freq = struct();
    %     cfg_freq.method = 'wavelet';
    %     cfg_freq.foi = 1:40; % frequencies of interest
    %     %cfg_freq.foilim = [1,40]; %frequency band of interest
    %     cfg_freq.toi = 0:0.01:1; % the times on which the analysis windows should be centered (in seconds)
    %     cfg_freq.width = 3; %(default = 7)
    % %   cfg.gwidth     = determines the length of the used wavelets in standard
    % %                    deviations of the implicit Gaussian kernel and should
    % %                    be choosen >= 3; (default = 3)
    %     cfg_freq.output = 'pow';
    %     cfg_freq.channel = 'all';
    %     cfg_freq.trials = 'all';
    %     cfg_freq.keeptrials = 'no';
    %     cfg_freq.keeptapers = 'no';
    %     cfg_freq.pad = 'nextpow2';
    %     cfg_freq.polyremoval = -1;
    %     data_b_freq = ft_freqanalysis(cfg_freq, data_b);
    
    %figure
    %imagesc(squeeze(data_b_freq.powspctrm(50,:,:)));

    params = struct();
    params.chans2remove = chans2remove;
    params.fixed_ori = study_config.recon.fixed_ori;
    params.patch_space = study_config.recon.patch_space;
    params.cholesky = false;
    params.normalizeDepth = true;
    params.normalizeDepthParam = 0.5;
    params.weightWithPatchSize = false;
    [cfg.sourcemodel, cov, cov_chol, dipsROI] = load_sourcemodel(N, {EEG_b.chanlocs.labels}, params);
    if ~params.cholesky
        cfg.mne.sourcecov = cov;
        %cfg.mne.sourcecov = eye(size(cov));
    end
    
    %% Simulate activity based on ROI information
    stim_strength = 1e-2;
    fields = fieldnames(dipsROI);
    raw_cross_talk = zeros(numel(fields));
    corr_coeff = zeros(numel(fields));
    rel_error = zeros(numel(fields));
    for fs = 1:numel(fields)
        source = fields{fs};
        stimulation = zeros(numel(cfg.sourcemodel.leadfield), size(data_b_time.avg,2));
        % Constant stimulation
        stimulation(dipsROI.(source),:) = stim_strength;
        % Sine stimulation
        %stimulation(dipsROI.(source),:) = sin(data_b_time.time);
        % Add noise
        %stimulation = stimulation + stim_strength * randn(size(stimulation))/10;
        % Merge ROIs from both sides ?
        simulation = cell2mat(cfg.sourcemodel.leadfield)*stimulation;
        
        data_b_time.avg = simulation;
        data_b_time.var = zeros(size(data_b_time.var));
        % Computation of source activity from the simulation
        stc = ft_sourceanalysis(cfg, data_b_time);
        
        % Transform mom into activations
        activations = nan(size(stc.avg.pow));
        for d=1:size(activations,1)
            if ~isempty(stc.avg.mom{d})
                activations(d,:) = stc.avg.mom{d};
            end
        end
        
        for ft = 1:numel(fields)
            target = fields{ft};
            % Simple mean over the ROI
            raw_cross_talk(fs,ft) = mean(activations(dipsROI.(target),:),'all', 'omitnan');
            % Metrics from Babiloni et al. 2007
            ROIstimulation = mean(stimulation(dipsROI.(target),:),1);
            ROIactivation = mean(activations(dipsROI.(target),:),1);
            % Correlation coefficient
            corr_coeff(fs,ft) = dot(ROIstimulation, ROIactivation)/(norm(ROIstimulation)*norm(ROIactivation));
            % Relative error
            rel_error(fs,ft) = norm(ROIstimulation-ROIactivation)/norm(ROIstimulation);
        end
    end
    
    % Normalize CT matrix
    cross_talk_abs = abs(raw_cross_talk)./stim_strength;
    cross_talk_rel = zeros(numel(fields));
    for f = 1:numel(fields)
        cross_talk_rel(f,:) = abs(raw_cross_talk(f,:))./max(raw_cross_talk(f,:));
    end
    
    %     figure
    %     subplot(1,2,1)
    %     hold on
    %     imagesc(1:numel(fields),1:numel(fields),cross_talk_abs, [0,1]);
    %     colormap('gray');
    %     xline(3.5,'--r', 'LineWidth', 2);
    %     yline(3.5,'--r', 'LineWidth', 2);
    %     axis image ij
    %     xticklabels(fields)
    %     xlabel('Receiving area')
    %     yticklabels(fields)
    %     ylabel('Seed area')
    %     colorbar;
    %     title('Cross-talk matrix scaled to original stimulation')
    %
    %     subplot(1,2,2)
    %     hold on
    %     imagesc(1:numel(fields),1:numel(fields),cross_talk_rel, [0,1]);
    %     colormap('gray');
    %     xline(3.5,'--r', 'LineWidth', 2);
    %     yline(3.5,'--r', 'LineWidth', 2);
    %     axis image ij
    %     xticklabels(fields)
    %     xlabel('Receiving area')
    %     yticklabels(fields)
    %     ylabel('Seed area')
    %     colorbar;
    %     title('Cross-talk matrix scaled to maximum activity recovered')
    %
    %     suptitle(sprintf('Subject %s - lambda=%.2f - depthnorm=%.2f', subject, cfg.mne.lambda, params.normalizeDepthParam))
    
    figure
    hold on
    imagesc(1:numel(fields),1:numel(fields),cross_talk_rel, [0,1]);
    colormap('gray');
    xline(3.5,'--r', 'LineWidth', 2);
    yline(3.5,'--r', 'LineWidth', 2);
    axis image ij
    xticklabels(fields)
    xlabel('Receiving area')
    yticklabels(fields)
    ylabel('Seed area')
    colorbar;
    title({'Cross-talk matrix scaled to maximum activity recovered',...
        sprintf('Subject %s - lambda=%.2f - depthnorm=%.2f', subject, cfg.mne.lambda, params.normalizeDepthParam)})
    
    % If option 1 was taken multiply solution with corr_chol
    %     if study_config.recon.fixed_ori
    %         load(fullfile(N.bemFolder, N.corrCholFile), 'corr_chol');
    %         % Wrong below: multiply mom instead of pow
    %         stc.avg.pow = corr_chol * stc.avg.pow;
    %         clear corr_chol
    %     else
    %         error('Not implemented')
    %     end
    
    
    %% Plot:
    %     if iscell(stc.avg.mom)
    %         stc.avg.mom = cell2mat(stc.avg.mom);
    %     end
    %     cfg_plot = struct();
    %     cfg_plot.funparameter = 'mom';
    %     %cfg_plot.maskparameter = 'mom';
    %     cfg_plot.method = 'surface';
    %     cfg_plot.latency = 'all';
    %     cfg_plot.avgovertime = 'yes';
    %     ft_sourceplot(cfg_plot, stc);
end