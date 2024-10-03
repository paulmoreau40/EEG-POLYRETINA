%% Start %%
%clear all;
config4SC_Alex;
addpath(genpath('C:\Users\Alexandre\Desktop\code\mne-matlab'))

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
%cfg.mne.lambda = 5e5; % regularization parameter set later in the loop
% lambda is estimated from snr if not specified.
%cfg.mne.snr = 3;
cfg.mne.prewhiten = 'no';
cfg.mne.scalesourcecov = 'no';

recompute = false;
if ~recompute
    fname = 'optimizationReport.mat';
    load(fullfile([study_config.figures_folder, 'SourceModelOptim'], fname))
end

fields_report = {'Sid','optimDepth','optimReg','optimDev','optimCT',...
    'optimSDepth','optimSReg','optimSDev','optimSCT',...
    'depthValues','regValues','regBounds','Tolerance',...
    'allGeneralDeviations','allDevSpread','allDevPerField','allCrossTalks'};
if ~exist('optDepthWithLambda', 'var')
    optDepthWithLambda = cell2struct(cell(numel(fields_report),1), fields_report', 1);
    lastLine = 0;
else
    lastLine = size(optDepthWithLambda,2);
end

for subject_ind = subject_inds
    if ~exist('ALLEEG','var')
        launchEEGLAB;
    end
    
    %subject_ind = 9;
    subject = study_config.subjects(subject_ind).id;
    disp(['Subject ' subject]);
    study_config.current_subject = subject_ind;
    
    if any(strcmp({optDepthWithLambda.Sid}, subject))
        continue
    else
        lastLine = size(optDepthWithLambda,2);
    end
    
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
    cfg_time.latency = [0 0.1];
    cfg_time.covariance = 'yes';
    cfg_time.covariancewindow = 'all';
    cfg_time.keeptrials = 'no';
    cfg_time.removemean = 'no'; % Remove mean for covariance computation
    data_b_time = ft_timelockanalysis(cfg_time, data_b);
    data_b_time.cov = eye(size(data_b_time.cov));
    %clear EEG_b data_b
    
    %% Params struct for loading the source model
    params = struct();
    params.chans2remove = chans2remove;
    params.fixed_ori = study_config.recon.fixed_ori;
    params.patch_space = study_config.recon.patch_space;
    params.cholesky = false;
    params.normalizeDepth = true;
    params.weightWithPatchSize = false;
    
    % Optimization of depth parameter
    % grid search between 0.1 and 1
    depthValues = 0.5:0.05:1;
    
    % Optimization of regularization parameter
    % grid search between 0.01 and 1000 (log scale)
    regValues = 10.^(-4:2);
    
    %% Compute cross-talk matrices
    stim_strength = 1e-2;
    chanLabels = {EEG_b.chanlocs.labels};
    % Load dipsROI
    params.normalizeDepthParam = depthValues(1);
    [~, ~, ~, dipsROI] = load_sourcemodel(N, chanLabels, params);
    fields = fieldnames(dipsROI);
    
    CTValues = nan(numel(fields),numel(fields),length(depthValues),length(regValues));
    devValues = nan(length(depthValues),length(regValues));
    devPerFieldValues = nan(numel(fields),length(depthValues),length(regValues));
    devSpreadValues = nan(length(depthValues),length(regValues));
    regBounds = nan(2,length(depthValues));
    
    poolobj = gcp('nocreate');
    if isempty(poolobj)
        parpool('local',6);
    end
    ppm = ParforProgressbar(length(depthValues), 'showWorkerProgress', true,...
        'title', sprintf('%s - Finding optimal values for Cross-Talk matrix', subject));
    parfor d = 1:length(depthValues)
        params_loop = params;
        cfg_loop = cfg;
        cfg_loop.mne.feedback = 'no';
        data_loop = data_b_time;
        
        CTValues_loop = nan(numel(fields),numel(fields),length(regValues));
        devValues_loop = nan(1,length(regValues));
        devPerFieldValues_loop = nan(numel(fields),length(regValues));
        devSpreadValues_loop = nan(1,length(regValues));
        regBounds_loop = nan(2,1);
        
        %fprintf('Computing for depth param = %.3f\n',depthValues(it));
        % Load source model
        params_loop.normalizeDepthParam = depthValues(d);
        [cfg_loop.sourcemodel, cov, ~, ~] = load_sourcemodel(N, chanLabels, params_loop);
        if ~params_loop.cholesky
            cfg_loop.mne.sourcecov = cov;
        end
        
        G = cell2mat(cfg_loop.sourcemodel.leadfield);
        G_non0 = G(:,all(G,1));
        [~,S,~] = csvd(G_non0);
        regBounds_loop(1) = min(S);
        regBounds_loop(2) = max(S);
        
        for r = 1:length(regValues)
            cfg_loop.mne.lambda = regValues(r);
            
            % Simulate activity based on ROI information
            raw_cross_talk = zeros(numel(fields));
            for fs = 1:numel(fields)
                source = fields{fs};
                stimulation = zeros(numel(cfg_loop.sourcemodel.leadfield), size(data_loop.avg,2));
                % Constant stimulation
                stimulation(dipsROI.(source),:) = stim_strength;
                simulation = cell2mat(cfg_loop.sourcemodel.leadfield)*stimulation;
                
                data_loop.avg = simulation;
                data_loop.var = zeros(size(data_loop.var));
                % Computation of source activity from the simulation
                stc = ft_sourceanalysis(cfg_loop, data_loop);
                
                % Transform mom into activations
                activations = nan(size(stc.avg.pow));
                for dip=1:size(activations,1)
                    if ~isempty(stc.avg.mom{dip})
                        activations(dip,:) = stc.avg.mom{dip};
                    end
                end
                
                for ft = 1:numel(fields)
                    target = fields{ft};
                    % Simple mean over the ROI
                    raw_cross_talk(fs,ft) = mean(activations(dipsROI.(target),:),'all', 'omitnan');
                end
            end
            
            % Normalize CT matrix
            cross_talk_rel = zeros(numel(fields));
            for f = 1:numel(fields)
                cross_talk_rel(f,:) = abs(raw_cross_talk(f,:))./max(raw_cross_talk(f,:));
            end
            CTValues_loop(:,:,r) = cross_talk_rel;
            
            % Quantity to minimize in the optimization
            deviation = 10*(numel(fields) - trace(cross_talk_rel))... % Compare diagonal to ideal diagonal
                + sum(cross_talk_rel(:)) - trace(cross_talk_rel); % Sum of non-diag elements
            devValues_loop(r) = deviation;
            
            % Details about the deviation per field
            dev_perfield = zeros(numel(fields),1);
            for f = 1:numel(fields)
                dev_perfield(f) = 10*(1-cross_talk_rel(f,f))+...
                    sum(cross_talk_rel(f,:)) + sum(cross_talk_rel(:,f)) - 2*cross_talk_rel(f,f);
            end
            devPerFieldValues_loop(:,r) = dev_perfield;
            devSpreadValues_loop(r) = mad(dev_perfield,0);
            
            % Figure for visualization of the cross-talk
            %             plotCrossTalk(cross_talk_rel, fields);
            %             title({sprintf('Cross-talk matrix for subject %s', subject),...
            %                 sprintf('Depth = %.3f - Regularization = %.0e - Deviation = %.3f',...
            %                 depthValues(d), depthValues(r), deviation)},...
            %                 'Fontsize', 14);
        end
        
        CTValues(:,:,d,:) = reshape(CTValues_loop,[numel(fields),numel(fields),1,length(regValues)]);
        devValues(d,:) = devValues_loop;
        devPerFieldValues(:,d,:) = reshape(devPerFieldValues_loop,[numel(fields),1,length(regValues)]);
        devSpreadValues(d,:) = devSpreadValues_loop;
        regBounds(:,d) = regBounds_loop;
        ppm.increment();
    end
    delete(ppm);
    
    %% Compute summary values
    [X,Y] = meshgrid(regValues,depthValues);
    aboveMinBound = X > repmat(regBounds(1,:)',1,length(regValues));
    belowMaxBound = X < repmat(regBounds(2,:)',1,length(regValues));
    inBounds = aboveMinBound & belowMaxBound;
    dev_optim = min(devValues(inBounds));
    [d_opt,r_opt] = find(devValues==dev_optim);
    
    tol = 0.05;
    tolerance_dev_values = devValues < dev_optim + tol*(max(devValues(:)) - min(devValues(:)));
    acceptable_dev_values = tolerance_dev_values & inBounds;
    dev_spread_optim = min(devSpreadValues(acceptable_dev_values));
    [ds_opt,rs_opt] = find(devSpreadValues==dev_spread_optim);
    
    %% Save info in struct
    values = {subject,depthValues(d_opt),regValues(r_opt),dev_optim,CTValues(:,:,d_opt,r_opt),...
        depthValues(ds_opt),regValues(rs_opt),devValues(ds_opt,rs_opt),CTValues(:,:,ds_opt,rs_opt),...
        depthValues, regValues, regBounds, tol,...
        devValues, devSpreadValues, devPerFieldValues, CTValues};
    
    optDepthWithLambda(lastLine + 1) = cell2struct(values', fields_report', 1);
    
    fname = 'optimizationReport.mat';
    save(fullfile([study_config.figures_folder, 'SourceModelOptim'], fname),'optDepthWithLambda')
    
    %% Figures
    figure;
    hold on;
    surf(X,Y,devValues);
    fill3([regBounds(1,:),flip(regBounds(1,:))],...
        [depthValues, flip(depthValues)],...
        reshape(repmat([0;20],1,length(depthValues))',[1,2*length(depthValues)]),...
        'r', 'FaceAlpha', 0.5);
    fill3([regBounds(2,:),flip(regBounds(2,:))],...
        [depthValues, flip(depthValues)],...
        reshape(repmat([0;20],1,length(depthValues))',[1,2*length(depthValues)]),...
        'g', 'FaceAlpha', 0.5);
    set(gca,'XScale','log')
    xlabel('Regularization parameter')
    ylabel('Depth parameter')
    zlabel('Deviation value')
    title(sprintf('Subject %s', subject))
    saveCurrentFig([study_config.figures_folder, 'SourceModelOptim', filesep],...
        sprintf('%s_depth-regOptim_DevValues_surf', subject), {'fig'}, [])
    
    figure;
    hold on;
    surf(X,Y,devSpreadValues);
    fill3([regBounds(1,:),flip(regBounds(1,:))],...
        [depthValues, flip(depthValues)],...
        reshape(repmat([0.5;1],1,length(depthValues))',[1,2*length(depthValues)]),...
        'r', 'FaceAlpha', 0.5);
    fill3([regBounds(2,:),flip(regBounds(2,:))],...
        [depthValues, flip(depthValues)],...
        reshape(repmat([0.5;1],1,length(depthValues))',[1,2*length(depthValues)]),...
        'g', 'FaceAlpha', 0.5);
    set(gca,'XScale','log')
    xlabel('Regularization parameter')
    ylabel('Depth parameter')
    zlabel('Deviation spread value')
    title(sprintf('Subject %s', subject))
    saveCurrentFig([study_config.figures_folder, 'SourceModelOptim', filesep],...
        sprintf('%s_depth-regOptim_DevSpreadValues_surf', subject), {'fig'}, [])
    
    plotCrossTalk(CTValues(:,:,d_opt,r_opt), fields);
    title({sprintf('Optimized cross-talk matrix for subject %s', subject),...
        sprintf('Opt depth = %.3f - Opt regularization = %.0e - Opt deviation = %.3f',...
        depthValues(d_opt),regValues(r_opt),dev_optim)},...
        'Fontsize', 14);
    saveCurrentFig([study_config.figures_folder, 'SourceModelOptim', filesep],...
        sprintf('%s_depth-regOptim_CToptim', subject), {'png'}, [800,800]);
    
    plotCrossTalk(CTValues(:,:,ds_opt,rs_opt), fields);
    title({sprintf('Optimized (with spread) cross-talk matrix for subject %s', subject),...
        sprintf('Opt depth = %.3f - Opt regularization = %.0e - Opt deviation = %.3f',...
        depthValues(ds_opt),regValues(rs_opt),devValues(ds_opt,rs_opt))},...
        'Fontsize', 14);
    saveCurrentFig([study_config.figures_folder, 'SourceModelOptim', filesep],...
        sprintf('%s_depth-regOptim_CToptimWithSpread', subject), {'png'}, [800,800]);
end