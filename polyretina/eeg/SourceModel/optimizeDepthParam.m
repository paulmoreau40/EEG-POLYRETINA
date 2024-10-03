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
cfg.mne.lambda = 5e5; % regularization parameter
% lambda is estimated from snr if not specified.
%cfg.mne.snr = 3;
cfg.mne.prewhiten = 'no';
cfg.mne.scalesourcecov = 'no';

depth_optim = nan(1,length(study_config.subjects));
CT_optim = nan(6,6,length(study_config.subjects));
dev_optim = nan(1,length(study_config.subjects));
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
    
    % 2 loop rounds
    for loop = 1:2
        if loop == 1
            % Parameters for the first optimization loop
            % grid search between 0.1 and 1 with 10 values
            lowerbound = 0.1;
            upperbound = 1;
            nb_values = 10; % including lower and upper bound
            
            % Define values to test
            step = (upperbound - lowerbound)/(nb_values-1);
            depthValues = lowerbound:step:upperbound;
            
            % No previous values to take into account
            depthValues_old = [];
            CTValues_old = [];
            devValues_old = [];
        else
            % Save previous values
            depthValues_old = depthValues;
            CTValues_old = CTValues;
            devValues_old = devValues;
            % Search for the minimum
            [~,m] = min(devValues);
            % Define new values to test
            if m > 1
                lowerbound = depthValues(m-1);
            else
                warning('Minimum value found for lowest parameter value tested: Consider enlarging the range of tested values');
                lowerbound = depthValues(m);
            end
            if m < nb_values
                upperbound = depthValues(m+1);
            else
                warning('Maximum value found for highest parameter value tested: Consider enlarging the range of tested values');
                upperbound = depthValues(m);
            end
            nb_values = 9;
            step = (upperbound - lowerbound)/(nb_values-1);
            depthValues = lowerbound:step:upperbound;
        end
        CTValues = nan(6,6,nb_values);
        devValues = nan(1,nb_values);
        
        %% Compute cross-talk matrices
        for it = 1:nb_values
            if ~isempty(depthValues_old) && any(depthValues_old == depthValues(it))
                CTValues(:,:,it) = CTValues_old(:,:,depthValues_old == depthValues(it));
                devValues(it) = devValues_old(depthValues_old == depthValues(it));
            end
            
            if isnan(devValues(it))
                fprintf('Computing for depth param = %.3f\n',depthValues(it));
                % Load source model
                params.normalizeDepthParam = depthValues(it);
                [cfg.sourcemodel, cov, dipsROI] = load_sourcemodel(N, {EEG_b.chanlocs.labels}, params);
                if ~params.cholesky
                    cfg.mne.sourcecov = cov;
                end
                
                % Simulate activity based on ROI information
                stim_strength = 1e-2;
                fields = fieldnames(dipsROI);
                raw_cross_talk = zeros(numel(fields));
                for fs = 1:numel(fields)
                    source = fields{fs};
                    stimulation = zeros(numel(cfg.sourcemodel.leadfield), size(data_b_time.avg,2));
                    % Constant stimulation
                    stimulation(dipsROI.(source),:) = stim_strength;
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
                    end
                end
                
                % Normalize CT matrix
                cross_talk_rel = zeros(numel(fields));
                for f = 1:numel(fields)
                    cross_talk_rel(f,:) = abs(raw_cross_talk(f,:))./max(raw_cross_talk(f,:));
                end
                CTValues(:,:,it) = cross_talk_rel;
                
                % Quantity to minimize in the optimization
                deviation = 10*(numel(fields) - trace(cross_talk_rel))... % Compare diagonal to ideal diagonal
                    + sum(cross_talk_rel(:)) - trace(cross_talk_rel); % Sum of non-diag elements
                devValues(it) = deviation;
                
                % Figure for visualization of the cross-talk
%                                 figure
%                                 hold on
%                                 imagesc(1:numel(fields),1:numel(fields),cross_talk_rel, [0,1]);
%                                 colormap('gray');
%                                 xline(3.5,'--r', 'LineWidth', 2);
%                                 yline(3.5,'--r', 'LineWidth', 2);
%                                 axis image ij
%                                 xticklabels(fields)
%                                 xlabel('Receiving area')
%                                 yticklabels(fields)
%                                 ylabel('Seed area')
%                                 colorbar;
%                                 title({'Cross-talk matrix scaled to maximum activity recovered',...
%                                     sprintf('Subject %s - lambda=%.2f - depthnorm=%.3f', subject, cfg.mne.lambda, params.normalizeDepthParam)})
            end
        end
        
        % Merge all values computed for this subject
        [depthValues, inds, ~] = unique([depthValues, depthValues_old],'sorted');
        CTValues = cat(3,CTValues,CTValues_old);
        CTValues = CTValues(:,:,inds);
        devValues = [devValues,devValues_old];
        devValues = devValues(inds);
        
        if loop == 2
            figure
            plot(depthValues, devValues, '-o');
            xlim([0 1])
            xlabel('Depth parameter')
            ylabel('Deviation value')
            title(sprintf('Subject %s', subject))
            saveCurrentFig([study_config.figures_folder, 'SourceModelOptim', filesep],...
                sprintf('%s_depthOptim_curve', subject), {'png'}, [])
        end
    end
    
    [dev_optim(subject_ind),m] = min(devValues);
    depth_optim(subject_ind) = depthValues(m);
    CT_optim(:,:,subject_ind) = CTValues(:,:,m);
    
    figure
    hold on
    imagesc(1:numel(fields),1:numel(fields),CT_optim(:,:,subject_ind), [0,1]);
    colormap('gray');
    xline(3.5,'--r', 'LineWidth', 2);
    yline(3.5,'--r', 'LineWidth', 2);
    axis image ij
    xticklabels(fields)
    xlabel('Receiving area')
    yticklabels(fields)
    ylabel('Seed area')
    colorbar;
    title({sprintf('Optimized cross-talk matrix for subject %s', subject),...
        sprintf('Optimum depth Param = %.3f - Optimum deviation = %.3f', depth_optim(subject_ind), dev_optim(subject_ind))})
    saveCurrentFig([study_config.figures_folder, 'SourceModelOptim', filesep],...
        sprintf('%s_depthOptim_CT', subject), {'png'}, [])
end