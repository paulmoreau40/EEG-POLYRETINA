function cfg = importConfig(study_config, func)
cfg = struct();
switch func
    case 'ft_sourceanalysis'
        % Operations on the forward model (only if the leadfield is to be computed on the fly)
        cfg.reducerank = 'no'; % Reduce rank of dipole orientations if the input leadfield is free.
        cfg.backproject = 'no'; % Back project the rank reduction of the leadfield
        % Depth normalization ('no' because done directly when loading the sourcemodel)
        cfg.normalize = 'no';
        
        cfg.method = study_config.recon.inv.method;
        switch cfg.method
            case 'mne'
                cfg.mne.noiselambda = study_config.recon.inv.noiselambda;
                
                regparam = 0;
                if isfield(study_config.recon.inv,'lambda')
                    cfg.mne.lambda = study_config.recon.inv.lambda;
                    regparam = regparam+1;
                end
                if isfield(study_config.recon.inv,'snr')
                    cfg.mne.lambda = study_config.recon.inv.snr;
                    regparam = regparam+1;
                end
                % Check that only 1 parameter was specified;
                if regparam==2
                    error('You cannot specify lambda and snr at the same time');
                elseif regparam == 0
                    warning('No regularization parameter set by study_config');
                end
                
                cfg.mne.prewhiten = study_config.recon.inv.prewhiten;
                cfg.mne.scalesourcecov = study_config.recon.inv.scalesourcecov;
            otherwise
                error('Not implemented')
        end
        
        cfg.keeptrials = 'no'; % Only useful to keep trials info in the output struct
        cfg.keepleadfield = 'no'; % Only useful to keep leadfield info in the output struct
        
    case 'load_sourcemodel'
        cfg.fixed_ori = study_config.recon.fixed_ori;
        cfg.patch_space = study_config.recon.patch_space;
        cfg.cholesky = study_config.recon.fwd.cholesky;
        cfg.normalizeDepth = study_config.recon.fwd.normalizeDepth;
        cfg.normalizeDepthParam = study_config.recon.fwd.normalizeDepthParam;
        cfg.weightWithPatchSize = study_config.recon.fwd.weightWithPatchSize;
        
    case {'ft_freqanalysis-main', 'ft_freqanalysis-base'}
        cfg.method = study_config.ersps.method;
        cfg.foi = study_config.ersps.FoI; % frequencies of interest
        
        % the times on which the analysis windows should be centered (in seconds)
        if strcmp(func,'ft_freqanalysis-main')
            switch study_config.epochs.window
                case 'fixed'
                    cfg.toi = study_config.epochs.limits_wdw(1):...
                        study_config.ersps.timeRes:...
                        study_config.epochs.limits_wdw(2);
                otherwise
                    error('Not implemented')
            end
        else
            cfg.toi = 0:study_config.ersps.timeRes:1;
        end
        
        switch cfg.method
            case 'superlet'
                cfg.superlet.basewidth = study_config.ersps.superlet.basewidth;
                cfg.superlet.gwidth = study_config.ersps.superlet.gwidth;
                cfg.superlet.combine = study_config.ersps.superlet.combine;
                cfg.superlet.order = study_config.ersps.superlet.order;
            otherwise
                error('Not implemented')
        end
        cfg.output = 'fourier';
        cfg.keeptrials = 'yes';
        cfg.keeptapers = 'yes';
        cfg.pad = 'nextpow2';
        cfg.polyremoval = -1;
        
    case {'ft_timelockanalysis-main', 'ft_timelockanalysis-base'}
        if strcmp(func,'ft_timelockanalysis-main')
            switch study_config.epochs.window
                case 'fixed'
                    cfg.latency = [study_config.epochs.limits_wdw(1),...
                        study_config.epochs.limits_wdw(2)];
                otherwise
                    error('Not implemented')
            end
            cfg.keeptrials = 'yes';
            cfg.removemean = 'no';
        else
            cfg.latency = [0,1];
            cfg.keeptrials = 'no';
            cfg.removemean = 'yes'; % Remove mean for covariance computation
        end
        
        cfg.covariance = 'yes';
        cfg.covariancewindow = 'all';
        
    otherwise
        error('Not implemented')
end
end