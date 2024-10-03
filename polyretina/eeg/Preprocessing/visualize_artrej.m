function visualize_artrej(subj, cfg, step)

switch step
    case 'ASR'
        switch lower(cfg.globalArchitecture)
            case 'simple'
                befFolder = [cfg.study_folder cfg.preprocessing_folder lower(cfg.globalArchitecture) filesep 'ASR'];
                aftFolder = befFolder;
                befFile = [subj '_' cfg.ASRin_filename];
                aftFile = [subj '_' cfg.ASRout_filename];
            case 'bemobil'
                befFolder = [cfg.study_folder cfg.preprocessing_folder lower(cfg.globalArchitecture)];
                aftFolder = [cfg.study_folder cfg.preprocessing_folder lower(cfg.globalArchitecture) filesep 'ASR'];
                befFile = [subj '_' cfg.BadChansRemoved_filename];
                aftFile = [subj '_' cfg.ASRout_filename];
        end
        
        EEG_bef = pop_loadset('filename',befFile,'filepath',befFolder);
        if strcmpi(cfg.globalArchitecture,'bemobil')
            % HP filter eeg bef for visualization
            lowcutoff = cfg.low_cut_off;
            highcutoff = [];
            [EEG_bef] = custom_filter(EEG_bef, lowcutoff, highcutoff);
        end
        
        EEG_aft = pop_loadset('filename',aftFile,'filepath',aftFolder);
        
        vis_artifacts(EEG_aft,EEG_bef);
        
    case 'BadChans' % Only for bemobil
        if ~strcmpi(cfg.globalArchitecture,'bemobil')
            error('Function not configured for that case')
        end
        befFolder = [cfg.study_folder cfg.preprocessing_folder];
        badChansFolder = [cfg.study_folder cfg.preprocessing_folder lower(cfg.globalArchitecture)];
        befFile = [subj '_' cfg.prepared_filename];
        badChansFile = [subj '_' cfg.BadChansRemoved_filename];
        
        EEG_bef = pop_loadset('filename',befFile,'filepath',befFolder);
        EEG_badChans = pop_loadset('filename',badChansFile,'filepath',badChansFolder);
        
        % HP filter eeg bef for visualization
        lowcutoff = cfg.low_cut_off;
        highcutoff = [];
        [EEG_bef] = custom_filter(EEG_bef, lowcutoff, highcutoff);
        
        EEG_aft = pop_select(EEG_bef,'nochannel',EEG_badChans.etc.noisyChannelsDetection.noisyChannels.all);
        EEG_aft.etc.clean_channel_mask = true(1,EEG_bef.nbchan);
        EEG_aft.etc.clean_channel_mask(EEG_badChans.etc.noisyChannelsDetection.noisyChannels.all) = false;
        
        figure
        topoplot(EEG_badChans.etc.noisyChannelsDetection.noisyChannels.all, EEG_badChans.chanlocs,...
            'style', 'blank', 'electrodes', 'ptslabels', 'plotdisk', 'on')
        
        vis_artifacts(EEG_aft,EEG_bef);
end
end