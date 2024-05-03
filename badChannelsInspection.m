configEEGAFF_Alex;

for subject_ind = subject_inds
    if ~exist('ALLEEG','var')
        launchEEGLAB;
    end
    
    subject_ind = 6;
    subject = study_config.subjects(subject_ind).id;
    disp(['Subject ' subject]);
    study_config.current_subject = subject_ind;
    
    %% Folders and Files names:
    N = makeFolderFileNames(study_config, subject);
    
    EEG_inter_avRef = pop_loadset('filename', N.nobadchansFile, 'filepath', N.searchFolder_2arch);
    [ALLEEG, EEG_inter_avRef, CURRENTSET] = eeg_store(ALLEEG, EEG_inter_avRef, CURRENTSET);
    noisyChansDetect = EEG_inter_avRef.etc.noisyChannelsDetection;
    threshCorr = noisyChansDetect.maximumCorrelations < noisyChansDetect.correlationThreshold;
    fractionBad = mean(threshCorr, 2);
    badChansFromCorr = find(fractionBad > noisyChansDetect.badTimeThreshold);
    
    fractionBadPerTime = mean(threshCorr(badChansFromCorr,:), 1);
    plot(fractionBadPerTime)
    
    threshCorr = noisyChansDetect.ransacCorrelations < noisyChansDetect.ransacCorrelationThreshold;
    fractionBadPerTime = mean(threshCorr(noisyChansDetect.ransacBadWindowFraction>noisyChansDetect.ransacUnbrokenTime,:),1);
    plot(fractionBadPerTime)
    
    plot_distribution(fractionBad, badChansFromCorr,...
        'Channel', 'PREP', noisyChansDetect.badTimeThreshold)
    if max(fractionBad) < 0.1
        xlim([0,0.1])
    end
    xlabel('Fraction of windows with insufficient correlation')
    title({[subject ' - Distribution of the correlation criterion'],...
        'for the channels inspected by PREP'})
    saveCurrentFig([study_config.figures_folder 'PREP_distributions' filesep],...
        [subject '_channels_correlation_fractionBad'], {'png'}, [600 500]);
end