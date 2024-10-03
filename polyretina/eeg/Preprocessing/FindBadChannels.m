function [ noisyOut ] = FindBadChannels(EEG, cfg)
% Search for noisy channels
% Includes some preparatory processing steps to improve the performance of the detection:
% 1. HP filter the data
% 2. Remove Line Noise
% 3. Re-reference the channels in case of APP pipeline

subject = cfg.subjects(cfg.current_subject).id;

% filter
lowcutoff = cfg.filterPreProc.low_cut_off;
highcutoff = cfg.filterPreProc.high_cut_off;
if ~isempty(lowcutoff)
    fprintf('Highpass Filtering (%.1f Hz) for automatic bad channel detection...\n', lowcutoff)
end
if ~isempty(highcutoff)
    fprintf('Lowpass Filtering (%.1f Hz) for automatic bad channel detection...\n', highcutoff)
end
[EEG_HP] = custom_filter(EEG, lowcutoff, highcutoff);

% Add path to prepPipeline subdirectories if not in the list
tmp = which('getPipelineDefaults');
if isempty(tmp)
    myPath = fileparts(which('prepPipeline'));
    addpath(genpath(myPath));
end

% Remove Line Noise with PREP pipeline functions
disp('Removing Line Noise...')
[EEG_HP_noLN, lineNoiseOut] = removeLineNoise_custom(EEG_HP, cfg.lineNoiseRemoval_method, false);
clear EEG_HP
% If you want to save the filteredEEG_noLN struct with the LineNoiseRemoval information for later:
%EEG_HP_noLN.etc.lineNoiseRemoval = lineNoiseOut;

% Remove out-of-interest data segments
% disp('Cropping out-of-interest data...')
% [OoI_segments_index] = select_data_of_interest(EEG_HP_noLN, 'fulldata');
% EEG_ready4BadChans = eeg_eegrej(EEG_HP_noLN, OoI_segments_index);
% clear filteredEEG_noLN

% use prep pipeline subfunction to determine noisy channels
disp('Detecting bad channels automatically')
% Get the parameters (default for now)
noisyIn=SetNoisyChansDetectionStruct(EEG_HP_noLN, cfg.globalArchitecture);

switch cfg.globalArchitecture
    case 'bemobil'
        noisyOut = findNoisyChannels(EEG_HP_noLN, noisyIn);
        
        if ~isempty(cfg.subjects(cfg.current_subject).badElectrodes)
            badElectrodes = cfg.subjects(cfg.current_subject).badElectrodes;
            for ch = 1:numel(badElectrodes)
                i = find(strcmp({noisyOut.channelLocations.labels}, badElectrodes{ch}));
                if find(noisyOut.noisyChannels.badChannelsFromNoData == i)
                    fprintf('%s already flagged as bad by PREP pipeline\n', badElectrodes{ch})
                else
                    fprintf('Adding %s to NoData channels\n', badElectrodes{ch})
                    noisyOut.noisyChannels.badChannelsFromNoData(end+1) = i;
                    noisyOut.noisyChannels.badChannelsFromNoData = sort(noisyOut.noisyChannels.badChannelsFromNoData);
                    noisyOut.noisyChannels.all(end+1) = i;
                    noisyOut.noisyChannels.all = unique(noisyOut.noisyChannels.all); % unique sorts the channels too
                end
            end
        end
        
        % Do some plots:
        plot_distribution(noisyOut.robustChannelDeviation, noisyOut.noisyChannels.badChannelsFromDeviation,...
            'Channel', 'PREP', noisyOut.robustDeviationThreshold)
        xlabel(['Robust zscore for robust standard deviation'])
        title({[subject ' - Distribution of the deviation criterion'],...
            'for the channels inspected by PREP'})
        saveCurrentFig([cfg.figures_folder 'PREP_distributions' filesep],...
            [subject '_channels_deviation'], {'png'}, [600 500]);
        
        plot_distribution(noisyOut.zscoreHFNoise, noisyOut.noisyChannels.badChannelsFromHFNoise,...
            'Channel', 'PREP', noisyOut.highFrequencyNoiseThreshold)
        xlabel(['Robust zscore for robust estimate of the power ratio HF/LF'])
        title({[subject ' - Distribution of the HF-noisiness criterion'],...
            'for the channels inspected by PREP'})
        saveCurrentFig([cfg.figures_folder 'PREP_distributions' filesep],...
            [subject '_channels_HF_noisy'], {'png'}, [600 500]);
        
        plot_distribution(noisyOut.medianMaxCorrelation, noisyOut.noisyChannels.badChannelsFromCorrelation,...
            'Channel', 'PREP')
        xlabel(['Median maximum correlation over ' num2str(noisyOut.correlationWindowSeconds) 's windows'])
        title({[subject ' - Distribution of the correlation criterion'],...
            'for the channels inspected by PREP'})
        saveCurrentFig([cfg.figures_folder 'PREP_distributions' filesep],...
            [subject '_channels_correlation'], {'png'}, [600 500]);
        
        plot_distribution(noisyOut.ransacBadWindowFraction, noisyOut.noisyChannels.badChannelsFromRansac,...
            'Channel', 'PREP', noisyOut.ransacUnbrokenTime)
        xlabel(['Fraction of windows (' num2str(noisyOut.ransacWindowSeconds) 's) failing to correlate with the RANSAC prediction'])
        title({[subject ' - Distribution of the RANSAC criterion'],...
            'for the channels inspected by PREP'})
        saveCurrentFig([cfg.figures_folder 'PREP_distributions' filesep],...
            [subject '_channels_RANSAC'], {'png'}, [600 500]);
        
        
    case 'APP'
        error('APP pipeline not ready in this function')
        % Review following code
        % + Add channels from bad cfg.badElectrodes !!
        %% Re-reference to the biweight mean of the channels
        % From APP comments: Because according to Biosemi the CMS electrode does not provide the full 80 dB CMRR
        
        % Avoid considering EOG for the computation of the biweight
        all_chans = 1:EEG_ready4BadChans.nbchan;
        EEG_chans=setdiff(all_chans, all_chans(strcmp({EEG_ready4BadChans.chanlocs.type},'EOG')));
        [avg_ch,~] = Biweight_custom(EEG_ready4BadChans.data(EEG_chans,:)',...
            noisyIn.samples, cfg.censorBiweight);
        
        EEG_ready4BadChans.data = EEG_ready4BadChans.data - repmat(avg_ch',EEG_ready4BadChans.nbchan,1);
        clear avg_ch
        
        %% Apply the APP pipeline measurements
        noisyOut = findNoisyChannelsAPP(EEG_ready4BadChans, noisyIn, cfg);
end
end

