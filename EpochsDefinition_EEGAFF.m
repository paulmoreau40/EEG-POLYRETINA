%clear all;
configEEGAFF_Ainhoa;

%% Epoch the data for further analyses
% Options for this script
recompute = false;
do_baseline = true;
rm_base = false;

N = makeFolderFileNames(study_config, 'CHE02');
fname = 'AllObservationsReport.mat';
load(fullfile(N.searchFolder_2, fname))
obsDurations_ms = [allObservations(:).duration_ms];

if ~exist('ALLEEG','var')
    launchEEGLAB;
end

%% Goal epoch structure
% Define a common epoch window for all:
switch lower(study_config.epochs.window)
    case 'full'
        % Definition depending on the durations (assumes timewarping later)
        obsDurations_ms_usedTrials = obsDurations_ms([allObservations.completeEEG]);
        
        % Take a 1s margin with respect to min and max values
        epoch_wdw = [-1,1].*study_config.epochs.bumper +...
            [-max(obsDurations_ms_usedTrials(:,1)), max(sum(obsDurations_ms_usedTrials(:,2:end),2))]/1000;
        
    case 'fixed'
        % Definition with no timewarping implied
        epoch_wdw = [-1,1].*study_config.epochs.bumper + study_config.epochs.limits_wdw;
    otherwise
        error('Unknown window type')
end

%% Load the data, filter the data, epoch the data, keep only the right epochs, save the data
subjects = {study_config.subjects(subject_inds).id};
n_sbjs = length(subjects);
for s = 1:n_sbjs
    disp(['Subject ' subjects{s}]);
    N = makeFolderFileNames(study_config, subjects{s});
    output_filepath = N.searchFolder_3arch_rej;
    if ~exist(output_filepath, 'dir')
        mkdir(output_filepath);
    end
    
    if ~recompute && exist(fullfile(output_filepath,N.epochedFile),'file')...
            && (exist(fullfile(output_filepath,N.baselineEpochedFile),'file') || ~do_baseline)
        continue
    end
    
    EEG = pop_loadset('filename', N.postLabelingFile, 'filepath', N.searchFolder_2arch_rej);
    
    % Filter the dataset with the desired upper and lower frequencies
    % (definitive changes before epoching)
    lowcutoff = study_config.filterAnalysis.low_cut_off;
    highcutoff = study_config.filterAnalysis.high_cut_off;
    fprintf('Filtering between %.1f Hz and %.1f Hz...\n', lowcutoff, highcutoff)
    [EEG] = custom_filter(EEG, lowcutoff, highcutoff);
    % No need for line Noise removal here (LP filtering)
    
    % Take the trials of interest for this subject only
    Obsv_subj = allObservations(strcmp({allObservations.Sid},subjects{s}) & ...
        [allObservations.completeEEG]);
    obsDurations_ms_subj = obsDurations_ms(strcmp({allObservations.Sid},subjects{s}) &...
        [allObservations.completeEEG]);
        Base_subj = allObservations(strcmp({allObservations.Sid},subjects{s}) & ...
        [allObservations.completeEEG_base]);
    
    %% Epoch Data:
    boundaryEvents_mask = strcmp({EEG.event.type}, 'boundary');
    % delete boundary events
    EEG2 = pop_editeventvals(EEG, 'delete', find(boundaryEvents_mask));
    clear EEG
    
    % add a bunch of zeros at the end of the data to avoid reject the last
    % epochs because of insufficient length
    EEG2.pnts = EEG2.pnts + round(epoch_wdw(2)*EEG2.srate);
    EEG2.xmax = EEG2.xmax + epoch_wdw(2)*1000;
    EEG2.times = EEG2.xmin:(1000/EEG2.srate):EEG2.xmax;
    EEG2.data = cat(2, EEG2.data, zeros(EEG2.nbchan, round(epoch_wdw(2)*EEG2.srate)));
    EEG2.icaact = [];
    EEG2 = eeg_checkset(EEG2);
    
    % Add a channel corresponding to bad samples
    EEG2.nbchan = EEG2.nbchan + 1;
    % Get rejected samples
    switch study_config.badSampsRejection
        case 'app'
            rejSamps = EEG2.etc.APP.rejectedSamples;
        case 'asr'
            rejSamps = EEG2.etc.ASR.rejectedSamples;
        case 'autoMoBI'
            rejSamps = EEG2.etc.autoMoBI.rejectedSamples;
    end
    EEG2.data = cat(1, EEG2.data, [rejSamps, ones(1,round(epoch_wdw(2)*EEG2.srate))]);
    EEG2.chanlocs(EEG2.nbchan) = struct('labels', 'rejSamps', 'type', 'MASK', 'ref', '', 'urchan', [],...
        'X', [], 'Y', [], 'Z', [], 'sph_theta', [], 'sph_phi', [], 'sph_radius', [], 'theta', [], 'radius', []);
    
    % Epoch the dataset according to observation events:
    epochEvents = find(logical(sum([Obsv_subj.urevent] == [EEG2.event.urevent]',2)));
    EEG_epoched = pop_epoch(EEG2, {}, epoch_wdw, 'eventindices', epochEvents, 'epochinfo', 'yes'); % time is indicated in seconds
    
    EEG_epoched.etc.epochInfo.Observations = Obsv_subj;
    EEG_epoched.etc.epochInfo.Durations_ms = obsDurations_ms_subj;
    % copy bad samples information and remove the corresponding channel
    EEG_epoched.etc.epochInfo.rejectedSamples = logical(squeeze(EEG_epoched.data(end,:,:)))';
    EEG_epoched.nbchan = EEG_epoched.nbchan - 1;
    EEG_epoched.data = EEG_epoched.data(1:end-1,:,:);
    EEG_epoched.chanlocs = EEG_epoched.chanlocs(1:end-1);
    EEG_epoched.icaact = [];
    EEG_epoched = eeg_checkset(EEG_epoched);
    
    pop_saveset(EEG_epoched, 'filename', N.epochedFile,'filepath', output_filepath);
    clear EEG_epoched
    
    if do_baseline
        % Epoch the dataset according to baseline events:
        epochEvents = find(logical(sum([Base_subj.urevent_base] == [EEG2.event.urevent]',2)));
        % Dark period is always timed to 1s in the experiment
        epoch_wdw_base = [-1,1].*study_config.epochs.bumper + [0,1];
        EEG_epoched_base = pop_epoch(EEG2, {}, epoch_wdw_base, 'eventindices', epochEvents, 'epochinfo', 'yes'); % time is indicated in seconds
        
        EEG_epoched_base.etc.epochInfo.Baselines = Base_subj;
        % copy bad samples information and remove the corresponding channel
        EEG_epoched_base.etc.epochInfo.rejectedSamples = logical(squeeze(EEG_epoched_base.data(end,:,:)))';
        EEG_epoched_base.nbchan = EEG_epoched_base.nbchan - 1;
        EEG_epoched_base.data = EEG_epoched_base.data(1:end-1,:,:);
        EEG_epoched_base.chanlocs = EEG_epoched_base.chanlocs(1:end-1);
        EEG_epoched_base.icaact = [];
        EEG_epoched_base = eeg_checkset(EEG_epoched_base);
        
        pop_saveset(EEG_epoched_base, 'filename', N.baselineEpochedFile, 'filepath', output_filepath);
        clear EEG_epoched_base
    end
    clear EEG2
end