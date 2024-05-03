%clear all;
configEEGAFF_Ainhoa;

%% Epoch the data for further analyses
% Options for this script
recompute = false;
mocap_needed = false; %boolean to load MOCAP data or not
switch study_config.badSampsRejection
    case 'app'
        badSampsMethod = 'APP';
    case 'asr'
        badSampsMethod = 'ASR';
    case 'autoMoBI'
        badSampsMethod = 'autoMoBI';
end

%% Get the folders names right:
pipe_name = study_config.globalArchitecture;
fig_path = fullfile(study_config.figures_folder, 'Epochs');
if ~exist(fig_path, 'dir')
    mkdir(fig_path)
end

if ~recompute
    N = makeFolderFileNames(study_config, 'CHE02');
    fname = 'AllObservationsReport.mat';
    load(fullfile(N.searchFolder_2, fname))
end

for subject_ind = subject_inds
    if ~exist('ALLEEG','var')
        launchEEGLAB;
    end
    
    % clear RAM
    %STUDY = []; CURRENTSTUDY = 0; ALLEEG = []; EEG=[]; CURRENTSET=[];
    
    %subject_ind = 41;
    subject = study_config.subjects(subject_ind).id;
    study_config.current_subject = subject_ind;
    
    if strcmp(study_config.subjects(subject_ind).excluded, 'Yes')
        continue
    elseif ~recompute && sum(strcmp({allObservations.Sid}, subject))>0
        continue
    else
        N = makeFolderFileNames(study_config, subject);
        
        EEG_prepared = pop_loadset('filename', N.preparedFile, 'filepath', N.searchFolder_2, 'loadmode', 'info');
        %[ALLEEG, EEG_prepared, CURRENTSET] = eeg_store(ALLEEG, EEG_prepared, CURRENTSET);
        events = EEG_prepared.event;
        urevents = EEG_prepared.urevent;
        ObsStarts = events(strcmp({events.type},'ObsStart'));
        ObsEnds = events(strcmp({events.type},'ObsEnd'));
        BaseStarts = events(strcmp({events.type},'BaseStart'));
        BaseEnds = events(strcmp({events.type},'BaseEnd'));
        Bounds = events(strcmp({events.type},'boundary'));
        noBounds = events(~strcmp({events.type},'boundary'));
        
        Trials = EEG_prepared.etc.CompAnalysis;
        
        EEG_preproc = pop_loadset('filename', N.preICAFile, 'filepath', N.searchFolder_2arch_rej, 'loadmode', 'info');
        badSamples = EEG_preproc.etc.(badSampsMethod).rejectedSamples;
        clear EEG_preproc
        
        fields_obs = {'Sid', 'Gender', 'AgeGroup',...
            'Trial', 'Pair', 'Type', 'LeftDiff', 'RightDiff', 'AffLevel', 'AbsAff',...
            'ObsvAnswer', 'ObstAttempts', 'ObstSuccess', 'Elimination',...
            'urevent', 'duration_ms', 'completeEEG', 'cleanEEG', 'percCleanEEG',...
            'urevent_base', 'duration_ms_base', 'completeEEG_base', 'cleanEEG_base', 'percCleanEEG_base'};
        if ~exist('allObservations', 'var')
            allObservations = cell2struct(cell(numel(fields_obs),1), fields_obs', 1);
            lastLine_obs = 0;
        else
            lastLine_obs = size(allObservations,2);
        end
        N_newObs = size(Trials,1);
        
        for i = 1:N_newObs
            trl = Trials.trial(i);
            dur_ms = 1000*(ObsEnds([ObsEnds.TrialIndex] == trl).latency - ObsStarts([ObsStarts.TrialIndex] == trl).latency)/EEG_prepared.srate;
            dur_ms_b = 1000*(BaseEnds([BaseEnds.TrialIndex] == trl).latency - BaseStarts([BaseStarts.TrialIndex] == trl).latency)/EEG_prepared.srate;
            
            %% EEG completeness
            startObs = ObsStarts([ObsStarts.TrialIndex] == trl);
            stopObs = ObsEnds([ObsEnds.TrialIndex] == trl);
            if ~isempty(startObs) && ~isempty(stopObs)
                obsSamps = ceil(startObs.latency):floor(stopObs.latency);
            elseif ~isempty(startObs)
                % missing stop
                nextBnd = Bounds(find([Bounds.latency]>startObs.latency,1));
                if isempty(nextBnd)
                    % reaching the end of the recording
                    obsSamps = ceil(startObs.latency):EEG_prepared.pnts;
                else
                    obsSamps = ceil(startObs.latency):floor(nextBnd.latency);
                end
            elseif ~isempty(stopObs)
                % missing start
                prevBnd = Bounds(find([Bounds.latency]<stopObs.latency,1,'last'));
                if isempty(prevBnd)
                    % beginning of the recording
                    obsSamps = 1:floor(stopObs.latency);
                else
                    obsSamps = ceil(prevBnd.latency):floor(stopObs.latency);
                end
            else
                obsSamps = [];
            end
            
            if ~isempty(obsSamps)
                percCleanEEG = 100*(1-(sum(badSamples(obsSamps))/length(badSamples(obsSamps))));
            end
            
            complete = EEG_prepared.etc.TrialsInspection.Observations_perc(i) == 100;
            clean = complete && (percCleanEEG == 100);
            
            startBase = BaseStarts([BaseStarts.TrialIndex] == trl);
            stopBase = BaseEnds([BaseEnds.TrialIndex] == trl);
            if ~isempty(startBase) && ~isempty(stopBase)
                baseSamps = ceil(startBase.latency):floor(stopBase.latency);
            elseif ~isempty(startBase)
                % missing stop
                nextBnd = Bounds(find([Bounds.latency]>startBase.latency,1));
                if isempty(nextBnd)
                    % reaching the end of the recording
                    baseSamps = ceil(startBase.latency):EEG_prepared.pnts;
                else
                    baseSamps = ceil(startBase.latency):floor(nextBnd.latency);
                end
            elseif ~isempty(stopBase)
                % missing start
                prevBnd = Bounds(find([Bounds.latency]<stopBase.latency,1,'last'));
                if isempty(prevBnd)
                    % beginning of the recording
                    baseSamps = 1:floor(stopBase.latency);
                else
                    baseSamps = ceil(prevBnd.latency):floor(stopBase.latency);
                end
            else
                baseSamps = [];
            end
            
            if ~isempty(baseSamps)
                percCleanEEG_b = 100*(1-(sum(badSamples(baseSamps))/length(badSamples(baseSamps))));
            end
            
            complete_b = EEG_prepared.etc.TrialsInspection.BlackBaselines_perc(i) == 100;
            clean_b = complete_b && (percCleanEEG_b == 100);
            
            values = {subject, study_config.subjects(subject_ind).gender, study_config.subjects(subject_ind).group,...
                trl, Trials.pair(i), Trials.type{i}, Trials.left_difficulty{i},...
                Trials.right_difficulty{i}, Trials.AffLevel{i}, Trials.AbsAff(i),...
                Trials.Answer(i), Trials.Attempted(i), Trials.Success(i), Trials.Elimination{i},...
                ObsStarts([ObsStarts.TrialIndex] == trl).urevent, dur_ms, complete, clean, percCleanEEG...
                BaseStarts([BaseStarts.TrialIndex] == trl).urevent, dur_ms_b, complete_b, clean_b, percCleanEEG_b};
            
            allObservations(lastLine_obs + i) = cell2struct(values', fields_obs', 1);
        end
        
        % Get rid of large structures:
        clear events urevents noBounds
        clear EEG_prepared
    end
end

if recompute
    fname = 'AllObservationsReport.mat';
    save(fullfile(N.searchFolder_2, fname),'allObservations')
end

%% Description of the data
youngObsv = strcmp({allObservations.AgeGroup}, 'Y');
maleObsv = strcmp({allObservations.Gender}, 'M');
GoObsv = strcmp({allObservations.Type}, 'Go');
baseDurations_ms = [allObservations(:).duration_ms_base]';
obsDurations_ms = [allObservations(:).duration_ms]';

figure
subplot(1,2,1)
[~, edges] = histcounts(baseDurations_ms);
hold on
h1 = histogram(baseDurations_ms(~GoObsv),edges);
h2 = histogram(baseDurations_ms(GoObsv),edges);
xlabel('Duration (ms)')
title('Baseline duration')
legend([h1,h2], {'NoGo trials', 'Go trials'}, 'Location', 'Best')

subplot(1,2,2)
[~, edges] = histcounts(obsDurations_ms);
hold on
h1 = histogram(obsDurations_ms(~GoObsv),edges);
h2 = histogram(obsDurations_ms(GoObsv),edges);
xlabel('Duration (ms)')
title('Observation duration')
legend([h1,h2], {'NoGo trials', 'Go trials'}, 'Location', 'Best')

%% Analyze the data in the table to get the trials you want:
completeObsPrep = [allObservations.completeEEG];
completeObsPreproc = [allObservations.cleanEEG];
completeObsPrepBase = [allObservations.completeEEG_base];
completeObsPreprocBase = [allObservations.cleanEEG_base];

figure
subjects = unique({allObservations.Sid},'stable');
n_sbjs = length(subjects);
countSubjObs = zeros(n_sbjs+1,3);
countSubjBase = zeros(n_sbjs+1,3);
countSubjObsBase = zeros(n_sbjs+1,3);
sep = false;
% n_sbjs = n_sbjs+1 %if you want to have a separation between age groups
for s = 1:n_sbjs
    s_ind = find(strcmp({study_config.subjects.id}, subjects(s)),1);
    if s>1 && ~sep && ...
            ~strcmp(study_config.subjects(s_ind).group, study_config.subjects(s_ind-1).group)
        subjects(s:end+1) = {'',subjects{s:end}};
        sep = true;
        continue
    else
        subjObs = strcmp({allObservations.Sid}, subjects(s));
        countSubjObs(s,1) = sum(subjObs & completeObsPreproc);
        countSubjObs(s,2) = sum(subjObs & completeObsPrep) - countSubjObs(s,1);
        countSubjObs(s,3) = sum(subjObs) - (countSubjObs(s,1) + countSubjObs(s,2));
        countSubjObs(s,:) = 100*countSubjObs(s,:)./sum(subjObs);
        
        countSubjBase(s,1) = sum(subjObs & completeObsPreprocBase);
        countSubjBase(s,2) = sum(subjObs & completeObsPrepBase) - countSubjBase(s,1);
        countSubjBase(s,3) = sum(subjObs) - (countSubjBase(s,1) + countSubjBase(s,2));
        countSubjBase(s,:) = 100*countSubjBase(s,:)./sum(subjObs);
        
        countSubjObsBase(s,1) = sum(subjObs & completeObsPreprocBase & completeObsPreproc);
        countSubjObsBase(s,2) = sum(subjObs & ~completeObsPreprocBase & completeObsPreproc);
        countSubjObsBase(s,3) = sum(subjObs & completeObsPreprocBase & ~completeObsPreproc);
        countSubjObsBase(s,4) = sum(subjObs) - sum(countSubjObsBase(s,1:3));
        countSubjObsBase(s,:) = 100*countSubjObsBase(s,:)./sum(subjObs);
    end
end
subplot(3,1,1)
bar(countSubjBase,'stacked')
ylabel('% of observations')
ylim([0,105])
xticks(1:n_sbjs)
xticklabels(subjects)
xtickangle(45)
legend({'Clean Baseline', 'Artifacts in Baseline', 'Missing Data'}, 'Location', 'SouthEast')

subplot(3,1,2)
bar(countSubjObs,'stacked')
ylabel('% of observations')
ylim([0,105])
xticks(1:n_sbjs)
xticklabels(subjects)
xtickangle(45)
legend({'Clean EEG', 'Artifacts in EEG', 'Missing Data'}, 'Location', 'SouthEast')

subplot(3,1,3)
bar(countSubjObsBase,'stacked')
ylabel('% of observations')
ylim([0,105])
xticks(1:n_sbjs)
xticklabels(subjects)
xtickangle(45)
legend({'Clean EEG & Baseline', 'Clean EEG only', 'Clean Baseline only', 'Else'}, 'Location', 'SouthEast')

