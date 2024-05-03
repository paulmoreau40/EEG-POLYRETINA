function EEG = events_check_EEGPOL(EEG, cfg)

%% Specifics of the experiment
n_Trials = 210;
n_Blocks = 70;
blockLength = floor(n_Trials/n_Blocks);

%% Events

events = EEG.event;
evts_noBounds = events(~strcmp({events.type},'boundary'));
%missing_base_start = false;

%% Prepare columns of the dataset
TrialInd = nan(n_Trials,1);
BlockInd = nan(n_Trials,1);
TrialType = cell(n_Trials,1);

Trial_urevent_seq = nan(n_Trials,5);

% EEG Completeness
Trial_perc = nan(n_Trials,1);
MaxBufferPre = nan(n_Trials,1);
MaxBufferPost = nan(n_Trials,1);
WalkingBaseline_perc = nan(n_Trials,1);
BlackBaseline_perc = nan(n_Trials,1);
Observation_perc = nan(n_Trials,1);
Question_perc = nan(n_Trials,1);
Exploration_perc = nan(n_Trials,1);

%% Create useful variables
startTrialEvents = contains({evts_noBounds.type},'BaseStart');
obsStartEvents = strcmp({evts_noBounds.type},'ObsStart');
questionStartEvents = strcmp({evts_noBounds.type},'QuestionStart');
exploStartEvents = strcmp({evts_noBounds.type},'ExploStart');
%endTrialEvents = strcmp({evts_noBounds.SoundName},'Close');
% walkInstructEvents = strcmp({evts_noBounds.SoundName},'Walk');

%% Loop over trials
for tr = 1:n_Trials
    % useful variables
    bl = ceil(tr/blockLength);
    tr_inBl = tr - (bl-1)*blockLength;
    trialEvents = [evts_noBounds.BlockIndex] == bl & [evts_noBounds.TrialIndex] == tr_inBl;
    
    % Fill columns
    TrialInd(tr) = tr_inBl;
    BlockInd(tr) = bl;
    
    % Check whether we are missing the trial
    if ~any(trialEvents & startTrialEvents)
        warning('Could not find trial start for Trial %d in Block %d.',tr_inBl, bl);
        if any(trialEvents)
            % Some events exist for this trial
            disp('Choosing a replacement for the start event...');
            if tr_inBl < 3 && any(trialEvents & walkInstructEvents)
                disp('Using walking instruction');
                trialStartEvent = trialEvents & walkInstructEvents;
            else
                % To do
                disp('...');
            end
        else
            disp('This trial was not recorded...');
            if tr_inBl < 3
                TrialType{tr} = 'Baseline';
            else
                error('Cannot infer trial type...')
            end
            continue
        end
    else
        trialStartEvent = trialEvents & startTrialEvents;
    end
    TrialType{tr} = evts_noBounds(trialStartEvent).TrialType;
    
    %% Fill urevent sequence
    Trial_urevent_seq(tr,1) = evts_noBounds(trialStartEvent).urevent;
    
    if ~strcmp(TrialType{tr}, 'Baseline')
        if any(trialEvents & obsStartEvents)
            Trial_urevent_seq(tr,2) = evts_noBounds(trialEvents & obsStartEvents).urevent;
        else
            warning('Could not find obs start for Trial %d in Block %d.', tr_inBl, bl);
        end
        
        if any(trialEvents & questionStartEvents)
            Trial_urevent_seq(tr,3) = evts_noBounds(trialEvents & questionStartEvents).urevent;
        else
            warning('Could not find question start for Trial %d in Block %d.', tr_inBl, bl);
        end
        
        if any(trialEvents & exploStartEvents)
            Trial_urevent_seq(tr,4) = evts_noBounds(trialEvents & exploStartEvents).urevent;
        else
            warning('Could not find explo start for Trial %d in Block %d.', tr_inBl, bl);
        end
    end
    
    if any(trialEvents & endTrialEvents)
        Trial_urevent_seq(tr,end) = evts_noBounds(trialEvents & endTrialEvents).urevent;
    else
        warning('Could not find trial end for Trial %d in Block %d.', tr_inBl, bl);
    end 
end


%% Second loop for EEG completeness
subject = cfg.subjects(cfg.current_subject).id;
for tr = 1:n_Trials
    bl = ceil(tr/blockLength);
    tr_inBl = tr - (bl-1)*blockLength;
    
    % % Special cases (trials that should be rejected from experimental notes)
    % if (strcmp(subject, 'P009') && bl == 1 && tr_inBl == 1) 
    %     Trial_perc(tr) = 0;
    % end
    
    % Full trial first
    limits = [evts_noBounds([evts_noBounds.urevent] == Trial_urevent_seq(tr,1)).latency,...
        evts_noBounds([evts_noBounds.urevent] == Trial_urevent_seq(tr,end)).latency];
    
    if isempty(limits)
        warning('Not enough latency information for trial %d', tr);
        Trial_perc(tr) = 0;
    elseif limits(1) >= limits(2)
        error('Wrong order of events for trial %d', tr);
    else
        if isnan(Trial_perc(tr))
            Trial_perc(tr) = 100 - computeMissingEEGData(EEG,limits(1),limits(2));
            % Otherwise this has already been set by a special case
        end
        
        if strcmp(TrialType{tr}, 'Baseline')
            WalkingBaseline_perc(tr) = Trial_perc(tr);
        else
            if Trial_perc(tr) == 100
                BlackBaseline_perc(tr) = 100;
                Observation_perc(tr) = 100;
                Question_perc(tr) = 100;
                Exploration_perc(tr) = 100;
            elseif Trial_perc(tr) > 0
                for j = 1:size(Trial_urevent_seq,2)-1
                    limits = [evts_noBounds([evts_noBounds.urevent] == Trial_urevent_seq(tr,j)).latency,...
                        evts_noBounds([evts_noBounds.urevent] == Trial_urevent_seq(tr,j+1)).latency];
                    
                    switch j
                        case 1
                            if isempty(limits)
                                warning('Not enough latency information for trial %d, in Black Baseline', tr);
                                BlackBaseline_perc(tr) = 0;
                            elseif limits(1) >= limits(2)
                                error('Wrong order of events for trial %d, in Black Baseline', tr);
                            else
                                BlackBaseline_perc(tr) = 100 - computeMissingEEGData(EEG, limits(1),limits(2));
                            end
                        case 2
                            if isempty(limits)
                                warning('Not enough latency information for trial %d, in Observation', tr);
                                Observation_perc(tr) = 0;
                            elseif limits(1) >= limits(2)
                                error('Wrong order of events for trial %d, in Observation', tr);
                            else
                                Observation_perc(tr) = 100 - computeMissingEEGData(EEG, limits(1),limits(2));
                            end
                        case 3
                            if isempty(limits)
                                warning('Not enough latency information for trial %d, in Question', tr);
                                Question_perc(tr) = 0;
                            elseif limits(1) >= limits(2)
                                error('Wrong order of events for trial %d, in Question', tr);
                            else
                                Question_perc(tr) = 100 - computeMissingEEGData(EEG, limits(1),limits(2));
                            end
                        case 4
                            if isempty(limits)
                                warning('Not enough latency information for trial %d, in Exploration', tr);
                                Exploration_perc(tr) = 0;
                            elseif limits(1) >= limits(2)
                                error('Wrong order of events for trial %d, in Exploration', tr);
                            else
                                Exploration_perc(tr) = 100 - computeMissingEEGData(EEG, limits(1),limits(2));
                            end
                    end
                end
            else
                BlackBaseline_perc(tr) = 0;
                Observation_perc(tr) = 0;
                Question_perc(tr) = 0;
                Exploration_perc(tr) = 0;
            end
        end
    end
   
    
    if Trial_perc(tr) == 100
        start_ev = find([EEG.event.urevent] == Trial_urevent_seq(tr,1));
        MaxBufferPre(tr) = (EEG.times(limits(1)) - EEG.times(searchNextNaN(EEG, start_ev, 'backward')))/1000;
        stop_ev = find([EEG.event.urevent] == Trial_urevent_seq(tr,end));
        MaxBufferPost(tr) = (EEG.times(searchNextNaN(EEG, stop_ev, 'forward')) - EEG.times(limits(2)))/1000;
    end
end

EEGCompletenessSummary = table(BlockInd, TrialInd, TrialType, Trial_urevent_seq, Trial_perc, MaxBufferPre, MaxBufferPost,...
    WalkingBaseline_perc, BlackBaseline_perc, Observation_perc, Question_perc, Exploration_perc);

%% Plot Summary info:
figure
for b = 1:n_Blocks
    subplot(n_Blocks, 1, b)
    hold on
    walkingBase_bl = BlockInd == b & strcmp(TrialType,'Baseline');
    trials_bl = BlockInd == b & ~strcmp(TrialType,'Baseline');
    
    bar(TrialInd(trials_bl), Trial_perc(trials_bl));
    bar(TrialInd(walkingBase_bl), Trial_perc(walkingBase_bl));
    title(sprintf('Block %d', b))
    ylim([0,100])
    xlabel('Trial')
    ylabel('% Complete')
    if b == 1
        legend({'Trials', 'WalkingBaselines'}, 'Position', [0.825,0.925,0.1,0.05])
    end
end
subject = cfg.subjects(cfg.current_subject).id;
suptitle(subject);
saveCurrentFig([cfg.figures_folder 'MissingData' filesep],...
    [subject, '_TrialsInspectionResults'], {'png'}, [1400,700]);

figure
for b = 1:n_Blocks
    subplot(n_Blocks, 1, b)
    trials_bl = BlockInd == b & ~strcmp(TrialType,'Baseline');
    
    bar(TrialInd(trials_bl), table2array(EEGCompletenessSummary(trials_bl,9:12)))
    xticks(TrialInd(trials_bl));
    xticklabels(TrialType(trials_bl))
    xtickangle(45)
    title(sprintf('Block %d', b))
    ylim([0,100])
    ylabel('% Complete')
    if b == 1
        legend({'BlackBaselines', 'Observations', 'Questions', 'Explorations'}, 'Position', [0.825,0.93,0.1,0.05])
    end
end
subject = cfg.subjects(cfg.current_subject).id;
suptitle(subject);
saveCurrentFig([cfg.figures_folder 'MissingData' filesep],...
    [subject, '_TrialsInspectionResults_Details'], {'png'}, [1400,700]);

EEG.etc.TrialsInspection = EEGCompletenessSummary;

%% Helper Functions ------------------------------------------
    function  perc = computeMissingEEGData(EEG,lat_start,lat_stop)
        pnts = lat_start:lat_stop;
        data2inspect = EEG.data(:,pnts);
        NanInspection = isnan(data2inspect);
        NanInspection = sum(NanInspection,1)>0;
        perc = 100*sum(NanInspection)/length(NanInspection);
    end

    function lat = searchNextNaN(EEG, ref_ind, direction)
        switch direction
            case 'backward'
                data2inspect = EEG.data(:,1:(EEG.event(ref_ind).latency-1));
            case 'forward'
                data2inspect = EEG.data(:,(EEG.event(ref_ind).latency+1):end);
        end
        
        NanInspection = isnan(data2inspect);
        NanInspection = sum(NanInspection,1)>0;
        
        switch direction
            case 'backward'
                next_bound = find(strcmp({EEG.event(1:(ref_ind-1)).type}, 'boundary'),1,'last');
                if isempty(next_bound)
                    lat = find(NanInspection, 1,'last');
                    if isempty(lat)
                        % No NaNs found before the event latency, probably
                        % the start of the recording (with no boundary
                        % event starting the recording)
                        lat = 1;
                    end
                else
                    lat = max(EEG.event(next_bound).latency, find(NanInspection,1,'last'));
                end
            case 'forward'
                next_bound = find(strcmp({EEG.event((ref_ind+1):end).type}, 'boundary'),1);
                if isempty(next_bound)
                    lat = EEG.event(ref_ind).latency + find(NanInspection,1);
                    if isempty(lat)
                        % No NaNs found after the event latency, probably
                        % the end of the recording (with no boundary
                        % event ending the recording)
                        lat = EEG.pnts;
                    end
                else
                    lat = min(EEG.event(ref_ind+next_bound).latency, EEG.event(ref_ind).latency + find(NanInspection,1));
                end
        end
        % Make sure the output is an integer
        lat = round(lat);
    end
end