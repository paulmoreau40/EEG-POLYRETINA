function EEG = events_check(EEG, cfg)

%% Specifics of the experiment
n_Trials = 78;
n_Blocks = 3;
blockLength = floor(n_Trials/n_Blocks);

%% Events
if strcmp(cfg.subjects(cfg.current_subject).id, 'P033')
    % Special case for P033: Block1 Trial7 was recorded twice, removing first part...
    toRemove = [EEG.event.urevent] >= 208 & [EEG.event.urevent] <= 220;
    EEG.event = EEG.event(~toRemove);
end

events = EEG.event;
evts_noBounds = events(~strcmp({events.type},'boundary'));
%missing_base_start = false;
% struct markers with the BaseStart

%% Prepare columns of the dataset
TrialInd = nan(n_Trials,1);
BlockInd = nan(n_Trials,1);
TrialType = cell(n_Trials,1);
% Latencies
Trial_urevent_seq = nan(n_Trials,5);
% Trial_lats = nan(n_Trials,2);
% WalkingBaseline_lats = nan(n_Trials,2);
% BlackBaseline_lats = nan(n_Trials,2);
% Observation_lats = nan(n_Trials,2);
% Question_lats = nan(n_Trials,2);
% Exploration_lats = nan(n_Trials,2);
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
endTrialEvents = strcmp({evts_noBounds.SoundName},'Close');
walkInstructEvents = strcmp({evts_noBounds.SoundName},'Walk');

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
    
    %     %% Fill Latency data
    %     Trial_lats(tr,1) = evts_noBounds(trialStartEvent).latency;
    %     if any(trialEvents & endTrialEvents)
    %         Trial_lats(tr,2) = evts_noBounds(trialEvents & endTrialEvents).latency;
    %     else
    %         warning('Could not find trial end for Trial %d in Block %d.',tr_inBl, bl);
    %     end
    %
    %     if strcmp(TrialType{tr},'Baseline')
    %         WalkingBaseline_lats(tr,:) = Trial_lats(tr,:);
    %     else
    %         BlackBaseline_lats(tr,1) = Trial_lats(tr,1);
    %
    %         if any(trialEvents & obsStartEvents)
    %             Observation_lats(tr,1) = evts_noBounds(trialEvents & obsStartEvents).latency;
    %         else
    %             warning('Could not find obs start for Trial %d in Block %d.', tr_inBl, bl);
    %         end
    %
    %         if any(trialEvents & questionStartEvents)
    %             Question_lats(tr,1) = evts_noBounds(trialEvents & questionStartEvents).latency;
    %         else
    %             warning('Could not find question start for Trial %d in Block %d.', tr_inBl, bl);
    %         end
    %
    %         if any(trialEvents & exploStartEvents)
    %             Exploration_lats(tr,1) = evts_noBounds(trialEvents & exploStartEvents).latency;
    %         else
    %             warning('Could not find explo start for Trial %d in Block %d.', tr_inBl, bl);
    %         end
    %
    %         BlackBaseline_lats(tr,2) = Observation_lats(tr,1);
    %         Observation_lats(tr,2) = Question_lats(tr,1);
    %         Question_lats(tr,2) = Exploration_lats(tr,1);
    %         Exploration_lats(tr,2) = Trial_lats(tr,2);
    %     end
end

% LatenciesSummary = table(BlockInd, TrialInd, TrialType, Trial_lats,...
%     WalkingBaseline_lats, BlackBaseline_lats, Observation_lats, Question_lats, Exploration_lats);

%% Second loop for EEG completeness
subject = cfg.subjects(cfg.current_subject).id;
for tr = 1:n_Trials
    bl = ceil(tr/blockLength);
    tr_inBl = tr - (bl-1)*blockLength;
    
    % Special cases (trials that should be rejected from experimental notes)
    if (strcmp(subject, 'P009') && bl == 1 && tr_inBl == 1) || ...
            (strcmp(subject, 'P010') && bl == 3 && tr_inBl == 2) || ...
            (strcmp(subject, 'P011') && bl == 1 && tr_inBl == 3) || ...
            (strcmp(subject, 'P012') && bl == 1 && tr_inBl == 1) || ...
            (strcmp(subject, 'P015') && bl == 3 && tr_inBl == 1) || ...
            (strcmp(subject, 'P016') && ((bl == 1 && tr_inBl == 3) || (bl == 3 && tr_inBl == 2))) || ...
            (strcmp(subject, 'P018') && bl == 1 && tr_inBl == 3) || ...
            (strcmp(subject, 'P020') && ((bl == 2 && tr_inBl == 3) || (bl == 3 && tr_inBl == 2))) || ...
            (strcmp(subject, 'P021') && bl == 1 && tr_inBl == 3) || ...
            (strcmp(subject, 'P024') && bl == 2 && tr_inBl == 2) || ...
            (strcmp(subject, 'P025') && bl == 1 && tr_inBl == 1) || ...
            (strcmp(subject, 'P027') && bl == 1 && tr_inBl == 1) || ...
            (strcmp(subject, 'P029') && ((bl == 1 && tr_inBl == 3) || (bl == 2 && any(tr_inBl == [1,13:14])))) || ...
            (strcmp(subject, 'P030') && bl == 1 && tr_inBl == 2) || ...
            (strcmp(subject, 'P033') && bl == 1 && any(tr_inBl == [1:3,6:7])) || ...
            (strcmp(subject, 'P034') && bl == 1 && tr_inBl == 1) || ...
            (strcmp(subject, 'P035') && bl == 1 && any(tr_inBl == [3,17]))
        
        Trial_perc(tr) = 0;
    end
    
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
    
    %     cols2inspect = [4];
    %     if strcmp(TrialType{tr}, 'Baseline')
    %         cols2inspect = [cols2inspect, 5];
    %     else
    %         cols2inspect = [cols2inspect, 6:9];
    %     end
    %
    %     for c = cols2inspect
    %         limits = table2array(LatenciesSummary(tr,c));
    %         if any(isnan(limits))
    %             warning('Not enough latency information for trial %d, in %s', tr, LatenciesSummary.Properties.VariableNames{c});
    %             switch c
    %                 case 4
    %                     Trial_perc(tr) = 0;
    %                 case 5
    %                     WalkingBaseline_perc(tr) = 0;
    %                 case 6
    %                     BlackBaseline_perc(tr) = 0;
    %                 case 7
    %                     Observation_perc(tr) = 0;
    %                 case 8
    %                     Question_perc(tr) = 0;
    %                 case 9
    %                     Exploration_perc(tr) = 0;
    %             end
    %         elseif limits(1) >= limits(2)
    %             error('Wrong order of events for trial %d, in %s', tr, LatenciesSummary.Properties.VariableNames{c});
    %         else
    %             switch c
    %                 case 4
    %                     Trial_perc(tr) = 100 - computeMissingEEGData(EEG,limits(1),limits(2));
    %                 case 5
    %                     WalkingBaseline_perc(tr) = 100 - computeMissingEEGData(EEG, limits(1),limits(2));
    %                 case 6
    %                     BlackBaseline_perc(tr) = 100 - computeMissingEEGData(EEG, limits(1),limits(2));
    %                 case 7
    %                     Observation_perc(tr) = 100 - computeMissingEEGData(EEG, limits(1),limits(2));
    %                 case 8
    %                     Question_perc(tr) = 100 - computeMissingEEGData(EEG, limits(1),limits(2));
    %                 case 9
    %                     Exploration_perc(tr) = 100 - computeMissingEEGData(EEG, limits(1),limits(2));
    %             end
    %         end
    %     end
    
    if Trial_perc(tr) == 100
        start_ev = find([EEG.event.urevent] == Trial_urevent_seq(tr,1));
        MaxBufferPre(tr) = (EEG.times(limits(1)) - EEG.times(searchNextNaN(EEG, start_ev, 'backward')))/1000;
        stop_ev = find([EEG.event.urevent] == Trial_urevent_seq(tr,end));
        MaxBufferPost(tr) = (EEG.times(searchNextNaN(EEG, stop_ev, 'forward')) - EEG.times(limits(2)))/1000;
    end
end

EEGCompletenessSummary = table(BlockInd, TrialInd, TrialType, Trial_urevent_seq, Trial_perc, MaxBufferPre, MaxBufferPost,...
    WalkingBaseline_perc, BlackBaseline_perc, Observation_perc, Question_perc, Exploration_perc);

% % Select all the start trials
% start_markers = events;
% %ToKeep = false(1,length(markers));
% WalkBaseStart = strcmp({start_markers.type}, 'WalkingBaseStart');
% TrialStart = strcmp({start_markers.type}, 'BaseStart');
% ToKeep = or(WalkBaseStart,TrialStart);
% start_markers = start_markers(ToKeep);
%
% % Check which start trial is missing
% for bl = 1:n_Blocks
%     i=1;
%     for tr = 1:length(start_markers) % # all trials
%         if (start_markers(i).BlockIndex == bl)
%             if ~(start_markers(i).TrialIndex==tr)
%                 warning(strcat('Trial ', num2str(tr),' in Block ',num2str(bl),' is missing from event start.'));
%                 %tr=tr+1;
%                 %i=i+1;
%             else
%                 i=i+1;
%             end
%         end
%     end
% end
%
% % Check if baseline start is missing
% if sum(WalkBaseStart)<6
%     warning('Missing WalkingBaseStart event(s). Trying to use "Open" sound instead.');
%     missing_base_start = true;
%     start_markers = events;
%     %ToKeep = false(1,length(markers));
%     Open = strcmp({start_markers.SoundName}, 'Open');
%     Baseline = strcmp({start_markers.TrialType}, 'Baseline');
%     WalkBaseStart = and(Open,Baseline);
%     TrialStart = strcmp({start_markers.type}, 'BaseStart');
%     ToKeep = or(WalkBaseStart,TrialStart);
%     start_markers = start_markers(ToKeep);
% end
%
% TrialsInd_block = (1:blockLength)';
% TrialsInd = vertcat(TrialsInd_block,TrialsInd_block,TrialsInd_block);
% BlockInd = (1:n_Blocks)';
% TrialsType = cell(n_Trials,1);
% Trials_lat = nan(n_Trials,2);
% WalkingBaselines_lat = nan(n_Trials,2);
% BlackBaselines_lat = nan(n_Trials,2);
% Observations_lat = nan(n_Trials,2);
% Explorations_lat = nan(n_Trials,2);
%
% % First fill trial types:
% missing_trials = struct;
% missing_trials.Trials = [];
% missing_trials.Block = [];
% for bl = 1:n_Blocks
%     for tr = 1:blockLength
%         TrialEvent = find([start_markers(:).BlockIndex] == BlockInd(bl) & [start_markers(:).TrialIndex] == TrialsInd(tr),1);
%         % Check if trial is missing
%         if isempty(TrialEvent)
%             warning(strcat('Trial ',num2str(tr),' in Block ',num2str(bl), ' is missing from events recording', TrialsInd(tr)));
%             missing_trials.Block = [missing_trials.Block, bl];
%             missing_trials.Trials = [missing_trials.Trials, tr];
%         else
% %             if bl == 1
% %                 TrialsType{tr} = markers(TrialEvent).TrialType;
% %             elseif bl == 2
% %                 TrialsType{tr+26} = markers(TrialEvent).TrialType;
% %             else
% %                 TrialsType{tr+52} = markers(TrialEvent).TrialType;
% %             end
%             TrialsType{(bl-1)*blockLength+tr} = start_markers(TrialEvent).TrialType;
%         end
%     end
% end
%
%
% % just Keep the important phase of the trial
% Impo_Markers = events;
% %ToRemove = false(1,length(Impo_Markers));
%
% % add the WalkingBaseEnd event
% ZoneChange = strcmp({Impo_Markers.type}, 'ZoneChange');
% Ending = strcmp({Impo_Markers.Phase}, 'ending');
% ZoneChangeEnd = and(ZoneChange, Ending);
% Baseline = strcmp({Impo_Markers.TrialType}, 'Baseline');
% WalkingBaseEnd = and(ZoneChangeEnd, Baseline);
% if (sum(WalkingBaseEnd)<6 && ~isempty(missing_trials)) % Not sure about this part, isempty may not be selective enough for this purpose and it duplicates events unnecessarily
%     PlayingSound = strcmp({Impo_Markers.type}, 'PlayingSound');
%     Ending = strcmp({Impo_Markers.SoundName}, 'Close');
%     ZoneChangeEnd = and(PlayingSound, Ending);
%     Baseline = strcmp({Impo_Markers.TrialType}, 'Baseline');
%     WalkingBaseEnd = and(ZoneChangeEnd, Baseline);
% end
%
% WalkingBaseEnd_ur = [Impo_Markers(WalkingBaseEnd).urevent];
% % add to Impo_Markers the WalkingBaseEnd
% for ur=WalkingBaseEnd_ur
%     Impo_Markers(ur).type = 'WalkingBaseEnd';
% end
%
% if missing_base_start
%     WalkingBaseStart_ur = [Impo_Markers(WalkBaseStart).urevent];
%     for ur=WalkingBaseStart_ur
%         Impo_Markers(ur).type = 'WalkingBaseStart';
%     end
% end
%
% % remove type not useful for the latencies and phase analysis
% ZoneChange = strcmp({Impo_Markers.type}, 'ZoneChange');
% Sound = strcmp({Impo_Markers.type}, 'PlayingSound');
% ToRemove = or(ZoneChange,Sound);
% SpatialDiscrExit = strcmp({Impo_Markers.type}, 'SpatialDiscrExit');
% ToRemove = or(ToRemove,SpatialDiscrExit);
% SpatialDiscrEnter = strcmp({Impo_Markers.type}, 'SpatialDiscrEnter');
% ToRemove = or(ToRemove,SpatialDiscrEnter);
% Boundary = strcmp({Impo_Markers.type}, 'boundary');
%
% ToRemove = or(ToRemove,Boundary);
%
% Impo_Markers = Impo_Markers(~ToRemove);
%
% startEvents = find(contains({Impo_Markers(:).type}, 'Start'));
% for ev = 1:length(startEvents)
%     tr = Impo_Markers(startEvents(ev)).TrialIndex;
%     bl = Impo_Markers(startEvents(ev)).BlockIndex;
%
% %     if bl == 2
% %         tr = tr+26;
% %     elseif bl == 3
% %         tr = tr+52;
% %     end
% tr = (bl-1)*blockLength+tr;
%
%     switch Impo_Markers(startEvents(ev)).type
%         case 'WalkingBaseStart'
%             WalkingBaselines_lat(tr,1) = Impo_Markers(startEvents(ev)).latency;
%         case 'BaseStart'
%             BlackBaselines_lat(tr,1) = Impo_Markers(startEvents(ev)).latency;
%         case 'TrialStart'
%             Trials_lat(tr,1) = Impo_Markers(startEvents(ev)).latency;
%         case 'ObsStart'
%             Observations_lat(tr,1) = Impo_Markers(startEvents(ev)).latency;
%         case 'ExploStart'
%             Explorations_lat(tr,1) = Impo_Markers(startEvents(ev)).latency;
%         case 'QuestionStart'
%             Question_lat(tr,1) = Impo_Markers(startEvents(ev)).latency;
%     end
% end
%
%
% endEvents = find(contains({Impo_Markers(:).type}, 'End'));
% for ev = 1:length(endEvents)
%     tr = Impo_Markers(endEvents(ev)).TrialIndex;
%     bl = Impo_Markers(endEvents(ev)).BlockIndex;
%
% %     if bl == 2
% %         tr = tr+26;
% %     elseif bl == 3
% %         tr = tr+52;
% %     end
%     tr = (bl-1)*blockLength+tr;
%
%     switch Impo_Markers(endEvents(ev)).type
%         case 'WalkingBaseEnd'
%             WalkingBaselines_lat(tr,2) = Impo_Markers(endEvents(ev)).latency;
%         case 'BaseEnd'
%             BlackBaselines_lat(tr,2) = Impo_Markers(endEvents(ev)).latency;
%         case 'TrialEnd'
%             Trials_lat(tr,2) = Impo_Markers(endEvents(ev)).latency;
%         case 'ObsEnd'
%             Observations_lat(tr,2) = Impo_Markers(endEvents(ev)).latency;
%         case 'ExploEnd'
%             Explorations_lat(tr,2) = Impo_Markers(endEvents(ev)).latency;
%         case 'QuestionEnd'
%             Question_lat(tr,2) = Impo_Markers(endEvents(ev)).latency;
%     end
% end
%
% LatenciesSummary = table(TrialsInd, TrialsType, Trials_lat,...
%     WalkingBaselines_lat, BlackBaselines_lat, Observations_lat, Explorations_lat, Question_lat);
%
%
% %% Check trials and their completeness
% Trials_perc = nan(n_Trials,1);
% WalkingBaselines_perc = nan(n_Trials,1);
% BlackBaselines_perc = nan(n_Trials,1);
% Observations_perc = nan(n_Trials,1);
% Explorations_perc = nan(n_Trials,1);
% Question_perc = nan(n_Trials,1);
%
% % for the trials presents check if one of the 2 delimiters in the latency is missing
% for tr = setdiff(1:n_Trials, missing_trials.Trials)
%     % the 3 columns in the summary
%     cols2inspect = [];%3;
%     if strcmp(LatenciesSummary.TrialsType{tr}, 'Baseline')
%         cols2inspect = [cols2inspect, 4];
%     else
%         cols2inspect = [cols2inspect, [3,5,6,7]];
%     end
%
%     if sum(isnan(table2array(LatenciesSummary(tr,cols2inspect))))>0
%         warning('Missing delimiter in trial %d', tr);
%     end
%
%     for c = cols2inspect
%         limits = table2array(LatenciesSummary(tr,c));
%         if limits(1) >= limits(2)
%             warning('Wrong order of events for trial %d, in %s', tr-1, LatenciesSummary.Properties.VariableNames{c});
%         else
%
%             switch c
%                 case 3
%                     Trials_perc(tr) = 100 - computeMissingEEGData(EEG,limits(1),limits(2));
%                 case 4
%                     WalkingBaselines_perc(tr) = 100 - computeMissingEEGData(EEG, limits(1),limits(2));
%                 case 5
%                     BlackBaselines_perc(tr) = 100 - computeMissingEEGData(EEG, limits(1),limits(2));
%                 case 6
%                     Observations_perc(tr) = 100 - computeMissingEEGData(EEG, limits(1),limits(2));
%                 case 7
%                     Explorations_perc(tr) = 100 - computeMissingEEGData(EEG, limits(1),limits(2));
%             end
%         end
%     end
% end
%
% EEGCompletenessSummary = table(TrialsInd, TrialsType, Trials_perc,...
%     WalkingBaselines_perc, BlackBaselines_perc, Observations_perc, Explorations_perc);

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

% figure
% for b = 1:n_Blocks
%     subplot(n_Blocks, 1, b)
%     hold on
%     range = (1+blockLength*(b-1)):(blockLength*b);
%     bar(table2array(EEGCompletenessSummary(range,2)),...
%         table2array(EEGCompletenessSummary(range,4)))
%     bar(table2array(EEGCompletenessSummary(range(1),2))-1,...
%         table2array(EEGCompletenessSummary(range(1),5)))
%     title(sprintf('Block %d', b))
%     ylim([0,100])
%     xlabel('Trial')
%     ylabel('% Complete')
%     if b == 1
%         legend({'Trials', 'WalkingBaselines'}, 'Position', [0.825,0.925,0.1,0.05])
%     end
% end
% subject = cfg.subjects(cfg.current_subject).id;
% suptitle(subject);
% saveCurrentFig([cfg.figures_folder 'MissingData' filesep],...
%     [subject, '_TrialsInspectionResults'], {'png'}, [1400,700]);
%
% figure
% for b = 1:n_Blocks
%     subplot(n_Blocks, 1, b)
%     range = (1+blockLength*(b-1)):(blockLength*b);
%     bar(table2array(EEGCompletenessSummary(range,1)),...
%         table2array(EEGCompletenessSummary(range,5:7)))
%     xticks(table2array(EEGCompletenessSummary(range,1)));
%     xticklabels(table2cell(EEGCompletenessSummary(range,2)))
%     xtickangle(45)
%     title(sprintf('Block %d', b))
%     ylim([0,100])
%     ylabel('% Complete')
%     if b == 1
%         legend({'BlackBaselines', 'Observations', 'Explorations'}, 'Position', [0.825,0.93,0.1,0.05])
%     end
% end
% subject = cfg.subjects(cfg.current_subject).id;
% suptitle(subject);
% saveCurrentFig([cfg.figures_folder 'MissingData' filesep],...
%     [subject, '_TrialsInspectionResults_Details'], {'png'}, [1400,700]);

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