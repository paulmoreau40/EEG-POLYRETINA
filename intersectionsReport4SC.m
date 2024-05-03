function intersectionsReport4SC(EEG, cfg)
% Should be called after fillNaNs and events_check but before resampling or selection

events = EEG.event;
markers = events([events(:).duration] == 0 & ~(strcmp({events(:).type}, 'Blink') |...
    strcmp({events(:).type}, 'Saccade') | strcmp({events(:).type}, 'Fixation')));

stops = find(strcmp({markers(:).type}, 'StopMoving'));

Intersections = struct('urevent', [], 'Block', [], 'Condition', [], 'Phase', [], 'Trial', [], 'Name', [],...
    'Order', [], 'Direction', [], 'Duration_ms', [], 'EEGCompleteness',[],...
    'StaticDuration_ms', [], 'Decisions', [], 'ureventSequence', []);

i = 1;
for s = 1:length(stops)-1
    cnd = markers(stops(s)).Condition;
    trl = markers(stops(s)).Trial;
    
    markersInCnd = strcmp({markers(:).Condition}, cnd);
    markersInTrl= [markers(:).Trial] == trl;
    
    endTrl = find(strcmp({markers(:).type}, 'TrialEnd') & markersInCnd & markersInTrl);
    ending = findMarker(markers, stops(s), 'StopMoving', 'ahead', endTrl);
    
    if ~isempty(ending)
        urev = markers(stops(s)).urevent;
        blk = markers(stops(s)).Block;
        phs = markers(stops(s)).Phase;
        name = markers(stops(s)).EnvPosition;
        
        % Keep the following sequence of events:
        % 1: Beginning of the hallway
        % 2: Stop for observation at the intersection
        % 3: (Optional) First subject decision (not necessarily good)
        % 4: Movement Initiation
        % 5: (Optional) Start of the turn movement
        % 6: (Optional) Stop of the turn movement
        % 7: Next Stop (Goal reach or next observation)
        seq = NaN(1,7);
        seq(2) = urev;
        seq(7) = markers(ending).urevent;
        
        % Get the order
        if i>1 && Intersections(i-1).Trial == trl
            o = o+1;
        else
            o = 1;
        end
        
        % Get the turn direction
        turn = determineTurn(name, markers(ending).EnvPosition);
        if ~strcmp(turn, 'Front')
            startTurn = findMarker(markers, stops(s), 'StartRotating', 'ahead', ending);
            seq(5) = markers(startTurn).urevent;
            stopTurn = findMarker(markers, startTurn, 'StopRotating', 'ahead', ending);
            seq(6) = markers(stopTurn).urevent;
        end
        
        % Get the duration
        startTrl = find(strcmp({markers(:).type}, 'TrialStart') & markersInCnd & markersInTrl);
        if o == 1
            beginning = findMarker(markers, stops(s), 'StartMoving', 'before', startTrl);
        else
            beginning = findMarker(markers, stops(s), 'StopRotating', 'before', startTrl);
        end
        seq(1) = markers(beginning).urevent;
        dur = 1 + markers(ending).latency - markers(beginning).latency;
        
        % Find completeness
        cmp = 100*(1 - (sum(sum(isnan(EEG.data(:,markers(beginning).latency:markers(ending).latency)),1)>0)/dur));
        
        % Get the static duration
        movInit = findMarker(markers, stops(s), 'StartMoving', 'ahead', ending);
        seq(4) = markers(movInit).urevent;
        statDur = 1 + markers(movInit).latency - markers(stops(s)).latency;
        
        % Get the decisions
        if strcmp(phs,'TestPhase')
            decision = findMarker(markers, stops(s), 'SubjectDecision', 'ahead', movInit);
            seq(3) = markers(decision).urevent;
            decs = '';
            while ~isempty(decision)
                decs = [decs, markers(decision).KeyboardInput, '-'];
                decision = findMarker(markers, decision, 'SubjectDecision', 'ahead', movInit);
            end
            decs = decs(1:end-1);
        else
            decs = 'NA';
        end       
        
        Intersections(i) = struct('urevent', urev,...
            'Block', blk, 'Condition', cnd, 'Phase', phs, 'Trial', trl,...
            'Name', name,...
            'Order', o,...
            'Direction', turn,...
            'Duration_ms', dur,...
            'EEGCompleteness', cmp,...
            'StaticDuration_ms', statDur,...
            'Decisions',decs,...
            'ureventSequence',seq);
        
        i = i+1;
    end
end

%% Save the output structure
subject = cfg.subjects(cfg.current_subject).id;
dirname = [cfg.study_folder cfg.preprocessing_folder];
fname = [subject '_IntersectionsReport.mat'];
save(fullfile(dirname, fname), 'Intersections')

%% Helper functions
    function ind = findMarker(evts, refInd, eName, direction, limitInd)
        % refInd: reference index in the events structure
        % eName: event name to find
        % direction: direction of search in the event structure: 'ahead' or 'before'
        % limitInd: do not search past this index
        
        switch direction
            case 'ahead'
                next = find(strcmp({evts(refInd+1:limitInd).type},eName),1);
                if ~isempty(next)
                    ind = refInd + next;
                else
                    ind = [];
                end
            case 'before'
                prev = find(strcmp({evts(limitInd:refInd-1).type},eName),1, 'last');
                if ~isempty(prev)
                    ind = limitInd + prev - 1;
                else
                    ind = [];
                end
        end
    end

    function turn = determineTurn(pointA, pointB)
        % pointA: intersection observation point - Stop1 or Stop2
        % pointB: intersection arrival point - Goal1 or Goal2 or Stop
        
        if contains(pointA, 'Stop1')
            if contains(pointB, 'Stop')
                turn = 'Right';
            elseif contains(pointB, 'Goal1')
                turn = 'Front';
            elseif contains(pointB, 'Goal2')
                turn = 'Left';
            else
                error('Unknown case')
            end
            
        elseif contains(pointA, 'Stop2')
            if contains(pointB, 'Stop')
                turn = 'Left';
            elseif contains(pointB, 'Goal1')
                turn = 'Right';
            elseif contains(pointB, 'Goal2')
                turn = 'Front';
            else
                error('Unknown case')
            end
            
        else
            error('Unknown case')
        end
    end
end