function events = export_events_EEGAFF2_build1(AllStreams, Event_streams, times, bl, acq)
% Interpret event XDF streams for 4SC experiment to put them in the EEGLAB set
% Inputs:
%   - AllStreams :          streams loaded with load_xdf.
%   - Event_streams :       Indices of the streams to consider for
%                               exporting events.
%   - times :               Complete time points vector of the recording
%                               (to compute latencies)
%
% Outputs:
%   - events :              Structure containing all events loaded.

events = [];
required_fields = {'type', 'latency', 'duration'};
required_types = {'str', 'num', 'num'};
required_inds = [];

for ev=1:length(Event_streams)
    stream = AllStreams{Event_streams(ev)};
    %channels = stream.info.desc.channels.channel;
    
    %if strcmp(stream.info.type,'Markers')
    if strcmp(stream.info.type,'Tag') %PAUL
        required_fields = [required_fields, 'BlockIndex', 'TrialIndex', 'TrialType', 'Phase', 'SoundName', 'Distance', 'Location'];
        required_types = [required_types, 'num', 'num', 'str', 'str', 'str', 'num', 'str'];
        
        %         for col = 1:numel(channels)
        %             if strcmp(channels{ch}.label,'EventName')
        %                 required_inds(1,ev)=ch;
        %             elseif strcmp(channels{ch}.label,'BlockNumber')
        %                 required_inds(strcmp(required_fields,'Block'),ev)=ch;
        %             elseif strcmp(channels{ch}.label,'ConditionName')
        %                 required_inds(strcmp(required_fields,'Condition'),ev)=ch;
        %             elseif strcmp(channels{ch}.label,'PhaseName')
        %                 required_inds(strcmp(required_fields,'Phase'),ev)=ch;
        %             elseif strcmp(channels{ch}.label,'TrialNumber')
        %                 required_inds(strcmp(required_fields,'Trial'),ev)=ch;
        %             elseif strcmp(channels{ch}.label,'Location')
        %                 required_inds(strcmp(required_fields,'EnvPosition'),ev)=ch;
        %             elseif strcmp(channels{ch}.label,'KeyboardInput')
        %                 required_inds(strcmp(required_fields,'KeyboardInput'),ev)=ch;
        %             end
        %         end
    end
end

for ev=1:length(Event_streams)
    stream = AllStreams{Event_streams(ev)};
    
    %if strcmp(stream.info.type,'Markers')
    if strcmp(stream.info.type,'Tag') %PAUL
        markers = SplitEventFields_EEGAFF(stream);
        
        start = 1;
        stop = length(markers.time);
        events_count = length(markers.time);
    end
    
    command = '';
    for f = 1:numel(required_fields)
        if strcmp(required_fields{f}, 'duration')
            command = [command,'''',required_fields{f},''',num2cell(ones(1, events_count)),'];
        elseif strcmp(required_types{f}, 'str')
            command = [command,'''',required_fields{f},''','''','];
        elseif strcmp(required_types{f}, 'num')
            command = [command,'''',required_fields{f},''',[],'];
        end
    end
    % Remove the last ','
    command = command(1:end-1);
    eval(['s_events = struct(',command,');']);
    
    for t=start:stop
        %if strcmp(stream.info.type,'Markers')
        if strcmp(stream.info.type,'Tag') %PAUL
            [~,s_events(t).latency] = min(abs(times - markers.time(t)));
            s_events(t).type = markers.event{t};
            s_events(t).duration = markers.duration(t);
            s_events(t).BlockIndex = markers.block(t);
            s_events(t).TrialIndex = markers.trial(t);
            s_events(t).TrialType = markers.type{t};
            s_events(t).Phase = markers.phase{t};
            s_events(t).SoundName = markers.name{t};
            s_events(t).Distance = markers.distance(t);
            s_events(t).Location = markers.location{t};
        end
    end
    
    events = [events, s_events];
    ToRemove = false(1,length(events));
    
    % Remove events not in the block of interest
    ToRemove = ToRemove | [events.BlockIndex] ~= bl-1;
    
    % Remove multiple block starts
    BlockStart = strcmp({events.type}, 'BlockStart');
    BlockStart(find(BlockStart,1)) = 0;
    ToRemove = ToRemove | BlockStart;
    
    % Remove duplicated events in the baseline trials
    for i=-2:-1
        trialBaseline = [events.TrialIndex] == i & [events.BlockIndex] == bl-1;
        
        if any(trialBaseline)
            BaseOpen = trialBaseline & strcmp({events.SoundName}, 'Open');
            BaseOpen(find(BaseOpen,1)) = 0;
            ToRemove = ToRemove | BaseOpen;
            
            BaseWalk = trialBaseline & strcmp({events.SoundName}, 'Walk');
            BaseWalk(find(BaseWalk,1)) = 0;
            ToRemove = ToRemove | BaseWalk;
            
            BaseNoneWalk = trialBaseline & strcmp({events.SoundName}, 'none') & strcmp({events.Phase}, 'walkingBaseline');
            BaseNoneWalk(find(BaseNoneWalk,1)) = 0;
            ToRemove = ToRemove | BaseNoneWalk;
            
            BaseClose = trialBaseline & strcmp({events.SoundName}, 'Close');
            BaseClose(find(BaseClose,1)) = 0;
            ToRemove = ToRemove | BaseClose;
            
            BaseNoneEnd = trialBaseline & strcmp({events.SoundName}, 'none') & strcmp({events.Phase}, 'ending');
            BaseNoneEnd(find(BaseNoneEnd,1)) = 0;
            ToRemove = ToRemove | BaseNoneEnd;
        end
    end
    
    events = events(~ToRemove);
    
    % Change number of blocks in 1,2,3
    % Change number of trials from 1 to 26
    for t=1:length(events)
        events(t).BlockIndex = events(t).BlockIndex + 1;
        events(t).TrialIndex = events(t).TrialIndex + 3; 
        if events(t).TrialIndex < 3
            events(t).TrialType = 'Baseline';
        end
    end
    
    % Assign the distance measure to all events in the same trial
    ToRemove = false(1,length(events));
    distEvents = strcmp({events.type}, 'DistMeasure');
    for tr=1:26
        trialEvents = [events.TrialIndex] == tr;
        if any(trialEvents)
            % At least one event in this trial
            if any(distEvents & trialEvents)
                % At least one distance measurement for this trial
                if sum(distEvents & trialEvents) == 1
                    % Only one distance measurement
                    distance = events(distEvents & trialEvents).Distance;
                else
                    % More than one distance measurement, use the one
                    % corresponding to the start of the trial (a bit different if the trial is a baseline)
                    if strcmp(events(find(trialEvents,1)).TrialType,'Baseline')
                        openEvents = strcmp({events.SoundName}, 'Open');
                        refLat = events(openEvents & trialEvents).latency;
                    else
                        baseStartEvents = strcmp({events.type}, 'BaseStart');
                        refLat = events(baseStartEvents & trialEvents).latency;
                    end
                    
                    distance = events(distEvents & trialEvents & [events.latency] == refLat).Distance;
                    ToRemove = ToRemove | distEvents & trialEvents & [events.latency] ~= refLat;
                    
                    if isempty(distance)
                        error('Distance not found')
                    end
                end
                
                trialInds = find(trialEvents);
                for i = trialInds
                    events(i).Distance = distance;
                end
            else
                warning('Trial %d does not have a distance measurement', tr)
            end
        end
    end
    
    events = events(~ToRemove);
    
    % Look for missing delimiters
    BlockStart = strcmp({events.type}, 'BlockStart');
    if ~any(BlockStart)
        FirstTrialStart = strcmp({events.SoundName}, 'Open') & [events.TrialIndex] == 1;
        if any(FirstTrialStart)
            startEvent = events(FirstTrialStart);
            startEvent.type = 'BlockStart';
            startEvent.latency = events(1).latency - 1;
            startEvent.SoundName = 'none';
            
            events = [startEvent, events];
        elseif acq == 1
            startEvent = events(1);
            startEvent.type = 'BlockStart';
            startEvent.latency = events(1).latency - 1;
            startEvent.duration = NaN;
            startEvent.SoundName = 'none';
            startEvent.Location = 'none';
            
            events = [startEvent, events];
        else
            warning('Couldn''t find the start of the block in this file')
        end
    end
    
    BlockEnd = strcmp({events.type}, 'BlockEnd');
    if ~any(BlockEnd)
        LastTrialEnd = strcmp({events.type}, 'TrialEnd') & [events.TrialIndex] == 26;
        if any(LastTrialEnd)
            events(end+1) = events(LastTrialEnd);
            events(end).type = 'BlockEnd';
            events(end).latency = events(end-1).latency + 1;
        else
            warning('Couldn''t find the end of the block in this file')
        end
    end
end
end