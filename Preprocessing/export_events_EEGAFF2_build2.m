function events = export_events_EEGAFF2_build2(AllStreams, Event_streams, times)
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
    
    if strcmp(stream.info.type,'Markers')        
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
    
    if strcmp(stream.info.type,'Markers')
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
        if strcmp(stream.info.type,'Markers')
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
end
end