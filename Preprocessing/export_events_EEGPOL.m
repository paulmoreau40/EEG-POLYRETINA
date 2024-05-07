function events = export_events_EEGPOL(AllStreams, Event_streams, times, block_ind, acq_ind)
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

% PAUL : only one stream for the Events, so NO LOOP -> CHANGE LATER IF
% NEEDED
stream = AllStreams{Event_streams};
%markers = SplitEventFields_EEGPOL(stream, times);


markers.events = stream.time_series';
markers.time = stream.time_stamps';
markers.index = find(markers.time); % simpler way to get the indices, but will be updated after event deletion ! ...
events_count = length(markers.events);

required_fields = {'Start fam_empty', 'End fam_empty', 'Start fam_symbol', 'End fam_symbol', 'Start baseline', 'Start training', 'End training', ...
    'Start c_baseline', 'End c_baseline', 'Start trial', 'End trial'}';
new_fields = {'FamEmptyStart', 'FamEmptyEnd', 'FamSymbolStart', 'FamSymbolEnd', 'BaselineStart', 'TrainStart', 'TrainEnd', ...
    'BaselineCoarseStart', 'BaselineCoarseEnd', 'TrialStart', 'TrialEnd'}';
dic = dictionary(required_fields, new_fields);

tbl = table(markers.index, markers.time, markers.events(:,1), repmat({NaN}, events_count,1), repmat({NaN}, events_count,1), ...
    cellfun(@str2double, markers.events(:,2)), cellfun(@str2double, markers.events(:,3)), NaN(events_count,1), zeros(events_count,1), ...
    zeros(events_count,1), 'VariableNames', {'Index', 'Time', 'RawName', 'Name', 'TrialType', 'Value', 'Duration','Latency','BlockIndex','TrialIndex'});

startTrialIndexing = 0;
iBlock = 1;
iTrial = 1;

for i=1:events_count
    row = tbl(i,:);

    if isKey(dic,row.RawName) % Rename the raw name with cleaner name (see new fields above)
        row.Name = dic(row.RawName);
    else
        error('Unkown name of event ?')
    end
    
    % Trial type
    if ismember(row.Name, {'FamEmptyStart','FamSymbolStart','TrainStart','TrialStart'})
        trialType = {"Angle"+num2str(row.Value)};
        row.TrialType = trialType;
    elseif ismember(row.Name, {'FamEmptyEnd','FamSymbolEnd','TrainEnd','TrialEnd'})
        % do nothing because is the same as the previous trial (start)
         row.TrialType = trialType;
    elseif ismember(row.Name, {'BaselineStart','BaselineCoarseStart','BaselineCoarseEnd'})
         row.TrialType = {"Baseline"};
    else
        error("Event not recognised")
    end

    % COUNT TRIALS AND BLOCKS -> BlockIndex, TrialIndex
    if startTrialIndexing && ~ismember(row.Name,'BaselineCoarseStart')
        row.BlockIndex = iBlock;
        row.TrialIndex = iTrial;
        
        % New trial
        if ismember(row.Name,{'BaselineStart','TrialEnd'})
            iTrial = iTrial + 1;
        end

        % New block
        if iTrial > 3
            iTrial = 1;
            iBlock = iBlock + 1;
        end

    else % if still indexing but reaches the last Coarse Baseline
        startTrialIndexing = 0;
    end


    
    % Boolean to start trial and block indexing (starting after first
    % Coarse Baseline)
    % 2nd condition to prevent restarting indexing again at the end
    if ismember(row.Name,'BaselineCoarseEnd') && iTrial*iBlock == 1
        startTrialIndexing = 1;
    end
        
    [~,row.Latency] = min(abs(times - row.Time));
    tbl(i,:) = row;
end

tbl.DurationComputed = [diff(tbl.Time); 0];
tbl.TimeZeroed = tbl.Time -  tbl.Time(1);



% COUNT OCCURENCE OF TRIALS WITH 20° AND 45° ANGLE AND BASELINE
TrialBaselineCount = sum(tbl(tbl.BlockIndex ~= 0, :).Value == -1);
Trial20Count = sum(tbl(tbl.BlockIndex ~= 0, :).Value == 20);
Trial45Count = sum(tbl(tbl.BlockIndex ~= 0, :).Value == 45);

% Display the count of occurrences
disp("There are " + num2str(iBlock-1) + " blocks")
disp("There are " + TrialBaselineCount + " baselines within the trials")
disp("There are " + Trial20Count + " trials of a 20-angle degree")
disp("There are " + Trial45Count + " trials of a 45-angle degree")


% STORE IN THE RIGHT STRUCTURE AND WAY

events = struct('type', {}, 'latency', {}, 'duration', {});

for t=1:events_count
    events(t).type = tbl.Name{t};
    events(t).duration = tbl.DurationComputed(t);
    events(t).latency = tbl.Latency(t);
    events(t).TrialType = tbl.TrialType{t};
    events(t).BlockIndex = tbl.BlockIndex(t);
    events(t).TrialIndex = tbl.TrialIndex(t);
end

return

  % 
  % for t=start:stop
  %       %if strcmp(stream.info.type,'Markers')
  %       if strcmp(stream.info.type,'Tag') %PAUL
  %           [~,s_events(t).latency] = min(abs(times - markers.time(t)));
  %           s_events(t).type = markers.event{t};
  %           s_events(t).duration = markers.duration(t);
  %           s_events(t).BlockIndex = markers.block(t);
  %           s_events(t).TrialIndex = markers.trial(t);
  %           s_events(t).TrialType = markers.type{t};
  %           s_events(t).Phase = markers.phase{t};
  %           s_events(t).SoundName = markers.name{t};
  %           s_events(t).Distance = markers.distance(t);
  %           s_events(t).Location = markers.location{t};
  %       end
  %   end
    






% 
% 
% 
% 
% 
% for ev=1:length(Event_streams)
%     stream = AllStreams{Event_streams(ev)};
% 
%     if strcmp(stream.info.type,'Markers')
%         markers = SplitEventFields_EEGAFF(stream);
% 
%         start = 1;
%         stop = length(markers.time);
%         events_count = length(markers.time);
%     end
% 
%     command = '';
%     for f = 1:numel(required_fields)
%         if strcmp(required_fields{f}, 'duration')
%             command = [command,'''',required_fields{f},''',num2cell(ones(1, events_count)),'];
%         elseif strcmp(required_types{f}, 'str')
%             command = [command,'''',required_fields{f},''','''','];
%         elseif strcmp(required_types{f}, 'num')
%             command = [command,'''',required_fields{f},''',[],'];
%         end
%     end
%     % Remove the last ','
%     command = command(1:end-1);
%     eval(['s_events = struct(',command,');']);
% 
%     for t=start:stop
%         if strcmp(stream.info.type,'Markers')
%             [~,s_events(t).latency] = min(abs(times - markers.time(t)));
%             s_events(t).type = markers.event{t};
%             s_events(t).duration = markers.duration(t);
%             s_events(t).BlockIndex = markers.block(t);
%             s_events(t).TrialIndex = markers.trial(t);
%             s_events(t).TrialType = markers.type{t};
%             s_events(t).Phase = markers.phase{t};
%             s_events(t).SoundName = markers.name{t};
%             s_events(t).Distance = markers.distance(t);
%             s_events(t).Location = markers.location{t};
%         end
%     end
%     events = [events, s_events];
% end
% end