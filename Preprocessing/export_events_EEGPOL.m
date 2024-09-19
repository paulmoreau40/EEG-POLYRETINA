function events = export_events_EEGPOL(AllStreams, Event_streams, times, block_ind, acq_ind, cfg)
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

if ~isscalar(Event_streams)
    error('More than 1 event streams ?')
end

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

% tbl = table(markers.index, markers.time, markers.events(:,1), repmat({NaN}, events_count,1), repmat({NaN}, events_count,1), ...
%     cellfun(@str2double, markers.events(:,2)), cellfun(@str2double, markers.events(:,3)), NaN(events_count,1), zeros(events_count,1), ...
%     zeros(events_count,1), 'VariableNames', {'Index', 'Time', 'RawName', 'Name', 'TrialType', 'Value', 'Duration','Latency','BlockIndex','TrialIndex'});

% % Try removing duration field
% tbl = table(markers.index, markers.time, markers.events(:,1), repmat({NaN}, events_count,1), repmat({NaN}, events_count,1), ...
%     cellfun(@str2double, markers.events(:,2)), NaN(events_count,1), zeros(events_count,1), ...
%     zeros(events_count,1), 'VariableNames', {'Index', 'Time', 'RawName', 'Name', 'TrialType', 'Value','Latency','BlockIndex','TrialIndex'}); %'Duration'

% Try removing Index
tbl = table(markers.time, markers.events(:,1), repmat({NaN}, events_count,1), repmat({NaN}, events_count,1), ...
    cellfun(@str2double, markers.events(:,2)), cellfun(@str2double, markers.events(:,3)), NaN(events_count,1), zeros(events_count,1), ...
    zeros(events_count,1), 'VariableNames', {'Time', 'RawName', 'Name', 'TrialType', 'Value', 'Duration','Latency','BlockIndex','TrialIndex'});


startTrialIndexing = 0;
iBlock = 1;
iTrial = 1;



% handle subject P006 missing lot of events at beginning
subject = cfg.subjects(cfg.current_subject).id;

% if strcmp(subject, 'P005')
%     trialType = {"Angle20"};
% end
% 
% % TEST TO ONLY KEEP ESSAI 1 of P009, bug après ?
if strcmp(subject, 'P009')
    events_count = 291;
end



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



% ADDING A BaselineEnd EVENT AFTER EACH BASELINE

baseline_start_indices = find(strcmp(tbl.Name, 'BaselineStart'));
new_rows = table(); % Initialize an empty table to hold new rows
for i = 1:length(baseline_start_indices)
    current_index = baseline_start_indices(i);
    if current_index < height(tbl)
        next_index = current_index + 1;
        
        % Create the new row
        new_row = tbl(current_index, :);
        new_row.Name = {'BaselineEnd'};
        new_row.RawName = {'End baseline'};
        new_row.Latency = tbl.Latency(next_index) - 1;
        new_row.Time = tbl.Time(next_index) - 0.01;
        %new_row.BlockIndex = tbl.BlockIndex(current_index);
        %new_row.TrialIndex = tbl.TrialIndex(current_index);
        
        new_rows = [new_rows; new_row];
    end
end
for i = 1:height(new_rows)
    insert_index = baseline_start_indices(i) + i - 1; % Adjust the index to account for inserted rows
    tbl = [tbl(1:insert_index, :); new_rows(i, :); tbl(insert_index + 1:end, :)];
end


tbl.DurationComputed = [diff(tbl.Time); 0];
% tbl.TimeZeroed = tbl.Time -  tbl.Time(1);



if strcmp(subject, 'P009')
    tbl(349:end,:) = [];
end



% LOOP TO CHECK FOR ABNORMAL EVENTS' DURATIONS
sample_rate = 500; % HARDCODING !! careful if changes in the future
for i=1:height(tbl)
    row = tbl(i,:);
    if ismember(row.Name,'BaselineStart') && row.DurationComputed > 6
        next_row = tbl(i + 1,:);
        next_row.Time = row.Time + 5; % 5 seconds
        next_row.Latency = row.Latency + 5*sample_rate;
        tbl(i + 1,:) = next_row;
    end

    if ismember(row.Name,'TrialStart') && row.DurationComputed > 16
        next_row = tbl(i + 1,:);
        next_row.Time = row.Time + 16; % 5 seconds
        next_row.Latency = row.Latency + 16*sample_rate;
        tbl(i + 1,:) = next_row;
    end 
end
tbl.NewDurationComputed = [diff(tbl.Time); 0]; % just to visually verify





% ADD THE EDGE COARSE BASELINE AFTER EACH BLACK COARSE BASELINE
BaseCoarse_idx = find(strcmp(tbl.Name,'BaselineCoarseStart'));
old_tbl = tbl;

for i=length(BaseCoarse_idx):-1:1
    edgeStart  = old_tbl(BaseCoarse_idx(i)+1,:);
    edgeStart.Time = edgeStart.Time + 0.01;
    edgeStart.RawName = {'Start c_baseline'};
    edgeStart.Name = {'BaselineCoarseStart'};
    edgeStart.TrialType = {['Baseline']};
    edgeStart.Latency = edgeStart.Latency + 1;
    
    edgeEnd = old_tbl(BaseCoarse_idx(i)+2,:);
    edgeEnd.Time = edgeEnd.Time - 0.01;
    edgeEnd.RawName = {'End c_baseline'};
    edgeEnd.Name = {'BaselineCoarseEnd'};
    edgeEnd.TrialType = {['Baseline']};
    edgeEnd.Latency = edgeEnd.Latency - 1;
    edgeEnd.BlockIndex = 0;
    edgeEnd.TrialIndex = 0;

    tbl = [tbl(1:BaseCoarse_idx(i)+1, :); edgeStart ; edgeEnd; tbl(BaseCoarse_idx(i)+2:end, :)];
end







% BIG CHECK UP TO VERIFY
% - 1 baseline per block
% - 2 trials of the same type per block (Angle45 or Angle20)
% - Time and Latencies are perfectly sorted, no mis-swapped rows

block_error = false;
baseline_error = false;
trial_error = false;
latency_sorted_error = false;
time_sorted_error = false;
count_error = false;

block_indices = unique(tbl.BlockIndex);
block_indices(block_indices == 0) = []; % to ignore familiarisation, training ...

for i = 1:length(block_indices)
    block_index = block_indices(i);
    block_rows = tbl(tbl.BlockIndex == block_index, :);
    
    % Check for one baseline (BaselineStart followed by BaselineEnd)
    baseline_start_indices = find(strcmp(block_rows.Name, 'BaselineStart'));
    baseline_end_indices = find(strcmp(block_rows.Name, 'BaselineEnd'));
    
    if length(baseline_start_indices) ~= 1 || length(baseline_end_indices) ~= 1 || ...
       baseline_start_indices + 1 ~= baseline_end_indices
        block_error = true;
        fprintf('Error: Block %d does not have exactly one baseline with a start and end in correct order.\n', block_index);
    end
    
    % Check for two trials of the same TrialType
    trial_start_indices = find(strcmp(block_rows.Name, 'TrialStart'));
    trial_end_indices = find(strcmp(block_rows.Name, 'TrialEnd'));
    
    if length(trial_start_indices) ~= 2 || length(trial_end_indices) ~= 2 || ...
       any(trial_start_indices + 1 ~= trial_end_indices)
        block_error = true;
        fprintf('Error: Block %d does not have exactly two trials with starts and ends in correct order.\n', block_index);
    end
    
    % Check that the two trials have the same TrialType
    trial_types = block_rows.TrialType(trial_start_indices);
    
    if (length(trial_types) ~= 2) && (trial_types{1} ~= trial_types{2})
        block_error = true;
        fprintf('Error: Block %d has trials of different types.\n', block_index);
    end
end

% 'BaselineStart' / 'TrialStart' followed by 'BaselineEnd' / 'TrialEnd'
for i = 1:height(tbl) - 1
    if strcmp(tbl.Name{i}, 'BaselineStart') && ~strcmp(tbl.Name{i+1}, 'BaselineEnd')
        baseline_error = true;
        fprintf('Error: BaselineStart at row %d is not followed by BaselineEnd.\n', i);
    end
    
    if strcmp(tbl.Name{i}, 'TrialStart') && ~strcmp(tbl.Name{i+1}, 'TrialEnd')
        trial_error = true;
        fprintf('Error: TrialStart at row %d is not followed by TrialEnd.\n', i);
    end
end

% Check that latencies are sorted in ascending order
if any(diff(tbl.Latency) <= 0)
    latency_sorted_error = true;
    fprintf('Error: Latencies are not sorted in ascending order.\n');
end

% Check that time values are sorted in ascending order
if any(diff(tbl.Time) <= 0)
    time_sorted_error = true;
    fprintf('Error: Time values are not sorted in ascending order.\n');
end

% COUNT OCCURRENCES OF TRIALS WITH 20° AND 45° ANGLE AND BASELINE
TrialBaselineCount = sum(tbl(tbl.BlockIndex ~= 0, :).Value == -1) / 2;
Trial20Count = sum(tbl(tbl.BlockIndex ~= 0, :).Value == 20);
Trial45Count = sum(tbl(tbl.BlockIndex ~= 0, :).Value == 45);
fprintf('There are %d blocks\n', length(block_indices));
fprintf('There are %d baselines within the trials\n', TrialBaselineCount);
fprintf('There are %d trials of a 20-angle degree\n', Trial20Count);
fprintf('There are %d trials of a 45-angle degree\n', Trial45Count);

if TrialBaselineCount == 0 || Trial20Count == 0 || Trial45Count == 0
    count_error = true;
    fprintf('Error: Incorrect count of baselines or trials.\n');
end

% Summary of errors
if ~baseline_error && ~trial_error && ~latency_sorted_error ...
    && ~time_sorted_error && ~block_error && ~count_error
    fprintf('All checks passed successfully.\n');
else
    error('Some checks failed. Please see the errors above.\n');
end





% events = struct('type', {}, 'latency', {}, 'duration', {});
events = struct('type', {}, 'latency', {});

for t=1:height(tbl)
    events(t).type = tbl.Name{t};
    %events(t).duration = tbl.DurationComputed(t);
    events(t).latency = tbl.Latency(t);
    events(t).TrialType = tbl.TrialType{t};
    events(t).BlockIndex = tbl.BlockIndex(t);
    events(t).TrialIndex = tbl.TrialIndex(t);
end



% Look for missing delimiters
BlockStart = strcmp({events.type}, 'BlockStart'); % Very first start ?
if ~any(BlockStart)
    if acq_ind == 1
        startEvent = events(1);
        startEvent.type = 'BlockStart';
        startEvent.latency = events(1).latency - 1;
        %startEvent.duration = NaN;
        startEvent.TrialType = NaN;
        startEvent.BlockIndex = NaN;
        startEvent.TrialIndex = NaN;
        
        events = [startEvent, events];
    else
        warning('Couldn''t find the start of the block in this file')
    end
end



return








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