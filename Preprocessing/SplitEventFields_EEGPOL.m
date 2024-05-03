function [ markers ] = SplitEventFields_EEGPOL(Stream, times)
% Split event's tag in the desired fields
%   each field is separated thanks to ';' and stored in the format key:value
%   Returns the markers struct with each field in a separate field of the
%   struct

EEGLABB_struct.type = [];
EEGLABB_struct.latency = [];
EEGLABB_struct.urevent = [];
EEGLABB_struct.urevent = [];





markers.events = Stream.time_series';
markers.time = Stream.time_stamps';
%markers.time= markers.time * (1000/EEG.srate); % conversion in ms
%markers.urevent=[EEG.event.urevent]; % WARNING: urevent field is not edited after event deletion to conserve a trace of previous indexing.
markers.index = find(markers.time); % simpler way to get the indices, but will be updated after event deletion ! ...


required_fields = {'Start fam_empty', 'End fam_empty', 'Start fam_symbol', 'End fam_symbol', 'Start baseline', 'Start training', 'End training', ...
    'Start c_baseline', 'End c_baseline', 'Start trial', 'End trial'}';

new_fields = {'FamEmptyStart', 'FamEmptyEnd', 'FamSymbolStart', 'FamSymbolEnd', 'BaselineStart', 'TrainStart', 'TrainEnd', ...
    'BaselineCooarseStart', 'BaselineCoarseEnd', 'TrialStart', 'TrialEnd'}';

dic = dictionary(required_fields, new_fields);


tbl = table(markers.index, markers.time, markers.events(:,1), repmat({NaN}, length(markers.time),1), cellfun(@str2double, markers.events(:,2)), ...
    cellfun(@str2double, markers.events(:,3)), 'VariableNames', {'Index', 'Time', 'RawName', 'Name', 'Value', 'Duration'});

for i=1:height(tbl)
    row = tbl(i,:);

    if isKey(dic,row.RawName) % Rename the raw name with cleaner name (see new fields above)
        row.Name = dic(row.RawName);
    else
        error('Unkown name of event ?')
    end
    
    if ~ismember(row.Value, [20, 45]) % Only keep 20 and 45 values
        row.Value = NaN;
    end
    tbl(i,:) = row;
end

tbl.DurationComputed = [diff(tbl.Time); 0];
tbl.TimeZeroed = tbl.Time -  tbl.Time(1);

markers.eventsFull = tbl;
clear tbl



EEGLABB_struct.type = [];
EEGLABB_struct.duration = [];
EEGLABB_struct.BlockIndex = [];
EEGLABB_struct.TrialIndex = [];
EEGLABB_struct.TrialType = [];

EEGLABB_struct.urevent = [];

%IN EXPORT
EEGLABB_struct.latency = [];


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





% 
% %% Extract EEG + Markers according to nomenclature:
% keys = {'block', 'trial', 'type', 'phase', 'event', 'duration', 'name', 'location', 'distance'};
% markers.block = NaN(length(markers.time),1);
% markers.trial = NaN(length(markers.time),1);
% markers.type = cell(length(markers.time),1);
% markers.phase = cell(length(markers.time),1);
% markers.event = cell(length(markers.time),1);
% markers.duration = NaN(length(markers.time),1);
% markers.name = cell(length(markers.time),1);
% markers.location = cell(length(markers.time),1);
% markers.distance = NaN(length(markers.time),1);
% 
% for m=1:length(markers.time)
%     pairs=strsplit(markers.fullName{m}, ';');
%     %disp(m);
%     for k=1:length(keys)
%         if (sum(contains(pairs, keys{k}))==1)
%             pair = pairs(contains(pairs, keys{k}));
%             if contains(pair{1},':')
%                 splitPair = strsplit(pair{1}, ':');
%                 value = splitPair{2};
%             else
%                 % Sometimes the ':' was forgotten
%                 value = pair{1}(length(keys{k})+1:end);
%             end
% 
%             if any(strcmp({'block', 'trial', 'duration', 'distance'}, keys{k}))
%                 % Numerical value to interpret
%                 markers.(keys{k})(m)= str2num(value);
%             else
%                 markers.(keys{k}){m}= value;
%             end
%         else
%             if any(strcmp({'type', 'phase', 'event', 'name', 'location'}, keys{k}))
%                 % Default value for strings
%                 markers.(keys{k}){m}= 'none';
%             end
%         end
%     end
% end
% 
% %markers = struct2table(markers);
% end

