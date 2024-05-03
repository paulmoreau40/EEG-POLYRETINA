function [ markers ] = SplitEventFields_EEGAFF(Stream)
% Split event's tag in the desired fields
%   each field is separated thanks to ';' and stored in the format key:value
%   Returns the markers struct with each field in a separate field of the
%   struct

markers.fullName = Stream.time_series';
markers.time = Stream.time_stamps';
%markers.time= markers.time * (1000/EEG.srate); % conversion in ms
%markers.urevent=[EEG.event.urevent]; % WARNING: urevent field is not edited after event deletion to conserve a trace of previous indexing.
markers.index = find(markers.time); % simpler way to get the indices, but will be updated after event deletion ! ...

%% Extract EEG + Markers according to nomenclature:
keys = {'block', 'trial', 'type', 'phase', 'event', 'duration', 'name', 'location', 'distance'};
markers.block = NaN(length(markers.time),1);
markers.trial = NaN(length(markers.time),1);
markers.type = cell(length(markers.time),1);
markers.phase = cell(length(markers.time),1);
markers.event = cell(length(markers.time),1);
markers.duration = NaN(length(markers.time),1);
markers.name = cell(length(markers.time),1);
markers.location = cell(length(markers.time),1);
markers.distance = NaN(length(markers.time),1);

for m=1:length(markers.time)
    pairs=strsplit(markers.fullName{m}, ';');
    %disp(m);
    for k=1:length(keys)
        if (sum(contains(pairs, keys{k}))==1)
            pair = pairs(contains(pairs, keys{k}));
            if contains(pair{1},':')
                splitPair = strsplit(pair{1}, ':');
                value = splitPair{2};
            else
                % Sometimes the ':' was forgotten
                value = pair{1}(length(keys{k})+1:end);
            end
            
            if any(strcmp({'block', 'trial', 'duration', 'distance'}, keys{k}))
                % Numerical value to interpret
                markers.(keys{k})(m)= str2num(value);
            else
                markers.(keys{k}){m}= value;
            end
        else
            if any(strcmp({'type', 'phase', 'event', 'name', 'location'}, keys{k}))
                % Default value for strings
                markers.(keys{k}){m}= 'none';
            end
        end
    end
end

%markers = struct2table(markers);
end

