function [gaze_events] = export_gaze_events4Polyretina(streams, p, events)
% Formats the gaze data accordingly in order to do a gaze analysis for the
% EEG experiment

gaze_events = struct('time_stamps', [], 'x_position', [], 'y_position', [], 'diameter', [], 'field_of_view', [], 'trial_id', []);

% Looping over the streams to get which one is the gaze data
for stream_number = 1:size(streams,2)
    
    if strcmp(streams{stream_number}.info.name,'Sim-EyeGaze')
        stream_gaze_data = streams{stream_number};        
    end

end

% For each time stamp, define which field of vision participant was under:
corresponding_FoV = [];
corresponding_trial_id = [];
time_stamp_list = [];
x_position_list = [];
y_position_list = [];
diameter_list = [];

parfor time = 1:size(stream_gaze_data.time_stamps,2)
    
    disp(['Computed: ' num2str(time/length(stream_gaze_data.time_stamps)*100) '% for participant ' num2str(p)])
    
    % Filling up systematically the data
    time_stamp_list = [time_stamp_list; stream_gaze_data.time_stamps(time)];
    x_position_list = [x_position_list; stream_gaze_data.time_series(1,time)];
    y_position_list = [y_position_list; stream_gaze_data.time_series(2, time)];
    diameter_list = [diameter_list; stream_gaze_data.time_series(3, time)];
    
    if stream_gaze_data.time_stamps(time) < events(1).time_stamp
        % Gaze data is before any events of interest
        corresponding_FoV = [corresponding_FoV; -1];
        corresponding_trial_id = [corresponding_trial_id; -1];
        continue
    elseif stream_gaze_data.time_stamps(time) > events(end).time_stamp
        % Gaze data is after any events of interest
        corresponding_FoV = [corresponding_FoV; -1];
        corresponding_trial_id = [corresponding_trial_id; -1];
        continue
    else
        % We are within the events of interest to get FoV and trial_id
        
        % Turning events times into a table to take substraction and then
        % get index of interst
        events_table = struct2table(events);
        array_events_times = table2array(events_table(:,2));
        [~, closest_index] = min(abs(stream_gaze_data.time_stamps(time) - array_events_times));
        corresponding_FoV = [corresponding_FoV; events(closest_index).field_of_view];
        corresponding_trial_id = [corresponding_trial_id; events(closest_index).trial_id];
        
    end
    
    
end

% Inserting the data into the structure of interest
gaze_events.time_stamps = time_stamp_list;
gaze_events.x_position = x_position_list;
gaze_events.y_position = y_position_list;
gaze_events.diameter = diameter_list;
gaze_events.field_of_view = corresponding_FoV;
gaze_events.trial_id = corresponding_trial_id;

    

end

