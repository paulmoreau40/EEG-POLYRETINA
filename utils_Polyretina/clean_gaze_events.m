function [results_participant_gaze] = clean_gaze_events(results_participant_gaze,p)
% Cleaning gaze data that is before and after the experiment

% Step 1: Removing the indices which are out of interest for us
% Getting indices of trial data that is of no interst to us
OoI_indices = find(~ismember(results_participant_gaze.(['Results_participant_' num2str(p)]).gaze_events.field_of_view, [20,45,110]));

% Remove it from gaze data structure for this participant in question
results_participant_gaze.(['Results_participant_' num2str(p)]).gaze_events.time_stamps(OoI_indices) = [];
results_participant_gaze.(['Results_participant_' num2str(p)]).gaze_events.x_position(OoI_indices) = [];
results_participant_gaze.(['Results_participant_' num2str(p)]).gaze_events.y_position(OoI_indices) = [];
results_participant_gaze.(['Results_participant_' num2str(p)]).gaze_events.diameter(OoI_indices) = [];
results_participant_gaze.(['Results_participant_' num2str(p)]).gaze_events.field_of_view(OoI_indices) = [];
results_participant_gaze.(['Results_participant_' num2str(p)]).gaze_events.trial_id(OoI_indices) = [];

% Step 2: Replacing the values where the data is wrong with NaNs
results_participant_gaze.(['Results_participant_' num2str(p)]).gaze_events.x_position(results_participant_gaze.(['Results_participant_' num2str(p)]).gaze_events.x_position < -0.1) = NaN;
results_participant_gaze.(['Results_participant_' num2str(p)]).gaze_events.x_position(results_participant_gaze.(['Results_participant_' num2str(p)]).gaze_events.x_position > 0.1) = NaN;
results_participant_gaze.(['Results_participant_' num2str(p)]).gaze_events.y_position(results_participant_gaze.(['Results_participant_' num2str(p)]).gaze_events.y_position < -0.1) = NaN;
results_participant_gaze.(['Results_participant_' num2str(p)]).gaze_events.y_position(results_participant_gaze.(['Results_participant_' num2str(p)]).gaze_events.y_position > 0.1) = NaN;



end
