function [standard_deviation_20,standard_deviation_45, standard_deviation_110,summed_norm_20,summed_norm_45,summed_norm_110] = compute_gaze_std(results_participant_gaze, p)
% Computing the standard deviation for each field of view

standard_deviation_20 = [];
standard_deviation_45 = [];
standard_deviation_110 = [];

% Computing sum of norms (how much they watched)
summed_norm_20 = [];
summed_norm_45 = [];
summed_norm_110 = [];


% Retrieving how many trials are present in the dataset
number_of_trials = results_participant_gaze.(['Results_participant_' num2str(p)]).gaze_events.trial_id(end);

% Looping over every trial
for trial = 1:number_of_trials
    
    % Defining variable of interest for each trial
    norm_vector = [];
    
    % Find the indices linked to that trial
    indices_trial = find(ismember(results_participant_gaze.(['Results_participant_' num2str(p)]).gaze_events.trial_id, [trial]));
    
    % Looping over every time_stamp of that trial to get positions
    for current_index = 1:length(indices_trial)-1
        
        % Retrieving current positions to compute norm
        current_x = results_participant_gaze.(['Results_participant_' num2str(p)]).gaze_events.x_position(indices_trial(current_index));
        next_x = results_participant_gaze.(['Results_participant_' num2str(p)]).gaze_events.x_position(indices_trial(current_index+1));
        
        current_y = results_participant_gaze.(['Results_participant_' num2str(p)]).gaze_events.y_position(indices_trial(current_index));
        next_y = results_participant_gaze.(['Results_participant_' num2str(p)]).gaze_events.y_position(indices_trial(current_index+1));
        
        % Computing the norm
        current_norm = sqrt((current_x - next_x)^2 + (current_y - next_y)^2 );
        
        % Adding norm to vector
        norm_vector = [norm_vector; current_norm];
        
    end
    
    % Seeing for that trial which is the field of view of interest
    current_FoV = results_participant_gaze.(['Results_participant_' num2str(p)]).gaze_events.field_of_view(indices_trial(end));
    
    if current_FoV == 20
        standard_deviation_20 = [standard_deviation_20; std(norm_vector)];
        summed_norm_20 = [summed_norm_20; sum(norm_vector, 'all', 'omitnan')];
    elseif current_FoV == 45
        standard_deviation_45 = [standard_deviation_45; std(norm_vector)];
        summed_norm_45 = [summed_norm_45; sum(norm_vector, 'all', 'omitnan')];
    elseif current_FoV == 110
        standard_deviation_110 = [standard_deviation_110; std(norm_vector)];
        summed_norm_110 = [summed_norm_110; sum(norm_vector, 'all', 'omitnan')];
    end
    
end

end
