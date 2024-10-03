function events = export_events_4Polyertina_behavioral(streams, participant_id, path_raw_data)
% Interpret event XDF streams for 4SC experiment to put them in the EEGLAB set
% Inputs:
%   - AllStreams :          streams loaded with load_xdf.
%   - times :               Complete time points vector of the recording
%                               (to compute latencies)
%   - subject:       id of participant that is currently considered
%
% Outputs:
%   - events :              Structure containing all events loaded.

% 2. Restructuring the streams into events file for EEG data

events = [];

disp('Identifying streams...')

for stream_id=1:size(streams,2)
    % Looping over every stream and identifying which one is which (the
    % index may be subject to change across participants)
           
    if strcmp(streams{stream_id}.info.name,'ControlCenter-Tags')
        % This stream corresponds to the tags sent by the Control Center to
        % keep track of the events of the EEG Protocol (begin
        % familiarizaiton, end trial...)
        stream_id_control_center_tags = stream_id;
        stream_control_center_tags = streams{stream_id_control_center_tags};
                
    elseif strcmp(streams{stream_id}.info.name,'Sim-VisionConfirmation')
        % This stream corresponds to the vision of the participant defined
        % by the Simulation        
        stream_id_vision = stream_id;
        stream_vision = streams{stream_id_vision};

    elseif strcmp(streams{stream_id}.info.name,'Sim-ButtonPress')
        % This stream corresponds to the tags sent when the participant
        % presses on the VR controller        
        stream_id_button_press = stream_id;
        stream_button_press = streams{stream_id_button_press};
    end
end

disp('...completed')

% 3. Saving the timestamps into an event structure adapted for EEGLAB

disp('Loading event tags for the control center...')

% SAVING THE PROTOCOL TAGS AND TIMESTAMPS FROM CONTROL CENTER
start = 1;
stop = length(stream_control_center_tags.time_stamps);

for t = start:stop % loop will be stopped before the end given that the last tags are a repetition of 'End trial'...
    
    % Considering specificities of participants where last trials not
    % concluded
    if (participant_id == 3) && (t == 442)
        continue
    end
        
    % Saving the type of event of the EEG protocol as well as its
    % corresponding time stamp
    events(t).type = stream_control_center_tags.time_series(t);
    events(t).time_stamp = stream_control_center_tags.time_stamps(t);
    events(t).duration = 1;

    % Completing field of view column
    if (strcmp(stream_control_center_tags.time_series(t), 'EyeCalibration')) || ...
            (strcmp(stream_control_center_tags.time_series(t), ' Start baseline')) || ...
            (strcmp(stream_control_center_tags.time_series(t), ' End baseline'))
        % Setting the field of view for eye calibration and
        % baseline tags to 0
        events(t).field_of_view = 0;
        
    elseif (strcmp(stream_control_center_tags.time_series(t), ' Start trial <b><color=#ff6040ff>110</color></b>')) || ...
            (strcmp(stream_control_center_tags.time_series(t), ' Start trial <b><color=#ff6040ff>45</color></b>')) || ...
            (strcmp(stream_control_center_tags.time_series(t), ' Start trial <b><color=#ff6040ff>20</color></b>')) || ...
            (strcmp(stream_control_center_tags.time_series(t), ' Start training <b><color=#ff6040ff>110</color></b>')) || ...
            (strcmp(stream_control_center_tags.time_series(t), ' Start familiar_empty <b><color=#ff6040ff>110</color></b>')) || ...
            (strcmp(stream_control_center_tags.time_series(t), ' Start familiar_symbol <b><color=#ff6040ff>110</color></b>')) || ...
            (strcmp(stream_control_center_tags.time_series(t), ' Start training <b><color=#ff6040ff>45</color></b>')) || ...
            (strcmp(stream_control_center_tags.time_series(t), ' Start familiar_empty <b><color=#ff6040ff>45</color></b>')) || ...
            (strcmp(stream_control_center_tags.time_series(t), ' Start familiar_symbol <b><color=#ff6040ff>45</color></b>')) || ...
            (strcmp(stream_control_center_tags.time_series(t), ' Start training <b><color=#ff6040ff>20</color></b>')) || ...
            (strcmp(stream_control_center_tags.time_series(t), ' Start familiar_empty <b><color=#ff6040ff>20</color></b>')) || ...
            (strcmp(stream_control_center_tags.time_series(t), ' Start familiar_symbol <b><color=#ff6040ff>20</color></b>'))
        % We are finding the field of view condition corresponding to the
        % start trial, start training or familiarizations
        [~,closest_timestamp] = min(abs(stream_control_center_tags.time_stamps(t)-stream_vision.time_stamps));
        events(t).field_of_view = str2num(stream_vision.time_series{closest_timestamp});
        % Setting the 'End trial' also to the same field of view condition
        events(t+1).field_of_view = str2num(stream_vision.time_series{closest_timestamp});
        
    end
    
    if (strcmp(stream_control_center_tags.time_series(t), ' Start training <b><color=#ff6040ff>110</color></b>'))
        % Getting the time stamp of interest where the experiment actually
        % begins
        time_stamp_start_experiment = stream_control_center_tags.time_stamps(t);
    end
    if (participant_id ~= 3) % Participant 3 did not finish trial
        % For participants 1 to 5: the last tag 'End Trial' of the last
        % trial is repeated several times. This is used as the condition to
        % stop the loop
        if ismember(participant_id, [1,2,3,4,5])
            if (strcmp(stream_control_center_tags.time_series(t), stream_control_center_tags.time_series(t+1))) && (~strcmp(stream_control_center_tags.time_series(t), 'EyeCalibration'))
                % If the following trial has the same name as the current
                % trial, we are either in the ‘End Trial' repetition at the
                % end of the experiment (excluding the eventual
                % 'EyeCalibrations' duplicates)
                break
            end
        elseif (participant_id == 8) && (t == length(stream_control_center_tags.time_series))
            % For participant 8: no repetition of the last tag, so previous
            % condition does not work. Code adapted accordingly
            break
        end
    end
    
end

disp('...completed')

% LOOPING OVER THE BUTTON PRESSED AND INSERTING THEM AND THEIR TIME STAMP
% INTO THE STRUCTURE, THEN SORTING ACCORDING TO TIME_STAMP
start=1;
stop=length(stream_button_press.time_stamps);
length_events_only_cctags = size(events,2);
counter_button_press = 1;

disp('Adding events of button press...')

for t = start:stop
    
    % Disregarding all of the tags pressed before the familiarization:
    if stream_button_press.time_stamps(t) < time_stamp_start_experiment
        continue
    end
    
    % Incorporating the button press events in the code
    events(length_events_only_cctags+counter_button_press).time_stamp = stream_button_press.time_stamps(t);
    events(length_events_only_cctags+counter_button_press).type = 'Button_Press';
    events(length_events_only_cctags+counter_button_press).duration = 1;
    counter_button_press = counter_button_press + 1;
    
    
end
disp('... completed')

% 4. Sorting Tags according to time_stamp

% Sorting the time stamps of the event files so that the button press tags
% are correclty incorporated

disp('Adding the missing fields of view for the Button Press Tags, sorting the Events table chronologically, and adding trial ID...')


T = struct2table(events); % convert the struct array to a table
sortedT = sortrows(T, 'time_stamp'); % sort the table by 'time_stamp'
events = table2struct(sortedT);


start = 1;
stop = size(events,1);

trial_count = 1;

for t = start:stop
    
    % Considering Button Presses and Adding Corresponding Field of Views
    if (strcmp(events(t).type, ' Start trial <b><color=#ff6040ff>110</color></b>')) || ...
            (strcmp(events(t).type, ' Start trial <b><color=#ff6040ff>45</color></b>')) || ...
            (strcmp(events(t).type, ' Start trial <b><color=#ff6040ff>20</color></b>')) || ...
            (strcmp(events(t).type, ' Start training <b><color=#ff6040ff>110</color></b>')) || ...
            (strcmp(events(t).type, ' Start familiar_empty <b><color=#ff6040ff>110</color></b>')) || ...
            (strcmp(events(t).type, ' Start familiar_symbol <b><color=#ff6040ff>110</color></b>')) || ...
            (strcmp(events(t).type, ' Start training <b><color=#ff6040ff>45</color></b>')) || ...
            (strcmp(events(t).type, ' Start familiar_empty <b><color=#ff6040ff>45</color></b>')) || ...
            (strcmp(events(t).type, ' Start familiar_symbol <b><color=#ff6040ff>45</color></b>')) || ...
            (strcmp(events(t).type, ' Start training <b><color=#ff6040ff>20</color></b>')) || ...
            (strcmp(events(t).type, ' Start familiar_empty <b><color=#ff6040ff>20</color></b>')) || ...
            (strcmp(events(t).type, ' Start familiar_symbol <b><color=#ff6040ff>20</color></b>'))
        
        labeling_button_presses = 1;
        
        while ~(contains(events(t+labeling_button_presses).type, 'End'))
            events(t+labeling_button_presses).field_of_view = events(t).field_of_view;
            labeling_button_presses = labeling_button_presses + 1;
        end
    end
    
    % Counting the number of trials which was done by the participant
    % (normally 107 trials, but eventually less)
    if (strcmp(events(t).type, 'EyeCalibration')) || ...
            (strcmp(events(t).type, ' Start familiar_empty <b><color=#ff6040ff>110</color></b>')) || ...
            (strcmp(events(t).type, ' Start familiar_symbol <b><color=#ff6040ff>110</color></b>')) || ...
            (strcmp(events(t).type, ' Start familiar_empty <b><color=#ff6040ff>45</color></b>')) || ...
            (strcmp(events(t).type, ' Start familiar_symbol <b><color=#ff6040ff>45</color></b>')) || ...
            (strcmp(events(t).type, ' Start familiar_empty <b><color=#ff6040ff>20</color></b>')) || ...
            (strcmp(events(t).type, ' Start familiar_symbol <b><color=#ff6040ff>20</color></b>')) || ...
            (strcmp(events(t).type, ' End familiar_empty')) || ...
            (strcmp(events(t).type, ' End familiar_symbol'))
        events(t).trial_id = 0;
        
        
    elseif (contains(events(t).type, ' Start baseline'))
        events(t).trial_id = trial_count;
        
        labeling_every_event = 1;
        
        while ~(contains(events(t+labeling_every_event).type,' Start baseline')) && (t+labeling_every_event < stop)
            events(t+labeling_every_event).trial_id = trial_count;
            labeling_every_event = labeling_every_event + 1;
        end
        
        % Incrementing the trial count for participants
        trial_count = trial_count +1;
    end
    
    % Considering last trial separately (given issues related to indexing):
    if (t == stop)
        events(t).trial_id = trial_count-1;
    end
end
disp('... completed')

% 5. Adding the Auditory Cues
disp('Loading the data related to the auditory task...')

% Loading the csv with the auditory cues generated by unity
if participant_id~=2
    auditory_cues_csv = readtable(fullfile(path_raw_data,['P' num2str(participant_id) filesep 'P' num2str(participant_id) '-auditory_cues.csv']),  'Delimiter', ' ');%, 'PreserveVariableNames', true);
    auditory_answers_csv = readtable(fullfile(path_raw_data, ['P' num2str(participant_id) filesep 'P' num2str(participant_id) '-answers.csv']));
    % Defining global variable: the total number of trials
    total_number_trials = length(auditory_answers_csv.NumberTouched_);

elseif participant_id == 2
    auditory_cues_and_answers_csv = readtable(fullfile(path_raw_data, ['P' num2str(participant_id) filesep 'P' num2str(participant_id) '-answers.csv']));
    % Defining global variable: the total number of trials
    total_number_trials = length(auditory_cues_and_answers_csv.NumberTouched_);

end

% Creating structure which will contain auditory cues and response given by
% participant, this will then be incorporated into the events structure
% formated for the EEG
auditory_task = [];

for t = 1:total_number_trials % Given that first auditory cue is in the header
    
    if (participant_id ~=2)
        if t == 1
            % Récupérer le premier élément du csv:
            auditory_task(t).cue = get_value_into_number(auditory_cues_csv.Properties.VariableNames{1});
            auditory_task(t).answer = auditory_answers_csv.NumberTouched_(t); % first 8 lines are comments written down regarding participant's performance during familiarization, skipped here
        else
            auditory_task(t).cue = get_value_into_number(auditory_cues_csv{t-1,1}{1});
            auditory_task(t).answer = auditory_answers_csv.NumberTouched_(t); % first 8 lines are comments written down regarding participant's performance during familiarization, skipped here
        end
    elseif (participant_id == 2)
        % Récupérer le premier élément du csv:
        auditory_task(t).cue = auditory_cues_and_answers_csv.NumberAsked_(t);
        auditory_task(t).answer = auditory_cues_and_answers_csv.NumberTouched_(t); % first 8 lines are comments written down regarding participant's performance during familiarization, skipped here
    end

end
disp('...completed')

% 7. Incorporating the auditory task into the events table
disp('Incorporating the auditory task and answers in the event structure...')

start = 1;
stop = size(events,1);

for t = start:stop
    % Looping over every timestmap and adding the corresponding auditory
    % cues and answers
    
    
    % Counting the number of trials which was done by the participant
    % (normally 107 trials, but eventually less)
    if (strcmp(events(t).type, 'EyeCalibration')) || ...
            (strcmp(events(t).type, ' Start familiar_empty <b><color=#ff6040ff>110</color></b>')) || ...
            (strcmp(events(t).type, ' Start familiar_symbol <b><color=#ff6040ff>110</color></b>')) || ...
            (strcmp(events(t).type, ' Start familiar_empty <b><color=#ff6040ff>45</color></b>')) || ...
            (strcmp(events(t).type, ' Start familiar_symbol <b><color=#ff6040ff>45</color></b>')) || ...
            (strcmp(events(t).type, ' Start familiar_empty <b><color=#ff6040ff>20</color></b>')) || ...
            (strcmp(events(t).type, ' Start familiar_symbol <b><color=#ff6040ff>20</color></b>')) || ...
            (strcmp(events(t).type, ' End familiar_empty')) || ...
            (strcmp(events(t).type, ' End familiar_symbol'))
        % If we are in the eye calibration of the familiarization phases,
        % then there is no auditory cue nor answer (they were given out
        % loud but not analyzed here)
        events(t).auditory_cue = NaN;
        events(t).auditory_answer = NaN;
        
        
    elseif (contains(events(t).type, ' Start baseline'))
        % Identifying which trial we are at by retrieving the id
        current_trial = events(t).trial_id;
        
        if current_trial <= size(auditory_task,2)
            % If there is a cue was given for this trial (some cues missing
            % in the end, and the cues are therefore NaN for them)
            events(t).auditory_cue = auditory_task(current_trial).cue;
            events(t).auditory_answer = auditory_task(current_trial).answer;
        else 
            events(t).auditory_cue = NaN;
            events(t).auditory_answer = NaN;
        end
        
        labeling_auditory_task = 1;
        
        while ~(contains(events(t+labeling_auditory_task).type,' Start baseline')) && (t+labeling_auditory_task < stop)
            if current_trial <= size(auditory_task,2)
                % If there is a cue was given for this trial (some cues missing
                % in the end)
                events(t+labeling_auditory_task).auditory_cue = auditory_task(current_trial).cue;
                events(t+labeling_auditory_task).auditory_answer = auditory_task(current_trial).answer;
            else 
                events(t+labeling_auditory_task).auditory_cue = NaN;
                events(t+labeling_auditory_task).auditory_answer = NaN;
            end
            labeling_auditory_task = labeling_auditory_task + 1;
        end
        
    end
    
    % Considering last "End trial" separately (given issues related to indexing):
    if (t == stop)
        events(t).auditory_cue = events(t-1).auditory_cue;
        events(t).auditory_answer = events(t-1).auditory_answer;
    end
    
end

% Considering exceptions for following participants:
%   - Participant 1: Last trial was not completed (answer missing) and not
%   sure why
%   - Participant 3: Felt ill in last trials, so need to also discard them

for t=start:stop
    if participant_id == 1
        % Disregard last trial (fill auditory cue and answer to empty)
        if (events(t).trial_id == 108) % Ignoring last trial
            events(t).auditory_cue = [];
            events(t).auditory_answer = [];
        end
    elseif participant_id == 3
        % Disregard last couple of trials (fill auditory cue and answer to
        % empty)
        if (events(t).trial_id == 108) || ...
                (events(t).trial_id == 107)
            events(t).auditory_cue = [];
            events(t).auditory_answer = [];
        end
    end
end

disp('...completed')

end