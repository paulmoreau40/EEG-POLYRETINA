function EEG_trial_data = extract_segments_EEG_compute_spectrum(EEG_no_doublons, EEG_trial_data, segment_of_interest, participant_id, make_plot)
% Custom function to extract baselines of interest and return a structure with the segmented regions of interest and associated metadata
%
% Inputs:
% EEG_no_doublons           - EEG struct where Button Press doublons need to be removed
% segment_of_interst        - Defining what segment in time is of interest
%                             to us
% participant_id            - id of participant
% make_plot                 - whether or not to make spectoplot
%
% Outputs:
% EEG_trial_data            = Struct with EEG segments and their respective
%                             time trials and metadata


% Defining Start and End of segments of interest
switch segment_of_interest
    case 'detection_time'
        begin_segment = {' Start training <b><color=#ff6040ff>110</color></b>'; ...
                         ' Start training <b><color=#ff6040ff>45</color></b>'; ...
                         ' Start training <b><color=#ff6040ff>20</color></b>'; ...
                         ' Start trial <b><color=#ff6040ff>110</color></b>'; ...
                         ' Start trial <b><color=#ff6040ff>45</color></b>'; ...
                         ' Start trial <b><color=#ff6040ff>20</color></b>'};
        end_segment = {'Button_Press'};
    case 'press_time'
        begin_segment = {'Button_Press'};
        end_segment = {' End training'; ...
                       ' End trial'};
    case 'total_time'
        begin_segment = {' Start training <b><color=#ff6040ff>110</color></b>'; ...
                         ' Start training <b><color=#ff6040ff>45</color></b>'; ...
                         ' Start training <b><color=#ff6040ff>20</color></b>'; ...
                         ' Start trial <b><color=#ff6040ff>110</color></b>'; ...
                         ' Start trial <b><color=#ff6040ff>45</color></b>'; ...
                         ' Start trial <b><color=#ff6040ff>20</color></b>'};
        end_segment = {' End training'; ...
                       ' End trial'};
    otherwise
        error('Segment of interest not correct: please try "detection_time" , "press_time" , "total_time" , "baseline" ')
end

% Which line we are at currently for metaInfo
current_line_trial = size({EEG_trial_data.metaInfo.trial_id},2);

% Creating Structure
EEG_trial_data.(['P' participant_id]).srate = EEG_no_doublons.srate;
EEG_trial_data.(['P' participant_id]).chanlocs = EEG_no_doublons.chanlocs;

% Retrieving how many trials are present
num_trials = EEG_no_doublons.event(end).trial_id;

% Issue with P2 data, see if can be fixed
if strcmp(participant_id, '2')
    num_trials = num_trials - 1; % Not considering trial 18 which only consists of ' Start Baseline' ....
end

% Issue with P1 data:
if strcmp(participant_id, '1')
    % For participant 1: after removing bad segments, skips from 73 to 108
    num_trials = EEG_no_doublons.event(end-3).trial_id - 1;
end

% Issue with P3 data:
if strcmp(participant_id, '3')
    % For participant 3: the last trial (87) is not complete
    num_trials = EEG_no_doublons.event(end).trial_id - 1;
end

% Creating Final Structure:
spectrum = [];

switch segment_of_interest
    case 'detection_time'
        
        for trial = 1:num_trials
            % Resetting values to identify if missing tag:
            start_index = -1;
            end_index = -1;

            % Retrieve indices of trial
            intermediate_table = struct2table(EEG_no_doublons.event);
            trial_indices = find(table2array(intermediate_table(:,5)) == trial);
            
            % Retrieving information for metadata
            field_of_view = EEG_no_doublons.event(trial_indices(end)).field_of_view;
            auditory_cue = EEG_no_doublons.event(trial_indices(end)).auditory_cue;
            auditory_answer = EEG_no_doublons.event(trial_indices(end)).auditory_answer;
            
            % Alexandre's nicer code: 
%             starts = contains({EEG_no_doublons.event(:).type}, 'Start training') | contains({EEG_no_doublons.event(:).type}, 'Start trial');
%             start_trial = starts & EEG_no_doublons.event(:).trial_id == trial;
%             
%             if isempty(start_trial)
%                 % There are no tags indicating the beginning of a trial,
%                 usng the tag of the end of trial of previous one
%             end
            
            for i = 1:length(trial_indices)
                % Looping to find start and end 
                if (ismember(EEG_no_doublons.event(trial_indices(i)).type, begin_segment)) && (start_index == -1)
                    start_index = trial_indices(i);
                elseif ismember(EEG_no_doublons.event(trial_indices(i)).type, end_segment) && (end_index == -1)
                    end_index = trial_indices(i);
                else
                    continue
                end
            end

            % Check if tags are missing or not: if missing than take the tag that
            % precedes/follows but note it in metadata
            if (start_index == -1)
                % Get the index of the tag that precedes (ie end baseline)
                for i = 1:length(trial_indices)
                    if ismember(EEG_no_doublons.event(trial_indices(i)).type, {' End baseline'})
                        start_index = trial_indices(i);
                    end
                end
            elseif (end_index == -1)
                % Button was not pressed: so selecting end of trial as end
                % of segment of interest
                for i = 1:length(trial_indices)
                    if ismember(EEG_no_doublons.event(trial_indices(i)).type, {' End training';' End trial'})
                        end_index = trial_indices(i);
                    end
                end
            end
            
            % Filling up the table with segment of interest and metadata
            % 1. Find time stamps to retrieve data
            start_segment_time_index = round(EEG_no_doublons.event(start_index).latency);
            end_segment_time_index = round(EEG_no_doublons.event(end_index).latency);
            
            data = EEG_no_doublons.data(:,start_segment_time_index:end_segment_time_index);
            
            % 2. Retrieve MetaData          
            EEG_trial_data.metaInfo(current_line_trial + trial).participant_id = participant_id;
            EEG_trial_data.metaInfo(current_line_trial + trial).trial_id = trial;
            EEG_trial_data.metaInfo(current_line_trial + trial).field_of_view = field_of_view;
            EEG_trial_data.metaInfo(current_line_trial + trial).auditory_cue = auditory_cue;
            EEG_trial_data.metaInfo(current_line_trial + trial).auditory_answer = auditory_answer;
            

            if make_plot
                % Need to close all other windows open because spectopo
                % will panic otherwise
                close all;
                % Setting plot option to true
                plot_option = 'on';
            else
                plot_option = 'off';
            end
            
            % 3. Compute spectrum            
            [spectrum(:,:,trial), freqs, ~, ~, ~] =...
                spectopo(data, size(data,2), EEG_no_doublons.srate,...
                'freq', [10.0, 20.0], ... %freq_of_interest, ... %10,...                %'plotchan', [],...
                'chanlocs', EEG_no_doublons.chanlocs,...
                'freqfac', 2,...
                'winsize', 250, ... % fréquence d'échantillonnage
                'overlap', 125,... 
                'wintype','hamming',...
                'freqrange',[2 40], ...
                'plot', plot_option,...
                'plotmean', 'off',...
                'verbose','off');
            
            disp(['Finished computing spectrum for trial ' num2str(trial) '/' num2str(num_trials) ' for Participant ' participant_id])
        end
        
        % Saving frequencie and spectre for participant considered
        EEG_trial_data.(['P' participant_id]).freqs = freqs;
        disp('Converting spectrum into linear scale for each TRIAL');
        EEG_trial_data.(['P' participant_id]).spectrum = 10.^(spectrum./10);
        
    otherwise
        error('Case not yet considered, please code it')
end

end

