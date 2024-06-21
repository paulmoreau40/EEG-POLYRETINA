function EEG_trial_data = extract_segments_EEG_compute_spectrum(EEG, EEG_trial_data, segment_of_interest, participant_id, make_plot)
% Custom function to extract baselines of interest and return a structure with the segmented regions of interest and associated metadata
%
% Inputs:
% EEG           - EEG struct where Button Press doublons need to be removed
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
    % case 'detection_time'
    %     begin_segment = {' Start training <b><color=#ff6040ff>110</color></b>'; ...
    %                      ' Start training <b><color=#ff6040ff>45</color></b>'; ...
    %                      ' Start training <b><color=#ff6040ff>20</color></b>'; ...
    %                      ' Start trial <b><color=#ff6040ff>110</color></b>'; ...
    %                      ' Start trial <b><color=#ff6040ff>45</color></b>'; ...
    %                      ' Start trial <b><color=#ff6040ff>20</color></b>'};
    %     end_segment = {'Button_Press'};
    % case 'press_time'
    %     begin_segment = {'Button_Press'};
    %     end_segment = {' End training'; ...
    %                    ' End trial'};
    % case 'total_time'
    %     begin_segment = {' Start training <b><color=#ff6040ff>110</color></b>'; ...
    %                      ' Start training <b><color=#ff6040ff>45</color></b>'; ...
    %                      ' Start training <b><color=#ff6040ff>20</color></b>'; ...
    %                      ' Start trial <b><color=#ff6040ff>110</color></b>'; ...
    %                      ' Start trial <b><color=#ff6040ff>45</color></b>'; ...
    %                      ' Start trial <b><color=#ff6040ff>20</color></b>'};
    %     end_segment = {' End training'; ...
    %                    ' End trial'};
    case 'trialPOL'
        begin_segment = {'TrialStart'};
        end_segment = {' TrialEnd'};

    otherwise
        error('Segment of interest not correct: please try "detection_time" , "press_time" , "total_time" , "baseline" ')
end

% Which line we are at currently for metaInfo
current_line_trial = size({EEG_trial_data.metaInfo.TrialIndex},2);

% Creating Structure
EEG_trial_data.([participant_id]).srate = EEG.srate;
EEG_trial_data.([participant_id]).chanlocs = EEG.chanlocs;

% Retrieving how many trials are present
num_trials = EEG.event(end).TrialIndex * EEG.event(end).BlockIndex;

% % Issue with P2 data, see if can be fixed
% if strcmp(participant_id, '2')
%     num_trials = num_trials - 1; % Not considering trial 18 which only consists of ' Start Baseline' ....
% end
% 
% % Issue with P1 data:
% if strcmp(participant_id, '1')
%     % For participant 1: after removing bad segments, skips from 73 to 108
%     num_trials = EEG.event(end-3).trial_id - 1;
% end
% 
% % Issue with P3 data:
% if strcmp(participant_id, '3')
%     % For participant 3: the last trial (87) is not complete
%     num_trials = EEG.event(end).trial_id - 1;
% end

% Creating Final Structure:
spectrum = [];
blockLength = 3;
baseline_trial = 0;

switch segment_of_interest
    case 'detection_time'
    
    case 'trialPOL'
        for trial = 1:num_trials
            % Resetting values to identify if missing tag:
            start_idx = -1;
            end_idx = -1;
         
            % Retrieve indices of trial
            bl = ceil(trial/blockLength);
            tr_inBl = trial - (bl-1)*blockLength;
            trial_idx = find([EEG.event.BlockIndex] == bl & [EEG.event.TrialIndex] == tr_inBl);
            
            if trial_idx ~= 2
                error('This trial has ' + trial_idx + " events")
            end

            % Retrieving information for metadata
            field_of_view = str2double(regexp(EEG.event(trial_idx(end)).TrialType, '\d+', 'match'));

            % Find start and end idx (TrialStart, TrialEnd) for this trial
            for i = 1:length(trial_idx)
                if ismember(EEG.event(trial_idx(i)).TrialType, "Baseline")
                    baseline_trial = 1;
                    break
                else
                    % Retrieving information for metadata
                    field_of_view = str2double(regexp(EEG.event(trial_idx(end)).TrialType, '\d+', 'match'));

                    if (ismember(EEG.event(trial_idx(i)).type, begin_segment))
                        start_idx = trial_idx(i);
                    elseif ismember(EEG.event(trial_idx(i)).type, end_segment)
                        end_idx = trial_idx(i);
                    else
                        error("Didn't find any start or end event for this trial")
                    end
                end
            end

            % if baseline trial, go to next iteration
            if baseline_trial
                break
            end


            % Alexandre's nicer code: 
%             starts = contains({EEG.event(:).type}, 'Start training') | contains({EEG.event(:).type}, 'Start trial');
%             start_trial = starts & EEG.event(:).trial_id == trial;
%             
%             if isempty(start_trial)
%                 % There are no tags indicating the beginning of a trial,
%                 usng the tag of the end of trial of previous one
%             end
            
            % for i = 1:length(trial_idx)
            %     % Looping to find start and end 
            %     if (ismember(EEG.event(trial_idx(i)).type, begin_segment)) && (start_idx == -1)
            %         start_idx = trial_idx(i);
            %     elseif ismember(EEG.event(trial_idx(i)).type, end_segment) && (end_idx == -1)
            %         end_idx = trial_idx(i);
            %     else
            %         continue
            %     end
            % end

            % % Check if tags are missing or not: if missing than take the tag that
            % % precedes/follows but note it in metadata
            % if (start_idx == -1)
            %     % Get the index of the tag that precedes (ie end baseline)
            %     for i = 1:length(trial_idx)
            %         if ismember(EEG.event(trial_idx(i)).type, {' End baseline'})
            %             start_idx = trial_idx(i);
            %         end
            %     end
            % elseif (end_idx == -1)
            %     % Button was not pressed: so selecting end of trial as end
            %     % of segment of interest
            %     for i = 1:length(trial_idx)
            %         if ismember(EEG.event(trial_idx(i)).type, {' End training';' End trial'})
            %             end_idx = trial_idx(i);
            %         end
            %     end
            % end
            
            % Filling up the table with segment of interest and metadata
            % 1. Find time stamps to retrieve data
            start_segment_time_index = round(EEG.event(start_idx).latency);
            end_segment_time_index = round(EEG.event(end_idx).latency);
            
            data = EEG.data(:,start_segment_time_index:end_segment_time_index);
            
            % 2. Retrieve MetaData          
            EEG_trial_data.metaInfo(current_line_trial + trial).participant_id = participant_id;
            EEG_trial_data.metaInfo(current_line_trial + trial).TrialIndex = trial;
            EEG_trial_data.metaInfo(current_line_trial + trial).field_of_view = field_of_view;
            % EEG_trial_data.metaInfo(current_line_trial + trial).auditory_cue = auditory_cue;
            % EEG_trial_data.metaInfo(current_line_trial + trial).auditory_answer = auditory_answer;
            

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
                spectopo(data, size(data,2), EEG.srate,...
                'freq', [10.0, 20.0], ... %freq_of_interest, ... %10,...                %'plotchan', [],...
                'chanlocs', EEG.chanlocs,...
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
        %EEG_trial_data.(['P' participant_id]).freqs = freqs;
        EEG_trial_data.(['P' participant_id]).freqs = freqs;
        disp('Converting spectrum into linear scale for each TRIAL');
        %EEG_trial_data.(['P' participant_id]).spectrum = 10.^(spectrum./10);
        
    otherwise
        error('Case not yet considered, please code it')
end

end

