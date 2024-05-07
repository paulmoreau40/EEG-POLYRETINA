% bemobil_process_all_mobilab - wrapper function that incorporates all necessary processing steps from raw .xdf to .set
% files in EEGLAB. Data is being loaded into mobilab, rigidbody mocap streams are processed (filtered, transformed into
% euler angles, derived) and data and marker streams are then exported to EEGLAB containing EEG and all other kinds of
% data. The dataset has the suffix '_MoBI'. This dataset is then split into individual channel types (e.g. 'EEG',
% 'MOCAP', 'EYE', 'OTHER'), and subsequently all EEG files (from several raw .xdf files) will be merged into one large
% EEG file for this participant, which can then be used for further processing (e.g. with bemobil_process_all_AMICA)
%
% The intermediate files are stored on the disk.
%
% Usage:
%   >>  [ALLEEG, EEG_merged, CURRENTSET] = bemobil_process_all_mobilab(subject, bemobil_config, ALLEEG, CURRENTSET, mobilab)
%
% Inputs:
%   subject                   - subject number of the current subject (necessary for filepaths and storage)
%   bemobil_config            - configuration struct with all necessary information. See EEG_processing_example file
%                                that comes with this function!
%   ALLEEG                    - complete EEGLAB data set structure
%   CURRENTSET                - index of current EEGLAB EEG structure within ALLEEG
%	mobilab					  - container for the mobilab application. execute "runmobilab" before this script to get it
%
% Outputs:
%   ALLEEG                    - complete EEGLAB data set structure
%   EEG_merged				  - merged EEGLAB EEG structure, contains EEG datasets of all conditions.
%   CURRENTSET                - index of current EEGLAB EEG structure within ALLEEG
%
%   .set data file of EEGLAB EEG structures are stored on disk according to their names in the bemobil_config
%
% See also:
%   EEGLAB
%
% Authors: Marius Klug, 2019

function [ALLEEG, EEG_merged, CURRENTSET] = xdf2set(ALLEEG, CURRENTSET, subject_ind, study_config, overwrite)

user = study_config.user;
subject = study_config.subjects(subject_ind).id;
%disp(['Subject ' num2str(subject)]);

input_filepath = [study_config.study_folder study_config.raw_data_folder subject];
output_filepath = [study_config.study_folder study_config.raw_EEGLAB_data_folder subject];

file_prefix{1} = [subject '_AllStreams_'];

for i_filename = 1:length(study_config.filenames)
    %i_filename = 4;
    full_filename = ['sub-', subject, '_', study_config.filenames{i_filename}];
    
    %continue
    
    %% import
    input_file = fullfile(input_filepath, [full_filename '.xdf']);
    if isfile(input_file)
        if isfile([output_filepath, filesep, file_prefix{1}, study_config.filenames{i_filename}, '.set']) && ~overwrite
            disp([file_prefix{1}, study_config.filenames{i_filename}, '.set already exists : not importing XDF file']);
            continue
        else
            disp(['Importing file: "' full_filename '.xdf" ...']);
            % if any(strcmp({'P007_block001_2','P036_block002_1','P038_block001_3'}, full_filename))
            %     [AllStreams, ~] = load_xdf_AD(input_file, 'JitterBreakThresholdSeconds', 0.1);
            % elseif any(strcmp({'P007_block001_1'}, full_filename))
            %     [AllStreams, ~] = load_xdf_AD(input_file, 'JitterBreakThresholdSeconds', 0.05);                
            % else
            %     [AllStreams, ~] = load_xdf_AD(input_file);
            %     %[AllStreams, ~] = load_xdf_AD(input_file, 'Verbose', true);
            %     %[AllStreams, ~] = load_xdf_AD(input_file, 'Verbose', true, 'HandleClockResets', false);
            %     %[AllStreams, ~] = load_xdf_AD(input_file, 'Verbose', true, 'HandleJitterRemoval', false);
            % end
            [AllStreams, ~] = load_xdf_AD(input_file); % paul
        end
    else
        continue
    end
    
    % Detect Stream types
    EEG_stream_inds = zeros(1,study_config.stream_count(1));
    ET_stream_inds = zeros(1,study_config.stream_count(2));
    MOCAP_stream_inds = zeros(1,study_config.stream_count(3));
    Event_stream_inds = zeros(1,study_config.stream_count(4));
    for s = 1:numel(AllStreams)
        stream_type = AllStreams{s}.info.type;
        if (sum(contains(study_config.eeg_streams,stream_type))>0) % EEG
            ind2replace = find(contains(study_config.eeg_streams_order, AllStreams{s}.info.name));
            if EEG_stream_inds(ind2replace) == 0
                EEG_stream_inds(ind2replace)=s;
            else
                % Catch if the same stream is present more than once
                fprintf('Old stream length: %d sample(s).\n', size(AllStreams{EEG_stream_inds(ind2replace)}.time_series,2))
                fprintf('New stream length: %d sample(s).\n', size(AllStreams{s}.time_series,2))
                kept = input('Which one do I keep? [0=old, 1=new]: '); % enter 0 or 1
                if kept
                    EEG_stream_inds(ind2replace)=s;
                end
            end
        elseif (sum(contains(study_config.eye_tracker_streams,stream_type))>0) % ET
            curr_ind = find(ET_stream_inds,1,'last');
            if isempty(curr_ind)
                ET_stream_inds(1)=s;
            else
                ET_stream_inds(curr_ind+1)=s;
            end
        elseif (sum(contains(study_config.rigidbody_streams,stream_type))>0) % MOCAP
            curr_ind = find(MOCAP_stream_inds,1,'last');
            if isempty(curr_ind)
                MOCAP_stream_inds(1)=s;
            else
                MOCAP_stream_inds(curr_ind+1)=s;
            end
        elseif (sum(contains(study_config.event_streams,stream_type))>0) % EVENTS (Tag)
            % OLD code of JB to deal with weird Markers
            % if sum(contains(study_config.eeg_streams_order,AllStreams{s}.info.name))==0
               % To avoid loading EEG Marker streams (useless in ANT)
            
               % PAUL: new code to only keep the good events Tag(remove
               % 1st, 4th, 6th
            if sum(contains({'ControlCenter-Vision', 'ControlCenter-Tags','Sim-VisionConfirmation'},AllStreams{s}.info.name))==0
                curr_ind = find(Event_stream_inds,1,'last');
                if isempty(curr_ind)
                    Event_stream_inds(1)=s;
                else
                    Event_stream_inds(curr_ind+1)=s;
                end
            end
        elseif (sum(contains(study_config.bad_streams,stream_type))>0) % if type 'Markers' -> Bad for Paul
            continue
            
        else
            error(['Could not recognize the stream type, please make sure all stream types fall into',...
                ' a category defined in the config file']);
        end
    end
    
    %% Complete Missing data from eego file (if necessary & possible)
    if (sum(strcmp(study_config.filenames{i_filename}, study_config.subjects(subject_ind).missingData))>0)
        AllStreams = completeEEGfromEEGO_EEGAFF(AllStreams, EEG_stream_inds, study_config.filenames{i_filename}, study_config.subjects,...
            input_filepath, study_config.figures_folder);
    end
    
    %% Read Streams
    Data_all = eeg_emptyset;
    
    %% Get info necessary to create the data:
    %% EEG:
    n_chans_eeg = 0;
    t_min = +inf;
    t_max = -inf;
    srate_eeg = -1;
    segments2keep_eeg = cell(length(EEG_stream_inds),1);
    for ee = 1:length(EEG_stream_inds)
        EEG_stream = AllStreams{EEG_stream_inds(ee)};
        n_chans_eeg = n_chans_eeg + size(EEG_stream.time_series,1);
        
        srate_eeg = str2num(EEG_stream.info.nominal_srate);
        
        % Look at Segments to understand if there is weird data
        for seg=1:numel(EEG_stream.segments)
            if EEG_stream.segments(seg).duration == 0
                % Invalid segment
            else
                % Valid segment
                last_valid_seg = segments2keep_eeg{ee}(find(segments2keep_eeg{ee},1,'last'));
                if ~isempty(last_valid_seg)
                    % A previous valid segment has been found
                    if EEG_stream.segments(seg).t_begin < EEG_stream.segments(last_valid_seg).t_end
                        % An overlap exists --> Remove the overlap period
                        if EEG_stream.segments(seg).t_begin <= EEG_stream.segments(last_valid_seg).t_begin
                            % Remove the whole last valid seg
                            segments2keep_eeg{ee} = [segments2keep_eeg{ee}(1:end-1), seg];
                            warning('Overlap: removed a full segment of %.2f seconds from EEG data',...
                                EEG_stream.segments(last_valid_seg).duration);
                            
                            t_min = min(t_min,EEG_stream.segments(seg).t_begin);
                            t_max = max(t_max,EEG_stream.segments(seg).t_end);
                        elseif EEG_stream.segments(seg).t_end <= EEG_stream.segments(last_valid_seg).t_end
                            % Do not add the current seg
                            warning('Overlap: removed a full segment of %.2f seconds from EEG data',...
                                EEG_stream.segments(seg).duration);
                        else
                            t_discard_begin = EEG_stream.segments(seg).t_begin;
                            t_discard_stop = EEG_stream.segments(last_valid_seg).t_end;
                            
                            % Modify last_valid_seg:
                            EEG_stream.segments(last_valid_seg).t_end = t_discard_begin;
                            EEG_stream.segments(last_valid_seg).duration = t_discard_begin - EEG_stream.segments(last_valid_seg).t_begin;
                            EEG_stream.segments(last_valid_seg).index_range(2) = ...
                                find(EEG_stream.time_stamps(EEG_stream.segments(last_valid_seg).index_range(1):...
                                EEG_stream.segments(last_valid_seg).index_range(2)) <= t_discard_begin, 1, 'last')...
                                + EEG_stream.segments(last_valid_seg).index_range(1) - 1;
                            EEG_stream.segments(last_valid_seg).num_samples = diff(EEG_stream.segments(last_valid_seg).index_range)+1;
                            EEG_stream.segments(last_valid_seg).effective_srate = (EEG_stream.segments(last_valid_seg).num_samples-1)...
                                /EEG_stream.segments(last_valid_seg).duration;
                            
                            % Modify seg:
                            EEG_stream.segments(seg).t_begin = t_discard_stop;
                            EEG_stream.segments(seg).duration = EEG_stream.segments(seg).t_end - t_discard_stop;
                            EEG_stream.segments(seg).index_range(1) = ...
                                find(EEG_stream.time_stamps(EEG_stream.segments(seg).index_range(1):...
                                EEG_stream.segments(seg).index_range(2)) >= t_discard_stop, 1)...
                                + EEG_stream.segments(seg).index_range(1) - 1;
                            EEG_stream.segments(seg).num_samples = diff(EEG_stream.segments(seg).index_range)+1;
                            EEG_stream.segments(seg).effective_srate = (EEG_stream.segments(seg).num_samples-1)...
                                /EEG_stream.segments(seg).duration;
                            
                            segments2keep_eeg{ee} = [segments2keep_eeg{ee}, seg];
                            warning('Overlap: removed a partial segment of %.2f seconds from EEG data',...
                                t_discard_stop - t_discard_begin);
                            
                            t_min = min(t_min,EEG_stream.segments(seg).t_begin);
                            t_max = max(t_max,EEG_stream.segments(seg).t_end);
                        end
                        
                        %                         % An overlap exists --> Remove the shortest segment
                        %                         if (EEG_stream.segments(seg).duration > EEG_stream.segments(last_valid_seg).duration)
                        %                             segments2keep_eeg{ee} = [segments2keep_eeg{ee}, seg];
                        %                             warning('Overlap: removed a segment of %.2f seconds from EEG data',...
                        %                                 EEG_stream.segments(last_valid_seg).duration);
                        %
                        %                             t_min = min(t_min,EEG_stream.segments(seg).t_begin);
                        %                             t_max = max(t_max,EEG_stream.segments(seg).t_end);
                        %                         elseif (EEG_stream.segments(seg).duration < EEG_stream.segments(last_valid_seg).duration)
                        %                             segments2keep_eeg{ee} = segments2keep_eeg{ee}(1:end-1);
                        %                             warning('Overlap: removed a segment of %.2f seconds from EEG data',...
                        %                                 EEG_stream.segments(seg).duration);
                        %                         else
                        %                             error('Cannot decide which segment to remove, they have equal lengths')
                        %                         end
                    else
                        segments2keep_eeg{ee} = [segments2keep_eeg{ee}, seg];
                        t_min = min(t_min,EEG_stream.segments(seg).t_begin);
                        t_max = max(t_max,EEG_stream.segments(seg).t_end);
                    end
                else
                    segments2keep_eeg{ee} = [segments2keep_eeg{ee}, seg];
                    t_min = min(t_min,EEG_stream.segments(seg).t_begin);
                    t_max = max(t_max,EEG_stream.segments(seg).t_end);
                end
            end
        end
        
        %{
        % First check that nominal srate and effective srate are not too far appart:
        nom_srate = str2num(EEG_stream.info.nominal_srate);
        eff_srate = EEG_stream.info.effective_srate;
        if (eff_srate < nom_srate*0.95 || eff_srate > nom_srate*1.05)
            error('Effective srate too far from nominal srate')
        end
        
        % Then check that nominal srate are the same for all EEG recordings
        if (srate_eeg==-1)
            srate_eeg = nom_srate;
        elseif (srate_eeg ~= nom_srate)
            error('FATAL ERROR: Not the same srate between EEG recordings')
        end
        %}
        
        AllStreams{EEG_stream_inds(ee)} = EEG_stream;
    end
    
    % Remove channels if indicated:
    if ~isempty(study_config.channels_to_remove)
        n_chans_eeg = n_chans_eeg - sum(sum(~isnan(study_config.channels_to_remove)));
    end
    
    %% ET:
    n_chans_et = 0;
    srate_et = -1;
    if ~isempty(ET_stream_inds) % Only one stream
        t_min = min(t_min, AllStreams{ET_stream_inds(1)}.time_stamps(1));
        t_max = max(t_max, AllStreams{ET_stream_inds(1)}.time_stamps(end));
        
        % Two channels are redundant and contain time data
        n_chans_et = str2num(AllStreams{ET_stream_inds(1)}.info.channel_count) - 2;
        
        % Check that nominal srate and effective srate are not too far appart:
        nom_srate = str2num(AllStreams{ET_stream_inds(1)}.info.nominal_srate);
        if (nom_srate > 0)
            eff_srate = AllStreams{ET_stream_inds(1)}.info.effective_srate;
            if (eff_srate < nom_srate*0.95 || eff_srate > nom_srate*1.05)
                warning('Eye tracking data: Effective srate (%.2f Hz) far from nominal srate (%.2f Hz)', eff_srate, nom_srate)
            end
            srate_et = round(eff_srate);
        else
            % No sampling rate defined for ET in this experiment
            srate_et = nom_srate;
        end
    end
    
    %% MOCAP:
    n_chans_mocap = 0;
    srate_mocap = -1;
    if ~isempty(MOCAP_stream_inds)
        for moc = MOCAP_stream_inds
            t_min = min(t_min, AllStreams{moc}.time_stamps(1));
            t_max = max(t_max, AllStreams{moc}.time_stamps(end));
            
            n_chans_mocap = n_chans_mocap + str2num(AllStreams{moc}.info.channel_count);
        end
        % First check that nominal srate and effective srate are not too far appart:
        nom_srate = str2num(AllStreams{moc}.info.nominal_srate);
        if nom_srate > 0
            eff_srate = AllStreams{moc}.info.effective_srate;
            if (eff_srate < nom_srate*0.95 || eff_srate > nom_srate*1.05)
                error('Effective srate too far from nominal srate')
            end
            
            % Then check that nominal srate are the same for all Mocap recordings
            if (srate_mocap==-1)
                srate_mocap = nom_srate;
            elseif (srate_mocap ~= nom_srate)
                error('FATAL ERROR: Not the same srate between Mocap recordings')
            end
        else
            srate_mocap = 0;
            n_chans_mocap = n_chans_mocap+1; % Add a channel to indicate when this is an orginal sample
        end
    end
    
    interval = 1/srate_eeg; % Assuming eeg as the highest sampling rate
    times = round(t_min,4):interval:round(t_max,4)+interval;
    




    for typ = 1:length(study_config.stream_count)
        switch typ
            case 1 %EEG
                
                disp('Importing EEG...')
                %% Data
                data_eeg = nan(n_chans_eeg,length(times));
                for ee = 1:length(EEG_stream_inds)

                    % (DE)COMMENT TO IGNORE THE EEG LOADING, time consuming
                    continue

                    EEG_stream = AllStreams{EEG_stream_inds(ee)};
                    
                    % To compute channels span:
                    local_data = EEG_stream.time_series;
                    if ~isempty(study_config.channels_to_remove)
                        chans2keep = setdiff(1:size(local_data,1),study_config.channels_to_remove(ee,:));
                        local_data = local_data(chans2keep,:);
                    end
                    
                    if ee == 1
                        channels_span = 1:size(local_data,1);
                    else
                        next_chan = channels_span(end)+1;
                        channels_span = next_chan:next_chan+size(local_data,1)-1;
                    end
                    
                    for seg = segments2keep_eeg{ee}
                        eff_srate = EEG_stream.segments(seg).effective_srate;
                        time_span = EEG_stream.segments(seg).index_range(1):EEG_stream.segments(seg).index_range(2);
                        local_data = EEG_stream.time_series(chans2keep, time_span);
                        local_times = EEG_stream.time_stamps(time_span);
                        
                        local2global = zeros(1,length(local_times));
                        if strcmpi(user, 'ila')
                            fprintf('Finding time samples for EEG%d;Seg%d\n', ee, seg)
                            parfor t = 1:length(local_times)
                                ind = find(times>=local_times(t),1);
                                if isempty(ind)
                                    if t == length(local_times) && local_times(t)<times(end)+interval/2
                                        local2global(t) = length(times);
                                    else
                                        error('Seg %d: Sample %d could not be placed',seg,t);
                                    end
                                elseif ind==1
                                    if t == 1 && local_times(t)>times(1)-interval/2
                                        local2global(t) = 1;
                                    else
                                        error('Seg %d: Sample %d could not be placed',seg,t);
                                    end
                                else
                                    if round(local_times(t)-times(ind-1),4)<=interval/2
                                        local2global(t) = ind-1;
                                    elseif round(times(ind)-local_times(t),4)<=interval/2
                                        local2global(t) = ind;
                                    else
                                        error('Seg %d: Sample %d could not be placed',seg,t);
                                    end
                                end
                                
                                %                             ind = find(round(local_times(t)-interval/2,4)<=times & times<round(local_times(t)+interval/2,4),1);
                                %                             if isempty(ind)
                                %                                 if t == length(local_times) && local_times(t)<times(end)+interval
                                %                                     local2global(t) = length(times);
                                %                                 elseif ~isempty(find(times == round(local_times(t)+interval/2,4),1))
                                %                                     local2global(t) = find(times == round(local_times(t)+interval/2,4),1);
                                %                                 else
                                %                                     error('Seg %d: Sample %d could not be placed',seg,t);
                                %                                 end
                                %                             elseif length(ind)==1
                                %                                 local2global(t) = ind;
                                %                             end
                            end
                        else
                            ppm = ParforProgressbar(length(local_times), 'showWorkerProgress', true,...
                                'title', sprintf('Finding time samples for EEG%d;Seg%d', ee, seg));
                            parfor t = 1:length(local_times)
                                ind = find(times>=local_times(t),1);
                                if isempty(ind)
                                    if t == length(local_times) && local_times(t)<times(end)+interval/2
                                        local2global(t) = length(times);
                                    else
                                        error('Seg %d: Sample %d could not be placed',seg,t);
                                    end
                                elseif ind==1
                                    if t == 1 && local_times(t)>times(1)-interval/2
                                        local2global(t) = 1;
                                    else
                                        error('Seg %d: Sample %d could not be placed',seg,t);
                                    end
                                else
                                    if round(local_times(t)-times(ind-1),4)<=interval/2
                                        local2global(t) = ind-1;
                                    elseif round(times(ind)-local_times(t),4)<=interval/2
                                        local2global(t) = ind;
                                    else
                                        error('Seg %d: Sample %d could not be placed',seg,t);
                                    end
                                end
                                
                                %                             ind = find(round(local_times(t)-interval/2,4)<=times & times<round(local_times(t)+interval/2,4),1);
                                %                             if isempty(ind)
                                %                                 if t == length(local_times) && local_times(t)<times(end)+interval
                                %                                     local2global(t) = length(times);
                                %                                 elseif ~isempty(find(times == round(local_times(t)+interval/2,4),1))
                                %                                     local2global(t) = find(times == round(local_times(t)+interval/2,4),1);
                                %                                 else
                                %                                     error('Seg %d: Sample %d could not be placed',seg,t);
                                %                                 end
                                %                             elseif length(ind)==1
                                %                                 local2global(t) = ind;
                                %                             end
                                ppm.increment();
                            end
                            delete(ppm);
                        end                        
                        
                        if (length(local2global) ~= length(unique(local2global)))
                            % Only dealing with contiguous repetitions for now
                            rep_list = find(local2global == [0,local2global(1:end-1)]);
                            %[C, ia, ic] = unique(local2global);
                            %reps = setdiff(1:length(local2global),ia);
                            for r = length(rep_list):-1:1
                                if rep_list(r) > 2 && (eff_srate >= srate_eeg*0.95 && eff_srate <= srate_eeg*1.05)...
                                        && local2global(rep_list(r)-2)==local2global(rep_list(r))-2
                                    % Correct local2global knowing samples
                                    % should be regularly spaced given the
                                    % srate
                                    local2global(rep_list(r)-1)=local2global(rep_list(r))-1;
                                else
                                    % Remove the duplicate sample and
                                    % average the values in the data
                                    local2global = [local2global(1:rep_list(r)-1),local2global(rep_list(r)+1:end)];
                                    local_times = [local_times(1:rep_list(r)-1),local_times(rep_list(r)+1:end)];
                                    local_data = [local_data(:,1:rep_list(r)-2),...
                                        mean(local_data(:,rep_list(r)-1:rep_list(r)),2),...
                                        local_data(:,rep_list(r)+1:end)];
                                end
                            end
                            
                            if (length(local2global) ~= length(unique(local2global)))
                                error('Still repetitions to deal with')
                            end
                        end
                        
                        % Check that eeg srate and effective srate are not too far appart:
                        srate_scale = round(srate_eeg/eff_srate);
                        if (eff_srate*srate_scale < srate_eeg*0.95 || eff_srate*srate_scale > srate_eeg*1.05)
                            warning('Impossible to combine the segment srate with global srate in a satisfying way')
                            
                            % portion to visualize what's happening
                            instant_srate = 1./diff(EEG_stream.time_stamps(EEG_stream.segments(seg).index_range(1):...
                                EEG_stream.segments(seg).index_range(2)));
                            histogram(instant_srate);
                            saveCurrentFig([output_filepath filesep],...
                                sprintf('%s_SrateHistogram_%s_EEG%dSeg%d',subject,full_filename, ee, seg), {'fig'},[]);
                            
                            data_eeg(channels_span,local2global)=local_data;
                        else
                            % We can consider the srate are multiples
                            fprintf('EEG%d;Seg%d: nomSR= %d, effSR= %.2f\n', ee, seg, srate_eeg, eff_srate)
                            for j = 1:srate_scale
                                data_eeg(channels_span,local2global(1:end-j+1)+j-1)=local_data(:,1:end-j+1);
                            end
                        end
                    end
                end
                
                %% chanlocs
                chanlocs_eeg = struct();
                if ~isempty(study_config.channel_locations_filename)
                    chanlocs_eeg = readlocs(fullfile(input_filepath, study_config.channel_locations_filename));
                else
                    % TODO
                    % Read from stream.info.desc.channels
                    % Get inspired from eeg_load_xdf
                end
                
                fields = fieldnames(chanlocs_eeg);
                command = '';
                for f = 1:(numel(fields)+2)
                    if f==1
                        command = [command,'''type'',','cell(1, n_chans_eeg+n_chans_et+n_chans_mocap),'];
                    elseif f==2
                        command = [command,'''unit'',[],'];
                    else
                        command = [command,'''',fields{f-2},''',[],'];
                    end
                end
                command = command(1:end-1);
                eval(['chanlocs = struct(',command,');']);
                for ch = 1:n_chans_eeg
                    if ~isempty(study_config.eog_channels)
                        if sum(strcmp(study_config.eog_channels, chanlocs_eeg(ch).labels),1)>0
                            chanlocs(ch).type = 'EOG';
                        else
                            chanlocs(ch).type = 'EEG';
                        end
                    else
                        chanlocs(ch).type = 'EEG';
                    end
                    chanlocs(ch).unit = 'microVolt';
                    for f = 1:numel(fields)
                        chanlocs(ch).(fields{f}) = chanlocs_eeg(ch).(fields{f});
                    end
                end
                
                %% etc
                for ee = 1:length(EEG_stream_inds)
                    Data_all.etc.(['descEEG', num2str(ee)]) = AllStreams{EEG_stream_inds(ee)}.info.desc;
                    Data_all.etc.(['infoEEG', num2str(ee)]) = rmfield(AllStreams{EEG_stream_inds(ee)}.info,'desc');
                    Data_all.etc.(['infoEEG', num2str(ee)]).('segments') = AllStreams{EEG_stream_inds(ee)}.segments;
                    Data_all.etc.(['infoEEG', num2str(ee)]).('keptSegments') = segments2keep_eeg{ee};
                end
                disp ('...done')
                
            case 2 %ET
                %continue
                if ~isempty(ET_stream_inds) % Only one stream
                    disp('Importing ET...')
                    %% Data
                    local_data = AllStreams{ET_stream_inds(1)}.time_series(1:n_chans_et,:);
                    local_times = AllStreams{ET_stream_inds(1)}.time_stamps;
                    if (srate_et > 0)
                        eff_srate = AllStreams{ET_stream_inds(1)}.info.effective_srate;
                    end
                    
                    local2global = zeros(1,length(local_times));
                    if strcmpi(user, 'ila')
                        disp('Finding time samples for ET');
                        parfor t = 1:length(local_times)
                            ind = find(times>=local_times(t),1);
                            if isempty(ind)
                                if t == length(local_times) && local_times(t)<times(end)+interval/2
                                    local2global(t) = length(times);
                                else
                                    error('Sample %d could not be placed',t);
                                end
                            elseif ind==1
                                if t == 1 && local_times(t)>times(1)-interval/2
                                    local2global(t) = 1;
                                else
                                    error('Sample %d could not be placed',t);
                                end
                            else
                                if round(local_times(t)-times(ind-1),4)<=interval/2
                                    local2global(t) = ind-1;
                                elseif round(times(ind)-local_times(t),4)<=interval/2
                                    local2global(t) = ind;
                                else
                                    error('Sample %d could not be placed',t);
                                end
                            end
                            
                            %                         ind = find(local_times(t)-interval/2<=times & times<local_times(t)+interval/2,1);
                            %                         if isempty(ind)
                            %                             if t == length(local_times) && local_times(t)<times(end)+interval
                            %                                 local2global(t) = length(times);
                            %                             else
                            %                                 error('Sample %d could not be placed',t);
                            %                             end
                            %                         elseif length(ind)==1
                            %                             local2global(t) = ind;
                            %                         end
                        end
                    else
                        ppm = ParforProgressbar(length(local_times), 'showWorkerProgress', true,...
                            'title', 'Finding time samples for ET');
                        parfor t = 1:length(local_times)
                            ind = find(times>=local_times(t),1);
                            if isempty(ind)
                                if t == length(local_times) && local_times(t)<times(end)+interval/2
                                    local2global(t) = length(times);
                                else
                                    error('Sample %d could not be placed',t);
                                end
                            elseif ind==1
                                if t == 1 && local_times(t)>times(1)-interval/2
                                    local2global(t) = 1;
                                else
                                    error('Sample %d could not be placed',t);
                                end
                            else
                                if round(local_times(t)-times(ind-1),4)<=interval/2
                                    local2global(t) = ind-1;
                                elseif round(times(ind)-local_times(t),4)<=interval/2
                                    local2global(t) = ind;
                                else
                                    error('Sample %d could not be placed',t);
                                end
                            end
                            
                            %                         ind = find(local_times(t)-interval/2<=times & times<local_times(t)+interval/2,1);
                            %                         if isempty(ind)
                            %                             if t == length(local_times) && local_times(t)<times(end)+interval
                            %                                 local2global(t) = length(times);
                            %                             else
                            %                                 error('Sample %d could not be placed',t);
                            %                             end
                            %                         elseif length(ind)==1
                            %                             local2global(t) = ind;
                            %                         end
                            ppm.increment();
                        end
                        delete(ppm);
                    end
                    
                    
                    if (length(local2global) ~= length(unique(local2global)))
                        % Only dealing with contiguous repetitions for now
                        rep_list = find(local2global == [0,local2global(1:end-1)]);
                        %[C, ia, ic] = unique(local2global);
                        %reps = setdiff(1:length(local2global),ia);
                        for r = length(rep_list):-1:1
                            if srate_et > 0 && rep_list(r) > 2 && (eff_srate >= srate_et*0.95 && eff_srate <= srate_et*1.05)...
                                    && local2global(rep_list(r)-2)==local2global(rep_list(r))-2
                                % Correct local2global knowing samples
                                % should be regularly spaced given the
                                % srate
                                local2global(rep_list(r)-1)=local2global(rep_list(r))-1;
                            else
                                % Remove the duplicate sample and
                                % average the values in the data
                                local2global = [local2global(1:rep_list(r)-1),local2global(rep_list(r)+1:end)];
                                local_times = [local_times(1:rep_list(r)-1),local_times(rep_list(r)+1:end)];
                                local_data = [local_data(:,1:rep_list(r)-2),...
                                    mean(local_data(:,rep_list(r)-1:rep_list(r)),2),...
                                    local_data(:,rep_list(r)+1:end)];
                            end
                        end
                        
                        if (length(local2global) ~= length(unique(local2global)))
                            error('Still repetitions to deal with')
                        end
                    end
                    
                    if srate_et > 0
                        data_et = nan(n_chans_et,length(times));
                        srate_scale = srate_eeg/srate_et;
                        if floor(srate_scale) == srate_scale
                            % EEG srate is a multiple of ET srate
                            % 'Stairs' completing (To upsample, the value is copied 'srate_scale' times until the next new sample)
                            for j = 1:srate_scale
                                data_et(:,local2global(1:end-j+1)+j-1)=local_data(:,1:end-j+1);
                            end
                        else
                            % TO DO
                        end
                    else
                        % Create an additional data channel to indicate
                        % whether the given sample is original or interpolated
                        data_et = nan(n_chans_et+1,length(times));
                        data_et(end,:) = zeros(1,length(times));
                        data_et(end,local2global) = ones(1,length(local2global));
                        
                        % Interpolation
                        %data_et(1:end-1,:) = interp1(times(local2global),local_data',times,'linear')';
                        % No interpolation
                        data_et(1:end-1,local2global)=local_data;
                    end
                    
                    %% chanlocs
                    for ch = 1:n_chans_et
                        chanlocs(ch+n_chans_eeg).type = 'ET';
                        channel_desc = AllStreams{ET_stream_inds(1)}.info.desc.channels.channel{ch};
                        
                        if isfield(channel_desc, 'unit')
                            data_unit = AllStreams{ET_stream_inds(1)}.info.desc.channels.channel{ch}.unit;
                        else
                            data_unit = '';
                        end
                        
                        if ~strcmp(data_unit, '')
                            if contains(lower(channel_desc.label),'validity') || ...
                                    contains(lower(channel_desc.label),'blink')
                                chanlocs(ch+n_chans_eeg).unit = 'bool';
                            else
                                chanlocs(ch+n_chans_eeg).unit = data_unit;
                            end
                        else
                            if isfield(channel_desc, 'type')
                                data_type = AllStreams{ET_stream_inds(1)}.info.desc.channels.channel{ch}.type;
                            else
                                data_type = '';
                            end
                            if strcmp(data_type, 'position')
                                chanlocs(ch+n_chans_eeg).unit = 'pixel';
                            elseif strcmp(data_type, 'area')
                                chanlocs(ch+n_chans_eeg).unit = 'a.u.';
                            elseif strcmp(data_type, 'ppd')
                                chanlocs(ch+n_chans_eeg).unit = 'pixelsPerDegree';
                            end
                        end
                        chanlocs(ch+n_chans_eeg).labels = channel_desc.label;
                    end
                    
                    if ~(srate_et > 0)
                        % Additional channel definition
                        n_chans_et = n_chans_et+1;
                        chanlocs(n_chans_et+n_chans_eeg).type = 'ET';
                        chanlocs(n_chans_et+n_chans_eeg).labels = 'Original_samples_eyeTracker';
                        chanlocs(n_chans_et+n_chans_eeg).unit = 'bool';
                    end
                    
                    %% etc
                    Data_all.etc.('descET') = AllStreams{ET_stream_inds(1)}.info.desc;
                    Data_all.etc.('infoET') = rmfield(AllStreams{ET_stream_inds(1)}.info,'desc');
                    disp ('...done')
                end
                
            case 3 %MOCAP
                %continue
                if ~isempty(MOCAP_stream_inds)
                    disp('Importing MOCAP...')
                    %% Data
                    data_mocap = nan(n_chans_mocap,length(times));
                    for moc = MOCAP_stream_inds
                        local_data = AllStreams{moc}.time_series;
                        local_times = AllStreams{moc}.time_stamps;
                        
                        n_chans = size(local_data,1);
                        if srate_mocap == 0
                            n_chans = n_chans+1;
                        end
                        
                        if moc == MOCAP_stream_inds(1)
                            channels_span = 1:n_chans;
                        else
                            next_chan = channels_span(end)+1;
                            channels_span = next_chan:next_chan+n_chans-1;
                        end
                        
                        local2global = zeros(1,length(local_times));
                        if strcmpi(user, 'ila')
                            disp('Finding time samples for MOCAP');
                            parfor t = 1:length(local_times)
                                ind = find(times>=local_times(t),1);
                                if isempty(ind)
                                    if t == length(local_times) && local_times(t)<times(end)+interval/2
                                        local2global(t) = length(times);
                                    else
                                        error('Sample %d could not be placed',t);
                                    end
                                elseif ind==1
                                    if t == 1 && local_times(t)>times(1)-interval/2
                                        local2global(t) = 1;
                                    else
                                        error('Sample %d could not be placed',t);
                                    end
                                else
                                    if round(local_times(t)-times(ind-1),4)<=interval/2
                                        local2global(t) = ind-1;
                                    elseif round(times(ind)-local_times(t),4)<=interval/2
                                        local2global(t) = ind;
                                    else
                                        error('Sample %d could not be placed',t);
                                    end
                                end
                            end
                        else
                            ppm = ParforProgressbar(length(local_times), 'showWorkerProgress', true,...
                                'title', ['Finding time samples for MOCAP', num2str(moc)]);
                            parfor t = 1:length(local_times)
                                ind = find(times>=local_times(t),1);
                                if isempty(ind)
                                    if t == length(local_times) && local_times(t)<times(end)+interval/2
                                        local2global(t) = length(times);
                                    else
                                        error('Sample %d could not be placed',t);
                                    end
                                elseif ind==1
                                    if t == 1 && local_times(t)>times(1)-interval/2
                                        local2global(t) = 1;
                                    else
                                        error('Sample %d could not be placed',t);
                                    end
                                else
                                    if round(local_times(t)-times(ind-1),4)<=interval/2
                                        local2global(t) = ind-1;
                                    elseif round(times(ind)-local_times(t),4)<=interval/2
                                        local2global(t) = ind;
                                    else
                                        error('Sample %d could not be placed',t);
                                    end
                                end
                                ppm.increment();
                            end
                            delete(ppm);
                        end
                        
                        if (length(local2global) ~= length(unique(local2global)))
                            % Only dealing with contiguous repetitions for now
                            rep_list = find(local2global == [0,local2global(1:end-1)]);
                            for r = length(rep_list):-1:1
                                if srate_mocap > 0 && rep_list(r) > 2 && (eff_srate >= srate_mocap*0.95 && eff_srate <= srate_mocap*1.05)...
                                        && local2global(rep_list(r)-2)==local2global(rep_list(r))-2
                                    % Correct local2global knowing samples
                                    % should be regularly spaced given the
                                    % srate
                                    local2global(rep_list(r)-1)=local2global(rep_list(r))-1;
                                else
                                    % Remove the duplicate sample and
                                    % average the values in the data
                                    local2global = [local2global(1:rep_list(r)-1),local2global(rep_list(r)+1:end)];
                                    local_times = [local_times(1:rep_list(r)-1),local_times(rep_list(r)+1:end)];
                                    local_data = [local_data(:,1:rep_list(r)-2),...
                                        mean(local_data(:,rep_list(r)-1:rep_list(r)),2),...
                                        local_data(:,rep_list(r)+1:end)];
                                end
                            end
                            
                            if (length(local2global) ~= length(unique(local2global)))
                                error('Still repetitions to deal with')
                            end
                        end
                        
                        if srate_mocap > 0
                            srate_scale = srate_eeg/srate_mocap;
                            if floor(srate_scale) == srate_scale
                                % EEG srate is a multiple of MOCAP srate
                                % 'Stairs' completing (To upsample, the value is copied 'srate_scale' times until the next new sample)
                                for j = 1:srate_scale
                                    data_mocap(channels_span,local2global(1:end-j+1)+j-1)=local_data(:,1:end-j+1);
                                end
                            else
                                % 'Stairs' completing (To upsample, the value is copied 'srate_scale' times until the next new sample)
                                for t = 1:length(local_times)-1
                                    reps = local2global(t+1)-local2global(t);
                                    data_mocap(channels_span,local2global(t):local2global(t+1)-1)=repmat(local_data(:,t),1,reps);
                                end
                                data_mocap(channels_span,local2global(end))=local_data(:,end);
                            end
                        else
                            % Create an additional data channel to indicate
                            % whether the given sample is original or interpolated
                            data_mocap(channels_span(end),:) = zeros(1,length(times));
                            data_mocap(channels_span(end),local2global) = ones(1,length(local2global));
                            
                            % Interpolation
                            data_mocap(channels_span(1:end-1),:) = interp1(times(local2global),local_data',times, 'spline', NaN)';
                            % No interpolation
                            %data_mocap(channels_span(1:end-1),local2global)=local_data;
                        end
                        
                        %% chanlocs
                        for ch = 1:length(channels_span)
                            chanlocs(channels_span(ch)+n_chans_eeg+n_chans_et).type = 'MOCAP';
                            if ch < length(channels_span)
                                channel_desc = AllStreams{moc}.info.desc.channels.channel{ch};
                                chanlocs(channels_span(ch)+n_chans_eeg+n_chans_et).unit = channel_desc.unit;
                                chanlocs(channels_span(ch)+n_chans_eeg+n_chans_et).labels = channel_desc.label;
                            else
                                channel_desc = AllStreams{moc}.info.desc.channels.channel{1};
                                chanlocs(channels_span(ch)+n_chans_eeg+n_chans_et).unit = 'bool';
                                chanlocs(channels_span(ch)+n_chans_eeg+n_chans_et).labels = ['Original_samples_', channel_desc.label(7:end)];
                            end
                        end
                        
                        %% etc
                        Data_all.etc.(['descMOCAP', num2str(find(MOCAP_stream_inds==moc))]) = AllStreams{moc}.info.desc;
                        Data_all.etc.(['infoMOCAP', num2str(find(MOCAP_stream_inds==moc))]) = rmfield(AllStreams{moc}.info,'desc');
                        
                    end
                    disp ('...done')
                end
                
            case 4 %Events
                disp('Importing Events...')
                list = strsplit(full_filename,'_');
                block_str = list{contains(list,'block')};
                block_ind = str2num(block_str(end));
                acq_ind = str2num(list{end});
                all_events = export_events_EEGPOL(AllStreams, Event_stream_inds, times, block_ind, acq_ind, study_config);
                
                disp ('...done')
        end
    end
    
    data = data_eeg;
    if ~isempty(ET_stream_inds)
        data = cat(1, data, data_et);
    end
    if ~isempty(MOCAP_stream_inds)
        data = cat(1, data, data_mocap);
    end
    
    Data_all.data = data;
    [Data_all.nbchan,Data_all.pnts,Data_all.trials] = size(Data_all.data);
    [Data_all.filepath,fname,fext] = fileparts(input_file);
    Data_all.filename = [fname fext];
    Data_all.srate = srate_eeg;
    Data_all.xmin = 0;
    Data_all.xmax = (Data_all.pnts-1)/Data_all.srate;
    
    Data_all.chanlocs = chanlocs;
    
    Data_all.event = all_events;
    %Data_all.urevent = all_events;
    
    % Save file
    if ~isfolder(output_filepath)
        mkdir(output_filepath); % make sure that folder exists
    end
    pop_saveset(Data_all, 'filename',[file_prefix{1} study_config.filenames{i_filename}],'filepath', output_filepath);
    
    clear AllStreams data
    
    %% Split into unique channel types
    % get rid of memory mapped object storage
    %pop_editoptions( 'option_storedisk', 0, 'option_savetwofiles', 1, 'option_saveversion6', 0, 'option_single', 0, 'option_memmapdata', 0, 'option_eegobject', 0, 'option_computeica', 1, 'option_scaleicarms', 1, 'option_rememberfolder', 1, 'option_donotusetoolboxes', 0, 'option_checkversion', 1, 'option_chat', 1);
    
    disp('Splitting dataset into unique channel types...')
    Data_all = pop_loadset('filename',[subject '_AllStreams_' study_config.filenames{i_filename} '.set'],'filepath', output_filepath);
    
    unique_types = unique({Data_all.chanlocs.type});
    for unique_type = 1:length(unique_types)
        if strcmp(unique_types(unique_type),'EEG') && any(strcmp(unique_types,'EOG'))
            indices = find(strcmp({Data_all.chanlocs.type}, 'EEG') | strcmp({Data_all.chanlocs.type}, 'EOG'));
        elseif strcmp(unique_types(unique_type),'EOG')
            continue
        else
            indices = find(strcmp({Data_all.chanlocs.type}, unique_types(unique_type)));            
        end
        file_prefix{end+1} = [subject '_' unique_types{unique_type} '_'];
        fprintf('Type "%s": %d of %d channels. ',unique_types{unique_type}, length(indices), length(Data_all.chanlocs))
        
        Data_split(unique_type) = pop_select(Data_all,'channel',indices);
        Data_split(unique_type) = eeg_checkset(Data_split(unique_type));
        
        % new data set in EEGLAB
        [ALLEEG, ~, CURRENTSET] = pop_newset(ALLEEG, Data_split(unique_type), CURRENTSET, 'gui', 'off');
        % save set        
        pop_saveset(Data_split(unique_type), 'filename', [file_prefix{end} study_config.filenames{i_filename}], 'filepath', output_filepath);
    end
    disp('...done!')
    
    % clear RAM
    STUDY = []; CURRENTSTUDY = 0; ALLEEG = []; EEG=[]; CURRENTSET=[];
    clear Data_all Data_split
    disp(['XDF exporting done for ' full_filename]);
end

% This merges all EEG data files into one big file
input_filepath = output_filepath;

%file_prefix{1} = [subject '_AllStreams_'];

file_prefix = unique(file_prefix);
for f = 1:numel(file_prefix)
    % Only merge EEG file
    if contains(file_prefix{f},'EEG') || contains(file_prefix{f},'AllStreams')
        % get rid of memory mapped object storage
        %pop_editoptions( 'option_storedisk', 0, 'option_savetwofiles', 1, 'option_saveversion6', 0, 'option_single', 0, 'option_memmapdata', 0, 'option_eegobject', 0, 'option_computeica', 1, 'option_scaleicarms', 1, 'option_rememberfolder', 1, 'option_donotusetoolboxes', 0, 'option_checkversion', 1, 'option_chat', 1);
        
        % make sure EEGLAB has no files other than the ones to be merged
        STUDY = []; CURRENTSTUDY = 0; ALLEEG = []; EEG=[]; CURRENTSET=[];
        
        for i_filename = 1:length(study_config.filenames)
            full_filename = [file_prefix{f} study_config.filenames{i_filename} '.set'];
            if isfile([input_filepath filesep full_filename])
                EEG = pop_loadset('filename', full_filename, 'filepath', input_filepath);
                [ALLEEG, EEG, CURRENTSET] = pop_newset(ALLEEG, EEG, 0,'study',0);
            else
                continue
            end
        end
        
        % merges all files currently loaded in EEGLAB into one file and stores
        % the original filenames in EEG.etc.appended_files
        [ALLEEG, EEG_merged, CURRENTSET] = bemobil_merge(ALLEEG,EEG,CURRENTSET,1:length(ALLEEG),...
            [file_prefix{f} study_config.merged_filename], output_filepath);
        %disp('Entire mobilab loading, processing, and exporting done!')
        %disp('You can start the EEGLAB processing now, using the merged dataset.')
    end
end