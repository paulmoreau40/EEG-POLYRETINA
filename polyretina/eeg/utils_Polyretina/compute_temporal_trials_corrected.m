function [corrected_trials] = compute_temporal_trials_corrected(EEG_trial_data, EEG_baseline_data, electrodes_of_interest, bool_plot, bool_all_electrodes, bool_export)

% Baseline corrects EEG trial data, segments it by Field of View (FoV), and prepares data for visualization and export.
%
% Inputs:
% EEG_trial_data            - Struct containing trial EEG data for participants
% EEG_baseline_data         - Struct containing baseline EEG data for participants
% electrodes_of_interest    - List of electrodes to retain for analysis
% bool_plot                 - Whether to plot the corrected trial data
% bool_all_electrodes       - Whether to keep all electrodes or only the ones specified
%
% Outputs:
% corrected_trials          - Struct containing corrected EEG data, segmented by FoV, with statistical summaries (mean, std) across participants

corrected_trials = struct();
participant_ids = unique({EEG_baseline_data.metaInfo.participant_id});

for i = 1:length(participant_ids)
    participant_id = participant_ids{i};
    
    participant_baseline_data = EEG_baseline_data.(participant_id).data2sec;
    mean_baseline = mean(participant_baseline_data, 2);  % Average across time points (dimension 2)
    
    trial_data = EEG_trial_data.(participant_id).data2sec;
    corrected_data = zeros(size(trial_data));
    num_trials = size(trial_data, 3); 
    
    if ismember(participant_id,'P009')
        num_trials = 72;
        corrected_data(:,:,73:end) = [];
    end

    for trial = 1:num_trials
        baseline_index = ceil(trial / 2);  % Use the same baseline for every two consecutive trials
        corrected_data(:, :, trial) = trial_data(:, :, trial) - mean_baseline(:, baseline_index);
    end
    corrected_trials.(participant_id).corrected_data2sec = corrected_data;
end

corrected_trials.metaInfo = EEG_trial_data.metaInfo;
corrected_trials.metaInfo(strcmp({corrected_trials.metaInfo.participant_id}, 'P009') & [corrected_trials.metaInfo.BlockIndex] == 37) = [];

% SELECTING ELECTRODES OF INTEREST
if ~bool_all_electrodes
    disp('Keeping only electrodes of brain region of interest...')
    if ischar(electrodes_of_interest)
        if strcmp(electrodes_of_interest, 'all')
            disp('Keeping all of the electrodes...')
        else
            error('Please enter correct option')
        end
    else
        for p = 1:length(participant_ids)
            chanlocs = EEG_trial_data.(participant_ids{p}).chanlocs;
    
            % Retrieve the indices that correspond to electrodes that we don't want to keep
            idx = find(~ismember({chanlocs.labels}, electrodes_of_interest));
    
            % Removing all data from electrodes that are out of interst (OoI)
            corrected_trials.(participant_ids{p}).corrected_data2sec(idx,:,:) = [];
        end
    end
end



% SEPARATE EEG data into two structures based on the FoV 20 / 45
trials_20 = struct();
trials_45 = struct();
mean_trials_20 = [];
mean_trials_45 = [];

for p = 1:length(participant_ids)
    participant_id = participant_ids{p};
    

    participant_meta_indices = find(strcmp({corrected_trials.metaInfo.participant_id}, participant_id));
    participant_meta_info = corrected_trials.metaInfo(participant_meta_indices);
    %meta_info = corrected_trials.metaInfo;

    trial_data = corrected_trials.(participant_id).corrected_data2sec;
    
    trials_20_indices = find([participant_meta_info.FieldOfView] == 20);
    trials_45_indices = find([participant_meta_info.FieldOfView] == 45);
    
    trials_20_data = trial_data(:, :, trials_20_indices);
    trials_45_data = trial_data(:, :, trials_45_indices);

    trials_20.(participant_id) = squeeze(mean(trials_20_data, 1));  % Average over electrodes
    trials_45.(participant_id) = squeeze(mean(trials_45_data, 1));

    % STORE IN MEGA OUTPUT STRUCTURE
    corrected_trials.(participant_id).corrected_data2sec_20 = trials_20_data;
    corrected_trials.(participant_id).corrected_data2sec_45 = trials_45_data;
end

% COMPUTE MEAN AND STD
for p = 1:length(participant_ids)
    participant_id = participant_ids{p};
    
    mean_trials_20_participant = mean(trials_20.(participant_id), 2); % Mean across trials
    mean_trials_20 = cat(2, mean_trials_20, mean_trials_20_participant); % Stack participants
    
    mean_trials_45_participant = mean(trials_45.(participant_id), 2); % Mean across trials
    mean_trials_45 = cat(2, mean_trials_45, mean_trials_45_participant); % Stack participants
end

mean_20_all_participants = mean(mean_trials_20, 2); % Mean across participants
std_20_all_participants = std(mean_trials_20, 0, 2); % Std across participants
mean_45_all_participants = mean(mean_trials_45, 2); 
std_45_all_participants = std(mean_trials_45, 0, 2);

corrected_trials.mean_20_all_participants = mean_20_all_participants;
corrected_trials.std_20_all_participants = std_20_all_participants;
corrected_trials.mean_45_all_participants = mean_45_all_participants;
corrected_trials.std_45_all_participants = std_45_all_participants;







if bool_plot
    % PLOT MEAN 45 AND 20 WITH STD
    color_20 = [0.1, 0.6, 0.8]; 
    color_45 = [0.8, 0.4, 0.0];
    
    time_points = linspace(0, 2, size(mean_20_all_participants, 1)); 
    figure;
    hold on;
    
    h1 = plot(time_points, mean_20_all_participants, 'Color', color_20, 'LineWidth', 2); % Store handle for legend
    patch('XData', [time_points, fliplr(time_points)], ...
          'YData', [mean_20_all_participants - std_20_all_participants; ...
                    flipud(mean_20_all_participants + std_20_all_participants)], ...
          'FaceColor', color_20, 'EdgeColor', 'none', 'FaceAlpha', 0.2);
    
    h2 = plot(time_points, mean_45_all_participants, 'Color', color_45, 'LineWidth', 2); % Store handle for legend
    patch('XData', [time_points, fliplr(time_points)], ...
          'YData', [mean_45_all_participants - std_45_all_participants; ...
                    flipud(mean_45_all_participants + std_45_all_participants)], ...
          'FaceColor', color_45, 'EdgeColor', 'none', 'FaceAlpha', 0.2);
    
    xlabel('Time (s)');
    ylabel('Amplitude');
    legend([h1, h2], {'20° Trials', '45° Trials'}); % Use plot handles for correct legend
    title('Mean and STD of EEG Trials (20° vs 45°)');
    grid on;
    hold off;
    
    saveas(gcf, fullfile(fullfile(fileparts(fileparts(pwd)), 'figures', 'AnalysisPlots'), 'Last2Seconds_20vs45_participantsAveraged.png'));

    
    
    
    
    
    
    
    
    % PLOT TEMPORAL TRIALS AVERAGED FOR EACH PARTICIPANT (20° FOV)
    figure;
    hold on;
    color_order = lines(size(mean_trials_20, 2)); % Generate colors for each participant
    
    for p = 1:size(mean_trials_20, 2)
        plot(time_points, mean_trials_20(:, p), 'Color', color_order(p, :), 'LineWidth', 2);
    end
    
    xlabel('Time (s)');
    ylabel('Amplitude');
    legend(participant_ids, 'Location', 'best');
    title('Mean Trials (20°) for Each Participant');
    grid on;
    hold off;

    saveas(gcf, fullfile(fullfile(fileparts(fileparts(pwd)), 'figures', 'AnalysisPlots'), 'Last2Seconds_20_eachParticipant.png'));
    
    
    % PLOT TEMPORAL TRIALS AVERAGED FOR EACH PARTICIPANT (45° FOV)
    figure;
    hold on;
    for p = 1:size(mean_trials_45, 2)
        plot(time_points, mean_trials_45(:, p), 'Color', color_order(p, :), 'LineWidth', 2);
    end
    xlabel('Time (s)');
    ylabel('Amplitude');
    legend(participant_ids, 'Location', 'best');
    title('Mean Trials (45°) for Each Participant');
    grid on;
    hold off;
    
    saveas(gcf, fullfile(fullfile(fileparts(fileparts(pwd)), 'figures', 'AnalysisPlots'), 'Last2Seconds_45_eachParticipant.png'));
    
    
    
    
    
    
    % % PLOT MEAN TEMPORAL TRIALS FOR 20 AND 45°, AVERAGED ACROSS ALL TRIALS AND PARTICIPANTS
    all_trials_20 = [];
    all_trials_45 = [];
    
    for i = 1:length(participant_ids)
        participant_id = participant_ids{i};
        
        participant_data_20 = trials_20.(participant_id); 
        participant_data_45 = trials_45.(participant_id); 
        
        all_trials_20 = cat(2, all_trials_20, participant_data_20); 
        all_trials_45 = cat(2, all_trials_45, participant_data_45); 
    end
    
    mean_20 = mean(all_trials_20, 2); 
    std_20 = std(all_trials_20, [], 2);
    
    mean_45 = mean(all_trials_45, 2);
    std_45 = std(all_trials_45, [], 2);
    
    time_points = linspace(0, 2, size(mean_20, 1));
    color_20 = [0.1, 0.6, 0.8]; % Couleur pour les trials à 20°
    color_45 = [0.8, 0.4, 0.0]; % Couleur pour les trials à 45°
    
    figure;
    hold on;
    
    plot(time_points, mean_20, 'Color', color_20, 'LineWidth', 2);
    patch('XData', [time_points, fliplr(time_points)], ...
          'YData', [mean_20 - std_20; flipud(mean_20 + std_20)], ...
          'FaceColor', color_20, 'EdgeColor', 'none', 'FaceAlpha', 0.2);
    
    plot(time_points, mean_45, 'Color', color_45, 'LineWidth', 2);
    patch('XData', [time_points, fliplr(time_points)], ...
          'YData', [mean_45 - std_45; flipud(mean_45 + std_45)], ...
          'FaceColor', color_45, 'EdgeColor', 'none', 'FaceAlpha', 0.2);
    
    xlabel('Time (s)');
    ylabel('Amplitude');
    legend({'20° Trials', '20° Trials (Std)', '45° Trials', '45° Trials (Std)'});
    title('Mean and STD of Concatenated EEG Trials (20° vs 45°)');
    grid on;
    hold off;
    
    saveas(gcf, fullfile(fullfile(fileparts(fileparts(pwd)), 'figures', 'AnalysisPlots'), 'Last2Seconds_20vs45_trialsAveraged.png'));
    
    
    
    
    
    % SPECTROGRAM OF AVERAGED TRIALS
    
    srate = EEG_baseline_data.(participant_ids{1}).srate;
    n_samples = 500;
    mean_20_adjusted = mean_20_all_participants(1:n_samples);
    mean_45_adjusted = mean_45_all_participants(1:n_samples);
    
    window_size = round(0.25 * srate);
    overlap_size = round(0.125 * srate);
    
    [s_20, f_20, t_20] = spectrogram(mean_20_adjusted, hamming(window_size), overlap_size, [], srate, 'yaxis');
    S_20_dB = 10*log10(abs(s_20));
    
    [s_45, f_45, t_45] = spectrogram(mean_45_adjusted, hamming(window_size), overlap_size, [], srate, 'yaxis');
    S_45_dB = 10*log10(abs(s_45));
    
    diff_spectrogram = S_20_dB - S_45_dB;
    
    figure;
    
    % Spectrogram for 20° trials
    subplot(3,1,1);
    imagesc(t_20, f_20, S_20_dB);
    axis xy;
    ylim([0 45]);
    title('Spectrogram for 20° Trials');
    colorbar;
    xlabel('Time (s)');
    ylabel('Frequency (Hz)');
    
    % Spectrogram for 45° trials
    subplot(3,1,2);
    imagesc(t_45, f_45, S_45_dB);
    axis xy;
    ylim([0 45]);
    title('Spectrogram for 45° Trials');
    colorbar;
    xlabel('Time (s)');
    ylabel('Frequency (Hz)');
    
    % Difference between spectrogram
    subplot(3,1,3);
    imagesc(t_20, f_20, diff_spectrogram);
    axis xy;
    ylim([0 45]);
    title('Difference between 20° and 45° Trials');
    colorbar;
    xlabel('Time (s)');
    ylabel('Frequency (Hz)');
    
    saveas(gcf, fullfile(fullfile(fileparts(fileparts(pwd)), 'figures', 'AnalysisPlots'), 'Last2Seconds_20vs45_spectrogram_averagedTrials.png'));
    
    
    
    
    % SPECTROGRAM FOR EACH TRIAL, THEN  AVERAGE OF SPECTROGRAMS (instead of spectrogram of
    % averaged trials)
    
    S_sum_20 = [];
    S_sum_45 = [];
    
    % Individual spectrogram for each 20° trial
    for i = 1:size(all_trials_20, 2)
        [s_20, f_20, t_20] = spectrogram(all_trials_20(:, i), hamming(window_size), overlap_size, [], srate);
        S_20_dB = 10 * log10(abs(s_20));  % Convert in dB
        
        if isempty(S_sum_20)
            S_sum_20 = S_20_dB;
        else
            S_sum_20 = S_sum_20 + S_20_dB;
        end
    end
        
    % Individual spectrogram for each 45° trial
    for i = 1:size(all_trials_45, 2)
        [s_45, f_45, t_45] = spectrogram(all_trials_45(:, i), hamming(window_size), overlap_size, [], srate);
        S_45_dB = 10 * log10(abs(s_45));  % Convertir en dB
        
        if isempty(S_sum_45)
            S_sum_45 = S_45_dB;
        else
            S_sum_45 = S_sum_45 + S_45_dB;  % Accumuler
        end
    end
    
    S_mean_20_dB = S_sum_20 / size(all_trials_20, 2); % Averaged spectrogram (20°)
    S_mean_45_dB = S_sum_45 / size(all_trials_45, 2); % Averaged spectrogram (45°)
    
    % Diff between averaged spectrogram
    diff_spectrogram = S_mean_20_dB - S_mean_45_dB;
    
    freq_indices_20 = find(f_20 <= 45); 
    freq_indices_45 = find(f_45 <= 45);
    
    S_mean_20_limited = S_mean_20_dB(freq_indices_20, :);
    S_mean_45_limited = S_mean_45_dB(freq_indices_45, :);
    diff_spectrogram_limited = diff_spectrogram(freq_indices_20, :);
    
    % Find min and max limits of spectrogram for each subplot, in the studied frequency range
    clim_20 = [min(S_mean_20_limited(:)), max(S_mean_20_limited(:))];
    clim_45 = [min(S_mean_45_limited(:)), max(S_mean_45_limited(:))];
    clim_diff = [min(diff_spectrogram_limited(:)), max(diff_spectrogram_limited(:))];
    
    figure;
    
    % Spectrogram for 20° trials
    subplot(3,1,1);
    imagesc(t_20, f_20(freq_indices_20), S_mean_20_limited);
    axis xy;
    ylim([0 45]); 
    title('Mean Spectrogram for 20° Trials');
    colorbar;
    clim(clim_20);  % Limit the colorbar for the 20° spectrogram
    xlabel('Time (s)');
    ylabel('Frequency (Hz)');
    
    % Spectrogram for 45° trials
    subplot(3,1,2);
    imagesc(t_45, f_45(freq_indices_45), S_mean_45_limited);
    axis xy;
    ylim([0 45]); 
    title('Mean Spectrogram for 45° Trials');
    colorbar;
    clim(clim_45); 
    xlabel('Time (s)');
    ylabel('Frequency (Hz)');
    
    % Spectrogram difference
    subplot(3,1,3);
    imagesc(t_20, f_20(freq_indices_20), diff_spectrogram_limited);
    axis xy;
    ylim([0 45]);  
    title('Difference of Spectrograms (20° - 45°)');
    colorbar;
    clim(clim_diff);
    xlabel('Time (s)');
    ylabel('Frequency (Hz)');
    
    saveas(gcf, fullfile(fullfile(fileparts(fileparts(pwd)), 'figures', 'AnalysisPlots'), 'Last2Seconds_20vs45_averagedSpectrograms.png'));
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    % PLOT ALL TEMPORAL TRIALS
    figure;
    
    % Subplot for 20° trials
    subplot(2, 1, 1);
    hold on;
    for i = 1:size(all_trials_20, 2)
        plot(time_points, all_trials_20(:, i), 'Color', [0.1, 0.6, 0.8, 0.1]); % More transparency for individual trials
    end
    plot(time_points, mean_20, 'Color', [0.1, 0.6, 0.8], 'LineWidth', 2); % Plot the mean over all trials
    xlabel('Time (s)');
    ylabel('Amplitude');
    title('All EEG Trials with Mean (20°)');
    grid on;
    hold off;
    
    % Subplot for 45° trials
    subplot(2, 1, 2);
    hold on;
    for i = 1:size(all_trials_45, 2)
        plot(time_points, all_trials_45(:, i), 'Color', [0.8, 0.4, 0.0, 0.1]); % More transparency for individual trials
    end
    plot(time_points, mean_45, 'Color', [0.8, 0.4, 0.0], 'LineWidth', 2); % Plot the mean over all trials
    xlabel('Time (s)');
    ylabel('Amplitude');
    title('All EEG Trials with Mean (45°)');
    grid on;
    hold off;
    
    
    
    
    
    % BOXPLOT OF TEMPORAL DATA ALONG TIME
    figure;
    hold on;
    
    % Boxplot for 20° trials
    subplot(2,1,1);
    boxplot(all_trials_20', 'PlotStyle', 'compact', 'Colors', color_20);
    title('Boxplot of 20° Trials');
    xlabel('Time points');
    ylabel('Amplitude');
    
    % Boxplot for 45° trials
    subplot(2,1,2);
    boxplot(all_trials_45', 'PlotStyle', 'compact', 'Colors', color_45);
    title('Boxplot of 45° Trials');
    xlabel('Time points');
    ylabel('Amplitude');
    
    hold off;
end







%% EXPORT DATA INTO CSV -> try Deep learning like LSTM to classify FoV ?

if bool_export
    data = corrected_trials;
    IDs = unique({data.metaInfo.participant_id});
    meta = corrected_trials.metaInfo;
    
    output_data_20 = [];
    output_data_45 = [];
    
    for i = 1:length(IDs)
        ID = IDs{i};
        disp(['Process data of participant ', ID, ' for csv exports'])
        
        output_data_20 = [output_data_20; process_FoV(data, meta, ID, 20, 'corrected_data2sec_20')];
        
        output_data_45 = [output_data_45; process_FoV(data, meta, ID, 45, 'corrected_data2sec_45')];
    end
    
    output_data_all = [output_data_20; output_data_45];
    
    var_names = [{'ParticipantID', 'FieldOfView', 'BlockIndex', 'TrialIndex', 'Electrode'}, arrayfun(@(x) ['Time_' num2str(x)], 1:size(output_data_20, 2)-5, 'UniformOutput', false)];
    
    output_table_20 = cell2table(output_data_20, 'VariableNames', var_names);
    output_table_45 = cell2table(output_data_45, 'VariableNames', var_names);
    output_table_all = cell2table(output_data_all, 'VariableNames', var_names);
    
    disp("Exporting into csv data for FoV 20")
    writetable(output_table_20, 'outputs/corrected_data_FoV_20.csv');
    disp("Exporting into csv data for FoV 45")
    writetable(output_table_45, 'outputs/corrected_data_FoV_45.csv');
    disp("Exporting into csv data for all FoV")
    writetable(output_table_all, 'outputs/corrected_data_FoV_all.csv');
    writetable(struct2table(meta), 'outputs/corrected_metadata.csv');
end


















return


% Function to extract EEG and meta data based on the FoV
function output_data = process_FoV(data, meta, ID, FoV, field_name)
    output_data = [];
    meta_i = meta(strcmp({meta.participant_id}, ID) & [meta.FieldOfView] == FoV);
    
    if isfield(data.(ID), field_name)
        data_FoV = data.(ID).(field_name);
        [E, ~, Tr] = size(data_FoV);
        
        for trial = 1:Tr
            for electrode = 1:E
                time_series = squeeze(data_FoV(electrode, :, trial));
                output_row = [{ID, FoV, meta_i(trial).BlockIndex, meta_i(trial).TrialIndex, electrode}, num2cell(time_series)];
                output_data = [output_data; output_row];
            end
        end
    end
return
