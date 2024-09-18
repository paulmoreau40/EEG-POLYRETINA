function [corrected_trials] = compute_temporal_trials_corrected(EEG_trial_data, EEG_baseline_data, electrodes_of_interest, bool_corrected)



corrected_trials = struct();
participant_ids = unique({EEG_baseline_data.metaInfo.participant_id});

for i = 1:length(participant_ids)
    participant_id = participant_ids{i};
    
    participant_baseline_data = EEG_baseline_data.(participant_id).data2sec;
    mean_baseline = mean(participant_baseline_data, 2);  % Average across time points (dimension 2)
    
    trial_data = EEG_trial_data.(participant_id).data2sec;
    corrected_data = zeros(size(trial_data));
    num_trials = size(trial_data, 3);
    
    for trial = 1:num_trials
        baseline_index = ceil(trial / 2);  % Use the same baseline for every two consecutive trials
        corrected_data(:, :, trial) = trial_data(:, :, trial) - mean_baseline(:, baseline_index);
    end
    corrected_trials.(participant_id).corrected_data2sec = corrected_data;
end

corrected_trials.metaInfo = EEG_trial_data.metaInfo;


% SELECTING ELECTRODES OF INTEREST
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




% SPEARATE 20 AND 45, COMPUTE MEAN AND STD
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

    trials_20.(participant_id) = squeeze(mean(trials_20_data, 1));  % Average over electrodes (dimension 1)
    trials_45.(participant_id) = squeeze(mean(trials_45_data, 1));  % Average over electrodes (dimension 1)
end


for p = 1:length(participant_ids)
    participant_id = participant_ids{p};
    
    mean_trials_20_participant = mean(trials_20.(participant_id), 2); % Mean across trials (dimension 2)
    mean_trials_20 = cat(2, mean_trials_20, mean_trials_20_participant); % Stack participants
    
    mean_trials_45_participant = mean(trials_45.(participant_id), 2); % Mean across trials (dimension 2)
    mean_trials_45 = cat(2, mean_trials_45, mean_trials_45_participant); % Stack participants
end

mean_20_all_participants = mean(mean_trials_20, 2); % Mean across participants (dimension 2)
std_20_all_participants = std(mean_trials_20, 0, 2); % Std across participants (dimension 2)
mean_45_all_participants = mean(mean_trials_45, 2); 
std_45_all_participants = std(mean_trials_45, 0, 2);

corrected_trials.mean_20_all_participants = mean_20_all_participants;
corrected_trials.std_20_all_participants = std_20_all_participants;
corrected_trials.mean_45_all_participants = mean_45_all_participants;
corrected_trials.std_45_all_participants = std_45_all_participants;








% PLOT MEAN 45 AND 20 WITH STD
color_20 = [0.1, 0.6, 0.8]; % Color for 20° trials
color_45 = [0.8, 0.4, 0.0]; % Color for 45° trials

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










% PLOT FOR 20° TRIALS (5 participants)
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


% PLOT FOR 45° TRIALS (5 participants)
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







% MEAN AND STD OF ALL TRIALS ACROSS ALL PARTICIPANTS
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






% SPECTROGRAM ANALYSIS AND PLOT

srate = EEG_baseline_data.P001.srate;
n_samples = 500; % Nombre d'échantillons ajustés
mean_20_adjusted = mean_20_all_participants(1:n_samples);
mean_45_adjusted = mean_45_all_participants(1:n_samples);

window_size = round(0.25 * srate); % Taille de la fenêtre en échantillons
overlap_size = round(0.125 * srate); % Chevauchement de 50%

[s_20, f_20, t_20] = spectrogram(mean_20_adjusted, hamming(window_size), overlap_size, [], srate, 'yaxis');
S_20_dB = 10*log10(abs(s_20));

[s_45, f_45, t_45] = spectrogram(mean_45_adjusted, hamming(window_size), overlap_size, [], srate, 'yaxis');
S_45_dB = 10*log10(abs(s_45));

diff_spectrogram = S_20_dB - S_45_dB;

% Plot
figure;

% Spectrogramme pour 20° trials
subplot(3,1,1);
imagesc(t_20, f_20, S_20_dB);
axis xy;
ylim([0 45]); % Limite de 0 à 45 Hz
title('Spectrogram for 20° Trials');
colorbar;
xlabel('Time (s)');
ylabel('Frequency (Hz)');

% Spectrogramme pour 45° trials
subplot(3,1,2);
imagesc(t_45, f_45, S_45_dB);
axis xy;
ylim([0 45]); % Limite de 0 à 45 Hz
title('Spectrogram for 45° Trials');
colorbar;
xlabel('Time (s)');
ylabel('Frequency (Hz)');

% Différence entre les spectrogrammes
subplot(3,1,3);
imagesc(t_20, f_20, diff_spectrogram);
axis xy;
ylim([0 45]); % Limite de 0 à 45 Hz
title('Difference between 20° and 45° Trials');
colorbar;
xlabel('Time (s)');
ylabel('Frequency (Hz)');






% MEAN OF SPECTROGRAMS OF EACH TRIAL

% Paramètres
srate = EEG_baseline_data.P001.srate;
window_size = round(0.25 * srate);  % Fenêtre de 0.25 seconde
overlap_size = round(0.125 * srate);  % Chevauchement de 50%

% Initialiser les matrices de somme pour les spectrogrammes
S_sum_20 = [];
S_sum_45 = [];

% Calculer et accumuler les spectrogrammes pour les trials à 20°
for i = 1:size(all_trials_20, 2)
    [s_20, f_20, t_20] = spectrogram(all_trials_20(:, i), hamming(window_size), overlap_size, [], srate);
    S_20_dB = 10 * log10(abs(s_20));  % Convertir en dB
    
    if isempty(S_sum_20)
        S_sum_20 = S_20_dB;
    else
        S_sum_20 = S_sum_20 + S_20_dB;  % Accumuler
    end
end

% Calculer la moyenne des spectrogrammes pour les trials à 20°
S_mean_20_dB = S_sum_20 / size(all_trials_20, 2);

% Calculer et accumuler les spectrogrammes pour les trials à 45°
for i = 1:size(all_trials_45, 2)
    [s_45, f_45, t_45] = spectrogram(all_trials_45(:, i), hamming(window_size), overlap_size, [], srate);
    S_45_dB = 10 * log10(abs(s_45));  % Convertir en dB
    
    if isempty(S_sum_45)
        S_sum_45 = S_45_dB;
    else
        S_sum_45 = S_sum_45 + S_45_dB;  % Accumuler
    end
end

% Calculer la moyenne des spectrogrammes pour les trials à 45°
S_mean_45_dB = S_sum_45 / size(all_trials_45, 2);

% Calculer la différence entre les spectrogrammes moyens
diff_spectrogram = S_mean_20_dB - S_mean_45_dB;

% Limiter la gamme de fréquences (0 à 45 Hz)
freq_limit = 45;
freq_indices_20 = find(f_20 <= freq_limit);  % Indices des fréquences <= 45 Hz
freq_indices_45 = find(f_45 <= freq_limit);  % Indices des fréquences <= 45 Hz

% Extraire les données correspondant aux fréquences limitées pour chaque spectrogramme
S_mean_20_limited = S_mean_20_dB(freq_indices_20, :);
S_mean_45_limited = S_mean_45_dB(freq_indices_45, :);
diff_spectrogram_limited = diff_spectrogram(freq_indices_20, :);

% Trouver les limites min et max des spectrogrammes dans la gamme de fréquences limitée pour chaque subplot
clim_20 = [min(S_mean_20_limited(:)), max(S_mean_20_limited(:))];
clim_45 = [min(S_mean_45_limited(:)), max(S_mean_45_limited(:))];
clim_diff = [min(diff_spectrogram_limited(:)), max(diff_spectrogram_limited(:))];

% Tracé des spectrogrammes avec les échelles de couleurs spécifiques à chaque subplot
figure;

% Spectrogramme pour les trials à 20°
subplot(3,1,1);
imagesc(t_20, f_20(freq_indices_20), S_mean_20_limited);
axis xy;
ylim([0 freq_limit]);  % Limiter de 0 à 45 Hz
title('Mean Spectrogram for 20° Trials');
colorbar;
caxis(clim_20);  % Limiter la colorbar spécifiquement pour 20°
xlabel('Time (s)');
ylabel('Frequency (Hz)');

% Spectrogramme pour les trials à 45°
subplot(3,1,2);
imagesc(t_45, f_45(freq_indices_45), S_mean_45_limited);
axis xy;
ylim([0 freq_limit]);  % Limiter de 0 à 45 Hz
title('Mean Spectrogram for 45° Trials');
colorbar;
caxis(clim_45);  % Limiter la colorbar spécifiquement pour 45°
xlabel('Time (s)');
ylabel('Frequency (Hz)');

% Différence des spectrogrammes
subplot(3,1,3);
imagesc(t_20, f_20(freq_indices_20), diff_spectrogram_limited);
axis xy;
ylim([0 freq_limit]);  % Limiter de 0 à 45 Hz
title('Difference of Spectrograms (20° - 45°)');
colorbar;
caxis(clim_diff);  % Limiter la colorbar spécifiquement pour la différence
xlabel('Time (s)');
ylabel('Frequency (Hz)');














% OTHER VISUALISATION


% PLOT ALL TRIALS
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





% BOXPLOT OF TEMPORAL DATA ACROSS TIME
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










return