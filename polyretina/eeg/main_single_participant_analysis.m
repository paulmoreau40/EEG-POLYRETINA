% This script processes EEG data using the EEGLAB framework and custom functions.
% It computes power spectra for all trials, baselines, and performs statistical
% comparisons between different conditions (20° vs 45° Field of View).
% The results are visualized through various plots, including absolute spectrum
% and topographical heatmaps for different regions and frequency bands.

% Summary:
% 1. Initialisation and Setup
% 2. Processing and Extracting Power Spectra
% 3. Frequency and Region of Interest Extraction
% 4. Absolute Spectrum Plots
% 5. Analysing last 2 Seconds of Trials
% 6. Permutation Analysis and heatmaps:
%    - Condition vs Baseline
%    - 20° vs 45° Condition Comparison
%    - Baseline-corrected Condition Comparison
% 7. Additional Plots for All Regions
% Helper functions


configEEGPOL;

fontsize_labels = 14;
fontsize_title = 16;

addpath(genpath(pwd));
addpath(fullfile(pwd, '..', '..', 'toolboxes', 'eeglab2024.0'));
addpath(genpath(fullfile(pwd, '..', '..', 'toolboxes', 'ParforProgMon')));
addpath(genpath(fullfile(pwd, '..', '..', 'toolboxes', 'bemobil_pipeline0.2')));

output_filepath = fullfile(fileparts(fileparts(pwd)), 'data', 'analysis', '3_single-subject-analysis', 'bemobil', 'autoMoBI');

overwriteSpectraComputations = false; % to recompute Spectra if needed

if ~exist('ALLEEG','var')
    launchEEGLAB;
end

if ~exist(output_filepath, 'dir')
    mkdir(output_filepath);
end


%% 1. PROCESSING AND EXTRACTING THE POWER SPECTRA FOR ALL TRIALS, BASELINES AND COARSE BASELINE

subject_inds = [1, 2, 3, 6, 7, 8, 9];
% subject_inds = 9; % Overwrite subject for testing (COMMENT / DECOMMENT)

if (~exist(fullfile(output_filepath, 'EEG_trial_data.mat'),'file') || ~exist(fullfile(output_filepath, 'EEG_baseline_data.mat'),'file') || overwriteSpectraComputations)
    for subject_ind = subject_inds
        % Load  and filter the data
        subject = study_config.subjects(subject_ind).id;
        disp(['Subject ' subject]);
        study_config.current_subject = subject_ind;
        N = makeFolderFileNames(study_config, subject);
        EEG = pop_loadset('filename', N.postLabelingFile, 'filepath', N.searchFolder_2arch_rej_ICcats);
                
        % Filter the dataset with the desired upper and lower frequencies
        lowcutoff = study_config.filterAnalysis.low_cut_off;
        highcutoff = study_config.filterAnalysis.high_cut_off;
        fprintf('Filtering between %.1f Hz and %.1f Hz...\n', lowcutoff, highcutoff)
        [EEG] = custom_filter(EEG, lowcutoff, highcutoff);

        % "Epoch" Data: Removing events which are not of interest
        boundaryEvents_mask = strcmp({EEG.event.type}, 'boundary');
        EEG2 = pop_editeventvals(EEG, 'delete', find(boundaryEvents_mask));

        % Retrieving Segments of Interest, corresponding Baselines, and computing their Spectra
        if ~exist('EEG_trial_data','var') % creating structures if don't exist
            disp('Creating EEG_trial_data structure...')
            EEG_trial_data = [];
            EEG_trial_data.metaInfo = struct('participant_id', [],'BlockIndex', [],...
                'TrialIndex', [], 'FieldOfView', []);
        end
        if ~exist('EEG_baseline_data','var')
            disp('Creating EEG_baseline_data structure...')
            EEG_baseline_data = [];
            EEG_baseline_data.metaInfo = struct('participant_id', [],'BlockIndex', [],...
                'TrialIndex', [], 'FieldOfView', []);
        end
        if ~exist('EEG_coarse_data','var')
            disp('Creating EEG_coarse_data structure...')
            EEG_coarse_data = [];
            EEG_coarse_data.metaInfo = struct('participant_id', [],'type', [], 'idx', []);
        end

        [EEG_trial_data, EEG_baseline_data, EEG_coarse_data] = extract_segments_EEG_compute_spectrum(EEG2, EEG_trial_data, EEG_baseline_data, EEG_coarse_data, subject, false);
    end

    % Removing initial shift which is created when concatenating metaInfo
    EEG_trial_data.metaInfo(1) = [];
    EEG_baseline_data.metaInfo(1) = [];
    EEG_coarse_data.metaInfo(1) = [];

    disp('Saving spectra of trial, baseline and coarse data...')
    save(fullfile(output_filepath, 'EEG_trial_data.mat'),'-struct','EEG_trial_data');
    save(fullfile(output_filepath, 'EEG_baseline_data.mat'),'-struct','EEG_baseline_data');
    save(fullfile(output_filepath, 'EEG_coarse_data.mat'),'-struct','EEG_coarse_data');
    disp('Finished saving...');
else
    EEG_trial_data = load(fullfile(output_filepath, 'EEG_trial_data.mat'), '-mat');
    EEG_baseline_data = load(fullfile(output_filepath, 'EEG_baseline_data.mat'), '-mat');
    EEG_coarse_data = load(fullfile(output_filepath, 'EEG_coarse_data.mat'), '-mat');
    disp('Finished loading EEG_trial_data EEG_baseline_data and EEG_coarse_data structures');
end




%% 2. CHOOSING GLOBAL PARAMETERS FOR WHAT FOLLOWS:

% Defining which electrodes we are considering for the respective brain regions
frontal_electrodes ={'Z1','Z2','Z3','Z4','Z5','L1','L2','L3','L4','L5','LL1','LL2','LL3','LL4','LL5','R1','R2','R3','R4','R5','RR1','RR2','RR3','RR4','RR5'}; %{'LL3','LL4', 'L4', 'Z3', 'Z4', 'R4', 'RR3', 'RR4'};
%{'Z1','Z2','Z3','Z4','Z5','L1','L2','L3','L4','L5','LL1','LL2','LL3','LL4','LL5','R1','R2','R3','R4','R5','RR1','RR2','RR3','RR4','RR5'};
parietal_electrodes = {'LL6','LL7','LL8','RR6','RR7','RR8','RA3','RA4','LA3','LA4','L7','L8','R7','R8'}; %{'L8', 'R8', 'Z8', 'LL8', 'RR8'};
%{'LL6','LL7','LL8','RR6','RR7','RR8','RA3','RA4','LA3','LA4','L7','L8','R7','R8'};
occipital_electrodes = {'Z12', 'Z11','Z10','R11','L11','R12','L12'};

% Which Electrodes to consider for Brain Regions of Interest:
brain_region_chosen = frontal_electrodes; % 'occipital_electrodes', 'parietal_electrodes', 'frontal_electrodes'
brain_region_name = 'Frontal'; % 'Occipital', 'Parietal', 'Frontal'

file_electrode_positions = strcat(study_config.study_folder,study_config.raw_data_folder,'P001\',study_config.channel_locations_filename);



%% 3. EXTRACT FREQUENCIES AND REGIONS OF INTEREST, FOR EACH FOV (FOR PLOTS)

% 1. Dropping the computed amplitudes for frequencies which are not of interest
if ~exist('lowcutoff','var')
    lowcutoff = study_config.filterAnalysis.low_cut_off;
    highcutoff = study_config.filterAnalysis.high_cut_off;
end

% 2. Removing frequencies which are filtered out (TRIALS / BASELINES / COARSE BASELINES)
EEG_trial_data = drop_useless_frequencies(EEG_trial_data,lowcutoff, highcutoff, 'absolute');
EEG_baseline_data = drop_useless_frequencies(EEG_baseline_data,lowcutoff,highcutoff, 'absolute');
EEG_coarse_data = drop_useless_frequencies(EEG_coarse_data,lowcutoff,highcutoff, 'absolute');


% 3. Extract region of interest for each FoV
[EEG_selected_absolute_spectrum_20] = extract_trials_according_to_brainregion_and_frequency(EEG_trial_data, brain_region_chosen, 'all', 1, 20, 'absolute');
[EEG_selected_absolute_spectrum_45] = extract_trials_according_to_brainregion_and_frequency(EEG_trial_data, brain_region_chosen, 'all', 1, 45, 'absolute');
[EEG_selected_absolute_base_spectrum_20] = extract_trials_according_to_brainregion_and_frequency(EEG_baseline_data, brain_region_chosen, 'all', 1, 20, 'absolute');
[EEG_selected_absolute_base_spectrum_45] = extract_trials_according_to_brainregion_and_frequency(EEG_baseline_data, brain_region_chosen, 'all', 1, 45, 'absolute');
[EEG_selected_absolute_coarse_spectrum_black] = extract_trials_according_to_brainregion_and_frequency(EEG_coarse_data, brain_region_chosen, 'all', 0, "black", 'absolute');
[EEG_selected_absolute_coarse_spectrum_edge] = extract_trials_according_to_brainregion_and_frequency(EEG_coarse_data, brain_region_chosen, 'all', 0, "edge", 'absolute');

% 4. Concatenating all of the data across the trials to format them for the plots:
[EEG_selected_absolute_spectrum_20_all_trials] = format_for_plotting_spectra(EEG_selected_absolute_spectrum_20); % TRIALS
[EEG_selected_absolute_spectrum_45_all_trials] = format_for_plotting_spectra(EEG_selected_absolute_spectrum_45);
[EEG_selected_absolute_base_spectrum_20_all_trials] = format_for_plotting_spectra(EEG_selected_absolute_base_spectrum_20); % BASELINE
[EEG_selected_absolute_base_spectrum_45_all_trials] = format_for_plotting_spectra(EEG_selected_absolute_base_spectrum_45);
[EEG_selected_absolute_coarse_spectrum_black_all_trials] = format_for_plotting_spectra(EEG_selected_absolute_coarse_spectrum_black); % COARSE
[EEG_selected_absolute_coarse_spectrum_edge_all_trials] = format_for_plotting_spectra(EEG_selected_absolute_coarse_spectrum_edge);

% 5. Averating across electrodes to get brain regions of interest:
EEG_absolute_spectrum_20_all_trials_averaged = mean(EEG_selected_absolute_spectrum_20_all_trials.spectrum,1); % TRIALS
EEG_absolute_spectrum_45_all_trials_averaged = mean(EEG_selected_absolute_spectrum_45_all_trials.spectrum,1);
EEG_absolute_base_spectrum_20_all_trials_averaged = mean(EEG_selected_absolute_base_spectrum_20_all_trials.spectrum,1); % BASELINE
EEG_absolute_base_spectrum_45_all_trials_averaged = mean(EEG_selected_absolute_base_spectrum_45_all_trials.spectrum,1);
EEG_absolute_coarse_spectrum_black_all_trials_averaged = mean(EEG_selected_absolute_coarse_spectrum_black_all_trials.spectrum,1); % COARSE
EEG_absolute_coarse_spectrum_edge_all_trials_averaged = mean(EEG_selected_absolute_coarse_spectrum_edge_all_trials.spectrum,1);

% 6. Averaging across trials (retrieving mean and standard deviation)
EEG_absolute_spectrum_20_averaged = mean(EEG_absolute_spectrum_20_all_trials_averaged,3);
EEG_absolute_spectrum_20_std = std(EEG_absolute_spectrum_20_all_trials_averaged, [], 3)/sqrt(size(EEG_absolute_spectrum_20_all_trials_averaged,3));
EEG_absolute_spectrum_45_averaged = mean(EEG_absolute_spectrum_45_all_trials_averaged,3);
EEG_absolute_spectrum_45_std = std(EEG_absolute_spectrum_45_all_trials_averaged, [], 3)/sqrt(size(EEG_absolute_spectrum_45_all_trials_averaged,3));

EEG_absolute_base_spectrum_20_averaged = mean(EEG_absolute_base_spectrum_20_all_trials_averaged,3);
EEG_absolute_base_spectrum_20_std = std(EEG_absolute_base_spectrum_20_all_trials_averaged, [], 3)/sqrt(size(EEG_absolute_base_spectrum_20_all_trials_averaged,3));
EEG_absolute_base_spectrum_45_averaged = mean(EEG_absolute_base_spectrum_45_all_trials_averaged,3);
EEG_absolute_base_spectrum_45_std = std(EEG_absolute_base_spectrum_45_all_trials_averaged, [], 3)/sqrt(size(EEG_absolute_base_spectrum_45_all_trials_averaged,3));

EEG_absolute_coarse_spectrum_black_averaged = mean(EEG_absolute_coarse_spectrum_black_all_trials_averaged,3); %EEG_selected_absolute_coarse_spectrum_black_all_trials
EEG_absolute_coarse_spectrum_black_std = std(EEG_absolute_coarse_spectrum_black_all_trials_averaged, [], 3)/sqrt(size(EEG_absolute_coarse_spectrum_black_all_trials_averaged,3));
EEG_absolute_coarse_spectrum_edge_averaged = mean(EEG_absolute_coarse_spectrum_edge_all_trials_averaged,3);
EEG_absolute_coarse_spectrum_edge_std = std(EEG_absolute_coarse_spectrum_edge_all_trials_averaged, [], 3)/sqrt(size(EEG_absolute_coarse_spectrum_edge_all_trials_averaged,3));

% 7. Fill in function to plot std
EEG_absolute_spectrum_20_averaged_std_above = EEG_absolute_spectrum_20_averaged + EEG_absolute_spectrum_20_std./2;
EEG_absolute_spectrum_20_averaged_std_below = EEG_absolute_spectrum_20_averaged - EEG_absolute_spectrum_20_std./2;
inBetween_absolute_20 = [EEG_absolute_spectrum_20_averaged_std_below(:); flipud(EEG_absolute_spectrum_20_averaged_std_above(:))];

EEG_absolute_spectrum_45_averaged_std_above = EEG_absolute_spectrum_45_averaged + EEG_absolute_spectrum_45_std./2;
EEG_absolute_spectrum_45_averaged_std_below = EEG_absolute_spectrum_45_averaged - EEG_absolute_spectrum_45_std./2;
inBetween_absolute_45 = [EEG_absolute_spectrum_45_averaged_std_below(:); flipud(EEG_absolute_spectrum_45_averaged_std_above(:))];

% Baseline
EEG_absolute_base_spectrum_20_averaged_above = EEG_absolute_base_spectrum_20_averaged + EEG_absolute_base_spectrum_20_std./2;
EEG_absolute_base_spectrum_20_averaged_below = EEG_absolute_base_spectrum_20_averaged - EEG_absolute_base_spectrum_20_std./2;
inBetween_absolute_baseline_20 = [EEG_absolute_base_spectrum_20_averaged_below(:); flipud(EEG_absolute_base_spectrum_20_averaged_above(:))];

EEG_absolute_base_spectrum_45_averaged_above = EEG_absolute_base_spectrum_45_averaged + EEG_absolute_base_spectrum_45_std./2;
EEG_absolute_base_spectrum_45_averaged_below = EEG_absolute_base_spectrum_45_averaged - EEG_absolute_base_spectrum_45_std./2;
inBetween_absolute_baseline_45 = [EEG_absolute_base_spectrum_45_averaged_below(:); flipud(EEG_absolute_base_spectrum_45_averaged_above(:))];

% Coarse
EEG_absolute_coarse_spectrum_black_averaged_above = EEG_absolute_coarse_spectrum_black_averaged + EEG_absolute_coarse_spectrum_black_std./2;
EEG_absolute_coarse_spectrum_black_averaged_below = EEG_absolute_coarse_spectrum_black_averaged - EEG_absolute_coarse_spectrum_black_std./2;
inBetween_absolute_coarse_black = [EEG_absolute_coarse_spectrum_black_averaged_below(:); flipud(EEG_absolute_coarse_spectrum_black_averaged_above(:))];

EEG_absolute_coarse_spectrum_edge_averaged_above = EEG_absolute_coarse_spectrum_edge_averaged + EEG_absolute_coarse_spectrum_edge_std./2;
EEG_absolute_coarse_spectrum_edge_averaged_below = EEG_absolute_coarse_spectrum_edge_averaged - EEG_absolute_coarse_spectrum_edge_std./2;
inBetween_absolute_coarse_edge = [EEG_absolute_coarse_spectrum_edge_averaged_below(:); flipud(EEG_absolute_coarse_spectrum_edge_averaged_above(:))];


% 8. Converting into dB for the log plots
EEG_absolute_spectrum_20_averaged_dB = 10*log10(EEG_absolute_spectrum_20_averaged); % TRIALS
EEG_absolute_spectrum_45_averaged_dB = 10*log10(EEG_absolute_spectrum_45_averaged);
EEG_absolute_spectrum_20_std_dB = 10*log10(EEG_absolute_spectrum_20_std);
EEG_absolute_spectrum_45_std_dB = 10*log10(EEG_absolute_spectrum_45_std);
inBetween_absolute_20_dB = 10*log10(inBetween_absolute_20);
inBetween_absolute_45_dB = 10*log10(inBetween_absolute_45);

EEG_absolute_base_spectrum_20_averaged_dB = 10*log10(EEG_absolute_base_spectrum_20_averaged); % BASELINE
EEG_absolute_base_spectrum_45_averaged_dB = 10*log10(EEG_absolute_base_spectrum_45_averaged);
EEG_absolute_base_spectrum_20_std_dB = 10*log10(EEG_absolute_base_spectrum_20_std);
EEG_absolute_base_spectrum_45_std_dB = 10*log10(EEG_absolute_base_spectrum_45_std);
inBetween_absolute_baseline_20_dB = 10*log10(inBetween_absolute_baseline_20);
inBetween_absolute_baseline_45_dB = 10*log10(inBetween_absolute_baseline_45);

EEG_absolute_coarse_spectrum_black_averaged_dB = 10*log10(EEG_absolute_coarse_spectrum_black_averaged); % COARSE
EEG_absolute_coarse_spectrum_edge_averaged_dB = 10*log10(EEG_absolute_coarse_spectrum_edge_averaged);
EEG_absolute_coarse_spectrum_black_std_dB = 10*log10(EEG_absolute_coarse_spectrum_black_std);
EEG_absolute_coarse_spectrum_edge_std_dB = 10*log10(EEG_absolute_coarse_spectrum_edge_std);
inBetween_absolute_coarse_black_dB = 10*log10(inBetween_absolute_coarse_black);
inBetween_absolute_coarse_edge_dB = 10*log10(inBetween_absolute_coarse_edge);

freqs_of_interest = EEG_selected_absolute_spectrum_20_all_trials.freqs;



% 9. PLOT TEST FOR COARSE BASELINES

color_black = [0, 0, 0]; % Black color for black baseline
color_edge = [0.8, 0.4, 0.0]; % Edge color for edge baseline

% Plot absolute spectra for black and edge baselines:
figure;
plot(freqs_of_interest, squeeze(EEG_absolute_coarse_spectrum_black_averaged_dB), 'Color', color_black, 'LineWidth', 2);
hold on;
plot(freqs_of_interest, squeeze(EEG_absolute_coarse_spectrum_edge_averaged_dB), 'Color', color_edge, 'LineWidth', 2);
patch('XData', [freqs_of_interest(:); flipud(freqs_of_interest(:))], 'YData', inBetween_absolute_coarse_black_dB, ...
      'FaceColor', color_black, 'EdgeColor', color_black, 'FaceAlpha', 0.2);
patch('XData', [freqs_of_interest(:); flipud(freqs_of_interest(:))], 'YData', inBetween_absolute_coarse_edge_dB, ...
      'FaceColor', color_edge, 'EdgeColor', color_edge, 'FaceAlpha', 0.2);
hold off;
grid on;
xlabel('Frequencies [Hz]');
ylabel('Power [dB]');
legend('Black Baseline', 'Edge Baseline');
title("Absolute Spectrum for " + brain_region_name + " Electrodes Across FoV");

saveas(gcf, fullfile(fullfile(fileparts(fileparts(pwd)), 'figures', 'AnalysisPlots'), 'SpectrumCoarseBaselinesComparison.png'));




%% 4. PLOTS FOR ABSOLUTE SPECTRUM

freqs_of_interest = EEG_selected_absolute_spectrum_20_all_trials.freqs;
x2 = [freqs_of_interest(:); flipud(freqs_of_interest(:))]; % needed for std plot

color_20 = [0.83 0.14 0.14];
color_45 = [1.00 0.54 0.00];
color_baseline = 1/255 * [0, 104, 87];

% Plotting absolute spectra separately for each FoV, averaged over all trials and all electrodes of brain RoI
figure;
plot(freqs_of_interest, EEG_absolute_spectrum_20_averaged_dB', 'Color', color_20, 'LineWidth', 2);
hold on;
plot(freqs_of_interest, EEG_absolute_spectrum_45_averaged_dB', 'Color', color_45, 'LineWidth', 2);
patch('XData',x2,'YData',inBetween_absolute_20_dB,'FaceColor', color_20,'EdgeColor',color_20,'FaceAlpha', 0.2);
patch('XData',x2,'YData',inBetween_absolute_45_dB,'FaceColor', color_45,'EdgeColor',color_45,'FaceAlpha', 0.2);
hold off;
grid on;
xlabel('Frequencies [Hz]');
ylabel('Power [dB]');
legend('20°','45°');
title({['Absolute Spectrum for ' brain_region_name ' Electrodes' ],'Across FoV'});
saveas(gcf, fullfile(fullfile(fileparts(fileparts(pwd)), 'figures', 'AnalysisPlots'), 'SpectrumAveragedTrials.png'));


% Plot of every single trial line
figure;
plot(freqs_of_interest, 10*log10(EEG_absolute_spectrum_20_all_trials_averaged(:,:,1)), 'Color', color_20)
hold on;
for line = 2:size(EEG_absolute_spectrum_20_all_trials_averaged,3)
    plot(freqs_of_interest, 10*log10(EEG_absolute_spectrum_20_all_trials_averaged(:,:,line)), 'Color', color_20);
end
for line = 1:size(EEG_absolute_spectrum_45_all_trials_averaged,3)
    plot(freqs_of_interest, 10*log10(EEG_absolute_spectrum_45_all_trials_averaged(:,:,line)), 'Color', color_45);
end
hold off;
grid on;
xlabel('Frequencies [Hz]');
ylabel('Power [dB]');
title(['Absolute Spectrum for ' brain_region_name ' Electrodes (all trials)']);
saveas(gcf, fullfile(fullfile(fileparts(fileparts(pwd)), 'figures', 'AnalysisPlots'), 'SpectrumAllTrials.png'));

% Plotting absolute spectra OF BASELINE for each FoV, averaged over all trials and all electrodes of brain RoI
figure;
plot(freqs_of_interest, 10*log10(EEG_absolute_base_spectrum_20_all_trials_averaged(:,:,1)), 'Color', color_20);
hold on;
for line = 2:size(EEG_absolute_base_spectrum_20_all_trials_averaged,3)
    plot(freqs_of_interest, 10*log10(EEG_absolute_base_spectrum_20_all_trials_averaged(:,:,line)), 'Color', color_20);
end
for line = 1:size(EEG_absolute_base_spectrum_45_all_trials_averaged,3)
    plot(freqs_of_interest, 10*log10(EEG_absolute_base_spectrum_45_all_trials_averaged(:,:,line)), 'Color', color_45);
end
hold off;
grid on;
xlabel('Frequencies [Hz]');
ylabel('Power [dB]');
title({['Absolute Baseline Spectrum for ' brain_region_name ' Electrodes'],'(1 baseline / trial)'});
saveas(gcf, fullfile(fullfile(fileparts(fileparts(pwd)), 'figures', 'AnalysisPlots'), 'SpectrumAllBaselines.png'));



%% 5. EXTRACTING AND PLOTTING THE LAST 2 SECONDS OF TRIALS
% Instead of looking at the whole trials durations, we only extract the
% data from the last 2 seconds of which trials to see if different patterns
% show (as cognitive load could reach a maximum towards the end of the
% trials). Every trial is normalised by its corresponding baseline

bool_plot = 0; bool_all_electrodes = 1; bool_export = 0;

EEG_trial_data2sec = compute_temporal_trials_corrected(EEG_trial_data, EEG_baseline_data, brain_region_chosen, bool_plot, bool_all_electrodes, bool_export);








%% 6. Permutation Analysis -- Multi Subject Analysis

%% Part 1: Condition VS. Baseline

N_colors_standard = 512; % colormap
y = linspace(1, 42, 83); % for plot
organise_alphabetically_electrodes = 1;

participants = unique({EEG_trial_data.metaInfo(:).participant_id});

[spectrum_trial_20, spectrum_trial_45] = format_data_for_multisubject_stats(EEG_trial_data, [1 42]);
[spectrum_baseline_20, spectrum_baseline_45] = format_data_for_multisubject_stats(EEG_baseline_data, [1 42]);

% Define common options structure
commonOptions = struct();
commonOptions.model = 'classic';
commonOptions.style = 'chanXfreq';
commonOptions.Chans = {EEG_trial_data.(participants{end}).chanlocs(:).labels};
commonOptions.Freqs = 1:83; 
commonOptions.ElecFile = strcat(study_config.study_folder, study_config.raw_data_folder, 'P001\', study_config.channel_locations_filename);
commonOptions.MaxDeg = 20;
commonOptions.pairing = 'on';
commonOptions.N_reps = 256;
commonOptions.reusePerms = false;
commonOptions.removeSmallestClusters = false;

if ~exist('clustered_stats_table_20vbase', 'var')
    [clustered_stats_table_20vbase, statistical_clusters_20vbase, stats_surrog_20vbase, pairwise_stats_20vbase, permutations_20vbase] = ...
    compute_permutations('FoV20', [], spectrum_trial_20, spectrum_baseline_20, commonOptions);
end

if ~exist('clustered_stats_table_45vbase', 'var')
    [clustered_stats_table_45vbase, statistical_clusters_45vbase, stats_surrog_45vbase, pairwise_stats_45vbase, permutations_45vbase] = ...
    compute_permutations('FoV45', [], spectrum_trial_45, spectrum_baseline_45, commonOptions);
end

% Formating Data for Heatmaps
data_heatmap_20vbase = format_for_heatmap_conditionvbaseline(clustered_stats_table_20vbase, statistical_clusters_20vbase, spectrum_trial_20, spectrum_baseline_20);
data_heatmap_45vbase = format_for_heatmap_conditionvbaseline(clustered_stats_table_45vbase, statistical_clusters_45vbase, spectrum_trial_45, spectrum_baseline_45);

[data_heatmap_20vbase, ~] = organize_by_electrodes(data_heatmap_20vbase, {EEG_trial_data.(participants{end}).chanlocs(:).labels}, organise_alphabetically_electrodes);
[data_heatmap_45vbase, new_electrode_labels] = organize_by_electrodes(data_heatmap_45vbase, {EEG_trial_data.(participants{end}).chanlocs(:).labels}, organise_alphabetically_electrodes);

plot_heatmap_baseline_or_condition('20° FoV vs. Baseline (permutation testing)', data_heatmap_20vbase, new_electrode_labels, y, N_colors_standard, "20vsBaseline");
plot_heatmap_baseline_or_condition('45° FoV vs. Baseline (permutation testing)', data_heatmap_45vbase, new_electrode_labels, y, N_colors_standard, "45vsBaseline");



%% Part 2: Condition VS. Condition (20 vs 45)

if ~exist('clustered_stats_table_conditionvcondition','var')
    [clustered_stats_table_conditionvcondition, statistical_clusters_conditionvcondition, stats_surrog_conditionvcondition, pairwise_stats_conditionvcondition, permutations_conditionvcondition] = ...
    compute_permutations('FoV20', 'FoV45', spectrum_trial_20, spectrum_trial_45, commonOptions);
end

data_heatmap_conditionvcondition = format_for_heatmap_conditionvcondition_anova(clustered_stats_table_conditionvcondition, statistical_clusters_conditionvcondition);

[data_heatmap_conditionvcondition, new_electrode_labels] = organize_by_electrodes(data_heatmap_conditionvcondition, {EEG_trial_data.(participants{end}).chanlocs(:).labels}, 1);

plot_heatmap_baseline_or_condition('20° FoV vs. 45° FoV (permutation testing)', data_heatmap_conditionvcondition, new_electrode_labels, y, N_colors_standard, "20vs45");



%% Part 3: Baseline-corrected Condition VS. Condition (20 vs 45)

spectrum_trial_20_dB = 10*log10(spectrum_trial_20);
spectrum_trial_45_dB = 10*log10(spectrum_trial_45);

spectrum_trial_20_base_dB = 10*log10(spectrum_baseline_20);
spectrum_trial_45_base_dB = 10*log10(spectrum_baseline_45);

spectrum_trial_relative_20vbase_dB = spectrum_trial_20_dB - spectrum_trial_20_base_dB; % NORMALISATION
spectrum_trial_relative_45vbase_dB =  spectrum_trial_45_dB - spectrum_trial_45_base_dB;

if ~exist('clustered_stats_table_baseline_corrected','var')
    [clustered_stats_table_baseline_corrected, statistical_clusters_baselinecorrected, stats_surrog_baselinecorrected, pairwise_stats_baselinecorrected, permutations_baselinecorrected] = ...
    compute_permutations('FoV20vbase', 'FoV45vbase', spectrum_trial_relative_20vbase_dB, spectrum_trial_relative_45vbase_dB, commonOptions);
end

data_heatmap_baselinecorrected = format_for_heatmap_baselinecorrected_dB(clustered_stats_table_baseline_corrected, statistical_clusters_baselinecorrected, spectrum_trial_relative_20vbase_dB, spectrum_trial_relative_45vbase_dB);

[data_heatmap_baselinecorrected, new_electrode_labels] = organize_by_electrodes(data_heatmap_baselinecorrected, {EEG_trial_data.(participants{end}).chanlocs(:).labels}, organise_alphabetically_electrodes);

plot_heatmap_baseline_or_condition('Baseline-corrected 20° VS. baseline-corrected 45° (permutation testing)', data_heatmap_baselinecorrected, new_electrode_labels, ...
    y, N_colors_standard, "20vs45corrected");








%SALUT PAUL LA PIOCHE :) travaille esclave
%% 7. ADDITIONAL PLOTS FOR ALL REGIONS

% 7.1 SPECTRA : Condition v Baseline & Condition v Condition

plot_spectrum_all('parietal', 20, [], spectrum_trial_20, spectrum_baseline_20, participants, EEG_trial_data, color_20, color_baseline, freqs_of_interest, x2, 'trial_vs_baseline');
plot_spectrum_all('occipital', 20, [], spectrum_trial_20, spectrum_baseline_20, participants, EEG_trial_data, color_20, color_baseline, freqs_of_interest, x2, 'trial_vs_baseline');
plot_spectrum_all('frontal', 20, [], spectrum_trial_20, spectrum_baseline_20, participants, EEG_trial_data, color_20, color_baseline, freqs_of_interest, x2, 'trial_vs_baseline');

plot_spectrum_all('parietal', 45, [], spectrum_trial_45, spectrum_baseline_45, participants, EEG_trial_data, color_45, color_baseline, freqs_of_interest, x2, 'trial_vs_baseline');
plot_spectrum_all('occipital', 45, [], spectrum_trial_45, spectrum_baseline_45, participants, EEG_trial_data, color_45, color_baseline, freqs_of_interest, x2, 'trial_vs_baseline');
plot_spectrum_all('frontal', 45, [], spectrum_trial_45, spectrum_baseline_45, participants, EEG_trial_data, color_45, color_baseline, freqs_of_interest, x2, 'trial_vs_baseline');

plot_spectrum_all('parietal', 20, 45, spectrum_trial_20, spectrum_trial_45, participants, EEG_trial_data, color_20, color_45, freqs_of_interest, x2, 'trial_vs_trial');
plot_spectrum_all('occipital', 20, 45, spectrum_trial_20, spectrum_trial_45, participants, EEG_trial_data, color_20, color_45, freqs_of_interest, x2, 'trial_vs_trial');
plot_spectrum_all('frontal', 20, 45, spectrum_trial_20, spectrum_trial_45, participants, EEG_trial_data, color_20, color_45, freqs_of_interest, x2, 'trial_vs_trial');



% 7.2 TOPOPLOTS : Condition v Condition (and Condition normalised by baseline)

spectrum_trial_20_avg = mean(spectrum_trial_20, 3); % average across all subjects
spectrum_baseline_20_avg = mean(spectrum_baseline_20, 3);
spectrum_trial_45_avg = mean(spectrum_trial_45, 3);
spectrum_baseline_45_avg = mean(spectrum_baseline_45, 3);

% [20; 35] Hz
freq_20_35 = 43:55; % indices
spectrum_20_dB = 10*log10(mean(spectrum_trial_20_avg(:, freq_20_35), 2));
spectrum_45_dB = 10*log10(mean(spectrum_trial_45_avg(:, freq_20_35), 2));
relative_20_45_dB = spectrum_20_dB - spectrum_45_dB;

plot_topoplot(spectrum_20_dB, {'Heatmap for 20-35 Hz', 'FoV 20°'}, file_electrode_positions, []);
plot_topoplot(spectrum_45_dB, {'Heatmap for 20-35 Hz', 'FoV 45°'}, file_electrode_positions, []);
plot_topoplot(relative_20_45_dB, {'Heatmap for 20-35 Hz', 'Relative 20° vs 45°'}, file_electrode_positions, []);

% [10; 14] Hz (Normalised by Baseline)
freq_10_14 = 20:28; 
spectrum_20_dB = 10*log10(mean(spectrum_trial_20_avg(:, freq_10_14), 2));
baseline_20_dB = 10*log10(mean(spectrum_baseline_20_avg(:, freq_10_14), 2));

spectrum_45_dB = 10*log10(mean(spectrum_trial_45_avg(:, freq_10_14), 2));
baseline_45_dB = 10*log10(mean(spectrum_baseline_45_avg(:, freq_10_14), 2));

relative_20_dB = spectrum_20_dB - baseline_20_dB; % normalisation
relative_45_dB = spectrum_45_dB - baseline_45_dB;

set(0,'DefaultFigureColormap', asymColorMapWhiteZero([-5, 5], N_colors_standard));
plot_topoplot(relative_20_dB, {'Topoplot for 10-14 Hz', 'Relative FoV 20°'}, file_electrode_positions, [-5, 5]);
plot_topoplot(relative_45_dB, {'Topoplot for 10-14 Hz', 'Relative FoV 45°'}, file_electrode_positions, [-5, 5]);


% [7; 12] Hz
freq_7_12 = 13:23;
spectrum_20_dB = 10*log10(mean(spectrum_trial_20_avg(:, freq_7_12), 2));
spectrum_45_dB = 10*log10(mean(spectrum_trial_45_avg(:, freq_7_12), 2));
relative_20_45_dB = spectrum_20_dB - spectrum_45_dB;

set(0,'DefaultFigureColormap', asymColorMapWhiteZero([-0.9, 0.9], N_colors_standard));
plot_topoplot(relative_20_45_dB, {'Heatmap for 7-12 Hz', 'Relative 20° vs 45°'}, file_electrode_positions, [-0.9, 0.9]);


% [4; 6.5] Hz
freq_4_6_5 = 7:14;
spectrum_20_dB = 10*log10(mean(spectrum_trial_20_avg(:, freq_4_6_5), 2));
spectrum_45_dB = 10*log10(mean(spectrum_trial_45_avg(:, freq_4_6_5), 2));
relative_20_45_dB = spectrum_20_dB - spectrum_45_dB;

plot_topoplot(relative_20_45_dB, {'Heatmap for 4-6.5 Hz', 'Relative 20° vs 45°'}, file_electrode_positions, [-0.9, 0.9]);





%% Helper functions

function [clustered_stats_table, statistical_clusters, stats_surrog, pairwise_stats, permutations] = compute_permutations(fov1, fov2, spectrum_trial1, spectrum_trial2, commonOptions)
    % do permutations analysis to compare conditions and baselines

    spectra_conditions = struct();
    spectra_conditions.BaselineModel = '';
    spectra_conditions.Chans = commonOptions.Chans;
    spectra_conditions.Freqs = commonOptions.Freqs;
    options = commonOptions;
    
    if isempty(fov2) % compare condition to baseline
        spectra_conditions.(fov1) = spectrum_trial1;
        spectra_conditions.(['Baseline', fov1]) = spectrum_trial2;  % Assuming baseline is the second input
        options.fields = {fov1, ['Baseline', fov1]};
    else % compare condition to condition (20 vs 45)
        spectra_conditions.(fov1) = spectrum_trial1;
        spectra_conditions.(fov2) = spectrum_trial2;
        options.fields = {fov1, fov2};
    end

    [clustered_stats_table, statistical_clusters, stats_surrog, pairwise_stats, permutations] = NP_statTest(spectra_conditions, options);
end

function plot_heatmap_baseline_or_condition(title_str, data_heatmap, electrode_labels, y, N_colors_standard, label)
    
    % colorLimits = [-10, 1]; % used to be colorLimits = [-1, 550]; for 20vs45 ?
    if contains(title_str, 'corrected')
        colorLimits = [-1.5, 1.5]; % For baseline comparisons
    else
        colorLimits = [-10, 10];    % For general comparisons
    end

    figure;
    myCmap = asymColorMapWhiteZero(colorLimits, N_colors_standard);
    heatmap_plot = heatmap(electrode_labels, y, data_heatmap', 'Colormap', myCmap, 'ColorLimits', colorLimits, 'ColorbarVisible', 'on', 'XLabel', 'Electrodes', 'YLabel', 'Frequencies [Hz]');
    heatmap_plot.Title = title_str;  % Title based on the type of comparison
    
    kept_frequencies = {'1', '2', '4', '7.5', '12', '30'};
    CustomYLabels = string(y);
    CustomYLabels(~ismember(CustomYLabels, kept_frequencies)) = " ";
    heatmap_plot.YDisplayLabels = CustomYLabels;
    grid off;
    origState = warning('query', 'MATLAB:structOnObject');
    cleanup = onCleanup(@() warning(origState));
    warning('off', 'MATLAB:structOnObject');
    S = struct(heatmap_plot); 
    ax = S.Axes; 
    clear('cleanup');
    
    row = [7, 14, 23, 59];
    col = [50, 63, 77, 107];
    arrayfun(@(x) xline(ax, x, 'k-', 'Alpha', 0.3), col - 0.25);
    arrayfun(@(x) yline(ax, x, 'k-', 'Alpha', 0.3), row - 0.25);

    saveas(gcf, fullfile(fullfile(fileparts(fileparts(pwd)), 'figures', 'AnalysisPlots'), "Heatmap" + label + ".png"));
end


function plot_spectrum_all(region, angle1, angle2, spectrum1, spectrum2, participants, EEG_trial_data, color1, color2, freqs_of_interest, x2, comparison_type)
    % Retrieve region of interest
    spectrum1_region_of_interest = select_frequencies_OI(spectrum1, region, {EEG_trial_data.(participants{end}).chanlocs(:).labels});
    spectrum2_region_of_interest = select_frequencies_OI(spectrum2, region, {EEG_trial_data.(participants{end}).chanlocs(:).labels});

    % Average over electrodes
    spectrum1_region_of_interest_averaged_electrodes = mean(spectrum1_region_of_interest, 1);
    spectrum2_region_of_interest_averaged_electrodes = mean(spectrum2_region_of_interest, 1);

    % Average and standard deviation over participants
    spectrum1_region_of_interest_averaged_subjects = mean(spectrum1_region_of_interest_averaged_electrodes, 3);
    spectrum1_region_of_interest_std_subjects = std(spectrum1_region_of_interest_averaged_electrodes, [], 3) / sqrt(size(spectrum1_region_of_interest_averaged_electrodes, 3));
    spectrum2_region_of_interest_averaged_subjects = mean(spectrum2_region_of_interest_averaged_electrodes, 3);
    spectrum2_region_of_interest_std_subjects = std(spectrum2_region_of_interest_averaged_electrodes, [], 3) / sqrt(size(spectrum2_region_of_interest_averaged_electrodes, 3));

    % Lines of standard deviation
    spectrum1_RoI_above = spectrum1_region_of_interest_averaged_subjects + spectrum1_region_of_interest_std_subjects / 2;
    spectrum1_RoI_below = spectrum1_region_of_interest_averaged_subjects - spectrum1_region_of_interest_std_subjects / 2;
    spectrum2_RoI_above = spectrum2_region_of_interest_averaged_subjects + spectrum2_region_of_interest_std_subjects / 2;
    spectrum2_RoI_below = spectrum2_region_of_interest_averaged_subjects - spectrum2_region_of_interest_std_subjects / 2;

    % Create fill function
    inBetween_spectrum1_RoI = [spectrum1_RoI_below(:); flipud(spectrum1_RoI_above(:))];
    inBetween_spectrum2_RoI = [spectrum2_RoI_below(:); flipud(spectrum2_RoI_above(:))];

    % Convert to dB
    spectrum1_region_of_interest_averaged_subjects_dB = 10 * log10(spectrum1_region_of_interest_averaged_subjects);
    spectrum2_region_of_interest_averaged_subjects_dB = 10 * log10(spectrum2_region_of_interest_averaged_subjects);
    inBetween_spectrum1_RoI_dB = 10 * log10(inBetween_spectrum1_RoI);
    inBetween_spectrum2_RoI_dB = 10 * log10(inBetween_spectrum2_RoI);

    % Make figure
    figure;
    plot(freqs_of_interest, spectrum1_region_of_interest_averaged_subjects_dB, 'Color', color1);
    hold on;
    plot(freqs_of_interest, spectrum2_region_of_interest_averaged_subjects_dB, 'Color', color2);
    patch('XData', x2, 'YData', inBetween_spectrum1_RoI_dB, 'FaceColor', color1, 'EdgeColor', color1, 'FaceAlpha', 0.2);
    patch('XData', x2, 'YData', inBetween_spectrum2_RoI_dB, 'FaceColor', color2, 'EdgeColor', color2, 'FaceAlpha', 0.2);

    if strcmp(comparison_type, 'trial_vs_baseline')
        title(['Spectrum over ' region ' Region',' (' num2str(angle1) '° vs. Baseline of ' num2str(angle1) '°)']);
        filename = "Spectrum" + num2str(angle1) + "vsBaseline_" + region;
        legend([num2str(angle1) '°'], ['Baseline ' num2str(angle1) '°']);
    elseif strcmp(comparison_type, 'trial_vs_trial')
        title(['Spectrum over ' region ' Region',' (' num2str(angle1) '° vs. ' num2str(angle2) '°)']);
        legend([num2str(angle1) '°'], [num2str(angle2) '°']);
        filename = "Spectrum" + num2str(angle1) + "vs" + num2str(angle2) + region;
    end

    xlabel('Frequencies [Hz]');
    ylabel('Power [dB]');
    grid on;
    saveas(gcf, fullfile(fullfile(fileparts(fileparts(pwd)), 'figures', 'AnalysisPlots'), filename + ".png"));

end


% COMMON FUNCTION FOR PLOTTING TOPOPLOTS
function plot_topoplot(data, title_text, file_electrode_positions, colormap_limits)
    figure;
    colorbar;
    
    if isempty(colormap_limits)
        topoplot(data, file_electrode_positions);
    else
        topoplot(data, file_electrode_positions, 'maplimits', colormap_limits);
    end
    title(title_text);
end