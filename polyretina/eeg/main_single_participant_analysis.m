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

subject_inds = [1, 2, 3, 6, 7, 8, 9];
%subject_inds = 9; % Overwrite subject for testing (COMMENT / DECOMMENT)

if (~exist(fullfile(output_filepath, 'EEG_trial_data.mat'),'file') || ~exist(fullfile(output_filepath, 'EEG_baseline_data.mat'),'file') || overwriteSpectraComputations)
    for subject_ind = subject_inds

        % 1. Load the data, filter the data
        subject = study_config.subjects(subject_ind).id;
        disp(['Subject ' subject]);
        study_config.current_subject = subject_ind;
        N = makeFolderFileNames(study_config, subject);
        EEG = pop_loadset('filename', N.postLabelingFile, 'filepath', N.searchFolder_2arch_rej_ICcats);
                
        % 2. Filter the dataset with the desired upper and lower frequencies
        % (definitive changes before epoching)
        lowcutoff = study_config.filterAnalysis.low_cut_off;
        highcutoff = study_config.filterAnalysis.high_cut_off;
        fprintf('Filtering between %.1f Hz and %.1f Hz...\n', lowcutoff, highcutoff)
        [EEG] = custom_filter(EEG, lowcutoff, highcutoff);

        % 3. "Epoch" Data: Removing events which are not of interest
        boundaryEvents_mask = strcmp({EEG.event.type}, 'boundary');
        EEG2 = pop_editeventvals(EEG, 'delete', find(boundaryEvents_mask));


        % 4. Retrieving Segments of Interest and Corresponding Baselines, and computing their Spectra
        
        % Defining whether the structures of interest already exist or not
        if ~exist('EEG_trial_data','var')
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







%% CHOOSING GLOBAL PARAMETERS FOR WHAT FOLLOWS:

% Which baseline to choose
choice_of_baseline = 'black'; % Choose: 'black'
choice_of_black_baseline = 'one_per_trial'; % Choose: 'one_per_trial'

% Defining which electrodes we are considering for the respective brain regions
frontal_electrodes ={'Z1','Z2','Z3','Z4','Z5','L1','L2','L3','L4','L5','LL1','LL2','LL3','LL4','LL5','R1','R2','R3','R4','R5','RR1','RR2','RR3','RR4','RR5'}; %{'LL3','LL4', 'L4', 'Z3', 'Z4', 'R4', 'RR3', 'RR4'};
%{'Z1','Z2','Z3','Z4','Z5','L1','L2','L3','L4','L5','LL1','LL2','LL3','LL4','LL5','R1','R2','R3','R4','R5','RR1','RR2','RR3','RR4','RR5'};
parietal_electrodes = {'LL6','LL7','LL8','RR6','RR7','RR8','RA3','RA4','LA3','LA4','L7','L8','R7','R8'}; %{'L8', 'R8', 'Z8', 'LL8', 'RR8'};
%{'LL6','LL7','LL8','RR6','RR7','RR8','RA3','RA4','LA3','LA4','L7','L8','R7','R8'};
occipital_electrodes = {'Z12', 'Z11','Z10','R11','L11','R12','L12'};

% Which Electrodes to consider for Brain Regions of Interest:
brain_region_chosen = frontal_electrodes; % 'occipital_electrodes', 'parietal_electrodes', 'frontal_electrodes'
brain_region_name = 'Frontal'; % 'Occipital', 'Parietal', 'Frontal'

% Decide which plots to make:
plot_scale = 'dB'; % Choose: 'linear', 'dB'
plot_absolute_spectrum = true;
plot_absolute_baseline_spectrum = true;
plot_std = true;

bool_divide_by_FoV = 1;

file_electrode_positions = strcat(study_config.study_folder,study_config.raw_data_folder,'P001\',study_config.channel_locations_filename);



%% 6. Distinguish FoV Conditions, and eventually Brain RoI and frequencies of interest

% 6.0. Dropping the computed amplitudes for frequencies which are not of
% interest
if ~exist('lowcutoff','var')
    lowcutoff = study_config.filterAnalysis.low_cut_off;
    highcutoff = study_config.filterAnalysis.high_cut_off;
end

% Removing frequencies which are filtered out
EEG_trial_data = drop_useless_frequencies(EEG_trial_data,lowcutoff, highcutoff, 'absolute');
EEG_baseline_data = drop_useless_frequencies(EEG_baseline_data,lowcutoff,highcutoff, 'absolute');
EEG_coarse_data = drop_useless_frequencies(EEG_coarse_data,lowcutoff,highcutoff, 'absolute');


% 6.2. Extract region of interest for each FoV (TRIALS / BASELINE)
[EEG_selected_absolute_spectrum_20] = extract_trials_according_to_brainregion_and_frequency(EEG_trial_data, brain_region_chosen, 'all', 1, 20, 'absolute');
[EEG_selected_absolute_spectrum_45] = extract_trials_according_to_brainregion_and_frequency(EEG_trial_data, brain_region_chosen, 'all', 1, 45, 'absolute');
[EEG_selected_absolute_base_spectrum_20] = extract_trials_according_to_brainregion_and_frequency(EEG_baseline_data, brain_region_chosen, 'all', 1, 20, 'absolute');
[EEG_selected_absolute_base_spectrum_45] = extract_trials_according_to_brainregion_and_frequency(EEG_baseline_data, brain_region_chosen, 'all', 1, 45, 'absolute');
[EEG_selected_absolute_coarse_spectrum_black] = extract_trials_according_to_brainregion_and_frequency(EEG_coarse_data, brain_region_chosen, 'all', 0, "black", 'absolute');
[EEG_selected_absolute_coarse_spectrum_edge] = extract_trials_according_to_brainregion_and_frequency(EEG_coarse_data, brain_region_chosen, 'all', 0, "edge", 'absolute');

% 6.3. Concatenating all of the data across the trials to format them for the plots:
[EEG_selected_absolute_spectrum_20_all_trials] = format_for_plotting_spectra(EEG_selected_absolute_spectrum_20); % TRIALS
[EEG_selected_absolute_spectrum_45_all_trials] = format_for_plotting_spectra(EEG_selected_absolute_spectrum_45);
[EEG_selected_absolute_base_spectrum_20_all_trials] = format_for_plotting_spectra(EEG_selected_absolute_base_spectrum_20); % BASELINE
[EEG_selected_absolute_base_spectrum_45_all_trials] = format_for_plotting_spectra(EEG_selected_absolute_base_spectrum_45);
[EEG_selected_absolute_coarse_spectrum_black_all_trials] = format_for_plotting_spectra(EEG_selected_absolute_coarse_spectrum_black); % COARSE
[EEG_selected_absolute_coarse_spectrum_edge_all_trials] = format_for_plotting_spectra(EEG_selected_absolute_coarse_spectrum_edge);

% 6.4. Averating across electrodes to get brain regions of interest:
EEG_absolute_spectrum_20_all_trials_averaged = mean(EEG_selected_absolute_spectrum_20_all_trials.spectrum,1); % TRIALS
EEG_absolute_spectrum_45_all_trials_averaged = mean(EEG_selected_absolute_spectrum_45_all_trials.spectrum,1);
EEG_absolute_base_spectrum_20_all_trials_averaged = mean(EEG_selected_absolute_base_spectrum_20_all_trials.spectrum,1); % BASELINE
EEG_absolute_base_spectrum_45_all_trials_averaged = mean(EEG_selected_absolute_base_spectrum_45_all_trials.spectrum,1);
EEG_absolute_coarse_spectrum_black_all_trials_averaged = mean(EEG_selected_absolute_coarse_spectrum_black_all_trials.spectrum,1); % COARSE
EEG_absolute_coarse_spectrum_edge_all_trials_averaged = mean(EEG_selected_absolute_coarse_spectrum_edge_all_trials.spectrum,1);

% 6.5. Averaging across trials (retrieving mean and standard deviation)
% TRIALS / BASELINE
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

% Fill in function to plot std (TRIALS / BASELINE)
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


% Converting into dB for the log plots (TRIALS, BASELINE, COARSE)
EEG_absolute_spectrum_20_averaged_dB = 10*log10(EEG_absolute_spectrum_20_averaged);
EEG_absolute_spectrum_45_averaged_dB = 10*log10(EEG_absolute_spectrum_45_averaged);
EEG_absolute_spectrum_20_std_dB = 10*log10(EEG_absolute_spectrum_20_std);
EEG_absolute_spectrum_45_std_dB = 10*log10(EEG_absolute_spectrum_45_std);
inBetween_absolute_20_dB = 10*log10(inBetween_absolute_20);
inBetween_absolute_45_dB = 10*log10(inBetween_absolute_45);

EEG_absolute_base_spectrum_20_averaged_dB = 10*log10(EEG_absolute_base_spectrum_20_averaged);
EEG_absolute_base_spectrum_45_averaged_dB = 10*log10(EEG_absolute_base_spectrum_45_averaged);
EEG_absolute_base_spectrum_20_std_dB = 10*log10(EEG_absolute_base_spectrum_20_std);
EEG_absolute_base_spectrum_45_std_dB = 10*log10(EEG_absolute_base_spectrum_45_std);
inBetween_absolute_baseline_20_dB = 10*log10(inBetween_absolute_baseline_20);
inBetween_absolute_baseline_45_dB = 10*log10(inBetween_absolute_baseline_45);

EEG_absolute_coarse_spectrum_black_averaged_dB = 10*log10(EEG_absolute_coarse_spectrum_black_averaged);
EEG_absolute_coarse_spectrum_edge_averaged_dB = 10*log10(EEG_absolute_coarse_spectrum_edge_averaged);
EEG_absolute_coarse_spectrum_black_std_dB = 10*log10(EEG_absolute_coarse_spectrum_black_std);
EEG_absolute_coarse_spectrum_edge_std_dB = 10*log10(EEG_absolute_coarse_spectrum_edge_std);
inBetween_absolute_coarse_black_dB = 10*log10(inBetween_absolute_coarse_black);
inBetween_absolute_coarse_edge_dB = 10*log10(inBetween_absolute_coarse_edge);


freqs_of_interest = EEG_selected_absolute_spectrum_20_all_trials.freqs;



% PLOT TEST FOR COARSE BASELINES

% 4. Define colors:
color_black = [0, 0, 0]; % Black color for black baseline
color_edge = [0.8, 0.4, 0.0]; % Edge color for edge baseline

% 5. Plot absolute spectra for black and edge baselines:
figure;
plot(freqs_of_interest, squeeze(EEG_absolute_coarse_spectrum_black_averaged_dB), 'Color', color_black, 'LineWidth', 2);
hold on;
plot(freqs_of_interest, squeeze(EEG_absolute_coarse_spectrum_edge_averaged_dB), 'Color', color_edge, 'LineWidth', 2);

% 6. Add standard deviation (optional):
if plot_std
    % Assuming `inBetween_absolute_coarse_black` and `inBetween_absolute_coarse_edge` are defined
    patch('XData', [freqs_of_interest(:); flipud(freqs_of_interest(:))], 'YData', inBetween_absolute_coarse_black_dB, ...
          'FaceColor', color_black, 'EdgeColor', color_black, 'FaceAlpha', 0.2);
    patch('XData', [freqs_of_interest(:); flipud(freqs_of_interest(:))], 'YData', inBetween_absolute_coarse_edge_dB, ...
          'FaceColor', color_edge, 'EdgeColor', color_edge, 'FaceAlpha', 0.2);
end

% 7. Finalize plot:
hold off;
grid on;
xlabel('Frequencies [Hz]');
ylabel('Power [dB]');
legend('Black Baseline', 'Edge Baseline');
title('Absolute Spectrum for Brain Region Electrodes Across FoV');







%% 7. Making plots for Absolute Spectra

% 7.1. Defining variables of interest for plots:
freqs_of_interest = EEG_selected_absolute_spectrum_20_all_trials.freqs;
x2 = [freqs_of_interest(:); flipud(freqs_of_interest(:))]; % needed for std plot

color_20 = [0.83 0.14 0.14];
color_45 = [1.00 0.54 0.00];
color_baseline = 1/255 * [0, 104, 87];

if strcmp(plot_scale, 'dB')
    % 7.2. Plotting absolute spectra separately for each FoV, averaged over all trials and all electrodes of brain RoI
    if plot_absolute_spectrum
        figure;
        plot(freqs_of_interest, EEG_absolute_spectrum_20_averaged_dB', 'Color', color_20, 'LineWidth', 2);
        hold on;
        plot(freqs_of_interest, EEG_absolute_spectrum_45_averaged_dB', 'Color', color_45, 'LineWidth', 2);
        if plot_std
           patch('XData',x2,'YData',inBetween_absolute_20_dB,'FaceColor', color_20,'EdgeColor',color_20,'FaceAlpha', 0.2);
           patch('XData',x2,'YData',inBetween_absolute_45_dB,'FaceColor', color_45,'EdgeColor',color_45,'FaceAlpha', 0.2);
        end
        hold off;
        grid on;
        xlabel('Frequencies [Hz]');
        ylabel('Power [dB]');
        legend('20°','45°');
        title({['Absolute Spectrum for ' brain_region_name ' Electrodes' ],'Across FoV'});

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
    end

    % 7.3. Plotting absolute spectra OF BASELINE for each FoV, averaged over all trials and all electrodes of brain RoI
    if strcmp(choice_of_baseline, 'black')
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
    end
end
    


%% SHORT DURATION ANALYSIS (LAST 2 SECONDS OF TRIALS)
% Instead of looking at the whole trials durations, we only extract the
% data from the last 2 seconds of which trials to see if different patterns
% show (as cognitive load could reach a maximum towards the end of the
% trials). Every trial is normalised by its corresponding baseline

bool_all_electrodes = 1; bool_plot = 0; bool_export = 0;

EEG_trial_data2sec = compute_temporal_trials_corrected(EEG_trial_data, EEG_baseline_data, ...
    brain_region_chosen, bool_plot, bool_all_electrodes, bool_export);






















































%% 9. Permutation Analysis -- Multi Subject Analysis

%% Part 1: Condition VS. Baseline

% Defining Parameters for plots
N_colors_standard = 512; % colormap
y = linspace(1, 42, 83); % for plot
oragnize_alphabetically_electrodes = 1;

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

% 9.3.1 Formating Data for Heatmaps
data_heatmap_20vbase = format_for_heatmap_conditionvbaseline(clustered_stats_table_20vbase, statistical_clusters_20vbase, spectrum_trial_20, spectrum_baseline_20);
data_heatmap_45vbase = format_for_heatmap_conditionvbaseline(clustered_stats_table_45vbase, statistical_clusters_45vbase, spectrum_trial_45, spectrum_baseline_45);

[data_heatmap_20vbase, ~] = organize_by_electrodes(data_heatmap_20vbase, {EEG_trial_data.(participants{end}).chanlocs(:).labels}, oragnize_alphabetically_electrodes);
[data_heatmap_45vbase, new_electrode_labels] = organize_by_electrodes(data_heatmap_45vbase, {EEG_trial_data.(participants{end}).chanlocs(:).labels}, oragnize_alphabetically_electrodes);

plot_heatmap_baseline_or_condition('20° FoV vs. Baseline', data_heatmap_20vbase, new_electrode_labels, y, N_colors_standard);
plot_heatmap_baseline_or_condition('45° FoV vs. Baseline', data_heatmap_45vbase, new_electrode_labels, y, N_colors_standard);



%% Part 2: Condition VS. Condition (20 vs 45)

if ~exist('clustered_stats_table_conditionvcondition','var')
    [clustered_stats_table_conditionvcondition, statistical_clusters_conditionvcondition, stats_surrog_conditionvcondition, pairwise_stats_conditionvcondition, permutations_conditionvcondition] = ...
    compute_permutations('FoV20', 'FoV45', spectrum_trial_20, spectrum_trial_45, commonOptions);
end

data_heatmap_conditionvcondition = format_for_heatmap_conditionvcondition_anova(clustered_stats_table_conditionvcondition, statistical_clusters_conditionvcondition);

[data_heatmap_conditionvcondition, new_electrode_labels] = organize_by_electrodes(data_heatmap_conditionvcondition, {EEG_trial_data.(participants{end}).chanlocs(:).labels}, 1);

plot_heatmap_baseline_or_condition('20° FoV vs. 45° FoV', data_heatmap_conditionvcondition, new_electrode_labels, y, N_colors_standard);












%% 11. Comparing baseline-corrected signals

% 11.1. Correcting the signals
    
% Given that the spectrums are considered in decibel scale,
% the signals are substracted

% Converting initial spectrum into logarithmic scale
spectrum_trial_20_dB = 10*log10(spectrum_trial_20);
spectrum_trial_45_dB = 10*log10(spectrum_trial_45);

% PAUL
spectrum_trial_20_base_dB = 10*log10(spectrum_baseline_20);
spectrum_trial_45_base_dB = 10*log10(spectrum_baseline_45);

spectrum_trial_relative_20vbase_dB = spectrum_trial_20_dB - spectrum_trial_20_base_dB;
spectrum_trial_relative_45vbase_dB =  spectrum_trial_45_dB - spectrum_trial_45_base_dB;

if ~exist('clustered_stats_table_baseline_corrected','var')
    [clustered_stats_table_baseline_corrected, statistical_clusters_baselinecorrected, stats_surrog_baselinecorrected, pairwise_stats_baselinecorrected, permutations_baselinecorrected] = ...
    compute_permutations('FoV20vbase', 'FoV45vbase', spectrum_trial_relative_20vbase_dB, spectrum_trial_relative_45vbase_dB, commonOptions);
end

data_heatmap_baselinecorrected = format_for_heatmap_baselinecorrected_dB(clustered_stats_table_baseline_corrected, statistical_clusters_baselinecorrected, spectrum_trial_relative_20vbase_dB, spectrum_trial_relative_45vbase_dB);

[data_heatmap_baselinecorrected, new_electrode_labels] = organize_by_electrodes(data_heatmap_baselinecorrected, {EEG_trial_data.(participants{end}).chanlocs(:).labels}, oragnize_alphabetically_electrodes);

plot_heatmap_baseline_or_condition('20° corrected with 20° baseline VS. 45° corrected with 45° baseline', data_heatmap_baselinecorrected, new_electrode_labels, y, N_colors_standard);





















%% 12. Making Additional Plots to Accompany Heatmaps

% 11.1. Condition v Baseline

plot_spectrum_all('parietal', 20, [], spectrum_trial_20, spectrum_baseline_20, participants, EEG_trial_data, color_20, color_baseline, freqs_of_interest, x2, 'trial_vs_baseline');
plot_spectrum_all('occipital', 20, [], spectrum_trial_20, spectrum_baseline_20, participants, EEG_trial_data, color_20, color_baseline, freqs_of_interest, x2, 'trial_vs_baseline');
plot_spectrum_all('frontal', 20, [], spectrum_trial_20, spectrum_baseline_20, participants, EEG_trial_data, color_20, color_baseline, freqs_of_interest, x2, 'trial_vs_baseline');

plot_spectrum_all('parietal', 45, [], spectrum_trial_45, spectrum_baseline_45, participants, EEG_trial_data, color_45, color_baseline, freqs_of_interest, x2, 'trial_vs_baseline');
plot_spectrum_all('occipital', 45, [], spectrum_trial_45, spectrum_baseline_45, participants, EEG_trial_data, color_45, color_baseline, freqs_of_interest, x2, 'trial_vs_baseline');
plot_spectrum_all('frontal', 45, [], spectrum_trial_45, spectrum_baseline_45, participants, EEG_trial_data, color_45, color_baseline, freqs_of_interest, x2, 'trial_vs_baseline');

plot_spectrum_all('parietal', 20, 45, spectrum_trial_20, spectrum_trial_45, participants, EEG_trial_data, color_20, color_45, freqs_of_interest, x2, 'trial_vs_trial');
plot_spectrum_all('occipital', 20, 45, spectrum_trial_20, spectrum_trial_45, participants, EEG_trial_data, color_20, color_45, freqs_of_interest, x2, 'trial_vs_trial');
plot_spectrum_all('frontal', 20, 45, spectrum_trial_20, spectrum_trial_45, participants, EEG_trial_data, color_20, color_45, freqs_of_interest, x2, 'trial_vs_trial');






    
% Fig 1: Topoplot over Frequencies [20; 35] Hz
    % Average over the freqs of interest
% spectrum_trial_20_averaged_trials_and_freqs= mean(spectrum_trial_20_averaged_trials(:,43:55), 2);
% spectrum_trial_45_averaged_trials_and_freqs= mean(spectrum_trial_45_averaged_trials(:,43:55), 2);
%     % Convert in dB
% spectrum_trial_20_averaged_trials_and_freqs_dB = 10*log10(spectrum_trial_20_averaged_trials_and_freqs);
% spectrum_trial_45_averaged_trials_and_freqs_dB = 10*log10(spectrum_trial_45_averaged_trials_and_freqs);
% relative_spectrum_trialvbaseline_20v45_dB = spectrum_trial_20_averaged_trials_and_freqs_dB - spectrum_trial_45_averaged_trials_and_freqs_dB;
%     % Topoplot
% figure;
% topoplot(spectrum_trial_20_averaged_trials_and_freqs_dB, file_electrode_positions);
% colorbar;
% title({'Heatmap for 20-35 Hz','for trials under FoV 20°'});
% figure;
% colorbar;
% topoplot(spectrum_trial_45_averaged_trials_and_freqs_dB, file_electrode_positions);
% title({'Heatmap for 20-35 Hz','for trials under FoV 45°'});
% figure;
% colorbar;
% topoplot(relative_spectrum_trialvbaseline_20v45_dB, file_electrode_positions);
% title({'Heatmap for 20-35 Hz','for relative spectrum 20° vs 45°'});


plot_illustrative_conclusion_plots = true;

if plot_illustrative_conclusion_plots

    % PART 1: CONDITION VS BASELINE
    % Fig 1: Topoplots for 20vB20, 45vB45
    % frequencies
        % Average over all the subjects
    spectrum_trial_20_averaged_trials = mean(spectrum_trial_20, 3);
    spectrum_baseline_20_averaged_trials = mean(spectrum_baseline_20,3);
    spectrum_trial_45_averaged_trials = mean(spectrum_trial_45, 3);
    spectrum_baseline_45_averaged_trials = mean(spectrum_baseline_45,3);

        % Average over the freqs of interest
    spectrum_trial_20_averaged_trials_and_freqs= mean(spectrum_trial_20_averaged_trials(:,20:28), 2);
    spectrum_baseline_20_averaged_trials_and_freqs= mean(spectrum_baseline_20_averaged_trials(:,20:28), 2);
    spectrum_trial_45_averaged_trials_and_freqs= mean(spectrum_trial_45_averaged_trials(:,20:28), 2);
    spectrum_baseline_45_averaged_trials_and_freqs= mean(spectrum_baseline_45_averaged_trials(:,20:28), 2);

        % Convert in dB
    spectrum_trial_20_averaged_trials_and_freqs_dB = 10*log10(spectrum_trial_20_averaged_trials_and_freqs);
    spectrum_baseline_20_averaged_trials_and_freqs_dB = 10*log10(spectrum_baseline_20_averaged_trials_and_freqs);
    relative_spectrum_trialvbaseline_20_dB = spectrum_trial_20_averaged_trials_and_freqs_dB - spectrum_baseline_20_averaged_trials_and_freqs_dB;
    spectrum_trial_45_averaged_trials_and_freqs_dB = 10*log10(spectrum_trial_45_averaged_trials_and_freqs);
    spectrum_baseline_45_averaged_trials_and_freqs_dB = 10*log10(spectrum_baseline_45_averaged_trials_and_freqs);
    relative_spectrum_trialvbaseline_45_dB = spectrum_trial_45_averaged_trials_and_freqs_dB - spectrum_baseline_45_averaged_trials_and_freqs_dB;

        % Topoplot
    % Set color map to values of realtive spectrum
    myCmap = asymColorMapWhiteZero([-5,5], N_colors_standard);
    set(0,'DefaultFigureColormap',myCmap);
    
    figure;
    colorbar;
    topoplot(relative_spectrum_trialvbaseline_20_dB, file_electrode_positions, 'colormap', myCmap, 'maplimits', [-5;5]);
    title({'Topoplot for 10-14 Hz','for relative spectrum under FoV 20°'});
    figure;
    colorbar;
    topoplot(relative_spectrum_trialvbaseline_45_dB, file_electrode_positions, 'colormap', myCmap, 'maplimits', [-5;5]);
    title({'Topoplot for 10-14 Hz','for relative spectrum under FoV 45°'});

    % Fig 2: Comparing over occipital electrodes the spectrum
    % Spectrum over Occipital for 20 v B20
        % Retrieving over the brain Region of Interest: Occipital
    spectrum_trial_20_region_of_interest = select_frequencies_OI(spectrum_trial_20, 'occipital', {EEG_trial_data.(participants{end}).chanlocs(:).labels});
    spectrum_baseline_20_region_of_interest = select_frequencies_OI(spectrum_baseline_20, 'occipital', {EEG_trial_data.(participants{end}).chanlocs(:).labels});
        % Get Average over electrodes of interest
    spectrum_trial_20_region_of_interest_averaged_electrodes = mean(spectrum_trial_20_region_of_interest, 1);
    spectrum_baseline_20_region_of_interest_averaged_electrodes = mean(spectrum_baseline_20_region_of_interest,1);
        % Get Average and Standard Deviation over Participants
    spectrum_trial_20_region_of_interest_averaged_subjects = mean(spectrum_trial_20_region_of_interest_averaged_electrodes,3);
    spectrum_trial_20_region_of_interest_std_subjects = std(spectrum_trial_20_region_of_interest_averaged_electrodes,[],3)/sqrt(size(spectrum_trial_20_region_of_interest_averaged_electrodes,3));
    spectrum_baseline_20_region_of_interest_averaged_subjects = mean(spectrum_baseline_20_region_of_interest_averaged_electrodes,3);
    spectrum_baseline_20_region_of_interest_std_subjects = std(spectrum_baseline_20_region_of_interest_averaged_electrodes,[],3)/sqrt(size(spectrum_baseline_20_region_of_interest_averaged_electrodes,3));
        % Getting lines of standard deviation above and under for plot
            % Step 1: Create Lines
    spectrum_trial_20_RoI_above = spectrum_trial_20_region_of_interest_averaged_subjects + spectrum_trial_20_region_of_interest_std_subjects./2;
    spectrum_trial_20_RoI_below = spectrum_trial_20_region_of_interest_averaged_subjects - spectrum_trial_20_region_of_interest_std_subjects./2;
    spectrum_baseline_20_RoI_above = spectrum_baseline_20_region_of_interest_averaged_subjects + spectrum_baseline_20_region_of_interest_std_subjects./2;
    spectrum_baseline_20_RoI_below = spectrum_baseline_20_region_of_interest_averaged_subjects - spectrum_baseline_20_region_of_interest_std_subjects./2;
            % Step 2: Create fill function
    inBetween_spectrum_20_RoI = [spectrum_trial_20_RoI_below(:); flipud(spectrum_trial_20_RoI_above(:))];
    inBetween_spectrum_baseline_20_RoI = [spectrum_baseline_20_RoI_below(:); flipud(spectrum_baseline_20_RoI_above(:))];
        % Converting everything to dB
    spectrum_trial_20_region_of_interest_averaged_subjects_dB = 10*log10(spectrum_trial_20_region_of_interest_averaged_subjects);
    spectrum_baseline_20_region_of_interest_averaged_subjects_dB = 10*log10(spectrum_baseline_20_region_of_interest_averaged_subjects);
    inBetween_spectrum_20_RoI_dB = 10*log10(inBetween_spectrum_20_RoI);
    inBetween_spectrum_baseline_20_RoI_dB = 10*log10(inBetween_spectrum_baseline_20_RoI);
        % Making Figure
    figure;
    plot(freqs_of_interest, spectrum_trial_20_region_of_interest_averaged_subjects_dB,'Color', color_20);
    hold on;
    plot(freqs_of_interest, spectrum_baseline_20_region_of_interest_averaged_subjects_dB, 'Color', color_baseline);
    patch('XData',x2,'YData',inBetween_spectrum_20_RoI_dB,'FaceColor', color_20,'EdgeColor',color_20,'FaceAlpha', 0.2);
    patch('XData',x2,'YData',inBetween_spectrum_baseline_20_RoI_dB,'FaceColor', color_baseline,'EdgeColor',color_baseline,'FaceAlpha', 0.2);
    title({'Spectrum over Occipital Region','(20° vs. Baseline of 20°)'});
    legend('20°', 'Baseline 20°');
    xlabel('Frequencies [Hz]');
    ylabel('Power [dB]');
    grid on;
    
    % Spectrum over Occipital for 45 v B45
        % Retrieving over the brain Region of Interest: Occipital
    spectrum_trial_45_region_of_interest = select_frequencies_OI(spectrum_trial_45, 'occipital', {EEG_trial_data.(participants{end}).chanlocs(:).labels});
    spectrum_baseline_45_region_of_interest = select_frequencies_OI(spectrum_baseline_45, 'occipital', {EEG_trial_data.(participants{end}).chanlocs(:).labels});
        % Get Average over electrodes of interest
    spectrum_trial_45_region_of_interest_averaged_electrodes = mean(spectrum_trial_45_region_of_interest, 1);
    spectrum_baseline_45_region_of_interest_averaged_electrodes = mean(spectrum_baseline_45_region_of_interest,1);
        % Get Average and Standard Deviation over Participants
    spectrum_trial_45_region_of_interest_averaged_subjects = mean(spectrum_trial_45_region_of_interest_averaged_electrodes,3);
    spectrum_trial_45_region_of_interest_std_subjects = std(spectrum_trial_45_region_of_interest_averaged_electrodes,[],3)/sqrt(size(spectrum_trial_45_region_of_interest_averaged_electrodes,3));
    spectrum_baseline_45_region_of_interest_averaged_subjects = mean(spectrum_baseline_45_region_of_interest_averaged_electrodes,3);
    spectrum_baseline_45_region_of_interest_std_subjects = std(spectrum_baseline_45_region_of_interest_averaged_electrodes,[],3)/sqrt(size(spectrum_baseline_45_region_of_interest_averaged_electrodes,3));
        % Getting lines of standard deviation above and under for plot
            % Step 1: Create Lines
    spectrum_trial_45_RoI_above = spectrum_trial_45_region_of_interest_averaged_subjects + spectrum_trial_45_region_of_interest_std_subjects./2;
    spectrum_trial_45_RoI_below = spectrum_trial_45_region_of_interest_averaged_subjects - spectrum_trial_45_region_of_interest_std_subjects./2;
    spectrum_baseline_45_RoI_above = spectrum_baseline_45_region_of_interest_averaged_subjects + spectrum_baseline_45_region_of_interest_std_subjects./2;
    spectrum_baseline_45_RoI_below = spectrum_baseline_45_region_of_interest_averaged_subjects - spectrum_baseline_45_region_of_interest_std_subjects./2;
            % Step 2: Create fill function
    inBetween_spectrum_45_RoI = [spectrum_trial_45_RoI_below(:); flipud(spectrum_trial_45_RoI_above(:))];
    inBetween_spectrum_baseline_45_RoI = [spectrum_baseline_45_RoI_below(:); flipud(spectrum_baseline_45_RoI_above(:))];
        % Converting everything to dB
    spectrum_trial_45_region_of_interest_averaged_subjects_dB = 10*log10(spectrum_trial_45_region_of_interest_averaged_subjects);
    spectrum_baseline_45_region_of_interest_averaged_subjects_dB = 10*log10(spectrum_baseline_45_region_of_interest_averaged_subjects);
    inBetween_spectrum_45_RoI_dB = 10*log10(inBetween_spectrum_45_RoI);
    inBetween_spectrum_baseline_45_RoI_dB = 10*log10(inBetween_spectrum_baseline_45_RoI);
        % Making Figure
    figure;
    plot(freqs_of_interest, spectrum_trial_45_region_of_interest_averaged_subjects_dB,'Color', color_45);
    hold on;
    plot(freqs_of_interest, spectrum_baseline_45_region_of_interest_averaged_subjects_dB, 'Color', color_baseline);
    patch('XData',x2,'YData',inBetween_spectrum_45_RoI_dB,'FaceColor', color_45,'EdgeColor',color_45,'FaceAlpha', 0.2);
    patch('XData',x2,'YData',inBetween_spectrum_baseline_45_RoI_dB,'FaceColor', color_baseline,'EdgeColor',color_baseline,'FaceAlpha', 0.2);
    title({'Spectrum over Occipital Region','(45° vs. Baseline of 45°)'});
    legend('45°', 'Baseline 45°');
    xlabel('Frequencies [Hz]');
    ylabel('Power [dB]');
    grid on;

    % For 20 v 45
        % Retrieving over the brain Region of Interest: Frontal
    spectrum_trial_20_region_of_interest = select_frequencies_OI(spectrum_trial_20, 'specific', {EEG_trial_data.(participants{end}).chanlocs(:).labels}, specific_frontal_electrodes);
    spectrum_trial_45_region_of_interest = select_frequencies_OI(spectrum_trial_45, 'specific', {EEG_trial_data.(participants{end}).chanlocs(:).labels}, specific_frontal_electrodes);
        % Get Average over electrodes of interest
    spectrum_trial_20_region_of_interest_averaged_electrodes = mean(spectrum_trial_20_region_of_interest, 1);
    spectrum_trial_45_region_of_interest_averaged_electrodes = mean(spectrum_trial_45_region_of_interest, 1);
        % Get Average and Standard Deviation over Participants
    spectrum_trial_20_region_of_interest_averaged_subjects = mean(spectrum_trial_20_region_of_interest_averaged_electrodes,3);
    spectrum_trial_20_region_of_interest_std_subjects = std(spectrum_trial_20_region_of_interest_averaged_electrodes,[],3)/sqrt(size(spectrum_trial_20_region_of_interest_averaged_electrodes,3));
    spectrum_trial_45_region_of_interest_averaged_subjects = mean(spectrum_trial_45_region_of_interest_averaged_electrodes,3);
    spectrum_trial_45_region_of_interest_std_subjects = std(spectrum_trial_45_region_of_interest_averaged_electrodes,[],3)/sqrt(size(spectrum_trial_45_region_of_interest_averaged_electrodes,3));
        % Getting lines of standard deviation above and under for plot
            % Step 1: Create Lines
    spectrum_trial_20_RoI_above = spectrum_trial_20_region_of_interest_averaged_subjects + spectrum_trial_20_region_of_interest_std_subjects./2;
    spectrum_trial_20_RoI_below = spectrum_trial_20_region_of_interest_averaged_subjects - spectrum_trial_20_region_of_interest_std_subjects./2;
    spectrum_trial_45_RoI_above = spectrum_trial_45_region_of_interest_averaged_subjects + spectrum_trial_45_region_of_interest_std_subjects./2;
    spectrum_trial_45_RoI_below = spectrum_trial_45_region_of_interest_averaged_subjects - spectrum_trial_45_region_of_interest_std_subjects./2;
            % Step 2: Create fill function
    inBetween_spectrum_20_RoI = [spectrum_trial_20_RoI_below(:); flipud(spectrum_trial_20_RoI_above(:))];
    inBetween_spectrum_45_RoI = [spectrum_trial_45_RoI_below(:); flipud(spectrum_trial_45_RoI_above(:))];
        % Converting everything to dB
    spectrum_trial_20_region_of_interest_averaged_subjects_dB = 10*log10(spectrum_trial_20_region_of_interest_averaged_subjects);
    spectrum_trial_45_region_of_interest_averaged_subjects_dB = 10*log10(spectrum_trial_45_region_of_interest_averaged_subjects);
    inBetween_spectrum_20_RoI_dB = 10*log10(inBetween_spectrum_20_RoI);
    inBetween_spectrum_45_RoI_dB = 10*log10(inBetween_spectrum_45_RoI);
        % Making Figure
    figure;
    plot(freqs_of_interest, spectrum_trial_20_region_of_interest_averaged_subjects_dB,'Color', color_20);
    hold on;
    plot(freqs_of_interest, spectrum_trial_45_region_of_interest_averaged_subjects_dB, 'Color', color_45);
    patch('XData',x2,'YData',inBetween_spectrum_20_RoI_dB,'FaceColor', color_20,'EdgeColor',color_20,'FaceAlpha', 0.2);
    patch('XData',x2,'YData',inBetween_spectrum_45_RoI_dB,'FaceColor', color_45,'EdgeColor',color_45,'FaceAlpha', 0.2);
    title({'Spectrum over Frontal Region','(20° vs. 45°)'});
    legend('20°', '45°');
    xlabel('Frequencies [Hz]');
    ylabel('Power [dB]');
    grid on;
    
    % Fig 4: Spectrum for Parietal Electrodes
    specific_parietal_electrodes = {'R10', 'L10', 'R9', 'RR8', 'Z9'};

    % For 20 v 45
        % Retrieving over the brain Region of Interest: Parietal
    spectrum_trial_20_region_of_interest = select_frequencies_OI(spectrum_trial_20, 'specific', {EEG_trial_data.(participants{end}).chanlocs(:).labels}, specific_parietal_electrodes);
    spectrum_trial_45_region_of_interest = select_frequencies_OI(spectrum_trial_45, 'specific', {EEG_trial_data.(participants{end}).chanlocs(:).labels}, specific_parietal_electrodes);
        % Get Average over electrodes of interest
    spectrum_trial_20_region_of_interest_averaged_electrodes = mean(spectrum_trial_20_region_of_interest, 1);
    spectrum_trial_45_region_of_interest_averaged_electrodes = mean(spectrum_trial_45_region_of_interest, 1);
        % Get Average and Standard Deviation over Participants
    spectrum_trial_20_region_of_interest_averaged_subjects = mean(spectrum_trial_20_region_of_interest_averaged_electrodes,3);
    spectrum_trial_20_region_of_interest_std_subjects = std(spectrum_trial_20_region_of_interest_averaged_electrodes,[],3)/sqrt(size(spectrum_trial_20_region_of_interest_averaged_electrodes,3));
    spectrum_trial_45_region_of_interest_averaged_subjects = mean(spectrum_trial_45_region_of_interest_averaged_electrodes,3);
    spectrum_trial_45_region_of_interest_std_subjects = std(spectrum_trial_45_region_of_interest_averaged_electrodes,[],3)/sqrt(size(spectrum_trial_45_region_of_interest_averaged_electrodes,3));
        % Getting lines of standard deviation above and under for plot
            % Step 1: Create Lines
    spectrum_trial_20_RoI_above = spectrum_trial_20_region_of_interest_averaged_subjects + spectrum_trial_20_region_of_interest_std_subjects./2;
    spectrum_trial_20_RoI_below = spectrum_trial_20_region_of_interest_averaged_subjects - spectrum_trial_20_region_of_interest_std_subjects./2;
    spectrum_trial_45_RoI_above = spectrum_trial_45_region_of_interest_averaged_subjects + spectrum_trial_45_region_of_interest_std_subjects./2;
    spectrum_trial_45_RoI_below = spectrum_trial_45_region_of_interest_averaged_subjects - spectrum_trial_45_region_of_interest_std_subjects./2;
            % Step 2: Create fill function
    inBetween_spectrum_20_RoI = [spectrum_trial_20_RoI_below(:); flipud(spectrum_trial_20_RoI_above(:))];
    inBetween_spectrum_45_RoI = [spectrum_trial_45_RoI_below(:); flipud(spectrum_trial_45_RoI_above(:))];
        % Converting everything to dB
    spectrum_trial_20_region_of_interest_averaged_subjects_dB = 10*log10(spectrum_trial_20_region_of_interest_averaged_subjects);
    spectrum_trial_45_region_of_interest_averaged_subjects_dB = 10*log10(spectrum_trial_45_region_of_interest_averaged_subjects);
    inBetween_spectrum_20_RoI_dB = 10*log10(inBetween_spectrum_20_RoI);
    inBetween_spectrum_45_RoI_dB = 10*log10(inBetween_spectrum_45_RoI);
        % Making Figure
    figure;
    plot(freqs_of_interest, spectrum_trial_20_region_of_interest_averaged_subjects_dB,'Color', color_20);
    hold on;
    plot(freqs_of_interest, spectrum_trial_45_region_of_interest_averaged_subjects_dB, 'Color', color_45);
    patch('XData',x2,'YData',inBetween_spectrum_20_RoI_dB,'FaceColor', color_20,'EdgeColor',color_20,'FaceAlpha', 0.2);
    patch('XData',x2,'YData',inBetween_spectrum_45_RoI_dB,'FaceColor', color_45,'EdgeColor',color_45,'FaceAlpha', 0.2);
    title({'Spectrum over Parietal Region','(20° vs. 45°)'});
    legend('20°', '45°');
    xlabel('Frequencies [Hz]');
    ylabel('Power [dB]');
    grid on;
    
    % Fig 5: Spectrum for Occipital Electrodes

    % For 20 v 45
        % Retrieving over the brain Region of Interest: Occipital
    spectrum_trial_20_region_of_interest = select_frequencies_OI(spectrum_trial_20, 'occipital', {EEG_trial_data.(participants{end}).chanlocs(:).labels});
    spectrum_trial_45_region_of_interest = select_frequencies_OI(spectrum_trial_45, 'occipital', {EEG_trial_data.(participants{end}).chanlocs(:).labels});
        % Get Average over electrodes of interest
    spectrum_trial_20_region_of_interest_averaged_electrodes = mean(spectrum_trial_20_region_of_interest, 1);
    spectrum_trial_45_region_of_interest_averaged_electrodes = mean(spectrum_trial_45_region_of_interest, 1);
        % Get Average and Standard Deviation over Participants
    spectrum_trial_20_region_of_interest_averaged_subjects = mean(spectrum_trial_20_region_of_interest_averaged_electrodes,3);
    spectrum_trial_20_region_of_interest_std_subjects = std(spectrum_trial_20_region_of_interest_averaged_electrodes,[],3)/sqrt(size(spectrum_trial_20_region_of_interest_averaged_electrodes,3));
    spectrum_trial_45_region_of_interest_averaged_subjects = mean(spectrum_trial_45_region_of_interest_averaged_electrodes,3);
    spectrum_trial_45_region_of_interest_std_subjects = std(spectrum_trial_45_region_of_interest_averaged_electrodes,[],3)/sqrt(size(spectrum_trial_45_region_of_interest_averaged_electrodes,3));
        % Getting lines of standard deviation above and under for plot
            % Step 1: Create Lines
    spectrum_trial_20_RoI_above = spectrum_trial_20_region_of_interest_averaged_subjects + spectrum_trial_20_region_of_interest_std_subjects./2;
    spectrum_trial_20_RoI_below = spectrum_trial_20_region_of_interest_averaged_subjects - spectrum_trial_20_region_of_interest_std_subjects./2;
    spectrum_trial_45_RoI_above = spectrum_trial_45_region_of_interest_averaged_subjects + spectrum_trial_45_region_of_interest_std_subjects./2;
    spectrum_trial_45_RoI_below = spectrum_trial_45_region_of_interest_averaged_subjects - spectrum_trial_45_region_of_interest_std_subjects./2;
            % Step 2: Create fill function
    inBetween_spectrum_20_RoI = [spectrum_trial_20_RoI_below(:); flipud(spectrum_trial_20_RoI_above(:))];
    inBetween_spectrum_45_RoI = [spectrum_trial_45_RoI_below(:); flipud(spectrum_trial_45_RoI_above(:))];
        % Converting everything to dB
    spectrum_trial_20_region_of_interest_averaged_subjects_dB = 10*log10(spectrum_trial_20_region_of_interest_averaged_subjects);
    spectrum_trial_45_region_of_interest_averaged_subjects_dB = 10*log10(spectrum_trial_45_region_of_interest_averaged_subjects);
    inBetween_spectrum_20_RoI_dB = 10*log10(inBetween_spectrum_20_RoI);
    inBetween_spectrum_45_RoI_dB = 10*log10(inBetween_spectrum_45_RoI);
        % Making Figure
    figure;
    plot(freqs_of_interest, spectrum_trial_20_region_of_interest_averaged_subjects_dB,'Color', color_20);
    hold on;
    plot(freqs_of_interest, spectrum_trial_45_region_of_interest_averaged_subjects_dB, 'Color', color_45);
    patch('XData',x2,'YData',inBetween_spectrum_20_RoI_dB,'FaceColor', color_20,'EdgeColor',color_20,'FaceAlpha', 0.2);
    patch('XData',x2,'YData',inBetween_spectrum_45_RoI_dB,'FaceColor', color_45,'EdgeColor',color_45,'FaceAlpha', 0.2);
    title({'Spectrum over Occipital Region','(20° vs. 45°)'});
    legend('20°', '45°');
    xlabel('Frequencies [Hz]');
    ylabel('Power [dB]');
    grid on;

    % Fig 6: Topoplot for frequencies [7; 12]

    % 20v45
        % Average over the freqs of interest
    spectrum_trial_20_averaged_trials_and_freqs= mean(spectrum_trial_20_averaged_trials(:,13:23), 2);
    spectrum_trial_45_averaged_trials_and_freqs= mean(spectrum_trial_45_averaged_trials(:,13:23), 2);
        % Convert in dB
    spectrum_trial_20_averaged_trials_and_freqs_dB = 10*log10(spectrum_trial_20_averaged_trials_and_freqs);
    spectrum_trial_45_averaged_trials_and_freqs_dB = 10*log10(spectrum_trial_45_averaged_trials_and_freqs);
    relative_spectrum_trialvbaseline_20v45_dB = spectrum_trial_20_averaged_trials_and_freqs_dB - spectrum_trial_45_averaged_trials_and_freqs_dB;
        % Topoplot
    myCmap = asymColorMapWhiteZero([-0.9,0.9], N_colors_standard);
    set(0,'DefaultFigureColormap',myCmap);
    figure;
    colorbar;
    topoplot(relative_spectrum_trialvbaseline_20v45_dB, file_electrode_positions, 'colormap', myCmap, 'maplimits', [-0.9;0.9]);
    title({'Heatmap for 7-12 Hz','for relative spectrum 20° vs 45°'});
    
    % Fig 7: Topoplot for frequencies [4; 6.5]

    % 20v45
        % Average over the freqs of interest
    spectrum_trial_20_averaged_trials_and_freqs= mean(spectrum_trial_20_averaged_trials(:,7:14), 2);
    spectrum_trial_45_averaged_trials_and_freqs= mean(spectrum_trial_45_averaged_trials(:,7:14), 2);
        % Convert in dB
    spectrum_trial_20_averaged_trials_and_freqs_dB = 10*log10(spectrum_trial_20_averaged_trials_and_freqs);
    spectrum_trial_45_averaged_trials_and_freqs_dB = 10*log10(spectrum_trial_45_averaged_trials_and_freqs);
    relative_spectrum_trialvbaseline_20v45_dB = spectrum_trial_20_averaged_trials_and_freqs_dB - spectrum_trial_45_averaged_trials_and_freqs_dB;
        % Topoplot
    myCmap = asymColorMapWhiteZero([-0.9,0.9], N_colors_standard);
    set(0,'DefaultFigureColormap',myCmap);
    figure;
    colorbar;
    topoplot(relative_spectrum_trialvbaseline_20v45_dB, file_electrode_positions, 'colormap', myCmap, 'maplimits', [-0.9;0.9]);
    title({'Heatmap for 4-7.5 Hz','for relative spectrum 20° vs 45°'});
end







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

function plot_heatmap_baseline_or_condition(title_str, data_heatmap, electrode_labels, y, N_colors_standard)
    
    % colorLimits = [-10, 1]; % used to be colorLimits = [-1, 550]; for 20vs45 ?
    if contains(title_str, 'corrected')
        colorLimits = [-1.5, 1.5]; % For baseline comparisons
    else
        colorLimits = [-10, 1];     % For general comparisons
    end


    figure;
    myCmap = asymColorMapWhiteZero(colorLimits, N_colors_standard);
    heatmap_plot = heatmap(electrode_labels, y, data_heatmap', 'Colormap', myCmap, 'ColorLimits', colorLimits, 'ColorbarVisible', 'on', 'XLabel', 'Electrodes', 'YLabel', 'Frequencies [Hz]');
    heatmap_plot.Title = title_str;  % Title based on the type of comparison
    
    % Showing only the ticks of interest
    kept_frequencies = {'1', '2', '4', '7.5', '12', '30'};
    CustomYLabels = string(y);
    CustomYLabels(~ismember(CustomYLabels, kept_frequencies)) = " ";
    heatmap_plot.YDisplayLabels = CustomYLabels;
    
    % Remove grid lines and place custom lines
    grid off;
    origState = warning('query', 'MATLAB:structOnObject');
    cleanup = onCleanup(@() warning(origState));
    warning('off', 'MATLAB:structOnObject');
    S = struct(heatmap_plot); 
    ax = S.Axes; 
    clear('cleanup');
    
    % Place lines around selected columns and rows
    row = [7, 14, 23, 59];
    col = [50, 63, 77, 107];
    arrayfun(@(x) xline(ax, x, 'k-', 'Alpha', 0.3), col - 0.25);
    arrayfun(@(x) yline(ax, x, 'k-', 'Alpha', 0.3), row - 0.25);
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
        legend([num2str(angle1) '°'], ['Baseline ' num2str(angle1) '°']);
    elseif strcmp(comparison_type, 'trial_vs_trial')
        title(['Spectrum over ' region ' Region',' (' num2str(angle1) '° vs. ' num2str(angle2) '°)']);
        legend([num2str(angle1) '°'], [num2str(angle2) '°']);
    end

    xlabel('Frequencies [Hz]');
    ylabel('Power [dB]');
    grid on;
end


