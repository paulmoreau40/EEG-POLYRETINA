%% ONLY CHANGE THESE PARTS!
user = 'paul';

%% General foldernames and filenames
[study_config.study_folder, study_config.raw_data_folder, study_config.electrodes_folder,...
    study_config.raw_EEGLAB_data_folder, study_config.preprocessing_folder,...
    study_config.single_subject_analysis_folder, study_config.multi_subject_analysis_folder,...
    study_config.figures_folder, sourceFile, tab, study_config.recon.path] = getMainFoldersNames(user);

study_config.user = user;

%% Subjects definition:
study_config.subjects = table2struct(readtable(sourceFile, 'Sheet', tab));
for s = 1:length(study_config.subjects)
    if isempty(study_config.subjects(s).badElectrodes)
        study_config.subjects(s).badElectrodes = {};
    else
        study_config.subjects(s).badElectrodes = split(study_config.subjects(s).badElectrodes,',');
    end
    
    if isempty(study_config.subjects(s).missingData)
        study_config.subjects(s).missingData = {};
    else
        study_config.subjects(s).missingData = split(study_config.subjects(s).missingData,',');
    end
    
%     if isempty(study_config.subjects(s).badLocElectrodes)
%         study_config.subjects(s).badLocElectrodes = {};
%     else
%         study_config.subjects(s).badLocElectrodes = split(study_config.subjects(s).badLocElectrodes,',');
%     end
end

included = find(strcmp({study_config.subjects.excluded},'No'));
with_elecs = [];
%%% Subjects with all ICs labeled:
iclabeled = [];
%%% Subjects with forward model ready:
forward_complete = [];
%%% Subjects with ROI reversed to subject space:
roi_reversed = [];
%%% Subjects with mapping between ROI and dipoles ready:
roi_mapped = intersect(roi_reversed,[]);

%subject_inds = setdiff(included,[10,35]); % COMMENTED BY PAUL 
subject_inds = included;
%subject_inds = intersect(forward_complete,roi_reversed);
subject_ind = subject_inds(1);
study_config.current_subject = subject_ind;

%study_config.figures_folder = 'C:\Users\Alexandre\Desktop\results\figures\';
study_config.filenames = {'block001_1','block001_2','block001_3','block001_4','block002_1','block002_2','block003_1','block003_2'};

% Enter the number of streams expected for each category:
% Categories: {EEG, ET, MOCAP, Events}
study_config.stream_count = [2,0,0,1];
% Enter the TYPES to put in each category
study_config.eeg_streams = {'EEG'};
study_config.eye_tracker_streams = {};
study_config.rigidbody_streams = {'rigidBody'};
%study_config.event_streams = {'Markers'};
study_config.event_streams = {'Tag'}; % PAUL
study_config.bad_streams = {'Markers'}; %PAUL

% Specific to ANT: 128ch creates 2 EEG streams (one per amplifier) and it
% is important to know which one should be imported first to assign the
% right channel labels (indicate stream NAMES here)
study_config.eeg_streams_order = {'EE225-000000-000430', 'EE225-000000-000433'};
study_config.recording_unit = 'Volt';

% Enter channels that you did not use at all:
%study_config.channels_to_remove = {'N29' 'N30' 'N31'};
%study_config.channels_to_remove = [];
% For ANT data, the last 2 channels are 'trigger' and 'counter' and shoud be removed
study_config.channels_to_remove = [65,66;65,66];

% Enter EOG channel names here:
study_config.eog_channels  = {'VEOG'};
%study_config.eog_channels  = {'G16' 'G32'};
%study_config.eog_channels  = {'EOG2'};
study_config.rename_channels = [];

% Leave this empty if you have standard channel names that should use standard locations:
%study_config.channel_locations_filename = 'CA-213_waveguard128_duke.elc'; % = [];
%study_config.channel_locations_filename = 'CA-213_waveguard128_duke_modif.elc'; % = []; % PAUL: elc file modified (because of error with loadtxt() in readeetraklocs
study_config.channel_locations_filename = 'CA-213_EOG.elc'; % PAUL: "EOG" instead of noEOG, to have 128 channels like in the EEG data (otherwise, 127 if removes VEOG)
study_config.indiv_channel_locations_filename = '-get_chanlocs.txt';

%% Definition of the global preprocessing architecture:
% 'bemobil' = architecture developped in Gramann lab
% 'simple' = simple pipeline for quick results
study_config.globalArchitecture = 'bemobil';

%% Preprocessing
% Electrodes that should be excluded right away (known recording issue)
% One line per subject concerned
study_config.moveElecInwards = 5; % in mm
%study_config.mocap_lowpass = 6;
%study_config.rigidbody_derivatives = 2;
% Time buffer to include before and after each trial selected for analysis
study_config.trialBuffer = 3; % in seconds
% Maximal length of NaN series to interpolate (fillNaNs function) 
study_config.maxNans2Replace = 1; % in samples @1kHz
% Resampling frequency
study_config.resample_freq = 250; % in Hz

% Filters depending on the step:
study_config.filterPreProc.low_cut_off = 1.5; % in Hz, put [] for no HP filter
study_config.filterPreProc.high_cut_off = []; % in Hz, put [] for no LP filter
study_config.filterICLabel.low_cut_off = 0.5; % in Hz, put [] for no HP filter
study_config.filterICLabel.high_cut_off = []; % in Hz, put [] for no LP filter
study_config.filterAnalysis.low_cut_off = 1; % in Hz, put [] for no HP filter
study_config.filterAnalysis.high_cut_off = 42; % in Hz, put [] for no LP filter

%% Parameters for bemobil pipeline:
%%%% Line noise removal technique:
% 'cleanLine' = use the Cleanline plugin (as in APP)
% 'cleanLinePREP' = use the PREP plugin, CleanLineNoise function
study_config.lineNoiseRemoval_method = 'cleanLinePREP';

%% Definition of the bad samps rejection method
study_config.badSampsRejection = 'autoMoBI'; % 'manual', 'app', 'asr', 'autoMoBI'

%% Parameters for 'autoMoBI'
% 'automatic' = method used in Berlin with automatic detection of bad segments, preceeded by a first ICA to remove eye movements
study_config.autoMoBI.DoI4AMICAeye= {'baselines'}; % POL: use every block baselines as baselines for Eye ICA

study_config.autoMoBI.cleaned_data_type='sensor data'; % ICA not implemented yet; usually bad segments found on sensor level are also fine for IC later on

% channels that should be considered for cleaning
% [] use all available channels for cleaning, else specify [1 2 ...]; currently same channels for all subjects alike
study_config.autoMoBI.selected_sensor_channels_for_cleaning=[];

% channel(s) for plotting before vs. after
study_config.autoMoBI.chan_select_plot_before_after=[89];  % [] use all available channels for cleaning

% define frequency band (band-pass filter) only for cleaning
% (No filter because already taken care of at this stage)
study_config.autoMoBI.band_artifact_cleaning=[];  % in Hz; [] for no filter
study_config.autoMoBI.band_stop_artifact_cleaning=[]; % [] for nothing; else e.g. [48 52] for removal of line artifacts (notch filter)
study_config.autoMoBI.band_filtorder=2;               % for IIR Butterworth filter; since filtfilt is used, the effective order is double (here 4)

%%% further settings for the cleaning algorithm
study_config.autoMoBI.analyzed_files_N=1;  % so far: analyze only 1 file! (no appended different conditions)
study_config.autoMoBI.crit_all=0.85; % e.g., 0.9=90% keep amount of epochs; 1 value if only 1 file (no appended recordings); else indicate, e.g., 4 values for 4 appended files
study_config.autoMoBI.wind_ms=1000;    % in ms, epochs for finding large artifacts
study_config.autoMoBI.crit_keep_sec=[]; % in seconds; for NOT removing any additional data, set: automatic_cleaning_settings.crit_keep_sec=automatic_cleaning_settings.wind_ms/1000; else: value should be multiple of "wind_ms" for additionally removing data, i.e., keep uninterrupted "good" data segments not shorter than this value
study_config.autoMoBI.crit_percent_sample_epoch=0.2;  % [] for nothing; e.g., 0.1 = 10%; remove epochs if they have more than x % of samples with zero or NaN ("flat line")
study_config.autoMoBI.weighting_factor_epoch_cleaning_methods=[1 1 2]; % method I mean of epochs, method II = channel heterogeneity --> SD across channels of mean_epochs; method III = channel heterogeneity --> Mahal. distance of mean_epochs across channels; recommended: put method I not at zero, because Mahal. might not yield results if data set is extremely short
study_config.autoMoBI.visual_inspection_mode=false;  % =false if visual threshold rejection after automatic cleaning should not be applied; in this case, bad segments from previous automatic artifact rejection are taken
if ~study_config.autoMoBI.visual_inspection_mode
    study_config.autoMoBI.threshold_visual_reject=zeros(1,study_config.autoMoBI.analyzed_files_N);
end

study_config.autoMoBI.buffer_length = 0.2; % in percent of the epoch length used for automatic cleaning
study_config.autoMoBI.output_filepath = fullfile(study_config.figures_folder,...
    'AutoBadSampsCleaning');% path for saving figures

%% Parameters for ASR:
study_config.ASR_use = 'reject'; % 'rewrite' or 'reject'
study_config.burst_crit = 25; % Burst criterion for ASR. Put [] for default

%% Parameters for the APP pipeline:
study_config.APP.censorBiweight = 6; % censor value, it corresponds to a certain number of std for a normal distribution
% c=6 is 4 std; c=7.5 is 5 std; c=9 is 6 std.
study_config.APP.z_criterion = 3.5; % rejection criterion according to literature
study_config.APP.inner_fence = 1.5; % inner_fence criterion for extreme outliers

% from M. Hubert, E. Vandervieren, An adjusted boxplot for skewed distributions,
% Computational Statistics & Data Analysis, Volume 52, Issue 12, 2008, Pages 5186-5201, ISSN 0167-9473,
% https://doi.org/10.1016/j.csda.2007.11.008.
%study_config.APP.skew_side_limit = 3;
%study_config.APP.opp_side_limit = 4;

% from Janir:
%study_config.APP.skew_side_limit = 4;
%study_config.APP.opp_side_limit = 3.5;

study_config.APP.skew_side_limit = 2;
study_config.APP.opp_side_limit = 4;

%% Different ICA methods:
% 'runica' 
% 'amica'
study_config.ICAmethod = 'amica';
study_config.max_threads = 4;
study_config.num_models = 1;
study_config.max_iter = 2000;

%% Dipole fitting:
study_config.doDipoleFitting = true;
study_config.residualVariance_threshold = 100;
study_config.do_remove_outside_head = 'off';
study_config.number_of_dipoles = 1;

%% IC_label
% classes: {'Brain'  'Muscle'  'Eye'  'Heart'  'Line Noise'  'Channel Noise'  'Other'}
% The thresholds are extracted for the latest released paper about IClabel
% They are extracted from Table 3 to optimize accuracy of the classifier (Test)

% Different ICLabel options:
% 'lite' ; [0.53,0.17,0.06,0.10,0.42,0.15,0.29];
% 'default'; [0.35,0.30,0.04,0.03,0.84,0.05,0.26];
study_config.iclabel = 'default';
study_config.ICdetect_thresholds = [0.35,0.30,0.04,0.03,0.84,0.05,0.26];
study_config.cats2keep = {'Brain','BrainWithNoise'};

%% Epoching
study_config.epochs.event = 'Observation';
study_config.epochs.window = 'fixed'; % 'full' = based on durations, 'fixed' = using limits_wdw
study_config.epochs.limits_wdw = [0,3]; % in seconds
study_config.epochs.bumper = 1; % in seconds (additional time to include before and after the trial
                                % to prevent artifacts in ERSP computation = padding)
                                
%% Computing ERSPs
study_config.ersps.method = 'superlet';
study_config.ersps.FoI = 2:40; % frequencies of interest, in Hz, make sure there is a sufficient margin with respect to filters
study_config.ersps.timeRes = 10/study_config.resample_freq; % Temporal resolution of the ERSPs (in seconds)  
% Parameters for superlet:
study_config.ersps.superlet.basewidth = 3; % 'width', or number of cycles, of the base wavelet (default = 3)
study_config.ersps.superlet.gwidth = 3; % determines the length of the used wavelets in standard deviations of the implicit Gaussian kernel and should be choosen >= 3; (default = 3)
study_config.ersps.superlet.combine = 'multiplicative'; % 'additive', 'multiplicative' (default = 'additive') determines if cycle numbers of wavelets comprising a superlet are chosen additively or multiplicatively
% Cycles factors mimicking at best implementation in newtimef:
cycles_factor = linspace(1, 0.5*study_config.ersps.FoI(end)/study_config.ersps.FoI(1), length(study_config.ersps.FoI));
study_config.ersps.superlet.order = ceil(cycles_factor); % vector 1 x numfoi, superlet order, i.e. number of combined wavelets, for individual frequencies of interest.

%% Source reconstruction and ROIs
study_config.recon.path = 'D:\SourceReconstruction\EEG-AFF\FS_FinalSeg';
% Importation of fwd model
study_config.recon.fixed_ori = true;
study_config.recon.patch_space = true;
study_config.recon.src_spacing = 'ico5';
study_config.recon.downsamp_bem = 'None';
% Working with rois
study_config.recon.hemispheres = {'lh', 'rh'};
study_config.recon.rois = {'PPA', 'RSC', 'OPA'};
study_config.recon.roi_creation = 'withGM_sph10mm_60vox';
study_config.recon.resolution = 1; % resolution at which the ROI mask was exported (in mm)
study_config.recon.res_vect = [-study_config.recon.resolution, study_config.recon.resolution, study_config.recon.resolution]/2;
study_config.recon.dist_limit_mm = [2]; % in mm
study_config.recon.infl_levels = [0.5,0.25]; % One value per neighbouring order considered
study_config.recon.infl_th = 0.1; % In case influence matrix is computed with distance
study_config.recon.slight_pos_val = 1e-5;
% Direct model adaptations
study_config.recon.fwd.cholesky = false;
study_config.recon.fwd.normalizeDepth = true;
study_config.recon.fwd.normalizeDepthParam = 0.9;
study_config.recon.fwd.weightWithPatchSize = false;
% Inverse model
study_config.recon.inv.method = 'mne';
% Parameters for mne inverse modelling
study_config.recon.inv.noiselambda = 0; % Regularization parameter for noise covariance - only used for pre-whitening.
%study_config.recon.inv.lambda = 1e-3; % regularization parameter (typical value but will be optimized per subject)
% lambda is estimated from snr if not specified.
%study_config.recon.inv.snr = 3;
study_config.recon.inv.prewhiten = 'no';
study_config.recon.inv.scalesourcecov = 'no';
%}

%% Filenames suffix
study_config.merged_filename = 'allBlocks.set'; % In all pipelines
study_config.prepared_filename = 'prepared.set'; % In all pipelines
study_config.ASRin_filename = 'inputASR.set'; % In simple pipeline
study_config.ASRout_filename = 'outputASR.set'; % In simple pipeline
study_config.BadChansRemoved_filename = 'inter_avRef.set'; % In bemobil pipeline
study_config.EyeDetection_filename = 'eyeCompsRemoval.set'; % In autoMoBI pipeline
study_config.beforeICA_filename = 'cleanedForICA.set'; % In all pipelines
study_config.icaOutput_filename = 'postICA.set'; % In all pipelines
study_config.dipolesFitted_filename = 'postICA_withdips.set'; % In all pipelines
study_config.icLabelled_filename = 'automaticIClabels.set'; % In all pipelines
study_config.icaSelect_filename = 'cleaned_with_ICA.set'; % In all pipelines
study_config.postICA_tempRej_filename = 'final.set';
study_config.epoched_filename = 'epoched.set'; % In all pipelines
study_config.base_epoched_filename = 'epoched_baselines.set'; % In all pipelines

study_config.interpolated_avRef_filename = 'interpolated_avRef.set'; %
study_config.filtered_filename = 'filtered.set';
study_config.without_bad_temp_segments = 'badTempSegmentsRemoved.set';
study_config.without_eyes = 'no_eyes.set';
study_config.ica_filename_output = 'postICA.set';
study_config.warped_dipfitted_filename = 'warped_dipfitted.set';
study_config.copy_weights_interpolate_avRef_filename = 'interp_avRef_ICA.set';
study_config.categorizedICs = 'ICLabel_computed.set';
study_config.single_subject_cleaned_ICA_filename = 'cleaned_with_ICA.set';

%{
%% AMICA
% on some PCs AMICA may crash before the first iteration if the number of
% threads and the amount the data does not suit the algorithm. Jason Palmer
% has been informed, but no fix so far. just roll with it. if you see the
% first iteration working there won't be any further crashes. in this case
% just press "close program" or the like and the bemobil_spatial_filter
% algorithm will AUTOMATICALLY reduce the number of threads and start AMICA
% again. this way you will always have the maximum number
% of threads that should be used for AMICA. check in the
% task manager how many threads you have theoretically available and think
% how much computing power you want to devote for AMICA. on the bpn-s1
% server, 12 is half of the capacity and can be used. be sure to check with
% either Ole or your supervisor and also check the CPU usage in the task
% manager before!

% 4 threads are most effective for single subject speed, more threads don't
% really shorten the calculation time much. best efficiency is using just 1
% thread and have as many matlab instances open as possible (limited by the
% CPU usage). Remember your RAM limit in this case.

study_config.filter_lowCutoffFreqAMICA = 1;
study_config.filter_highCutoffFreqAMICA = [];
study_config.max_threads = 4;
study_config.num_models = 1;
study_config.max_iter = 2000;

% warp electrodemontage and run dipfit
%bemobil_config.warping_channel_names = [];
study_config.residualVariance_threshold = 100;
study_config.do_remove_outside_head = 'off';
study_config.number_of_dipoles = 1;

% FHs cleaning
study_config.buffer_length = 0.2; % percent of the epoch length used for automatic cleaning
study_config.automatic_cleaning_threshold_to_keep = 0.85;

%% IC_label
% classes: {'Brain'  'Muscle'  'Eye'  'Heart'  'Line Noise'  'Channel Noise'  'Other'}
% The thresholds are extracted for the latest released paper about IClabel
% They are extracted from Table 3 to optimize accuracy of the classifier (Test)
study_config.ICdetection_thresholds=[0.35 0.30 0.04 0.03 0.84 0.05 0.26];
% old parameters by Marius:
%bemobil_config.eye_threshold = 0.7;
%bemobil_config.brain_threshold = 0.4;

% Classes of interest to keep after main ICA:
study_config.classes2keep = [1];

%% Parameters for the final filtering step:
study_config.final_lowCutoffFreq = 1;
study_config.final_highCutoffFreq = 40;

% Function to get the right string for naming:
if round(study_config.final_lowCutoffFreq)~=study_config.final_lowCutoffFreq
    intPart = floor(study_config.final_lowCutoffFreq);
    firstDec = round((study_config.final_lowCutoffFreq - intPart)*10);
    lowFreqString = [num2str(intPart) '-' num2str(firstDec)];
else
    lowFreqString = num2str(study_config.final_lowCutoffFreq);
end

if round(study_config.final_highCutoffFreq)~=study_config.final_highCutoffFreq
    intPart = floor(study_config.final_highCutoffFreq);
    firstDec = round((study_config.final_highCutoffFreq - intPart)*10);
    highFreqString = [num2str(intPart) '-' num2str(firstDec)];
else
    highFreqString = num2str(study_config.final_highCutoffFreq);
end

study_config.final_filtering_filename = [lowFreqString '_' highFreqString '_Hz_' study_config.filtered_filename];
%}