function [EEG, automatic_cleaning]=autoClean_continuousEEG_custom(EEG, settings, subject)
%%% automatic epoch cleaning (including time-warp, if selected)
%%% requires specifications given in the "diary_automatic_cleaning_continuous_EEG.m"
%%% uses the function"EEG_find_artifact_indices.m"
%%% cleans and plots: EEG.data (sensor space)
%%% saves: results of artifact correction only (no EEG)
%%%
%%% USAGE:
%%%   >> [EEG]=wrapper_automatic_cleaning_continuous_EEG(datapath_specifications,filename_specifications,automatic_cleaning_settings)
%%%
%%% INPUTS:
%%%     EEG         The continuous EEG dataset to clean

%%%     automatic_cleaning_settings.datapath_save_files='';      %%% keep last \; path for saving updated EEG
%%%     automatic_cleaning_settings.datapath_save_figures='';    %%% keep last \;???path for saving figures of cleaning
%%%     automatic_cleaning_settings.cleaned_data_type='sensor';   %%% "sensor" or "ICA"; ICA not implemented yet; usually bad segments found on sensor level are also fine for IC later on
%%%     automatic_cleaning_settings.selected_sensor_channels_for_cleaning=[];  %%% [] use all available channels for cleaning, else specify [1 2 ...]; currently same channels for all subjects alike
%%%     automatic_cleaning_settings.chan_select_plot_before_after=[];          %%% [] use all available channels for cleaning; else specify digits
%%%     automatic_cleaning_settings.band_artifact_cleaning=[];                 %%% in Hz; [] for no filter
%%%     automatic_cleaning_settings.band_stop_artifact_cleaning=[];  %%% [] for nothing; else e.g. [48 52] for removal of line artifacts (notch filter)
%%%     automatic_cleaning_settings.band_filtorder=2;                %%% for IIR Butterworth filter; since filtfilt is used, the effective order is double (here 4)
%%%     automatic_cleaning_settings.analyzed_files_N=1;              %%% at least 1 file
%%%     automatic_cleaning_settings.crit_all=[];                     %%% e.g., 0.9=90% keep amount of epochs; 1 value if only 1 file (no appended recordings); else indicate, e.g., 4 values for 4 appended files
%%%     automatic_cleaning_settings.wind_ms=1000;                    %%% in ms, epochs for finding large artifacts
%%%     automatic_cleaning_settings.crit_keep_sec=10;                %%% in seconds; value should be multiple of "wind_ms"; additionally remove data, i.e., keep uninterrupted "good" data segments not shorter than this
%%%     automatic_cleaning_settings.crit_percent_sample_epoch=0.1;   %%% [] for nothing; e.g., 0.1 = 10%; remove epochs if they have more than x % of samples with zero or NaN ("flat line")
%%%     automatic_cleaning_settings.weighting_factor_epoch_cleaning_methods=[1 1 3]; %%% method I mean of epochs, method II = channel heterogeneity --> SD across channels of mean_epochs; method III = channel heterogeneity --> Mahal. distance of mean_epochs across channels; recommended: put method I not at zero, because Mahal. might not yield results if data set is extremely short
%%%     automatic_cleaning_settings.visual_inspection_mode=false;    %%% =false if visual threshold rejection after automatic cleaning should not be applied; in this case, bad segments from previous automatic artifact rejection are taken
%%%     if ~automatic_cleaning_settings.visual_inspection_mode
%%%         automatic_cleaning_settings.threshold_visual_reject=zeros(1,automatic_cleaning_settings.analyzed_files_N);
%%%     end
%%%
%%%
%%% OUTPUTS:
%%% EEG                     EEG structure with unchanged data
%%%                             All artifact correction information should be found in EEG.etc.automatic_cleaning
%%% automatic_cleaning     Direct access to the resulting artifact correction information
%%%                             Copy of EEG.etc.automatic_cleaning
%%%
%%% so far: implemented ONLY FOR 1 CONDITION (no markers between appended files)
%%%         visual inspected mode not implemented here yet
%%%
%%% 0. Initialization
%%% 1. load data and do some checks
%%% 2. automatic continuous data cleaning, step I: find "bad" segments; save figures
%%% 3. step 2: indicate final invalid segments, plot selected channel(s) before and after cleaning, save figures
%%%
%%%
%%% internal use only ;)
%%% FH, 2018

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 0. Initialization
% get path to save figures
datapath_save=settings.output_filepath;

% Import settings from the struct
cleaned_data_type=settings.cleaned_data_type;
selected_chans4cleaning=settings.selected_sensor_channels_for_cleaning;
chan_plot_before_after=settings.chan_select_plot_before_after;
band_artifact_cleaning=settings.band_artifact_cleaning;
band_stop_artifact_cleaning=settings.band_stop_artifact_cleaning;
band_filtorder=settings.band_filtorder;
analyzed_files_N=settings.analyzed_files_N;
crit_all=settings.crit_all;
visual_inspect=settings.visual_inspection_mode;
if visual_inspect
    th_visual_reject=settings.threshold_visual_reject;
end

%%% if not existing already, then create folders
if ~exist(datapath_save, 'dir')
    mkdir(datapath_save);
end

%% 1. load data and do some checks
% copy original data for safety check later
EEG_original=EEG;
timesMin=EEG.times./(1000*60);

% security checks
if visual_inspect
    error('not implemented yet in this script; please contact FH for updates')
end
if analyzed_files_N > 1
    warning('"analyzed_files_N": cleaning of several appended files only implemented in the "spot rotation" study yet; please contact FH for updates')
end

if length(crit_all) ~= analyzed_files_N
    error('"crit_all": specify for each of the analyzed files ("analyzed_files_N")')
end
if contains(cleaned_data_type,'ICA')
    error('cleaning for ICA not implemented yet; please contact FH for updates')
end

%% 2. automatic continuous data cleaning, step I: find "bad" segments; save figures
% find indices of bad segments; EEG.data remains unchanged here
% band-pass filtering, Hilbert transform (amplitude envelope) for artifact detection
% indicates final "valid" data segments as "1" in EEG.invalid_segments_cleaning_start_stop_sample
% EEG.event remains unchanged
% saves figures in specified folder

%% 2.1 indicate selected channels for artifact cleaning
if isempty(selected_chans4cleaning)  % use all available channels
    sel_chans=1:size(EEG.data,1);
else % user-selected channels
    sel_chans=selected_chans4cleaning;
    error('Not using all channels: not safe to run (script not up to date)')
end
settings.selected_sensor_channels_for_cleaning = sel_chans;

% Filter the data as specified
tempdata4plot=double(EEG.data);
if ~isempty(band_artifact_cleaning)
    if length(band_artifact_cleaning)==2
        [b,a]=butter(band_filtorder, band_artifact_cleaning/(EEG.srate/2));
        for chan=1:size(tempdata4plot,1)
            tempdata4plot(chan,:)=filtfilt(b,a,tempdata4plot(chan,:));
        end
    else
        error('band not well defined for band-pass filtering; e.g., [1 40] Hz')
    end
end

if ~isempty(band_stop_artifact_cleaning)
    if length(band_stop_artifact_cleaning)==2
        [b,a]=butter(band_filtorder,band_stop_artifact_cleaning/(EEG.srate/2),'stop');
        for chan=1:size(tempdata4plot,1)
            tempdata4plot(chan,:)=filtfilt(b,a,tempdata4plot(chan,:));
        end
    else
        error('band not well defined for band-stop filtering; e.g., [48 52] Hz')
    end
end

% plot selected channels of the raw data (filtered in the selected band for cleaning)
ax=[timesMin(1) timesMin(end) -500 500];
figure('units','normalized','outerposition',[0 0 1 1])
colors=cool(length(sel_chans));
hold on
for chan=1:length(sel_chans)
    plot(timesMin, tempdata4plot(sel_chans(chan),:),'color',colors(chan,:));
end
axis(ax)
xlabel({'Time [min]',['Length of raw data set: ' num2str(EEG.xmax/60) ' min']})
ylabel('Potential [microV]')
if ~isempty(band_artifact_cleaning)
    title([subject, ' before cleaning - Bandpass ['...
        num2str(band_artifact_cleaning(1)) ';' num2str(band_artifact_cleaning(2)) '] Hz - all channels superimposed'])
else
    title([subject, ' before cleaning - Raw data - all channels superimposed'])
end

% save figure
saveCurrentFig(fullfile(datapath_save, filesep), sprintf('%s_continuous_data_autoclean_before',subject), {'png'}, []); % too large for .fig

%% 2.2 Find artifacts (not removed here yet)
EEG = EEG_find_artifact_indices(EEG, settings, subject);

%%% store info in "EEG.etc.automatic_cleaning"
EEG.etc.automatic_cleaning.data_length_original_Nsamples=size(EEG.data,2);
EEG.etc.automatic_cleaning.data_length_original_minutes=size(EEG.data,2)/EEG.srate/60;

if isfield(EEG.etc.automatic_cleaning,'invalid_segments_start_stop_sample')
    % all appended files: how much data would be removed by previous automatic break rejection only?
    invalid_segments_index=EEG.etc.automatic_cleaning.invalid_segments_start_stop_sample;
    EEG_rej=eeg_eegrej(EEG,invalid_segments_index);
    EEG.etc.automatic_cleaning.data_length_after_autobreakrejection_Nsamples=size(EEG_rej.data,2);
    EEG.etc.automatic_cleaning.data_length_after_autobreakrejection_minutes=size(EEG_rej.data,2)/EEG_rej.srate/60;
end

%% 3. Indicate final invalid segments, plot selected channel(s) before and after cleaning, save figures
% 3.1 indicate final invalid segments
invalid_segments_index=EEG.etc.automatic_cleaning.invalid_segments_final_start_stop_sample;
EEG_rej=eeg_eegrej(EEG,invalid_segments_index);

%%% auto-cleaning without previous automatic break detection
if ~isfield(EEG.etc.automatic_cleaning,'invalid_segments_start_stop_sample')...
        && isfield(EEG.etc.automatic_cleaning,'invalid_segments_final_start_stop_sample')
    EEG.etc.automatic_cleaning.data_length_after_autocleaning_Nsamples=size(EEG_rej.data,2);
    EEG.etc.automatic_cleaning.data_length_after_autocleaning_min=size(EEG_rej.data,2)/EEG_rej.srate/60;
end

%%% auto-cleaning with previous automatic break detection
if isfield(EEG.etc.automatic_cleaning,'invalid_segments_start_stop_sample')...
        && isfield(EEG.etc.automatic_cleaning,'invalid_segments_final_start_stop_sample')
    EEG.etc.automatic_cleaning.data_length_after_autobreakautocleaning_Nsamples=size(EEG_rej.data,2);
    EEG.etc.automatic_cleaning.data_length_after_autobreakautocleaning_min=size(EEG_rej.data,2)/EEG_rej.srate/60;
end

if isfield(EEG.etc.automatic_cleaning,'boundaries_between_files_raw_firstGoodSampleOfNewFile')
    for n=1:length(EEG.etc.automatic_cleaning.boundaries_between_files_raw_firstGoodSampleOfNewFile)
        removed_data_length_temp1=EEG.etc.automatic_cleaning.boundaries_between_files_raw_firstGoodSampleOfNewFile(n);
        eval(['removed_data_length_temp2=EEG.etc.automatic_cleaning.auto_epoch_rejection_appendedfile' num2str(n) '.data_length_final_after_autobreakautoclean_keepcleanseg_Nsamp;'])
        if n==1
            if ~isfield(EEG.etc.automatic_cleaning,'invalid_segments_start_stop_sample')...
                    && isfield(EEG.etc.automatic_cleaning,'invalid_segments_final_start_stop_sample')
                EEG.etc.automatic_cleaning.after_autocleaning_firstGoodSampleOfNewFile(n)=removed_data_length_temp2+1;
            end
            if isfield(EEG.etc.automatic_cleaning,'invalid_segments_start_stop_sample')...
                    && isfield(EEG.etc.automatic_cleaning,'invalid_segments_final_start_stop_sample')
                EEG.etc.automatic_cleaning.after_autobreakautocleaning_firstGoodSampleOfNewFile(n)=removed_data_length_temp2+1;
            end
        else
            if ~isfield(EEG.etc.automatic_cleaning,'invalid_segments_start_stop_sample')...
                    && isfield(EEG.etc.automatic_cleaning,'invalid_segments_final_start_stop_sample')
                EEG.etc.automatic_cleaning.after_autocleaning_firstGoodSampleOfNewFile(n)=EEG.etc.automatic_cleaning.after_autocleaning_firstGoodSampleOfNewFile(n-1)+removed_data_length_temp2;
            end
            if isfield(EEG.etc.automatic_cleaning,'invalid_segments_start_stop_sample')...
                    && isfield(EEG.etc.automatic_cleaning,'invalid_segments_final_start_stop_sample')
                EEG.etc.automatic_cleaning.after_autobreakautocleaning_firstGoodSampleOfNewFile(n)=EEG.etc.automatic_cleaning.after_autobreakautocleaning_firstGoodSampleOfNewFile(n-1)+removed_data_length_temp2;
            end
        end
    end
end

% 3.2 Plot selected channel(s) before and after cleaning, save figures
if isempty(chan_plot_before_after)
    tempdata4plot_befcleaning=tempdata4plot;
else
    tempdata4plot_befcleaning=tempdata4plot(chan_plot_before_after,:);
end
%clear tempdata4plot

colors=cool(2*length(chan_plot_before_after));

%%% before cleaning
figure ('units','normalized','outerposition',[0 0 1 1])
hold on
for chan=1:size(tempdata4plot_befcleaning,1)
    plot(timesMin,tempdata4plot_befcleaning(chan,:),'color',colors(chan,:));
end

%%% after cleaning
invalid_segments_index=EEG.etc.automatic_cleaning.invalid_segments_final_start_stop_sample;
tempdata4plot_aftercleaning=tempdata4plot_befcleaning;
% place NANs where the segments were removed
for seg=1:size(invalid_segments_index,1)
    tempdata4plot_aftercleaning(:,invalid_segments_index(seg,1):invalid_segments_index(seg,2))=NaN;
end

for chan=1:size(tempdata4plot_aftercleaning,1)
    plot(timesMin,tempdata4plot_aftercleaning(chan,:),'color',colors(length(chan_plot_before_after)+chan,:));
end
legend({'before cleaning', 'after cleaning'})
axis(ax)
xlabel({'Time [min]',['Length of raw data set: ' num2str(EEG.xmax/60) ' min'],...
    ['Length after cleaning: ' num2str(length(find(isnan(tempdata4plot_aftercleaning(1,:))==0))/(EEG.srate*60)) ' min']});
ylabel('Potential [??V]')
if ~isempty(band_artifact_cleaning)
    title({[subject ' comparison before/after cleaning - Bandpass ['...
        num2str(band_artifact_cleaning(1)) ';' num2str(band_artifact_cleaning(2)) '] Hz'],...
        ['Channel(s): ' num2str(chan_plot_before_after)]})
else
    title({[subject ' comparison before/after cleaning - Raw data'],...
        ['Channel(s): ' num2str(chan_plot_before_after)]})
end

% save figure
saveCurrentFig(fullfile(datapath_save, filesep), sprintf('%s_continuous_data_autoclean_final_before_after_selChans', subject), {'png'}, []); % too large for .fig

% finally: externalize results of artifact correction for saving; saving it with EEG takes too much disk space and takes too long!
EEG.etc.automatic_cleaning.settings=settings;
EEG.etc.automatic_cleaning.cleaning_continuousOrEpochs='continuous';

automatic_cleaning=EEG.etc.automatic_cleaning;

% display results
disp('To be removed segments are in: EEG.etc.automatic_cleaning.invalid_segments_final_start_stop_sample')

disp(['length ORIGINAL data set = ' num2str(size(EEG.data,2)/EEG.srate/60) ' min'])
disp(['length SHORTENED data set after cleaning = ' num2str(length(find(isnan(tempdata4plot_aftercleaning(1,:))==0))/EEG.srate/60) ' min'])
temp_removed_data=size(EEG.data,2)/EEG.srate/60-length(find(isnan(tempdata4plot_aftercleaning(1,:))==0))/EEG.srate/60;
disp(['length REMOVED data = ' num2str(temp_removed_data) ' min'])

% safety check: EEG.data shouldn't be changed
if sum(sum(EEG_original.data)) ~= sum(sum(EEG.data))
    error('EEG.data was erroneously changed')
end
end