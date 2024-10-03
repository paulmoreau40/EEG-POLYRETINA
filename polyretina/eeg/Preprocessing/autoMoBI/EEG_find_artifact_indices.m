function [EEG]=EEG_find_artifact_indices(EEG, settings, subject)
%%% cleaning of artifacts in continuous data before spatial filtering
%%% band-pass filtering, Hilbert transform (amplitude envelope) for artifact detection
%%% calls 'EEG_find_artifact_indices_epoch_rejection'
%%% indicates final "valid" data segments as "1" in EEG.invalid_segments_cleaning_start_stop_sample
%%% EEG.event remains unchanged
%%% saves figures in specified folder
%%%
%%% 1. Filtering
%%% 2. If previous automatic break rejection is present, indicate segments as NaN
%%% 3. Create epochs of the continuous data
%%% 4) check "bad" epochs (automatically)
%%% 5) keep only "good" data segments which are uninterrupted and longer than criterion in seconds
%%% 6) keep indices of bad data segments
%%% 7) plot all channels of the finally cleaned data, save figure
%%%
%%% internal use only ;)
%%% FH, 2017

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 0. Initialization of parameters
% get path to save figures
datapath_save_figures=settings.output_filepath;

% Import settings from the struct
selected_chans4cleaning=settings.selected_sensor_channels_for_cleaning;
band=settings.band_artifact_cleaning;
band_stop=settings.band_stop_artifact_cleaning;
band_filtorder=settings.band_filtorder;
crit_all=settings.crit_all;
wind_ms=settings.wind_ms;
crit_keep_sec=settings.crit_keep_sec;
crit_percent_sample_epoch=settings.crit_percent_sample_epoch;
weights4methods=settings.weighting_factor_epoch_cleaning_methods;

%% 1. Bandpass filter, bandstop filter (if opted)
data=double(EEG.data(selected_chans4cleaning,:));
% 1.1 band-pass filter
if ~isempty(band)
    if length(band)==2
        [b,a]=butter(band_filtorder,band/(EEG.srate/2));
        disp('band-pass filtering, channel-wise...')
        for chan=1:size(data,1)
            data(chan,:)=filtfilt(b,a,data(chan,:));
        end
    else
        error('band not well defined for band-pass filtering; e.g., [1 40] Hz')
    end
end

% 1.2 band-stop (notch) filter
if ~isempty(band_stop)
    if length(band_stop)==2
        [b,a]=butter(band_filtorder,band_stop/(EEG.srate/2),'stop');
        disp('band-stop filtering, channel-wise...')
        for chan=1:size(data,1)
            data(chan,:)=filtfilt(b,a,data(chan,:));
        end
    else
        error('band not defined for band-stop filtering; e.g., [48 52] Hz')
    end
end

data_original4plot=data; % copy data for plotting uncleaned vs. cleaned in the end

% 1.3 take amplitude envelope by abs(Hilbert transform)
for chan=1:size(data,1)
    data(chan,:)=abs(hilbert(data(chan,:)));
end

%% 2. If previous automatic break rejection is present, indicate segments as NaN
%{
if isfield(EEG.etc.automatic_cleaning,'invalid_segments_start_stop_sample')
    for n=1:size(EEG.etc.automatic_cleaning.invalid_segments_start_stop_sample,1)
        data(:,(EEG.etc.automatic_cleaning.invalid_segments_start_stop_sample(n,1):EEG.etc.automatic_cleaning.invalid_segments_start_stop_sample(n,2)))=NaN;
    end
end
%}

%% 3. Create epochs of the continuous data
sampInWdw = EEG.srate*(wind_ms/1000);
if sampInWdw<1
    error('window too short')
end

windows=[0:sampInWdw:size(data,2)]';
windows(:,2) = [windows(2:end,1);size(data,2)];
windows(:,1) = windows(:,1)+1;
windows = windows(1:end-1,:);
if windows(end,2)<size(data,2)  %%% in case the very last data snippet was shorter than window_samples
    windows(end+1,1)=windows(end,2)+1;
    windows(end,2)=size(data,2);
end

%%% cut epochs
for n=1:size(windows,1)
    data_segmented{:,n}=data(:,windows(n,1):windows(n,2)); %%% cell, because epochs might have different length (last data snippet)
    data_segmented_raw{:,n}=EEG.data(selected_chans4cleaning,windows(n,1):windows(n,2));
end

%% 4. check "bad" epochs (automatically)
%%% either in whole continuous data (all appended files) = might be problematic if you have strong
%%% differences between conditions (e.g., baseline rest vs. strong movements)
%%% ... or look in each appended file (=condition) separately and finally merge the found epochs

% 4.1 Perform analysis in each single recording:
if length(crit_all) > 1
    error('The code for multiple files is deprecated: come back after updating it')
    disp('epoch rejection performed for each of the appended files')
    
    %     if length(crit_all) < length(EEG.etc.automatic_cleaning.boundaries_after_auto_breakrejection_firstGoodSampleOfNewFile)+1
    if length(crit_all) < length(EEG.etc.automatic_cleaning.boundaries_between_files_raw_firstGoodSampleOfNewFile)+1
        error('dimension mismatch; criterions should match the number of the appended files')
    end
    
    appended_file_member=zeros(1,size(windows,1));
    for n=1:length(EEG.etc.automatic_cleaning.boundaries_between_files_raw_firstGoodSampleOfNewFile)
        for nn=1:size(windows,1)
            if ismember(EEG.etc.automatic_cleaning.boundaries_between_files_raw_firstGoodSampleOfNewFile(n), windows(nn,1):windows(nn,2))
                appended_file_member(nn)=1; %%% check which epoch contains the boundary marker between appended files
            end
        end
    end
    appended_file_member=find(appended_file_member==1);
    
    if length(appended_file_member) ~= length(EEG.etc.automatic_cleaning.boundaries_between_files_raw_firstGoodSampleOfNewFile)
        error('error in detected boundary markers between files')
    end
    
    
    for n=1:length(appended_file_member)+1
        
        clear window_samples_vector
        
        if n==1
            epoch_index_appended_file=1:appended_file_member(n);
            window_samples_vector=windows(epoch_index_appended_file,:);  %%% start until the epoch containing the boundary
        end
        if n>1 && n<=length(appended_file_member)
            epoch_index_appended_file=[appended_file_member(n-1)+1:appended_file_member(n)];
            window_samples_vector=windows(epoch_index_appended_file,:);
        end
        if n==length(appended_file_member)+1
            epoch_index_appended_file=appended_file_member(n-1)+1:size(windows,1);
            window_samples_vector=windows(epoch_index_appended_file,:);  %%% epoch AFTER the epoch containing the boundary until the end
        end
        
        %%% reject epochs, see script for details, plots figure
        crit=crit_all(n);
        data_segmented_current=data_segmented(epoch_index_appended_file);
        data_segmented_raw_current=data_segmented_raw(epoch_index_appended_file);
        sampInWdw=(wind_ms/1000)*EEG.srate;
        
        [auto_epoch_inspection]=EEG_find_artifact_indices_epoch_rejection(data_segmented_current,data_segmented_raw_current,window_samples_vector,crit,crit_percent_sample_epoch,sampInWdw,EEG.srate,weighting_factor_epoch_cleaning_methods,selected_channels);
        
        auto_epoch_inspection.epochs_finalNote='final epoch indices';
        auto_epoch_inspection.good_epochs_final=auto_epoch_inspection.indices_good_epochs;
        auto_epoch_inspection.bad_epochs_final=auto_epoch_inspection.indices_bad_epochs;
        
        savefig([datapath_save_figures 'autocleaning_appended_file_' num2str(n) ''])
        close
        
        
        %%% mark data that should be removed as NaN for plotting
        window_samples_vector=auto_epoch_inspection.window_samples_vector_final_final;
        data_new=data;
        temp_index=find(window_samples_vector(:,3)==1);
        for epoch=1:length(temp_index)
            data_new(:,window_samples_vector(temp_index(epoch),1):window_samples_vector(temp_index(epoch),2))=NaN;
        end
        
        %%% how much data would be removed with respect to the raw data?
        if n==1
            boundary_temp=EEG.etc.automatic_cleaning.boundaries_between_files_raw_firstGoodSampleOfNewFile(n)-1;
            auto_epoch_inspection.data_length_raw_Nsamples=length(1:boundary_temp);
            auto_epoch_inspection.data_length_raw_min=length(1:boundary_temp)/EEG.srate/60;
            
            data_new_boundary=data_new(:,1:boundary_temp);  %%% data_new is created by the previous script "EEG_find_artifact_indices_epoch_rejection"
            data_new_good_length=length(find(isnan(data_new_boundary(1,:))==0));  %%% how many good samples?
            auto_epoch_inspection.data_length_after_autobreakautocleaning_Nsamples=data_new_good_length;
            auto_epoch_inspection.data_length_after_autobreakautocleaning_min=data_new_good_length/EEG.srate/60;
        end
        
        if n>1 && n<=length(appended_file_member)
            boundary_temp1=EEG.etc.automatic_cleaning.boundaries_between_files_raw_firstGoodSampleOfNewFile(n-1);
            boundary_temp2=EEG.etc.automatic_cleaning.boundaries_between_files_raw_firstGoodSampleOfNewFile(n)-1;
            auto_epoch_inspection.data_length_raw_Nsamples=length(boundary_temp1:boundary_temp2);
            auto_epoch_inspection.data_length_raw_min=length(boundary_temp1:boundary_temp2)/EEG.srate/60;
            
            data_new_boundary=data_new(:,boundary_temp1:boundary_temp2);
            data_new_good_length=length(find(isnan(data_new_boundary(1,:))==0));
            auto_epoch_inspection.data_length_after_autobreakautocleaning_Nsamples=data_new_good_length;
            auto_epoch_inspection.data_length_after_autobreakautocleaning_min=data_new_good_length/EEG.srate/60;
        end
        
        if n==length(appended_file_member)+1
            boundary_temp=EEG.etc.automatic_cleaning.boundaries_between_files_raw_firstGoodSampleOfNewFile(n-1);
            auto_epoch_inspection.data_length_raw_Nsamples=length(boundary_temp:size(EEG.data,2));
            auto_epoch_inspection.data_length_raw_min=length(boundary_temp:size(EEG.data,2))/EEG.srate/60;
            
            data_new_boundary=data_new(:,boundary_temp:size(data_new,2));
            data_new_good_length=length(find(isnan(data_new_boundary(1,:))==0));  %%% how many good samples?
            auto_epoch_inspection.data_length_after_autobreakautocleaning_Nsamples=data_new_good_length;
            auto_epoch_inspection.data_length_after_autobreakautocleaning_min=data_new_good_length/EEG.srate/60;
        end
        
        
        eval(['data_new' num2str(n) '=data_new;']); clear data_new
        eval(['auto_epoch_inspection' num2str(n) '=auto_epoch_inspection;']); clear auto_epoch_inspection
        
        eval(['window_samples_vector_file_' num2str(n) '=window_samples_vector;'])
    end
    
    %     %%% now merge bad segments from all appended files
    %     data_new=data;
    %     for n=1:length(appended_file_member)+1
    %         eval(['bad_segments_index=auto_epoch_inspection' num2str(n) '.window_samples_vector_final_final;'])
    %         bad_segments_index(find(bad_segments_index(:,3)==0),:)=[];
    %         for nn=1:size(bad_segments_index,1)
    %             data_new(:,bad_segments_index(nn,1):bad_segments_index(nn,2))=NaN;
    %         end
    %     end
    
else
    %% 4.2 look for "bad" epochs in a single file
    crit=crit_all;
    window_samples_vector=windows;
    epoch_index_appended_file=1:size(window_samples_vector,1);
    
    appended_file_member=zeros(1,size(windows,1));
    for nn=1:size(windows,1)
        if ismember(size(EEG.data,2), windows(nn,1):windows(nn,2))
            appended_file_member(nn)=1; %%% check which epoch contains the boundary marker between appended files
        end
    end
    appended_file_member=find(appended_file_member==1);
    
    n=1;
    
    % reject epochs, see script for details, plots figure
    crit=crit_all(n);
    sampInWdw=(wind_ms/1000)*EEG.srate;
    
    [auto_epoch_inspection] = EEG_find_artifact_indices_epoch_rejection(data_segmented,data_segmented_raw,...
        window_samples_vector, crit, crit_percent_sample_epoch, sampInWdw, EEG.srate,...
        weights4methods, selected_chans4cleaning, datapath_save_figures, subject);
    
    window_samples_vector=auto_epoch_inspection.window_samples_vector_final_final;
    windows=window_samples_vector;
    
    % saved directly in EEG_find_artifact_indices_epoch_rejection now (fragmented into multiple figures)
    %savefig([datapath_save_figures filesep 'autocleaning_appended_files_all'])
    %close
    
    %     %%% now merge bad segments from all appended files
    %     data_new=data;
    %     bad_segments_index=auto_epoch_inspection.window_samples_vector_final_final;
    %     bad_segments_index(find(bad_segments_index(:,3)==0),:)=[];
    %     for nn=1:size(bad_segments_index,1)
    %         data_new(:,bad_segments_index(nn,1):bad_segments_index(nn,2))=NaN;
    %     end
    
end


%% 5) keep only "good" data segments which are uninterrupted and longer than criterion in seconds

%%% thus avoid extremely segmented data (which would be suboptimal,
%%% e.g., for calculating long-range temporal correlations)

if length(crit_all) > 1
    eval(['window_samples_vector_all=auto_epoch_inspection' num2str(1) '.window_samples_vector_final_final;'])
    for n=2:length(crit_all)  %%% append all epochs back for continuous data
        eval(['window_samples_vector_all' num2str(n) '=auto_epoch_inspection' num2str(n) '.window_samples_vector_final_final;'])
        eval(['window_samples_vector_all=[window_samples_vector_all; window_samples_vector_all' num2str(n) '];'])
    end
end

crit_keep_samples=(crit_keep_sec*1000)/(1000/EEG.srate);
auto_epoch_inspection.keep_clean_segments_minimum_length_sec=crit_keep_sec;
if mod(crit_keep_samples/sampInWdw,2)>0 & mod(crit_keep_samples/sampInWdw,2)<1
    error('criterion should be multiple of window length')
end

%%% check whether after a bad epoch you have an uninterrupted segment of
%%% critical length; if not, mark as NaN until the next bad epoch
bad_epoch=find(windows(:,3)==1);
data_new=data;
bad_segments_index=windows;
bad_segments_index(find(bad_segments_index(:,3)==0),:)=[];
for nn=1:size(bad_segments_index,1)
    data_new(:,bad_segments_index(nn,1):bad_segments_index(nn,2))=NaN;
end

for epoch=1:length(bad_epoch)-1
    temp1=windows(bad_epoch(epoch),:);
    temp2=windows(bad_epoch(epoch+1),:);
    if diff([temp1(1) temp2(1)]) < crit_keep_samples   %%% mark as NaN in the data
        data_new(:,windows(bad_epoch(epoch),1):windows(bad_epoch(epoch+1),2))=NaN;
    end
end

%%% check if any data still available
if sum(isnan(data_new(1,:))) == size(data_new,2)
    error('all data removed; check removal criteria, especially "crit_keep_sec"')
end

%%% how much data would be removed with respect to the raw data?
%%% if you have different appended files, now you have also removed epochs
%%% which were very close, see step 6) above; thus you need to adjust
%%% information on the final data length in each of the appended files
if length(crit_all) > 1
    for n=1:length(appended_file_member)+1
        if n==1
            eval(['auto_epoch_inspection=auto_epoch_inspection' num2str(n) ';'])
            boundary_temp=EEG.etc.automatic_cleaning.boundaries_between_files_raw_firstGoodSampleOfNewFile(n)-1;
            
            data_new_boundary=data_new(:,1:boundary_temp);
            data_new_good_length=length(find(isnan(data_new_boundary(1,:))==0));  %%% how many good samples?
            auto_epoch_inspection.data_length_final_after_autobreakautoclean_keepcleanseg_Nsamp=data_new_good_length;
            auto_epoch_inspection.data_length_final_after_autobreakautoclean_keepcleanseg_min=data_new_good_length/EEG.srate/60;
            eval(['auto_epoch_inspection' num2str(n) '=auto_epoch_inspection;'])
        end
        
        if n>1 & n<=length(appended_file_member)
            eval(['auto_epoch_inspection=auto_epoch_inspection' num2str(n) ';'])
            boundary_temp1=EEG.etc.automatic_cleaning.boundaries_between_files_raw_firstGoodSampleOfNewFile(n-1);
            boundary_temp2=EEG.etc.automatic_cleaning.boundaries_between_files_raw_firstGoodSampleOfNewFile(n)-1;
            
            data_new_boundary=data_new(:,boundary_temp1:boundary_temp2);
            data_new_good_length=length(find(isnan(data_new_boundary(1,:))==0));
            auto_epoch_inspection.data_length_final_after_autobreakautoclean_keepcleanseg_Nsamp=data_new_good_length;
            auto_epoch_inspection.data_length_final_after_autobreakautoclean_keepcleanseg_min=data_new_good_length/EEG.srate/60;
            eval(['auto_epoch_inspection' num2str(n) '=auto_epoch_inspection;'])
        end
        
        if n==length(appended_file_member)+1
            eval(['auto_epoch_inspection=auto_epoch_inspection' num2str(n) ';'])
            boundary_temp=EEG.etc.automatic_cleaning.boundaries_between_files_raw_firstGoodSampleOfNewFile(n-1);
            auto_epoch_inspection.data_length_raw_Nsamples=length(boundary_temp:size(EEG.data,2));
            auto_epoch_inspection.data_length_raw_min=length(boundary_temp:size(EEG.data,2))/EEG.srate/60;
            
            data_new_boundary=data_new(:,boundary_temp:size(data_new,2));
            data_new_good_length=length(find(isnan(data_new_boundary(1,:))==0));  %%% how many good samples?
            auto_epoch_inspection.data_length_final_after_autobreakautoclean_keepcleanseg_Nsamp=data_new_good_length;
            auto_epoch_inspection.data_length_final_after_autobreakautoclean_keepcleanseg_min=data_new_good_length/EEG.srate/60;
            eval(['auto_epoch_inspection' num2str(n) '=auto_epoch_inspection;'])
        end
    end
end


%% 6) keep indices of bad data segments

if length(crit_all) > 1
    for n=1:length(crit_all)
        eval(['EEG.etc.automatic_cleaning.auto_epoch_rejection_appendedfile' num2str(n) '=auto_epoch_inspection' num2str(n) ';'])
    end
else
    EEG.etc.automatic_cleaning.auto_epoch_rejection_singlefile=auto_epoch_inspection;
end


EEG.etc.automatic_cleaning.auto_epoch_rejection=ones(1,size(EEG.data,2));
EEG.etc.automatic_cleaning.auto_epoch_rejection(find(isnan(data_new(1,:))==1))=0;  %%% take channel 1; it's the same for all channels anyway
if isfield(EEG.etc.automatic_cleaning,'event_latencyNote')
    EEG.etc.automatic_cleaning.auto_epoch_rejection_Note=EEG.etc.automatic_cleaning.event_latencyNote;
end
EEG.etc.automatic_cleaning.auto_epoch_rejection_invalid_segments_start_stop_sample=windows(find(windows(:,3)==1),1:2);

%%% find final bad segments
temp=EEG.etc.automatic_cleaning.auto_epoch_rejection;
temp_diff=diff(temp);
x_stop=find(temp_diff==1);
x_start=find(temp_diff==-1);
x_start=x_start+1;

if temp(1)==0  %%% bad segment starts directly at the first sample
    x_start=[1 x_start];
end

if temp(end)==0  %%% bad segment ends at last sample
    x_stop=[x_stop length(temp)];
end

remove_final=[x_start;x_stop]';
%%% contains all bad segments: of automatic break rejection and of continuous data cleaning


%%% security check if final start/stops are correct
%%% should be identical to those in EEG.etc.automatic_cleaning.auto_epoch_rejection
data1=EEG.data(1,:);
for n=1:size(remove_final,1)
    data1(remove_final(n,1):remove_final(n,2))=NaN;
end
temp2=find(isnan(data1)==1);
if length(unique([sum(ismember(find(EEG.etc.automatic_cleaning.auto_epoch_rejection==0),temp2)) length(find(EEG.etc.automatic_cleaning.auto_epoch_rejection==0)) length(temp2)])) > 1
    error ('artifact removal not consistent')
end


%%% these are the finally rejected segments from AUTO break rejection and
%%% epoch cleaning, as performed in the script here
EEG.etc.automatic_cleaning.invalid_segments_final_start_stop_sample=remove_final;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 7) plot all channels of the finally cleaned data, save figure

%%% mark invalid data
for nn=1:size(remove_final,1)
    data_original4plot(:,remove_final(nn,1):remove_final(nn,2))=NaN;
end

figure('units','normalized','outerposition',[0 0 1 1])
timesMin=EEG.times/(1000*60);
ax=[timesMin(1) timesMin(end) -500 500];
colorspec=cool(size(data_original4plot,1));
hold on
for chan=1:size(data_original4plot,1)
    plot(timesMin,data_original4plot(chan,:),'color',colorspec(chan,:));
end
axis(ax)
xlabel({'Time [min]',['Length of raw data set: ' num2str(EEG.xmax/60) ' min'],...
    ['Length of cleaned data set = ' num2str(sum(isnan(data_original4plot(1,:))==0)/(EEG.srate*60)) ' min']})
ylabel('Potential [??V]')
if ~isempty(band)
    title([subject, ' after cleaning - Bandpass ['...
        num2str(band(1)) ';' num2str(band(2)) '] Hz - all channels superimposed'])
else
    title([subject, ' after cleaning - Raw data - all channels superimposed'])
end

if isfield(EEG.etc.automatic_cleaning,'boundaries_between_files_raw_firstGoodSampleOfNewFile')
    plot([timesMin(EEG.etc.automatic_cleaning.boundaries_between_files_raw_firstGoodSampleOfNewFile(1)) timesMin(EEG.etc.automatic_cleaning.boundaries_between_files_raw_firstGoodSampleOfNewFile(1))],[ax(3) ax(4)],'r');
    for nb=1:length(EEG.etc.automatic_cleaning.boundaries_between_files_raw_firstGoodSampleOfNewFile)
        plot([timesMin(EEG.etc.automatic_cleaning.boundaries_between_files_raw_firstGoodSampleOfNewFile(nb)) timesMin(EEG.etc.automatic_cleaning.boundaries_between_files_raw_firstGoodSampleOfNewFile(nb))],[ax(3) ax(4)],'k');
        h=text('Position',[timesMin(EEG.etc.automatic_cleaning.boundaries_between_files_raw_firstGoodSampleOfNewFile(nb))  ax(4)-ax(4)*0.2],'string',EEG.etc.appended_files(nb+1)); set(h,'rotation',80);
    end
end

%% save figure
saveCurrentFig(fullfile(datapath_save_figures, filesep), sprintf('%s_continuous_data_autoclean_after', subject),{'png'},[])
% savefig([datapath_save_figures '\autocleaning_final'])
%print('-djpeg',[datapath_save_figures '\autocleaning_final'],'-r300') %%% too large file for saving as .fig
%close
end