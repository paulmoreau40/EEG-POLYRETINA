function[auto_epoch_cleaning]=EEG_find_artifact_indices_epoch_rejection(data, data_raw,...
    windows_samp, crit, crit_percent_sample_epoch, sampInwdw, srate,...
    weighting_factor, selected_channels, datapath_save_figures, subject)
%%% performs fully automatic epoch rejection on absolute values
%%% currently: 3 methods (mean, SD, Mahalanobian distance)
%%% creates a joint measure (by weighting factors) of a "bad" epoch (based on all methods)
%%% removes defined % of "worst" epochs (based on the joint measure)
%%% RECOMMENDED: apply amplitude envelope (e.g., Hilbert transform) of oscillatory data and then create epochs for cleaning
%%% creates figure
%%%
%%% USAGE:
%%%   >>  [auto_epoch_cleaning]=EEG_find_artifact_indices_epoch_rejection(data,data_raw,window_samples_vector,crit,crit_percent_sample_epoch,wind_samples,srate,weighting_factor,selected_channels)
%%%
%%% INPUTS:
%%% data        - epoched data (can be filtered; absolute values will be considered for automatic cleaning)
%%%               with the structure channels x samples x epochs; can be cell
%%% data_raw    - epoched data, raw (inspect raw data for flat lines); if not available, enter "data", nothing will happen
%%%               with the structure channels x samples x epochs; can be cell
%%% window_samples_vector - with start x end (samples) for each epoch as rows,
%%%                if you enter cell data from previously epoched continuous data (not done with EEGLAB epoching)
%%%                else: window_samples_vector=[] if data is already epoched
%%% crit        - keep amount of valid epochs; e.g. 0.9 for keeping 90 %
%%%               "valid" epochs means: epochs without NaN or "flat lines" (zeros)
%%%               of available valid epochs then, e.g., take 90 %
%%% crit_percent_sample_epoch -  [] for nothing; e.g., 0.5 for removing epochs if they have more than 50 % of samples with zero or NaN
%%% wind_samples - epoch length in samples; just used for plotting; if you have different epochs, enter average length
%%% srate        - sampling rate
%%% weighting_factor - weighting for each removal method; =0 for not applying indicated method
%%%                weighting_factor=[1 1 1] for equally weigthing all methods
%%%                put at zero for using only selected methods, e.g., weighting_factor=[0 1 2] or [1.5 0 0]
%%% selected_channels - [] for using all available channels in data; else specify channels that you want to use for cleaning
%%%
%%% OUTPUTS:
%%% auto_epoch_cleaning - information on the cleaned epochs
%%%
%%% most important output: epoch sample (start,end) x use (=1 remove, 0=keep epoch)
%%% auto_epoch_cleaning.window_samples_vector_methods_joined_final   <- if different epochs were found by all methods
%%% auto_epoch_cleaning.window_samples_vector_final                  <- if same epochs were found by all methods
%%%
%%%%%%%%%%%%%%% steps
%%% 0) basic settings
%%% 1) epoch evaluation method I: take MEAN in each epoch and take mean ACROSS channels
%%% 2) epoch evaluation method II: channel heterogeneity --> for each epoch_mean calculate SD across channels
%%% 3) epoch evaluation method III: channel heterogeneity --> calculate Mahalanobian distance across epoch_mean across channels
%%% 4) plot all epochs; NaN indicate previously rejected breaks
%%% 5) indicate bad epochs in the epoch vector
%%% 6) indicate bad epochs in the data as NaN and plot
%%%
%%%
%%% internal use only ;)
%%% FH, 2017

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 0. basic settings
if isempty(weighting_factor) || sum(weighting_factor)==0
    error('no weighting factors inserted')
end

if iscell(data)
    N_epochs=size(data,2);
else
    N_epochs=size(data,3);
end
disp(['number of epochs to examine = ' num2str(N_epochs) ''])

if isempty(windows_samp) %%% no vector, then created automatically
    windows_samp(1:N_epochs,1)=1;
    windows_samp(1:N_epochs,2)=sampInwdw;
end

wind_ms=(sampInwdw/srate)*1000;

auto_epoch_cleaning=[];
auto_epoch_cleaning.epochs_start_stop_samples=windows_samp;

%%% window_samples_vector is created by function "EEG_find_artifact_indices"
%%% rows = epochs, columns = start (samples) end (samples) of epoch
%%% column 3 indicates 1=remove this epoch, 0=keep epoch
%%% copy values for 3 used methods
window_samples_method1=windows_samp;  %%% for MEAN in epoch
window_samples_method2=windows_samp;  %%% for SD of epoch_MEAN across channels
window_samples_method3=windows_samp;  %%% for Mahal. distance of epoch_MEAN across channels
window_samples_methods_joined=windows_samp;  %%% for joining epochs of all methods

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 1. epoch evaluation I: take MEAN in each epoch and take mean ACROSS channels
%%%   mark epoch also as NaN if more than certain percentage of sample contains NaN or zeros (flat line)
auto_epoch_cleaning.removed_epochs_sorted_previous_artifact_rejection_applied=[];
auto_epoch_cleaning.removed_epochs_sorted_flat_line_detection=[];

x1=0; x2=0;
for epoch=1:N_epochs
    
    if iscell(data)  %%% if epochs are of different lengths
        temp=data{epoch};
        temp_raw=data_raw{epoch};
    else
        temp=squeeze(data(:,:,epoch));   %%% data should have the format channels x samples x epochs
        temp_raw=squeeze(data_raw(:,:,epoch));
    end
    
    
    if epoch==1
        selected_channels=sort(selected_channels);
        if isempty(selected_channels)
            warning('no channels selected for epoch cleaning...using all')
            selected_channels=1:size(temp_raw,1);
        else
            selected_channels=sort(selected_channels);
        end
    end
    
    %%% select channels that should be analyzed
    remove_channels_for_cleaning=setdiff(1:size(temp,1),selected_channels);
    temp(remove_channels_for_cleaning,:)=[];
    temp_raw(remove_channels_for_cleaning,:)=[];
    
    if sum([length(unique([size(temp,1) size(temp_raw,1)])) length(unique([size(temp,2) size(temp_raw,2)]))]) > 2
        error('channels do not match between data sets')
    end
    
auto_epoch_cleaning.used_channel_indices = selected_channels;
    
    for chan=1:size(temp,1)
        
        %%% for each channel, take mean across samples (absolute values)
        auto_epoch_cleaning.epoch_mean(chan,epoch)=mean(abs(temp(chan,:)),2);
        auto_epoch_cleaning.epoch_mean_raw(chan,epoch)=mean(abs(temp_raw(chan,:)),2);  %%% is not touched, for plotting later
        
        
        %         %%% alternatively: for each channel, square the log of the mean across samples (absolute values)
        %         auto_epoch_cleaning.epoch_mean(chan,epoch)=log(mean(abs(temp(chan,:)),2))^2;  %%% thus, amplify deviations
        %         auto_epoch_cleaning.epoch_mean_raw(chan,epoch)=log(mean(abs(temp_raw(chan,:)),2))^2;  %%% is not touched, for plotting later
        %
        
        if ~isempty(crit_percent_sample_epoch) %%% if criterion is set
            ttemp=temp(chan,:);  %%% check in each channel
            x1=x1+1;
            if length(find(isnan(ttemp)==1))/length(ttemp) > crit_percent_sample_epoch  %%% remove epoch, if NaN because data was removed from previous automatic artifact rejection (e.g., automatic break rejection)
                disp(['reject epoch, not enough data in epoch #' num2str(epoch)])
                auto_epoch_cleaning.removed_epochs_sorted_previous_artifact_rejection_applied(x1)=epoch;
                auto_epoch_cleaning.epoch_mean(chan,epoch)=NaN;  %%% take mean across epoch for each channel separately
            end
            
            %%% find epochs with "flat line": those are always removed, won't go into calculation of % of epochs according to criterion
            %%% TO DO: now it looks only for 0; in the future, look, e.g., for identical values, in case the flat line is not at zero
            %%% check in each channel; if flat line found, then epoch is removed from all channels alike
            ttemp=temp_raw(chan,:); %%% use RAW EEG here to find flat lines, because by filtering and Hilbert transform you might lose this info
            x2=x2+1;
            if length(find(ttemp==0))/length(ttemp) > crit_percent_sample_epoch
                disp(['flat line detected in epoch #' num2str(epoch)])
                auto_epoch_cleaning.epoch_mean(chan,epoch)=NaN;
                auto_epoch_cleaning.removed_epochs_sorted_flat_line_detection(x2)=epoch;
            end
        end
        auto_epoch_cleaning.removed_epochs_sorted_previous_artifact_rejection_applied=unique(auto_epoch_cleaning.removed_epochs_sorted_previous_artifact_rejection_applied);
        auto_epoch_cleaning.removed_epochs_sorted_flat_line_detection=unique(auto_epoch_cleaning.removed_epochs_sorted_flat_line_detection);
    end
end
if ~isempty(auto_epoch_cleaning.removed_epochs_sorted_previous_artifact_rejection_applied)
    auto_epoch_cleaning.removed_epochs_sorted_previous_artifact_rejection_applied=auto_epoch_cleaning.removed_epochs_sorted_previous_artifact_rejection_applied(find(auto_epoch_cleaning.removed_epochs_sorted_previous_artifact_rejection_applied>0));
end
if ~isempty(auto_epoch_cleaning.removed_epochs_sorted_flat_line_detection)
    auto_epoch_cleaning.removed_epochs_sorted_flat_line_detection=auto_epoch_cleaning.removed_epochs_sorted_flat_line_detection(find(auto_epoch_cleaning.removed_epochs_sorted_flat_line_detection>0));
end

%%% mean in epoch: take mean across channels
auto_epoch_cleaning.epoch_mean_across_channels=[];
for N_epoch=1:size(auto_epoch_cleaning.epoch_mean,2)
    ind=find(isnan(auto_epoch_cleaning.epoch_mean(:,N_epoch))==0);
    auto_epoch_cleaning.epoch_mean_across_channels(N_epoch)=mean(auto_epoch_cleaning.epoch_mean(ind,N_epoch));
end

%%% mean in epoch: take mean across epochs
auto_epoch_cleaning.epoch_mean_across_epochs=[];
for N_chan=1:size(auto_epoch_cleaning.epoch_mean,1)
    ind=find(isnan(auto_epoch_cleaning.epoch_mean(N_chan,:))==0);
    auto_epoch_cleaning.epoch_mean_across_epochs(N_chan)=mean(auto_epoch_cleaning.epoch_mean(N_chan,ind));
end


%%% now find epochs with largest mean value (across channels),
%%% based on % criterion
if size(auto_epoch_cleaning.epoch_mean,1)>1
    temp=auto_epoch_cleaning.epoch_mean_across_channels;  %%% mean across channels per epoch
else
    temp=auto_epoch_cleaning.epoch_mean;  %%% if only 1 channel was selected, mean cannot be taken
end
if length(find(isnan(temp)==1)) < length(temp)  %%% if any good epoch in the appended file
    [i1,i2]=sort(temp);  %%% sort in ascending order; i1 are values, i2 are epoch indices
    if ~isempty(find(isnan(i1)==1)) %%% in case you had removed segments from previous automatic break rejection
        i1=i1(find(isnan(i1)==0));
        i2=i2(find(isnan(i1)==0));
    end
    i1_method1=i1; %%% copy for plotting later
    i2_method1=i2;
    
    %%% define threshold, i.e., how many epochs should be kept
    %%% based on the number of VALID epochs (i.e., without NaN/zeros, latter defined by crit_percent_sample_epoch)
    all_valid_epochs=length(i2);   %%% number of all available epochs
    auto_epoch_cleaning.valid_epochs_number=all_valid_epochs;
    
    thrshold=round(length(i2)*crit);  %%% number of valid epochs that should be kept (according to % criterion), exluding flat lines (those are removed independently)
    bad_epoch_unsorted=i2(thrshold+1:end); bad_epoch=sort(bad_epoch_unsorted);
    auto_epoch_cleaning.removed_valid_epochs_number=all_valid_epochs-thrshold;
else
    disp('no valid epochs found')
    bad_epoch_unsorted=[];
    i1_method1=[];
    i2_method1=[];
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 2) epoch evaluation II: channel heterogeneity --> for each epoch_mean calculate SD ACROSS channels
%%%   i.e.,operationalize "channel heterogeneity". In each epoch you have
%%%   taken the mean value, separately for each channel. Now calculate in
%%%   each epoch the SD across channels, i.e., if among all channels there
%%%   are some channels which have a LARGE mean (most likely artifact) this will increase the SD

if size(auto_epoch_cleaning.epoch_mean,1)>1  %%% you have more than 1 channel
    for epoch=1:size(auto_epoch_cleaning.epoch_mean,2)
        auto_epoch_cleaning.epoch_mean_std_across_channels(epoch)=std(auto_epoch_cleaning.epoch_mean(:,epoch));
    end
    %%% now indicate epochs with largest SD across channels based on their MEAN within epoch
    %%% based on % criterion
    temp=auto_epoch_cleaning.epoch_mean_std_across_channels;
    [i1,i2]=sort(temp);  %%% i1 are values, i2 are epoch indices
    if ~isempty(find(isnan(i1)==1)) %%% in case you had removed segments from automatic break rejection
        i1=i1(find(isnan(i1)==0));
        i2=i2(find(isnan(i1)==0));
    end
    i1_method2=i1; %%% copy for plotting later
    i2_method2=i2;
    bad_epoch_2_unsorted=i2(thrshold+1:end); bad_epoch_2=sort(bad_epoch_2_unsorted);
else
    auto_epoch_cleaning.epoch_mean_std_across_channels_Note='mean SD across channels not calculated; only 1 channel available';
    bad_epoch_2_unsorted=[];
    bad_epoch_2=[];
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 3) epoch evaluation method III: channel heterogeneity --> calculate Mahalanobian distance across epoch_mean (single channels)

if size(auto_epoch_cleaning.epoch_mean,1)>1  %%% you have more than 1 channel
    
    temp=auto_epoch_cleaning.epoch_mean;
    %%% remove NaN for calculating Mahal., otherwise it wouldn't work
    %%% if there are any bad epochs from previous flat line and/or artifact detection
    all_epochs_temp=1:size(auto_epoch_cleaning.epoch_mean,2);
    remove_temp=[];
    x1=0;
    for epoch=1:size(temp,2)
        if sum(isnan(temp(:,epoch))) > 0   %%% if any channel in this epoch is bad, remove the whole epoch
            x1=x1+1;
            remove_temp(x1)=epoch;
        end
    end
    if ~isempty(remove_temp)
        temp(:,remove_temp)=[];
        keep_temp=setdiff(all_epochs_temp,remove_temp);  %%% those epoch indices refer to those calculated in M
    end
    
    temp_size=size(temp);
    
    if temp_size(2) > temp_size(1)
        [M]=mahal(temp',temp'); M=M';   %%% get an index of the Mahal. distance for each epoch
        if ~isempty(remove_temp)        %%% insert back NaN to be able to work with the full epoch set
            M2=NaN(1,length(all_epochs_temp)); M2(keep_temp)=M; M=M2; clear M2
        end
    else
        M=zeros(1,size(temp,2));  %%% Mahal. distance needs more columns (here epochs) than rows (here channels);
    end                           %%% if too few epochs, then insert 0
    
    
    auto_epoch_cleaning.epoch_mean_Mahal_distance_acrossChannels=M;
    
    if temp_size(2) > temp_size(1)
        auto_epoch_cleaning.epoch_mean_Mahal_distance_acrossChannels_Note=[];
        [i1,i2]=sort(M);  %%% i1 are values, i2 are epoch indices
        if ~isempty(find(isnan(i1)==1)) %%% in case you had removed segments from automatic break rejection
            i1=i1(find(isnan(i1)==0));
            i2=i2(find(isnan(i1)==0));
        end
        i1_method3=i1; %%% copy for plotting later
        i2_method3=i2;
        bad_epoch_3_unsorted=i2(thrshold+1:end); bad_epoch_3=sort(bad_epoch_3_unsorted);
    else
        auto_epoch_cleaning.epoch_mean_Mahal_distance_acrossChannels_Note={'too few epochs, therefore method not applied; 0 inserted'};
        bad_epoch_3_unsorted=[]; bad_epoch_3=[];
    end
else
    auto_epoch_cleaning.epoch_mean_Mahal_distance_acrossChannels_Note='MD across channels not calculated; only 1 channel available';
    bad_epoch_3_unsorted=[];
    bad_epoch_3=[];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 4) plot all epochs; NaN indicate previously rejected breaks
%%%  (and, if present, epochs with flat line defined % of samples)

figure
subplot(2,2,1)
plot(auto_epoch_cleaning.epoch_mean');
title({'Epoch mean', ['over ' num2str(length(auto_epoch_cleaning.used_channel_indices)) ' SINGLE chans']});
axis([0 size(auto_epoch_cleaning.epoch_mean,2) 0 inf]);
xlabel({['Epoch #; len: ' num2str(wind_ms) 'ms'],...
    ['checked epochs = ' num2str(size(auto_epoch_cleaning.epoch_mean,2)), ...
    ', valid for cleaning = ' num2str(num2str(auto_epoch_cleaning.valid_epochs_number))]});
ylabel('Mean (microV)');

subplot(2,2,2)
plot(auto_epoch_cleaning.epoch_mean_across_channels);
title({'Epoch mean', 'ACROSS channels'});
axis([0 length(auto_epoch_cleaning.epoch_mean) 0 inf]);
xlabel({['Epoch #; len: ' num2str(wind_ms) 'ms']});
ylabel('Mean (microV)')

subplot(2,2,3)
if ~isempty(bad_epoch_2)
    plot(auto_epoch_cleaning.epoch_mean_std_across_channels);
    title({'SD of single channel epoch mean', 'ACROSS channels'});
    axis([0 length(auto_epoch_cleaning.epoch_mean_std_across_channels) 0 inf]);
    xlabel({['Epoch #; len: ' num2str(wind_ms) 'ms']});
    ylabel('SD')
else
    title({'SD across channels N.A.', '(only 1 channel available)'});
end

subplot(2,2,4)
if ~isempty(bad_epoch_3)
    if isempty(auto_epoch_cleaning.epoch_mean_Mahal_distance_acrossChannels_Note) %%% if analysis was performed
        plot(auto_epoch_cleaning.epoch_mean_Mahal_distance_acrossChannels);
        title({'Mahalanobian distance across channels', 'based on single channel epoch mean'});
        axis([0 length(auto_epoch_cleaning.epoch_mean_Mahal_distance_acrossChannels) 0 inf]);
        xlabel({['Epoch #; len: ' num2str(wind_ms) 'ms']});
        ylabel('Mahalanobian distance')
    else
        title({'Mahalanobian distance not applicable','(too few epochs or only 1 channel available)'})
    end
else
    title({'Mahalanobian distance not applicable','(too few epochs or only 1 channel available)'});
end

sgtitle([subject, ' Measures report for bad segment removal'])
%suptitle([subject, ' Measures report for bad segment removal'])
saveCurrentFig(fullfile(datapath_save_figures, filesep),...
    sprintf('%s_autocleaning_measures_report',subject), {'fig','png'},[1100 800])
%savefig([datapath_save_figures filesep 'autocleaning_measures_report'])
%close

figure
if length(crit)>1
    subplot(2,2,1)
    title(['evaluated condition: ' EEG.etc.appended_files{n} ''])
end

subplot(2,2,2)
hold on
plot(i1_method1);
plot([thrshold thrshold],[0 i1_method1(end)],'r');
axis([0 length(i1_method1) 0 i1_method1(end)]);
xlabel('Sorted epochs');
if ~isempty(bad_epoch_2)
    ylabel({'Epoch mean', 'SINGLE channels (more N.A.)'});
else
    ylabel({'Epoch mean', 'ACROSS channels'});
end
legend({'Epoch distribution', ['threshold: ' num2str(crit*100) '% of epochs']})
title({'Method I:', 'Largest epoch mean'})

if ~isempty(bad_epoch_2)
    subplot(2,2,3)
    hold on
    plot(i1_method2);
    plot([thrshold thrshold],[0 i1_method2(end)],'r');
    axis([0 length(i1_method2) 0 i1_method2(end)])
    xlabel('Sorted epochs');
    ylabel('SD of single channel epoch mean');
    legend({'Epoch distribution', ['threshold: ' num2str(crit*100) '% of epochs']})
    title({'Method II:', 'Largest Standard Deviation'})
end

if ~isempty(bad_epoch_3)
    subplot(2,2,4)
    hold on
    plot(i1_method3);
    plot([thrshold thrshold],[0 i1_method3(end)],'r');
    axis([0 length(i1_method3) 0 i1_method3(end)])
    xlabel('Sorted epochs');
    ylabel('Mahal. distance across channels');
    legend({'Epoch distribution', ['threshold: ' num2str(crit*100) '% of epochs']})
    title({'Method III:', 'Mahal. distance across channels',...
        'based on single channel epoch mean'})
end

sgtitle([subject, ' Epochs distribution and selection depending on methods'])
%suptitle([subject, ' Epochs distribution and selection depending on methods'])
saveCurrentFig(fullfile(datapath_save_figures, filesep),...
    sprintf('%s_autocleaning_distributions_methods', subject), {'fig','png'},[1100 800])
%savefig([datapath_save_figures filesep 'autocleaning_distributions_methods'])
%close

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 5) indicate bad epochs in the epoch vector
%%%   1= bad epoch and as NaN in the data

if length(find(isnan(temp)==1)) < length(temp)  %%% if any good epoch in the appended file
    
    auto_epoch_cleaning.found_epochs_unsorted_method_I_mean=bad_epoch;
    auto_epoch_cleaning.found_epochs_unsorted_method_II_std=bad_epoch_2;
    auto_epoch_cleaning.found_epochs_unsorted_method_III_Mahal=bad_epoch_3;
    
    auto_epoch_cleaning.window_samples_vector_note='0 = keep epoch, 1 = remove epoch';
    
    
    %%% 5.1) indicate bad epochs in the final vector separately for each method
    window_samples_method1(bad_epoch,3)=1;  %%% =1 remove this epoch
    window_samples_method1(auto_epoch_cleaning.removed_epochs_sorted_previous_artifact_rejection_applied,3)=1;  %%% =1 remove this epoch
    window_samples_method1(auto_epoch_cleaning.removed_epochs_sorted_flat_line_detection,3)=1;
    auto_epoch_cleaning.window_samples_vector_method1=window_samples_method1;
    
    if ~isempty(bad_epoch_2)
        window_samples_method2(bad_epoch_2,3)=1;  %%% =1 remove this epoch
    end
    window_samples_method2(auto_epoch_cleaning.removed_epochs_sorted_previous_artifact_rejection_applied,3)=1;  %%% =1 remove this epoch
    window_samples_method2(auto_epoch_cleaning.removed_epochs_sorted_flat_line_detection,3)=1;
    auto_epoch_cleaning.window_samples_vector_method2=window_samples_method2;
    
    if ~isempty(bad_epoch_3)
        window_samples_method3(bad_epoch_3,3)=1;  %%% =1 remove this epoch
    end
    window_samples_method3(auto_epoch_cleaning.removed_epochs_sorted_previous_artifact_rejection_applied,3)=1;  %%% =1 remove this epoch
    window_samples_method3(auto_epoch_cleaning.removed_epochs_sorted_flat_line_detection,3)=1;
    auto_epoch_cleaning.window_samples_vector_method3=window_samples_method3;
    
    auto_epoch_cleaning.weighting_factor_method_1=weighting_factor(1);
    auto_epoch_cleaning.weighting_factor_method_2=weighting_factor(2);
    auto_epoch_cleaning.weighting_factor_method_3=weighting_factor(3);
    
    
    %%% 5.2) identical epochs were found by all methods, thus you can
    %%%      keep all epochs in order to match the requirement of % of epochs that should be kept (crit)
    if length(find(window_samples_method1(:,3)+window_samples_method2(:,3)+window_samples_method3(:,3)==3)) == size(window_samples_method1,1)
        disp('identical bad epochs found by all methods')
        window_samples_method1(bad_epoch,3)=1;
        auto_epoch_cleaning.window_samples_vector_final=window_samples_method1;
        auto_epoch_cleaning.removed_epochs_final_sorted=sort(find(window_samples_method1(:,3)==1));
        auto_epoch_cleaning.removed_epochs_method_note='all bad epochs found by all applied methods were removed';
        auto_epoch_cleaning.window_samples_vector_methods_joined_final=[];
        
        %%% 5.3) dissimilar epochs found by some/all methods
    else
        auto_epoch_cleaning.window_samples_vector_final=[];
        bad_epochs_all=[];
        
        disp('dissimilar bad epochs found by all methods, joining according to weighting factor')
        auto_epoch_cleaning.removed_epochs_method_note='bad epochs were removed by weighting epochs found by the applied methods';
        
        %%% 5.3.1) create variable containing epochs as rows and indication of
        %%%        being "bad" found by the different methods starting with column 3
        %%% (=1 epoch is good, =2 epoch is "bad")
        index_all_bad_epochs=sort(unique([bad_epoch bad_epoch_2 bad_epoch_3]));
        for epoch=1:length(index_all_bad_epochs)
            bad_epochs_all(epoch,1)=index_all_bad_epochs(epoch);  %%% index of found epoch by any applied method
        end
        
        %%% method I
        for epoch=1:length(bad_epoch_unsorted)  %%% give each bad epoch a value = largest for worst epoch; largest value is defined by length(removed_epochs_method_1)
            x1=find(bad_epochs_all(:,1)==bad_epoch_unsorted(epoch));
            bad_epochs_all(x1,2)=epoch;  %%% "weighting" of epoch
        end
        %%% e.g., by method I 100 epochs should be removed, their order, sorted according to largest mean,
        %%% could be [...10 5 8 9 98 4] where epoch 4 is the worst out of
        %%% 100 epochs; this epoch would get weigth=100 is worst epoch; epoch 98 get 99, epoch 9 get 98 and so on
        
        %%% method II
        if ~isempty(bad_epoch_2)
            for epoch=1:length(bad_epoch_2_unsorted)  %%% give each bad epoch a value = largest for worst epoch; largest value is defined by length(removed_epochs_method_1)
                x1=find(bad_epochs_all(:,1)==bad_epoch_2_unsorted(epoch));
                bad_epochs_all(x1,3)=epoch; %%% "weighting" of epoch
            end
            %%% same weighting idea as above
        else
            bad_epochs_all(:,3)=0;   %%% if method couldn't be performed (only 1 channel available, see above)
        end
        
        %%% method III
        if ~isempty(bad_epoch_3)
            for epoch=1:length(bad_epoch_3_unsorted)  %%% give each bad epoch a value = largest for worst epoch; largest value is defined by length(removed_epochs_method_1)
                x1=find(bad_epochs_all(:,1)==bad_epoch_3_unsorted(epoch));
                bad_epochs_all(x1,4)=epoch; %%% "weighting" of epoch
            end
            %%% same weighting idea as above
        else
            bad_epochs_all(:,4)=0;   %%% if method couldn't be performed (only 1 channel available or too few epochs, see above)
        end
        
        
        %%% 5.3.2) if selected, multiply "badness" of epoch by weighting factor
        %%%        for each method
        if ~isempty(weighting_factor(1))
            bad_epochs_all(:,2)=bad_epochs_all(:,2)*weighting_factor(1);
        end
        
        if ~isempty(weighting_factor(2))
            bad_epochs_all(:,3)=bad_epochs_all(:,3)*weighting_factor(2);
        end
        
        if ~isempty(weighting_factor(3))
            bad_epochs_all(:,4)=bad_epochs_all(:,4)*weighting_factor(3);
        end
        
        %%% 5.3.3) now sum up values from the methods for getting the final score
        %%%       (0 = good, largest is worst epoch)
        bad_epochs_all(:,5)=sum(bad_epochs_all(:,2:4),2);
        
        %%% 5.3.4) find final bad epochs and sort according to their final score
        %%%        (0 = good, largest is worst epoch)
        [i1,i2]=sort(bad_epochs_all(:,5));  %%% sort in ascending order
        %%% i1= weighting values, largest = worst
        %%% i2= epoch index
        %%% in total, number of epochs is larger then the defined %
        %%% criterion; thus, take now only number of epochs that should be
        %%% removed, take worst epochs found here; length(bad_epoch) and length(bad_epoch_2) are identical per definition
        bad_epoch_index_final=i2(length(i2)-length(bad_epoch)+1:end);
        bad_epoch_index_final=bad_epochs_all(bad_epoch_index_final,1);
        bad_epoch_weighting_final=i1(length(i1)-length(bad_epoch)+1:end);
        
        auto_epoch_cleaning.found_epochs_methods_joined=bad_epochs_all;
        auto_epoch_cleaning.found_epochs_unsorted_methods_joined_index=bad_epoch_index_final;
        auto_epoch_cleaning.found_epochs_unsorted_methods_joined_weighting=bad_epoch_weighting_final;
        
        window_samples_methods_joined(bad_epoch_index_final,3)=1;  %%% =1 remove this epoch
        window_samples_methods_joined(auto_epoch_cleaning.removed_epochs_sorted_previous_artifact_rejection_applied,3)=1;  %%% =1 remove this epoch
        window_samples_methods_joined(auto_epoch_cleaning.removed_epochs_sorted_flat_line_detection,3)=1;
        auto_epoch_cleaning.window_samples_vector_methods_joined_final=window_samples_methods_joined;
    end
    
else  %%% no valid epochs found at all
    auto_epoch_cleaning.removed_epochs_note='no valid epochs found';
end

%%% for convenience, indicate the final bad epochs, which shouldn't be used later
if isempty(auto_epoch_cleaning.window_samples_vector_final)
    auto_epoch_cleaning.window_samples_vector_final_final=auto_epoch_cleaning.window_samples_vector_methods_joined_final;
else
    auto_epoch_cleaning.window_samples_vector_final_final=auto_epoch_cleaning.window_samples_vector_final;
end


auto_epoch_cleaning.indices_all_epochs=1:size(auto_epoch_cleaning.epoch_mean,2);
auto_epoch_cleaning.indices_good_epochs=find(auto_epoch_cleaning.window_samples_vector_final_final(:,3)==0);
auto_epoch_cleaning.indices_bad_epochs=find(auto_epoch_cleaning.window_samples_vector_final_final(:,3)==1);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 6) indicate bad epochs in the data as NaN and plot

if ~isempty(auto_epoch_cleaning.window_samples_vector_final)  %%% if same epochs were found by all methods
    window_samples_vector_final=auto_epoch_cleaning.window_samples_vector_final;
end

if ~isempty(auto_epoch_cleaning.window_samples_vector_methods_joined_final)  %%% if different epochs were found, thus weight epochs and select % of the worst
    window_samples_vector_final=auto_epoch_cleaning.window_samples_vector_methods_joined_final;
end


%%% plot which epochs were removed
figure
subplot(2,2,1)
temp=auto_epoch_cleaning.epoch_mean;
temp(:,find(window_samples_vector_final(:,3)==1))=NaN;
temp(:,auto_epoch_cleaning.removed_epochs_sorted_previous_artifact_rejection_applied)=NaN;
temp(:,auto_epoch_cleaning.removed_epochs_sorted_flat_line_detection)=NaN;
hold on
plot(auto_epoch_cleaning.epoch_mean','r');
plot(temp','g');
%legend({'removed epochs', 'kept epochs'});
title({'Effect of joined methods removal',...
    ['weights: m1=' num2str(weighting_factor(1))...
    ', m2=' num2str(weighting_factor(2))...
    ', m3=' num2str(weighting_factor(3))]});
axis([0 size(auto_epoch_cleaning.epoch_mean,2) 0 inf]);
xlabel({['Epoch #; len: ' num2str(wind_ms) 'ms'],...
    ['removed epochs = ' num2str(length(auto_epoch_cleaning.indices_bad_epochs))...
    ', good epochs = ' num2str(length(auto_epoch_cleaning.indices_good_epochs))]});
ylabel('Epoch Mean')

subplot(2,2,2)
temp=mean(auto_epoch_cleaning.epoch_mean,1);
for chan=1:size(temp,1)
    temp(chan,find(auto_epoch_cleaning.window_samples_vector_method1(:,3)==1))=NaN;
end
hold on
plot(mean(auto_epoch_cleaning.epoch_mean,1),'r');
plot(temp','g');
legend({'removed epochs', 'kept epochs'});
title({'Effect of method I removal',...
    ['weight: m1=' num2str(weighting_factor(1))]});
axis([0 size(auto_epoch_cleaning.epoch_mean,2) 0 inf]);
xlabel({['Epoch #; len: ' num2str(wind_ms) 'ms']});
ylabel('Epoch Mean')

if ~isempty(bad_epoch_2)
    subplot(2,2,3)
    temp=auto_epoch_cleaning.epoch_mean_std_across_channels;
    for chan=1:size(temp,1)
        temp(chan,find(auto_epoch_cleaning.window_samples_vector_method2(:,3)==1))=NaN;
    end
    hold on
    plot(auto_epoch_cleaning.epoch_mean_std_across_channels,'r');
    plot(temp','g');
    legend({'removed epochs', 'kept epochs'});
    title({'Effect of method II removal',...
        ['weight: m2=' num2str(weighting_factor(2))]});
    axis([0 length(auto_epoch_cleaning.epoch_mean_std_across_channels) 0 inf]);
    xlabel({['Epoch #; len: ' num2str(wind_ms) 'ms']});
    ylabel('SD of single channel epoch mean');
end

if ~isempty(bad_epoch_3)
    subplot(2,2,4)
    temp=auto_epoch_cleaning.epoch_mean_Mahal_distance_acrossChannels;
    temp(find(auto_epoch_cleaning.window_samples_vector_method3(:,3)==1))=NaN;
    hold on
    plot(auto_epoch_cleaning.epoch_mean_Mahal_distance_acrossChannels,'r');
    plot(temp,'g');
    legend({'removed epochs', 'kept epochs'});
    title({'Effect of method III removal',...
        ['weight: m3=' num2str(weighting_factor(3))]});
    axis([0 length(auto_epoch_cleaning.epoch_mean_Mahal_distance_acrossChannels) 0 inf]);
    xlabel({['Epoch #; len: ' num2str(wind_ms) 'ms']});
    ylabel('Mahal. distance across channels');
end

sgtitle([subject, ' Effect of removal according to the used method'])
%suptitle([subject, ' Effect of removal according to the used method'])
saveCurrentFig(fullfile(datapath_save_figures, filesep),...
    sprintf('%s_autocleaning_removal_effect', subject), {'fig','png'},[1100 800])
%savefig([datapath_save_figures filesep 'autocleaning_removal_effect'])
%close
end

