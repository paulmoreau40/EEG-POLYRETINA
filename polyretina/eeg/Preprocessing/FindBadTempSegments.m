function [EEG, invalid_segments_index] = FindBadTempSegments(ALLEEG, EEG, CURRENTSET, subject, bemobil_config)
% Automatic Detection of bad temporal segments
% The EEG data should be HP filtered
% The data is not saved at this step

output_path = [bemobil_config.study_folder...
    bemobil_config.preprocessing_folder bemobil_config.preprocess_pipeline_folder];

output_filepath = [output_path bemobil_config.filename_prefix num2str(subject)];

switch bemobil_config.globalArchitecture
    case 'Gramann'
        switch bemobil_config.badSegmentsRemovalMethod
            case 'simple'
                %% ----------------SIMPLE PIPELINE--------------------------
                [invalid_segments_index] = select_data_of_interest(EEG, 'fulldata');
                EEG.etc.simple_cleaning.invalid_segments_index = invalid_segments_index;
            case 'automatic'
                %% ----------------BERLIN PIPELINE--------------------------
                % The EEG set is saved at different steps in the process
                %% Perform ICA on explo and long baselines only to detect eye components
                dataOfInterest = bemobil_config.dataOfInterest4AMICAeye;
                OoI_segments_index = select_data_of_interest(EEG, dataOfInterest);
                EEG_cleaned_for_ICAeye = eeg_eegrej(EEG, OoI_segments_index);
                
                % save the dataset (to have the data locally, better for AMICA)
                EEG_cleaned_for_ICAeye = pop_saveset(EEG_cleaned_for_ICAeye,...
                    'filename', [bemobil_config.filename_prefix num2str(subject) '_'...
                    bemobil_config.without_eyes], 'filepath', output_filepath);
                
                [ALLEEG, EEG_AMICAeye, CURRENTSET] =...
                    FindEyesWithICA(ALLEEG, EEG_cleaned_for_ICAeye, CURRENTSET, subject, bemobil_config);
                clear EEG_cleaned_for_ICAeye
                
                % copy the dataset and save it (overwrite the no_eyes file)
                [ALLEEG, EEG_prep4auto_cleaning, CURRENTSET] =...
                    bemobilCustom_copy_spatial_filter(EEG, ALLEEG, CURRENTSET, EEG_AMICAeye,...
                    [bemobil_config.filename_prefix num2str(subject) '_' bemobil_config.without_eyes], output_filepath);
                
                %% Define eyes and remove
                class_results= EEG_prep4auto_cleaning.etc.ic_classification.ICLabel.classifications;
                
                %% Visual inspection
                %{
                eye_comps = DefineEyeComps(class_results, bemobil_config.ICdetection_thresholds, 'all', [], []);
                freq_range = [1 100]; % in Hz
                spec_opt = {'freqrange', freq_range}; % cell array of options which are passed to spectopo()
                erp_opt = {}; % cell array of options which are passed to erpimage()
                pop_viewprops(EEG_prep4auto_cleaning, 0, eye_comps', spec_opt , erp_opt, 1, 'ICLabel');
                %}
                
                %% Only the unique ones                
                eye_comps = DefineEyeComps(class_results, bemobil_config.ICdetection_thresholds, 'unique',...
                    [bemobil_config.figures_folder, 'EyeICsAfterAMICA1/'],...
                    ['S',num2str(subject),'_Eye_ICs_distribution_AMICA1_' cell2str(dataOfInterest, '_')]);
                
                %% Subtract components:
                EEG_prep4auto_cleaning = pop_subcomp(EEG_prep4auto_cleaning, eye_comps);
                %{
                % project out eyes and save
                EEG_AMICA_no_eyes = pop_subcomp(EEG_AMICA_raw, eyes);
                % save the dataset without eye components, overwritting preceeding one
                pop_saveset(EEG_AMICA_no_eyes, 'filename', [bemobil_config.filename_prefix num2str(subject) '_'...
                    bemobil_config.without_eyes], 'filepath', output_filepath);
                %}
                
                %% determine automatic time domain cleaning boundaries on the channel level
                % Remove events to save RAM
                EEG_prep4auto_cleaning = pop_editeventvals(EEG_prep4auto_cleaning,'delete', 1:length(EEG_prep4auto_cleaning.event));
                
                automatic_cleaning_settings = AutomaticCleaningSettings(bemobil_config, dataOfInterest);
                
                disp('Determining continuous data cleaning boundaries...')
                [~, automatic_cleaning] = autoClean_continuousEEG_custom(...
                    EEG_prep4auto_cleaning, automatic_cleaning_settings, subject);
                
                EEG.etc.spatial_filter_eyerejection = EEG_prep4auto_cleaning.etc.spatial_filter;
                EEG.etc.ic_classification_eyerejection = EEG_prep4auto_cleaning.etc.ic_classification;
                EEG.etc.ic_classification_eyerejection.removed_comps = eye_comps;
                EEG.etc.automatic_cleaning = automatic_cleaning;
                clear EEG_prep4auto_cleaning
                
                % LUKAS adding: Add buffers to bad epochs
                disp('Adding buffers to data cleaning boundaries...')
                buffer_length_sample = (bemobil_config.buffer_length*automatic_cleaning_settings.wind_ms)...
                    *(EEG.srate/1000);
                EEG = add_buffers_continousCleaning(EEG, buffer_length_sample);            

                % get cleaning indices
                invalid_segments_index = EEG.etc.automatic_cleaning.invalid_segments_final_start_stop_sample;
        end
    case 'APP'
        %% -------------------APP PIPELINE--------------------------
        if bemobil_config.crop_OoI_data
            % Extract interesting data
            OoI_segments_index = select_data_of_interest(EEG, 'fulldata');
        else
            OoI_segments_index = [];
        end
        % Remove events (simpler for later)
        EEG = pop_editeventvals(EEG,'delete', 1:length(EEG.event));
        
        % Epoch the EEG
        disp('APP: Epoching data for Bad segments removal...')
        duration = EEG.xmax; %duration in s
        Nmax_wdws = floor(duration/1); % 1 indicates that we aim for 1s epochs
        epoch_wdw = duration/Nmax_wdws; % in s1
        
        EEG_epoched = eeg_regepochs(EEG, 'recurrence', epoch_wdw, 'limits', [0 epoch_wdw],...
            'rmbase', NaN, 'extractepochs', 'on');
        
        % Some additional events have appeared...
        events2del = find([EEG_epoched.event.latency]~=round([EEG_epoched.event.latency]));
        EEG_epoched = pop_editeventvals(EEG_epoched,'delete', events2del);
        
        % Remove out-of-interest epochs
        if ~isempty(OoI_segments_index)
            disp('APP: Cropping out-of-interest data...')
            epochs2rej= zeros(EEG_epoched.trials,1);
            for inter = 1:size(OoI_segments_index,1)
                if OoI_segments_index(inter,1)==0
                    startInd = 1;
                else
                    startInd = ceil(OoI_segments_index(inter,1)/EEG_epoched.pnts)+1;
                end
                
                if OoI_segments_index(inter,2)==EEG.pnts
                    endInd = EEG_epoched.trials;
                else
                    endInd = floor(OoI_segments_index(inter,2)/EEG_epoched.pnts);
                end
                
                epochs2rej(startInd:endInd)=1;
            end
            
            EEG_epoched = pop_rejepoch(EEG_epoched, epochs2rej, 0);
        end
        
        % Remove non boundary events
        %bound_events = find(strcmp({EEG_ready4Epoching.event.type}, 'boundary'));
        %non_bound_events = setdiff(1:length(EEG_ready4Epoching.event), bound_events);
        %EEG_ready4Epoching = pop_editeventvals(EEG_ready4Epoching,'delete', non_bound_events);
        
        % Max amplitude difference
        Max_amps = squeeze(max(EEG_epoched.data, [], 2));
        Min_amps = squeeze(min(EEG_epoched.data, [], 2));
        Amp_diffs = Max_amps - Min_amps;
        
        [Epoch_mean_amp_diff,~] = Biweight_custom(Amp_diffs', EEG_epoched.nbchan, bemobil_config.censorBiweight);
        disp('APP: Detecting bad epochs by Mean amplitude difference criterion')
        RejEpochs4AmpDiff = FindOutliers_custom(Epoch_mean_amp_diff,...
            bemobil_config.z_criterion, bemobil_config.inner_fence, 'right');
        plot_distribution(Epoch_mean_amp_diff, RejEpochs4AmpDiff, 'Epoch', 'APP_right')
        xlabel('Mean amplitude difference [??V]')
        title({['S' num2str(subject) '- Distribution of the mean amplitude difference'],...
            'for the epochs inspected by APP'})
        saveCurrentFig([bemobil_config.figures_folder 'APP_distributions/'],...
            ['S' num2str(subject) '_epochs_mean_amp_diff'], {'png'}, [600 500]);
        
        % Epoch variance or the mean GFP
        Epoch_GFP = mean(squeeze(std(EEG_epoched.data,0,1)));
        disp('APP: Detecting bad epochs by GFP criterion')
        RejEpochs4GFP = FindOutliers_custom(Epoch_GFP, bemobil_config.z_criterion,...
            bemobil_config.inner_fence, 'right');
        plot_distribution(Epoch_GFP, RejEpochs4GFP, 'Epoch', 'APP_right')
        xlabel('Global field potential [??V]')
        title({['S' num2str(subject) '- Distribution of the global field potential'],...
            'for the epochs inspected by APP'})
        saveCurrentFig([bemobil_config.figures_folder 'APP_distributions/'],...
            ['S' num2str(subject) '_epochs_GFP'], {'png'}, [600 500]);
        
        % Epoch's mean deviation from channel means.
        [Chans_mean,~] = Biweight_custom(EEG_epoched.data(:,:),...
            EEG_epoched.trials*EEG_epoched.pnts, bemobil_config.censorBiweight); % channel mean for all epochs
        Epoch_mean_dev = mean(abs(squeeze(mean(EEG_epoched.data,2)) - Chans_mean));
        disp('APP: Detecting bad epochs by deviation criterion')
        RejEpochs4MeanDev = FindOutliers_custom(Epoch_mean_dev,...
            bemobil_config.z_criterion, bemobil_config.inner_fence, 'right');
        plot_distribution(Epoch_mean_dev, RejEpochs4MeanDev, 'Epoch', 'APP_right')
        xlabel('Mean deviation from channel mean  [??V]')
        title({['S' num2str(subject) '- Distribution of the mean deviation'],...
            'for the epochs inspected by APP'})
        saveCurrentFig([bemobil_config.figures_folder 'APP_distributions/'],...
            ['S' num2str(subject) '_epochs_mean_dev'], {'png'}, [600 500]);
        
        RejEpochs = union(union(RejEpochs4AmpDiff, RejEpochs4GFP), RejEpochs4MeanDev);
        
        % Define the segment indices to reject:
        rejected_segments_index = zeros(length(RejEpochs), 2);
        for i = 1:length(RejEpochs)
            epoch_local_ind = find([EEG_epoched.event.epoch]==RejEpochs(i));
            epoch_global_ind = EEG_epoched.event(epoch_local_ind).urevent;
            rejected_segments_index(i,1) = EEG_epoched.urevent(epoch_global_ind).latency;
            rejected_segments_index(i,2) = rejected_segments_index(i,1) + EEG_epoched.pnts;
        end
        
        % Combining the rejected segments and the out-of-interest ones
        invalid_segments_index = [OoI_segments_index;rejected_segments_index];
        [~, sorting_vect] = sort(invalid_segments_index(:,1));
        invalid_segments_index = invalid_segments_index(sorting_vect,:);
        
        epoch_cleaning = struct('epoch_wdw', epoch_wdw,...
            'censorBiweigth', bemobil_config.censorBiweight,...
            'z_criterion', bemobil_config.z_criterion,...
            'inner_fence', bemobil_config.inner_fence,...
            'Rejected_epochs',RejEpochs,...
            'epochs_ampDiff', Epoch_mean_amp_diff,...
            'Rejected_epochs_ampDiff',RejEpochs4AmpDiff,...
            'epochs_GFP', Epoch_GFP,...
            'Rejected_epochs_GFP',RejEpochs4GFP,...
            'epochs_MeanDev', Epoch_mean_dev,...
            'Rejected_epochs_MeanDev',RejEpochs4MeanDev,...
            'Note', "Apart from the following, indices are not valid in the complete dataset (out-of-interest data was cropped)",...
            'cropped_segments_index',OoI_segments_index,...
            'rejected_byAPP_segment_index', rejected_segments_index,...
            'invalid_segments_index',invalid_segments_index);
        
        EEG.etc.APP_epoch_cleaning = epoch_cleaning;
end
end