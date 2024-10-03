function [ALLEEG, EEG_prepared, CURRENTSET] = prepareData( ALLEEG, EEG_to_prepare, CURRENTSET, subject, bemobil_config )
% Prepartion of data: enter chanlocs, remove unused channels, declare EOG, resample

output_filepath = [bemobil_config.study_folder bemobil_config.preprocessing_folder];

if ~isempty(bemobil_config.channel_locations_filename)
    channel_locations_filepath = [bemobil_config.study_folder bemobil_config.raw_data_folder...
        bemobil_config.filename_prefix num2str(subject) '/' bemobil_config.filename_prefix num2str(subject) '_'...
        bemobil_config.channel_locations_filename];
else
    channel_locations_filepath = [];
end

output_filename = [bemobil_config.filename_prefix num2str(subject) '_' bemobil_config.prepared_filename];

[ALLEEG, EEG_prepared, CURRENTSET] = bemobil_preprocess(ALLEEG, EEG_to_prepare, CURRENTSET,...
    channel_locations_filepath, bemobil_config.channels_to_remove, bemobil_config.eog_channels,...
    bemobil_config.resample_freq,...
    output_filename, output_filepath, bemobil_config.rename_channels);

disp('Data prepared for preprocessing !')
end

