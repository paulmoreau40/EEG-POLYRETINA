config4SC_Alex;

for subject_ind = subject_inds
    if ~exist('ALLEEG','var')
        launchEEGLAB;
    end
    
    % clear RAM
    STUDY = []; CURRENTSTUDY = 0; ALLEEG = []; EEG=[]; CURRENTSET=[];
    
    % Overwrite subject for testing
    %subject = 1;
    subject = study_config.subject_names{subject_ind};
    disp(['Subject ' subject]);
    study_config.current_subject = subject_ind;
    
    input_filepath = [study_config.study_folder study_config.raw_data_folder subject];
    list = dir(input_filepath);
    input_file = list(contains({list.name},'.cnt')).name;
    EEG = pop_loadeep_v4_custom(fullfile(input_filepath, input_file));
end