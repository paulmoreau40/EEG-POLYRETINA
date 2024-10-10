function [main, raw, elec, raw_EEGLAB, preproc, SS_ana, MS_ana, Figs, MetaFile, MetaTab, SR_path] =...
    getMainFoldersNames(user)
% Define here the folders access so there is only one place to change
% Input:
%   user        - permits to specify the user for the main path

switch lower(user)
    case 'default' % should work automatically
        base_path = fileparts(fileparts(pwd));
        main = [fullfile(base_path, 'data', 'analysis'), filesep];
        Figs = [fullfile(base_path, 'figures'), filesep];
        MetaFile = fullfile(main, 'Polyretina_meta.xlsx');
        MetaTab = 'SubjectsInfo';
        SR_path = '';

    % if default paths don't work, custom your personal paths :
    case 'alex'
        main = 'D:\Data_EEG-AFF\v2\analysis\';
        Figs = 'D:\Data_EEG-AFF\v2\figures\';
        MetaFile = fullfile('D:\Data_EEG-AFF\v2', 'TrackingEEGanalysis_EEGAFF.xlsx');
        MetaTab = 'SubjectsInfo';
        SR_path = '';
        
    case 'ilaria'
        main = '/Users/ilaria/Documents/MATLAB/data/EEG-AFF/analysis/';
        Figs = '/Users/ilaria/Documents/MATLAB/data/EEG-AFF/figures/';
        MetaFile = fullfile('/Users/ilaria/Documents/MATLAB/data/EEG-AFF/analysis', 'TrackingEEGanalysis_EEGAFF2.xlsx');
        MetaTab = 'SubjectsInfo';
        SR_path = '';

    case 'paulcoarse'
        main = 'F:\EEG-POL\data\analysis\';
        Figs = 'F:\EEG-POL\figures\';
        MetaFile = fullfile('F:\EEG-POL\data\analysis\', 'Polyretina_meta.xlsx');
        MetaTab = 'SubjectsInfo';
        SR_path = '';
    otherwise
        error('Unknown user');
end
raw = ['0_raw-data' filesep];
elec = ['0_electrodes' filesep];
raw_EEGLAB = ['1_raw-eeglab' filesep];
preproc = ['2_preprocessing' filesep];
SS_ana = ['3_single-subject-analysis' filesep];
MS_ana = ['4_multi-subject-analysis' filesep];
end

