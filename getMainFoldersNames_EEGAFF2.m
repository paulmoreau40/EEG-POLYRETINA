function [main, raw, elec, raw_EEGLAB, preproc, SS_ana, MS_ana, Figs, MetaFile, MetaTab, SR_path] =...
    getMainFoldersNames_EEGAFF2(user)
% Define here the folders access so there is only one place to change
% Input:
%   user        - permits to specify the user for the main path (Alex or JB)

switch lower(user)
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

    case 'paul'
        main = 'F:\PM_Polyretina\Data\analysis\';
        Figs = 'F:\PM_Polyretina\Data\figures\';
        MetaFile = fullfile('F:\PM_Polyretina\Data\analysis\', 'Polyretina_meta.xlsx');
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

