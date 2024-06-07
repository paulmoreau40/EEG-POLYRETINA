% bemobil_dipfit() - Prepares data for dipole fitting and runs the dipole fitting procedure
%
% Inputs:
%   EEG                     - current EEGLAB EEG structure
%   RV_threshold            - number percentage of residual variance accepted, default is '15'
%   remove_outside_head     - 'on' or 'off' to remove dipoles located outside the head
%   number_of_dipoles       - '1' or '2', 2 meaning bilateral dipole fitting
%
% Outputs:
%   EEG                     - current EEGLAB EEG structure
%
%   .set data file of current EEGLAB EEG structure stored on disk (OPTIONALLY)
%
% See also:
%   EEGLAB, coregister, pop_dipfit_settings, pop_multifit

function [EEG] = dipfit_EEGAFF(EEG, RV_threshold, remove_outside_head, number_of_dipoles)
% load some standard data for dipfit
dipfitdefs;
dipfit_path = fileparts(which('pop_multifit'));

% use standard BEM headmodel
% change relevant fiducial labels so that it matches standard 5/10 template
chaninfo = EEG.chaninfo;
% for f = 1:numel(chaninfo.nodatchans)
%     switch chaninfo.nodatchans(f).labels
%         case 'nas'
%             chaninfo.nodatchans(f).labels = 'Nz';
%         case 'lhj'
%             chaninfo.nodatchans(f).labels = 'LPA';
%         case 'rhj'
%             chaninfo.nodatchans(f).labels = 'RPA';
%         otherwise
%             error('Unknown fiducial name')
%     end
% end

disp('Coregistering electrodes to 10-5 template...')
% [~, transform] = coregister(EEG.chanlocs,...
%     fullfile(dipfit_path, 'standard_BEM', 'elec','standard_1005.elc'),'chaninfo1', chaninfo,...
%         'manual', 'off', 'alignfid', {'Nz', 'LPA','RPA'});
[~, transform] = coregister(EEG.chanlocs,...
   fullfile(dipfit_path, 'standard_BEM', 'elec','standard_1005.elc'), 'chaninfo1', chaninfo);

% Exclude EOG channels
% chans2select = [EEG.chanlocs.urchan]'; % POL CHANGED
chans2select = (1:1:128)';
EEG_indices = strcmp({EEG.chanlocs.type}, 'EEG');
chans2select = chans2select(EEG_indices);

% settings
EEG = pop_dipfit_settings(EEG, 'hdmfile',fullfile(dipfit_path, 'standard_BEM', 'standard_vol.mat'),...
    'coordformat','MNI',...
    'mrifile',fullfile(dipfit_path, 'standard_BEM', 'standard_mri.mat'),...
    'chanfile',fullfile(dipfit_path, 'standard_BEM', 'elec', 'standard_1005.elc'),...
    'coord_transform', transform, 'chansel', chans2select);

% do the dipole fitting
EEG = pop_multifit(EEG, [1:size(EEG.icaweights,1)], 'threshold', RV_threshold,...
    'dipoles', number_of_dipoles, 'rmout', remove_outside_head);

% save dipfit info in EEG.etc
EEG.etc.dipfit.transform = transform;
EEG.etc.dipfit.RV_threshold = RV_threshold;
EEG.etc.dipfit.remove_outside_head = remove_outside_head;
EEG.etc.dipfit.number_of_dipoles = number_of_dipoles;

EEG = eeg_checkset(EEG);