function [EEG] = add_buffers_continousCleaning(EEG, buff_len)
%%% Adds 'buffers' at the beginning and ending of each identified 'bad' epoch
% This function may yield out-of-bound rejection segments (before 0 and after EEG.pnts)
% as well as overlapping regions 
% but this poses no problem to the rejection function eeg_eegrej (corrected).

% The output EEG.etc.(cleaning_field) structure is modified accordingly.
% The buffered rejection segments can be found in:
% EEG.etc.(cleaning_field).invalid_segments_final_start_stop_sample
% The unbuffered boundaries are moved to the field:
% EEG.etc.(cleaning_field).invalid_segments_beforeBuff_start_stop_sample

% Inputs:
%   EEG           - EEG struct to modify. Should contain an
%                   epoch cleaning field in the etc field.
%   buff_len      - Buffer length in number of samples
%                   buff_len will be removed on both sides of the bad epoch
% Outputs:
%   EEG           - Modified EEG struct

% Get the field corresponding to the epoch cleaning data
fieldNames = fieldnames(EEG.etc);
if sum(contains(fieldNames,'autoMoBI'))==0
    error('No cleaning boundaries were found, did you forget something?');
end
cleaning_field = fieldNames{contains(fieldNames,'autoMoBI')};

% Copy from EEG struct
artifact_corr = EEG.etc.(cleaning_field);
original_bounds = artifact_corr.invalid_segments_final_start_stop_sample;
% Put the original bounds in another field
artifact_corr.invalid_segments_beforeBuff_start_stop_sample = original_bounds;

% Do the stuff
buffered_bounds = [squeeze(original_bounds(:,1) - buff_len) ...
    squeeze(original_bounds(:,2) + buff_len)];

artifact_corr.invalid_segments_final_start_stop_sample = buffered_bounds;

% boolean indicating rejected samples in the original dataset
rejSamps = zeros(1,EEG.pnts);
for i = 1:size(buffered_bounds,1)
    if buffered_bounds(i,1) < 1
        inter_rej = 1:buffered_bounds(i,2);
        rejSamps(1:buffered_bounds(i,2)) = ones(1,length(inter_rej));
    elseif buffered_bounds(i,2) > EEG.pnts
        inter_rej = buffered_bounds(i,1):EEG.pnts;
        rejSamps(buffered_bounds(i,1):EEG.pnts) = ones(1,length(inter_rej));
    else
        inter_rej = buffered_bounds(i,1):buffered_bounds(i,2);
        rejSamps(buffered_bounds(i,1):buffered_bounds(i,2)) = ones(1,length(inter_rej));
    end
end
artifact_corr.rejectedSamples = rejSamps;

% Final copy to the EEG struct
EEG.etc.(cleaning_field) = artifact_corr;
end

