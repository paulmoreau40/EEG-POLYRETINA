function baselineReport4SC(EEG, cfg)
% Should be called after fillNaNs and events_check but before resampling or selection

events = EEG.event;
markers = events([events(:).duration] == 0 & ~(strcmp({events(:).type}, 'Blink') |...
    strcmp({events(:).type}, 'Saccade') | strcmp({events(:).type}, 'Fixation')));

starts = find(strcmp({markers(:).type}, 'TrialStart'));

Baselines = struct('Block', [], 'Condition', [], 'Trial', [], 'Duration_ms', [],...
    'EEGCompleteness',[], 'urevStart', [], 'urevEnd', []);

for s = 1:length(starts)
    blk = markers(starts(s)).Block;
    cnd = markers(starts(s)).Condition;
    trl = markers(starts(s)).Trial;
    
    markersInCnd = strcmp({markers(:).Condition}, cnd);
    if trl == 1
        startBaseline = find(strcmp({markers(:).type}, 'ConditionStart') & markersInCnd);
    else
        markersInPrevTrl= [markers(:).Trial] == trl-1;
        startBaseline = find(strcmp({markers(:).type}, 'TrialEnd') & markersInCnd & markersInPrevTrl);
    end
    
    % Find duration
    dur = 1 + markers(starts(s)).latency - markers(startBaseline).latency;
    
    % Find completeness
    cmp = 100*(1 - (sum(sum(isnan(EEG.data(:,markers(startBaseline).latency:markers(starts(s)).latency)),1)>0)/dur));
    
    Baselines(s) = struct('Block', blk, 'Condition', cnd, 'Trial', trl,...
        'Duration_ms', dur, 'EEGCompleteness', cmp,...
        'urevStart', markers(startBaseline).urevent, 'urevEnd', markers(starts(s)).urevent);
end

%% Save the output structure
subject = cfg.subjects(cfg.current_subject).id;
dirname = [cfg.study_folder cfg.preprocessing_folder];
fname = [subject '_BaselinesReport.mat'];
save(fullfile(dirname, fname), 'Baselines')

end