function out_of_interest_intervals = select_data_of_interest_EEGAFF2(EEG, DoI)
events = EEG.event;
TrialsInspection = EEG.etc.TrialsInspection;

intervals_of_interest = [];
i = 1;
for d = 1:numel(DoI)
    switch lower(DoI{d})
        case 'walkingbaselines'
            walkBase_tr = find(TrialsInspection.WalkingBaseline_perc == 100);
            for tr = 1:length(walkBase_tr)                 
                start_ev = findInStructWithEmpties(events, 'urevent', TrialsInspection.Trial_urevent_seq(walkBase_tr(tr),1));
                if isempty(start_ev)
                    % It may happen that the trial has been rejected anyway
                    continue
                end
                stop_ev = findInStructWithEmpties(events, 'urevent', TrialsInspection.Trial_urevent_seq(walkBase_tr(tr),end));
                if isempty(stop_ev)
                    % It may happen that the trial has been rejected anyway
                    continue
                end
                intervals_of_interest(i,:) = [ceil(events(start_ev).latency), floor(events(stop_ev).latency)];
                i = i+1;
            end
            
        case 'blackbaselines'
            blackBase_tr = find(TrialsInspection.BlackBaseline_perc == 100);
            for tr = 1:length(blackBase_tr)
                start_ev = findInStructWithEmpties(events, 'urevent', TrialsInspection.Trial_urevent_seq(walkBase_tr(tr),1));
                if isempty(start_ev)
                    % It may happen that the trial has been rejected anyway
                    continue
                end
                stop_ev = findInStructWithEmpties(events, 'urevent', TrialsInspection.Trial_urevent_seq(walkBase_tr(tr),2));
                if isempty(stop_ev)
                    % It may happen that the trial has been rejected anyway
                    continue
                end
                intervals_of_interest(i,:) = [ceil(events(start_ev).latency), floor(events(stop_ev).latency)];
                i = i+1;
            end
    end
end

out_of_interest_intervals = zeros(i, 2);
for j = 1:i
    if j == 1
        out_of_interest_intervals(j,1) = 1;
    else
        out_of_interest_intervals(j,1) = intervals_of_interest(j-1,2)+1;
    end
    
    if j == i
        out_of_interest_intervals(j,2) = EEG.pnts;
    else
        out_of_interest_intervals(j,2) = intervals_of_interest(j,1)-1;
    end
end
end