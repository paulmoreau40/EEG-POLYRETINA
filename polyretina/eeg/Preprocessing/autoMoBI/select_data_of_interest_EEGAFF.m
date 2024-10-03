function out_of_interest_intervals = select_data_of_interest_EEGAFF(EEG, DoI)
events = EEG.event;

intervals_of_interest = [];
i = 0;
for d = 1:numel(DoI)
    switch lower(DoI{d})
        case 'walkingbaselines'
            walkBaseEvents = find(contains({events(:).type}, 'WalkingBase'));
            for e = 1:length(walkBaseEvents)
                if strcmp(events(walkBaseEvents(e)).type, 'WalkingBaseStart')
                    i = i+1;
                    intervals_of_interest(i,1) = ceil(events(walkBaseEvents(e)).latency);
                else
                    intervals_of_interest(i,2) = floor(events(walkBaseEvents(e)).latency);
                end
            end
        case 'blackbaselines'
            baseEvents = find(contains({events(:).type}, 'Base'));
            for e = 1:length(baseEvents)
                if strcmp(events(baseEvents(e)).type, 'BaseStart')
                    i = i+1;
                    intervals_of_interest(i,1) = ceil(events(baseEvents(e)).latency);
                elseif strcmp(events(baseEvents(e)).type, 'BaseEnd')
                    intervals_of_interest(i,2) = floor(events(baseEvents(e)).latency);
                end
            end
        case 'eyecalibration'
            % First category, eye calibration
            eyeCalibEvents = find(contains({events(:).type}, 'EyeCalib'));
            if ~isempty(eyeCalibEvents)
                
                if contains(EEG.filename, 'ANG06')
                    eyeCalibEvents = eyeCalibEvents(3:end);
                end
                
                for e = 1:length(eyeCalibEvents)
                    if strcmp(events(eyeCalibEvents(e)).type, 'StartEyeCalibration')
                        i = i+1;
                        intervals_of_interest(i,1) = ceil(events(eyeCalibEvents(e)).latency);
                    else
                        intervals_of_interest(i,2) = floor(events(eyeCalibEvents(e)).latency);
                    end
                end
            end
            
            % Second category, eye validation (single series of focus events)
            eyeValidEvents = find(contains({events(:).Location}, 'EyeCalib'));
            if ~isempty(eyeValidEvents)
                if length(eyeCalibEvents)==2
                    i = i+1;
                    intervals_of_interest(i,1) = ceil(events(eyeValidEvents(1)).latency);
                    intervals_of_interest(i,2) = floor(events(eyeValidEvents(end)).latency);
                else
                    error('Did not plan for more than 1 calibration');
                end
            end
    end
end

out_of_interest_intervals = zeros(i+1, 2);
for j = 1:i+1
    if j == 1
        out_of_interest_intervals(j,1) = 1;
    else
        out_of_interest_intervals(j,1) = intervals_of_interest(j-1,2)+1;
    end
    
    if j == i+1
        out_of_interest_intervals(j,2) = EEG.pnts;
    else
        out_of_interest_intervals(j,2) = intervals_of_interest(j,1)-1;
    end
end
end