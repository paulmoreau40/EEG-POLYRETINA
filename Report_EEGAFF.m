function EEG = Report_EEGAFF(EEG, cfg)
% Should be called after fillNaNs and events_check but before resampling or selection
% if all the trials are presents

events = EEG.event;
evts_noBounds = events(~strcmp({events.type},'boundary'));
% startTrialEvents = contains({evts_noBounds.type},'BaseStart');
% obsStartEvents = strcmp({evts_noBounds.type},'ObsStart');
% questionStartEvents = strcmp({evts_noBounds.type},'QuestionStart');
% exploStartEvents = strcmp({evts_noBounds.type},'ExploStart');
% endTrialEvents = strcmp({evts_noBounds.SoundName},'Close');
distMeasureEvents = strcmp({evts_noBounds.type},'DistMeasure');

Answers = ReadAnswers_EEGAFF(cfg);

DataBase = struct('ID', [], 'Gender', [], 'Age', [], 'Group', [], 'Block', [], 'Trial', [], 'Condition', [],...
    'Distance', [], 'Response', [], 'urevent_seq', [], 'EEGComplete_trial',[], 'EEGComplete_details',[]);

% Subject specific variables
subject = cfg.subjects(cfg.current_subject).id;
gender = cfg.subjects(cfg.current_subject).gender;
age = cfg.subjects(cfg.current_subject).age;
group = cfg.subjects(cfg.current_subject).group;

TrialInspection = EEG.etc.TrialsInspection;

i=1;
for tr = 1:size(TrialInspection,1)    
    cnd = TrialInspection.TrialType(tr);
    
    if strcmp(cnd,'Baseline')
        continue
    else
        bl = TrialInspection.BlockInd(tr);
        tr_inBl = TrialInspection.TrialInd(tr);
        trialEvents = [evts_noBounds.BlockIndex] == bl & [evts_noBounds.TrialIndex] == tr_inBl;        
        
%         % Get participant's answer
%         answer = Answers.(sprintf('Block%d',bl-1))(Answers.Essai==tr_inBl);
        
%         % Create sequence of events
%         seq = [evts_noBounds(trialEvents & startTrialEvents).urevent,...
%             evts_noBounds(trialEvents & obsStartEvents).urevent,...
%             evts_noBounds(trialEvents & questionStartEvents).urevent,...
%             evts_noBounds(trialEvents & exploStartEvents).urevent,...
%             evts_noBounds(trialEvents & endTrialEvents).urevent];
        
        % Concatenate EEG completness data
        eeg_perc = [TrialInspection.BlackBaseline_perc(tr),...
            TrialInspection.Observation_perc(tr),...
            TrialInspection.Question_perc(tr),...
            TrialInspection.Exploration_perc(tr)];
        
        %fprintf('Block %d, Trial %d.\n', bl, tr_inBl)
        DataBase(i) = struct('ID', subject, 'Gender', gender, 'Age', age, 'Group', group,...
            'Block', bl, 'Trial', tr_inBl, 'Condition', cnd,...
            'Distance', evts_noBounds(trialEvents & distMeasureEvents).Distance,...
            'Response', Answers.(sprintf('Block%d',bl-1))(Answers.Essai==tr_inBl),...
            'urevent_seq', TrialInspection.Trial_urevent_seq(tr,:),...
            'EEGComplete_trial', TrialInspection.Trial_perc(tr),...
            'EEGComplete_details', eeg_perc);
        i=i+1;
    end
end

DataBase = struct2table(DataBase);
EEG.etc.TrialData = DataBase;

% names = fields(Answers);
% id = names(1);
% id = char(id);
% num = id(end-1:end);
%      
% Answers(1).(sprintf(id))=0;
% events = EEG.event;
% markers = events;
% %missing_base_start = false;
% 
% block1 = markers([markers.BlockIndex]==1);
% n_trials_b1 = length(unique([block1.TrialIndex]));
% if n_trials_b1<26
%     warning(strcat('Block 1 miss some trial'));
% end
% block2 = markers([markers.BlockIndex]==2);
% n_trials_b2 = length(unique([block2.TrialIndex]));
% if n_trials_b2<26
%     warning(strcat('Block 2 miss some trial'));
% end
% block3 = markers([markers.BlockIndex]==3);
% n_trials_b3 = length(unique([block3.TrialIndex]));
% if n_trials_b3<26
%     warning(strcat('Block 3 miss some trial'));
% end
% 
% WalkingBaseEnd = strcmp({markers.type}, 'WalkingBaseEnd');
% WalkBaseStart = strcmp({markers.type}, 'WalkingBaseStart');
% 
% % add the WakingBaseEnd event
% % ZoneChange = strcmp({markers.type}, 'ZoneChange');
% % Ending = strcmp({markers.Phase}, 'ending');
% % ZoneChangeEnd = and(ZoneChange, Ending);
% % Baseline = strcmp({markers.TrialType}, 'Baseline');
% % WalkingBaseEnd = and(ZoneChangeEnd, Baseline);
% % WalkingBaseEnd_ur = [markers(WalkingBaseEnd).urevent];
% if sum(WalkingBaseEnd)<6
%     warning(strcat('One Baseline Trial is missing from event end.'));
% %     PlayingSound = strcmp({markers.type}, 'PlayingSound');
% %     Ending = strcmp({markers.SoundName}, 'Close');
% %     ZoneChangeEnd = and(PlayingSound, Ending);
% %     Baseline = strcmp({markers.TrialType}, 'Baseline');
% %     WalkingBaseEnd = and(ZoneChangeEnd, Baseline);
% end
% % add to Impo_Markers the WalkingBaseEnd
% % for ur=WalkingBaseEnd_ur
% %     markers(ur).type = 'WalkingBaseEnd';
% % end
% 
% 
% % if baseline start is missing 
% if sum(WalkBaseStart)<6
%     warning(strcat('One Baseline Trial is missing from event start.'));
% %     missing_base_start = true;
% %     markers = events;
% %     Open = strcmp({markers.SoundName}, 'Open');
% %     Baseline = strcmp({markers.TrialType}, 'Baseline');
% %     WalkBaseStart = and(Open,Baseline);
% end
% 
% % add WalkingBaseStart   
% % if missing_base_start
% %     WalkingBaseStart_ur = [events(WalkBaseStart).urevent];
% %     for ur=WalkingBaseStart_ur
% %         markers(ur).type = 'WalkingBaseStart';        
% %     end
% % end
% 
% 
% markers = markers(~(strcmp({markers(:).type}, 'PlayingSound') |...
%     strcmp({markers(:).type}, 'ZoneChange') | strcmp({markers(:).type}, 'SpatialDiscrExit') |...
%     strcmp({markers(:).type}, 'SpatialDiscrEnter') | strcmp({markers(:).type}, 'boundary') | strcmp({markers(:).type}, 'DistMeasure')));
% 
% stops_Baselines = find(strcmp({markers(:).type}, 'BaseEnd'));
% stops_WalkBaselines = find(strcmp({markers(:).type}, 'WalkingBaseEnd'));
% stops_obs = find(strcmp({markers(:).type}, 'ObsEnd'));
% stops_ques = find(strcmp({markers(:).type}, 'QuestionEnd'));
% stops_explo = find(strcmp({markers(:).type}, 'ExploEnd'));
% stops_Trials = find(strcmp({markers(:).type}, 'TrialEnd'));
% end_trials = markers(sort([stops_WalkBaselines,stops_Trials]));
% 
% DataBase = struct('ID', [], 'Gender', [], 'Age', [], 'Group', [], 'urevent', [], 'Block', [], 'Trial', [], 'Condition', [],...
%     'Distances', [], 'Responses', [], 'ErrUnderEst', [], 'Err_absolue', [], 'Mean_err_BlockCond', [], 'st_dev_BlockCond', [],...
%     'Duration_ms', [], 'EEGCompleteness',[]);
% 
% % inizialize variables
% seq=0;
% dur=0;
% cmp=0;
% mean_err=0;
% st_dev=0;
% 
% % ispection des Stops (walking baseline End and trial End) for all trials (without baselines) 
% i=1;
% for s = [end_trials.urevent]     
%      trial = end_trials(find([end_trials.urevent]==s)).TrialIndex;
%      cond = end_trials(find([end_trials.urevent]==s)).TrialType;
%      block = end_trials(find([end_trials.urevent]==s)).BlockIndex; 
%      dist = end_trials(find([end_trials.urevent]==s)).Distance; 
%      urev = end_trials(find([end_trials.urevent]==s)).urevent;
%      
%      if strcmp(cond, 'Baseline')
%          answer = 0;
%      else
%          if (block == 1)
%             answer = Answers(find([Answers.(sprintf(id))]==trial)).Var2;
%         elseif (block == 2)
%             answer = Answers(find([Answers.(sprintf(id))]==trial)).Var4;
%         elseif (block == 3)
%             answer = Answers(find([Answers.(sprintf(id))]==trial)).Var6;
%          end   
%      end
%     
%      err_abs = abs(dist-answer);
%      err = dist-answer;
%      
%      DataBase(i) = struct('ID', id, 'Gender',cfg.subjects(str2double(num)).gender, 'Age', cfg.subjects(str2double(num)).age, 'Group', cfg.subjects(str2double(num)).group,'urevent', urev,...
%          'Block', block, 'Trial', trial, 'Condition', cond,...
%          'Distances', dist, 'Responses', answer,'ErrUnderEst', err, 'Err_absolue', err_abs,'Mean_err_BlockCond', mean_err, 'st_dev_BlockCond', st_dev,...
%          'Duration_ms', dur,...
%          'EEGCompleteness', cmp);
%      
%      i=i+1;
%      
% end     
% 
% %% Add to the Database the mean response of each block, the std dev of each block and condition, 
% % the error of each response like (abs(distance-response))
% 
% % separate trials into blocks 
% ToRemove = strcmp({DataBase.Condition}, 'Baseline');
% DataBase = DataBase(~ToRemove);
% Block1 = DataBase([DataBase.Block]==1);
% Block2 = DataBase([DataBase.Block]==2);
% Block3 = DataBase([DataBase.Block]==3);
% 
% for c = 1:3 % blocks
%     mean_blocks = [];
%     if c==1 %short
%         b1 = Block1(strcmp({Block1.Condition}, 'Short'));
%         b2 = Block2(strcmp({Block2.Condition}, 'Short'));
%         b3 = Block3(strcmp({Block3.Condition}, 'Short'));
%         Data_ur = [DataBase(strcmp({DataBase.Condition}, 'Short')).urevent];
%     elseif c==2 %medium
%         b1 = Block1(strcmp({Block1.Condition}, 'Medium'));
%         b2 = Block2(strcmp({Block2.Condition}, 'Medium'));
%         b3 = Block3(strcmp({Block3.Condition}, 'Medium'));
%         Data_ur = [DataBase(strcmp({DataBase.Condition}, 'Medium')).urevent];
%     elseif c==3 %long
%         b1 = Block1(strcmp({Block1.Condition}, 'Long'));
%         b2 = Block2(strcmp({Block2.Condition}, 'Long'));
%         b3 = Block3(strcmp({Block3.Condition}, 'Long'));
%         Data_ur = [DataBase(strcmp({DataBase.Condition}, 'Long')).urevent];
%     end
% 
%     ans_1 = [b1.Responses]'; % resp B1
%     ans_2 = [b2.Responses]'; % resp B2
%     ans_3 = [b3.Responses]'; % resp B3
%     
%     meanErr = [mean([b1.Err_absolue]'),mean([b2.Err_absolue]'),mean([b3.Err_absolue]')];
%     
%     % create column with answers
%     answers = [ans_1;ans_2;ans_3];
%     mean_blocks = [mean(ans_1),mean(ans_2),mean(ans_3)];
%     st_dev = [std(ans_1),std(ans_2),std(ans_2)];    
%     
%     % assign std deviation (one par block, par condition)
%     for ur=Data_ur
%         [DataBase(([DataBase.Block]==1) & ([DataBase.urevent]==ur)).st_dev_BlockCond] = st_dev(1);
%         [DataBase(([DataBase.Block]==2) & ([DataBase.urevent]==ur)).st_dev_BlockCond] = st_dev(2);
%         [DataBase(([DataBase.Block]==3) & ([DataBase.urevent]==ur)).st_dev_BlockCond] = st_dev(3);
%         [DataBase(([DataBase.Block]==1) & ([DataBase.urevent]==ur)).Mean_err_BlockCond] = meanErr(1);
%         [DataBase(([DataBase.Block]==2) & ([DataBase.urevent]==ur)).Mean_err_BlockCond] = meanErr(2);
%         [DataBase(([DataBase.Block]==3) & ([DataBase.urevent]==ur)).Mean_err_BlockCond] = meanErr(3);       
%     end
%    
% end

%% Save the output structure
dirname = [cfg.study_folder cfg.preprocessing_folder];
fname = [subject '_DataBase.mat'];
save(fullfile(dirname, fname), 'DataBase')

%% Helper functions
    function ind = findMarker(evts, refInd, eName, direction, limitInd)
        % refInd: reference index in the events structure
        % eName: event name to find
        % direction: direction of search in the event structure: 'ahead' or 'before'
        % limitInd: do not search past this index
        
        switch direction
            case 'ahead'
                next = find(strcmp({evts(refInd+1:limitInd).type},eName),1);
                if ~isempty(next)
                    ind = refInd + next;
                else
                    ind = [];
                end
            case 'before'
                prev = find(strcmp({evts(limitInd:refInd-1).type},eName),1, 'last');
                if ~isempty(prev)
                    ind = limitInd + prev - 1;
                else
                    ind = [];
                end
        end
    end

end