%% Start %%
%clear all;
configEEGAFF_Ainhoa;
addpath(genpath('C:\Users\Ainhoa\Documents\M2\MATLAB\codeM2\mne-matlab'))

if ~exist('ALLEEG','var')
    launchEEGLAB;
end

%% Parameters
ROIs = {'OPAlh','OPArh','PPAlh','PPArh','RSClh','RSCrh'};
comparison = 'condition'; % condition, AbsAff, Level, (Only3Levels)
mode = 'subjbased'; % 'subjbased', 'trialbased'
ERSPsclims = [-5,5]; % in dBs
N_colors_standard = 512;
myCmap = asymColorMapWhiteZero(ERSPsclims, N_colors_standard);
set(0,'DefaultFigureColormap',myCmap)

% Folders and Files names:
N = makeFolderFileNames(study_config, study_config.subjects(1).id);
delims = strfind(N.searchFolder_3arch_rej,'\');
arch = N.searchFolder_3arch_rej(delims(end-2)+1:end);

load(fullfile(N.searchFolder_4arch_rej,'ERSPsources_ROIs'));
load(fullfile(N.searchFolder_4arch_rej,'ERSPsources_ROIs_baselines'));
subjects = unique(stc_freq_rois_all.trialinfo.Sid);

options.model = 'classic';
switch mode
    case 'trialbased'
        options.pairing = 'off'; % Not the same number of trials per subject
    case 'subjbased'
        options.pairing = 'on';
end
options.N_reps = 1000;
options.reusePerms = false;
options.removeSmallestClusters = false;

stat_params.apply_clusterCorr = false;
stat_params.correctionType = 'bonf';
stat_params.apply_pairwiseCorr = true;

ERSPs_comp.Subjects = subjects;
ERSPs_comp.Times = stc_freq_rois_all.time*1000;
ERSPs_comp.Freqs = stc_freq_rois_all.freq;
ERSPs_comp.BaselineModel = 'FullTBNorm_Gain_Log';

plot_params.saveFigFolder = fullfile(study_config.figures_folder,'ERSPs',arch);
plot_params.autoClimsCommon = true;
plot_params.principalFreqValues = [4,8,12,20,30];
plot_params.events.time_ms = [0];
plot_params.events.labels = {{'Observation','Start'}};

% Type of analysis
switch comparison
    case 'condition'
        options.fields = {'Go', 'NoGo'};
        plot_params.diff_inds = [1,2];
    case 'AbsAff'
        options.fields = {'Aff_0', 'Aff_1', 'Aff_2'};
        plot_params.diff_inds = [1,2;1,3;2,3];
    case 'Level'
        options.fields = {'level_1', 'level_2','level_3','level_4','level_5'};
        plot_params.diff_inds = [2,1;3,1;4,1;5,1;3,2;4,2;5,2;4,3;5,3;5,4];
%     case 'Only3Levels'
%         options.fields = {'level_1', 'level_3','level_5'};
%         plot_params.diff_inds = [2,1;3,1;3,2];
end

%% Subject-level
options.pairing = 'off';
for r = 1:numel(ROIs)
    plot_params.ROI = ROIs{r};
    ERSPs_all = stc_freq_rois_all.(sprintf('%s_pow',ROIs{r}));
    ERSPs_all_base = stc_freq_rois_all_b.(sprintf('%s_pow',ROIs{r}));
    ERSPs_base_surrog = nan(size(ERSPs_all));

    for sbj = 1:numel(subjects)
        fprintf('Computing statistics for %s - %s\n', subjects{sbj}, ROIs{r});
        plot_params.saveFigName_ave = sprintf('%s_%s_%sAve', subjects{sbj}, ROIs{r}, comparison);
        plot_params.saveFigName_diff = sprintf('%s_%s_%sDiff', subjects{sbj}, ROIs{r}, comparison);

        if isfile(fullfile(plot_params.saveFigFolder, plot_params.saveFigName_ave))
            continue
        end
        subj_trials = contains(stc_freq_rois_all.trialinfo.Sid, subjects{sbj});
        subj_trials_b = contains(stc_freq_rois_all_b.trialinfo.Sid, subjects{sbj});

        %% Prepare ERSPs
        % Baseline per subject and per condition
        for fld = 1:numel(options.fields)
            switch comparison
                case 'condition'
                    field_trials = strcmp(stc_freq_rois_all.trialinfo.Type, options.fields{fld});
                    field_trials_b = strcmp(stc_freq_rois_all_b.trialinfo.Type, options.fields{fld});
                case 'AbsAff'
                    field_trials = stc_freq_rois_all.trialinfo.AbsAff == str2num(options.fields{fld}(end));
                    field_trials_b = stc_freq_rois_all_b.trialinfo.AbsAff == str2num(options.fields{fld}(end));
%                 case 'Only3Levels'
%                     field_trials = strcmp(stc_freq_rois_all.trialinfo.AffLevel, options.fields{fld});
%                     field_trials_b = strcmp(stc_freq_rois_all_b.trialinfo.AffLevel, options.fields{fld});
            end

            % Create baseline surrogate trials
            %ERSPs_base_surrog(field_trials&subj_trials,:,base_inds) = ERSPs_all(field_trials&subj_trials,:,base_inds);
            for frq = 1:length(stc_freq_rois_all.freq)
                base_distrib = ERSPs_all_base(field_trials_b&subj_trials_b,frq,:);
                base_distrib = base_distrib(:);
                % Simply draw some samples from the distribution
                draw = ceil(length(base_distrib)*...
                    rand(sum(field_trials&subj_trials),size(ERSPs_base_surrog,3)));
                ERSPs_base_surrog(field_trials&subj_trials,frq,:) = base_distrib(draw);
                % Draw many samples and average them together to create new values
                %n_ave = 100;
                %draw = ceil(length(base_distrib)*...
                %    rand(sum(field_trials&subj_trials),length(base_inds),n_ave));
                %ERSPs_base_surrog(field_trials&subj_trials,frq,:) = mean(base_distrib(draw),3);
            end

            ERSPs_comp.(options.fields{fld}) = ERSPs_all(field_trials&subj_trials,:,:);
            ERSPs_comp.(sprintf('%s_base',options.fields{fld})) = ERSPs_base_surrog(field_trials&subj_trials,:,:);
            % Should have freq x time x rep structure
            ERSPs_comp.(options.fields{fld}) = permute(ERSPs_comp.(options.fields{fld}),[2,3,1]);
            ERSPs_comp.(sprintf('%s_base',options.fields{fld})) = permute(ERSPs_comp.(sprintf('%s_base',options.fields{fld})),[2,3,1]);
        end

        %% Statistical analysis
        stats = ERSPstats(ERSPs_comp, options, stat_params, numel(ROIs), plot_params.diff_inds);
        masks = struct();
        for fld = 1:numel(options.fields)
            masks.(options.fields{fld}) = stats.(options.fields{fld}).Masks;
        end
        masks.BETWEEN = stats.BETWEEN.Masks;
        plot_ERSPcomparisons3(ERSPs_comp, options.fields, plot_params, masks);
    end
end

%% Group-level
for r = 1:numel(ROIs)
    fprintf('Computing group-level statistics for %s (%s contrast)\n', ROIs{r}, comparison);
    plot_params.ROI = ROIs{r};
    plot_params.saveFigName_ave = sprintf('Allsubjects_%s_%sAve_%s', ROIs{r}, comparison, mode);
    plot_params.saveFigName_diff = sprintf('Allsubjects_%s_%sDiff_%s', ROIs{r}, comparison, mode);
    
    %% Prepare ERSPs
    ERSPs_all = stc_freq_rois_all.(sprintf('%s_pow',ROIs{r}));
    ERSPs_all_base = stc_freq_rois_all_b.(sprintf('%s_pow',ROIs{r}));
    ERSPs_base_surrog = nan(size(ERSPs_all));
    for fld = 1:numel(options.fields)
        switch comparison
            case 'condition'
                field_trials = strcmp(stc_freq_rois_all.trialinfo.Type, options.fields{fld});
                field_trials_b = strcmp(stc_freq_rois_all_b.trialinfo.Type, options.fields{fld});
            case 'AbsAff'
                field_trials = stc_freq_rois_all.trialinfo.AbsAff == str2num(options.fields{fld}(end));
                field_trials_b = stc_freq_rois_all_b.trialinfo.AbsAff == str2num(options.fields{fld}(end));
        end
        
        % Baseline per subject and per condition
        for sbj = 1:numel(subjects)
            subj_trials = contains(stc_freq_rois_all.trialinfo.Sid, subjects{sbj});
            subj_trials_b = contains(stc_freq_rois_all_b.trialinfo.Sid, subjects{sbj});
            
            % Create baseline surrogate trials
            %ERSPs_base_surrog(field_trials&subj_trials,:,base_inds) = ERSPs_all(field_trials&subj_trials,:,base_inds);
            for frq = 1:length(stc_freq_rois_all.freq)
                base_distrib = ERSPs_all_base(field_trials_b&subj_trials_b,frq,:);
                base_distrib = base_distrib(:);
                % Simply draw some samples from the distribution
                draw = ceil(length(base_distrib)*...
                    rand(sum(field_trials&subj_trials),size(ERSPs_base_surrog,3)));
                ERSPs_base_surrog(field_trials&subj_trials,frq,:) = base_distrib(draw);
                % Draw many samples and average them together to create new values
                %n_ave = 100;
                %draw = ceil(length(base_distrib)*...
                %    rand(sum(field_trials&subj_trials),length(base_inds),n_ave));
                %ERSPs_base_surrog(field_trials&subj_trials,frq,:) = mean(base_distrib(draw),3);
            end
            
            if strcmp(mode,'subjbased')
                % Average all trials from a single subject
                if sbj == 1
                    ERSPs_comp.(options.fields{fld}) = mean(ERSPs_all(field_trials&subj_trials,:,:),1);
                    ERSPs_comp.(sprintf('%s_base',options.fields{fld})) = mean(ERSPs_base_surrog(field_trials&subj_trials,:,:),1);
                else
                    ERSPs_comp.(options.fields{fld}) = cat(1, ERSPs_comp.(options.fields{fld}),...
                        mean(ERSPs_all(field_trials&subj_trials,:,:),1));
                    ERSPs_comp.(sprintf('%s_base',options.fields{fld})) = cat(1, ERSPs_comp.(sprintf('%s_base',options.fields{fld})),...
                        mean(ERSPs_base_surrog(field_trials&subj_trials,:,:),1));
                end
            end
        end
        if strcmp(mode,'trialbased')
            % Keep all trials
            ERSPs_comp.(options.fields{fld}) = ERSPs_all(field_trials,:,:);
            ERSPs_comp.(sprintf('%s_base',options.fields{fld})) = ERSPs_base_surrog(field_trials,:,:);
            ERSPs_comp.(sprintf('trials_%s',options.fields{fld})) = stc_freq_rois_all.trialinfo(field_trials,1:14);
            fprintf('%d trials in %s\n',size(ERSPs_comp.(options.fields{fld}),1),options.fields{fld})
        end
        % Should have freq x time x rep structure
        ERSPs_comp.(options.fields{fld}) = permute(ERSPs_comp.(options.fields{fld}),[2,3,1]);
        ERSPs_comp.(sprintf('%s_base',options.fields{fld})) = permute(ERSPs_comp.(sprintf('%s_base',options.fields{fld})),[2,3,1]);
        
        % load existing permutations to save time (computed with previous ROI)
        if r > 1 && exist('stats', 'var')
            options.permutations.(options.fields{fld}) = stats.(options.fields{fld}).clusteredStats.permutations;
        end
    end
    
    %% Statistical analysis
    %%% Comparison with baseline %%%
    % From newtimef
    %formula = 'mean(arg1,3);';
    %[ERSPboot, ~, ERSPboottrials] = bootstat(ESRPs_comp.(options.fields{1}), formula,...
    %    'boottype', 'shuffle', 'label', 'ERSP', 'bootside', 'both', 'naccu', options.N_reps,...
    %    'basevect', find(ESRPs_comp.Times<0), 'alpha', 0.5, 'dimaccu', 2);
    %res = 10*log10(ERSPboot);
    % From what I understood how to use bootstrap:
    % but what is the formula then?
    %formula = 'mean(arg1,3);';
    %[ERSPboot, ~, ERSPboottrials] = bootstat(ESRPs_comp.(options.fields{1}), formula,...
    %    'boottype', 'shuffle', 'suffledim', [2,3], 'shufflemode', 'regular',...
    %    'label', 'ERSP', 'bootside', 'both', 'naccu', options.N_reps,...
    %    'basevect', find(ESRPs_comp.Times<0), 'alpha', 0.5);
    
    % load existing permutations to save time (computed with previous ROI)
    if r > 1 && exist('stats', 'var')
        options.permutations.BETWEEN = stats.BETWEEN.clusteredStats.permutations;
        options.reusePerms = true;
    end
    stats = ERSPstats(ERSPs_comp, options, stat_params, numel(ROIs), plot_params.diff_inds);
    masks = struct();
    for fld = 1:numel(options.fields)
        masks.(options.fields{fld}) = stats.(options.fields{fld}).Masks;
    end
    masks.BETWEEN = stats.BETWEEN.Masks;
    plot_ERSPcomparisons3(ERSPs_comp, options.fields, plot_params, masks);
end

