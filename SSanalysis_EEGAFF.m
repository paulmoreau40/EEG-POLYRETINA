%% Start %%
%clear all;
configEEGAFF_Alex;
addpath(genpath('C:\Users\Alexandre\Desktop\code\mne-matlab'))

if ~exist('ALLEEG','var')
    launchEEGLAB;
end

%% Main parameters
merge4MSanalysis = true;
ROIs = {'PPAlh','PPArh','RSClh','RSCrh','OPAlh','OPArh'};
ERSPsclims = [-5,5]; % in dBs
N_colors_standard = 512;
myCmap = asymColorMapWhiteZero(ERSPsclims, N_colors_standard);
set(0,'DefaultFigureColormap',myCmap);

%% Parameters for statistics
comparison = 'condition'; % condition, phase, phase-detail, decision-time
mode = 'trialbased';

% Folders and Files names:
N = makeFolderFileNames(study_config, study_config.subjects(1).id);
delims = strfind(N.searchFolder_3arch_rej,'\');
arch = N.searchFolder_3arch_rej(delims(end-2)+1:end);

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

% ERSPs_comp.Subjects = subjects;
% ERSPs_comp.Times = stc_freq_rois_all.time*1000;
% ERSPs_comp.Freqs = stc_freq_rois_all.freq;
% ERSPs_comp.BaselineModel = 'FullTBNorm_Gain_Log';

plot_params.saveFigFolder = fullfile(study_config.figures_folder,'ERSPs',arch);
plot_params.autoClimsCommon = true;
plot_params.principalFreqValues = [4,8,12,20,30];
plot_params.events.time_ms = [0];
plot_params.events.labels = {{'Intersection','reached'}};

% Type of analysis
switch comparison
    case 'condition'
        options.fields = {'DOWN', 'UP', 'MIX'};
        plot_params.diff_inds = [2,1;3,2;3,1];
    case 'phase'
        options.fields = {'Encoding','Test'};
        plot_params.diff_inds = [2,1];
    case 'phase-detail'
        options.fields = {'Encoding','TestOld','TestNew'};
        plot_params.diff_inds = [2,1;3,2;3,1];
    case 'decision-time'
        options.fields = {'ShortDecTime','LongDecTime'};
        plot_params.diff_inds = [2,1];
end

for subject_ind = subject_inds(1:end)
    
    %subject_ind = 7;
    subject = study_config.subjects(subject_ind).id;
    disp(['Subject ' subject]);
    study_config.current_subject = subject_ind;
    
    % Folders and Files names:
    N = makeFolderFileNames(study_config, subject);
    delims = strfind(N.searchFolder_3arch_rej,'\');
    arch = N.searchFolder_3arch_rej(delims(end-2)+1:end);
    
    %% Baseline ERSPs
    load(fullfile(N.searchFolder_3arch_rej,sprintf('%s_baselines_ERSPsources_ROIs',subject)));
    
    operations = struct();
    operations.singleTrialNorm = true;
    operations.preStimBaseline = false;
    operations.baseline = [];
    
    plot_options = struct();
    plot_options.type = 'Baseline';
    plot_options.style = 'GrandAverage';
    plot_options.clims = ERSPsclims;
    plot_options.singlePlot = false;
    plot_options.nameEvents = false;
    plot_options.events.time_ms = [0];
    plot_options.events.labels = {{'Intersection','reached'}};
    
    figure;
    subplot(2,3,1);
    plot_options.xlabel = false;
    plot_options.ylabel = true;
    plotERSProi(stc_freq_rois_b, 'PPAlh', operations, plot_options)
    subplot(2,3,2);
    plot_options.ylabel = false;
    plotERSProi(stc_freq_rois_b, 'RSClh', operations, plot_options)
    subplot(2,3,3);
    plotERSProi(stc_freq_rois_b, 'OPAlh', operations, plot_options)
    subplot(2,3,4);
    plot_options.xlabel = true;
    plot_options.ylabel = true;
    plotERSProi(stc_freq_rois_b, 'PPArh', operations, plot_options)
    subplot(2,3,5);
    plot_options.ylabel = false;
    plotERSProi(stc_freq_rois_b, 'RSCrh', operations, plot_options)
    subplot(2,3,6);
    plotERSProi(stc_freq_rois_b, 'OPArh', operations, plot_options)
    sgtitle(sprintf('%s - Mean baseline activity for each ROI - %d trials', subject, size(stc_freq_rois_b.trialinfo,1)));
    % Custom color map
    c = colorbar;
    c.Position = [0.93 0.33, 0.015, 0.3];
    c.Label.String = 'Log power (dB)';
    c.Label.FontSize = 12;
    
    saveCurrentFig(fullfile(study_config.figures_folder,'ERSPs',arch),...
        sprintf('%s_GA-ERSProis_baseline', subject), {'png'}, [1600,800]);
    
    %% Main inter ERSPs
    load(fullfile(N.searchFolder_3arch_rej,sprintf('%s_ERSPsources_ROIs',N.epochedFile(1:end-12))));
    
    %% All trials
    % reuse already existing structures
    operations.preStimBaseline = true;
    plot_options.type = 'Intersection';
    plot_options.events.time_ms = [0];
    plot_options.events.labels = {{'Intersection','reached'}};
    
    figure;
    subplot(2,3,1);
    plot_options.xlabel = false;
    plot_options.ylabel = true;
    plotERSProi(stc_freq_rois, 'PPAlh', operations, plot_options)
    subplot(2,3,2);
    plot_options.ylabel = false;
    plotERSProi(stc_freq_rois, 'RSClh', operations, plot_options)
    subplot(2,3,3);
    plotERSProi(stc_freq_rois, 'OPAlh', operations, plot_options)
    subplot(2,3,4);
    plot_options.xlabel = true;
    plot_options.ylabel = true;
    plotERSProi(stc_freq_rois, 'PPArh', operations, plot_options)
    subplot(2,3,5);
    plot_options.ylabel = false;
    plotERSProi(stc_freq_rois, 'RSCrh', operations, plot_options)
    subplot(2,3,6);
    plotERSProi(stc_freq_rois, 'OPArh', operations, plot_options)
    sgtitle(sprintf('%s - Mean activity for each ROI - %d trials', subject, size(stc_freq_rois.trialinfo,1)));
    % Custom color map
    c = colorbar;
    c.Position = [0.93 0.33, 0.015, 0.3];
    c.Label.String = 'Log power wrt prestimulus baseline (dB)';
    c.Label.FontSize = 12;
    
    saveCurrentFig(fullfile(study_config.figures_folder,'ERSPs',arch),...
        sprintf('%s_GA-ERSProis_inter_all', subject), {'png'}, [1600,800]);
    
    %{
    %% Encoding trials
    % reuse already existing structures
    operations.preStimBaseline = true;
    plot_options.type = 'Intersection';
    plot_options.style = 'Encoding';
    plot_options.nameEvents = true;
    plot_options.events.time_ms = [0, 4000];
    plot_options.events.labels = {{'Intersection','reached'},{'Movement','Start'}};
    
    figure;
    subplot(2,3,1);
    plot_options.xlabel = false;
    plot_options.ylabel = true;
    plotERSProi(stc_freq_rois, 'PPAlh', operations, plot_options)
    subplot(2,3,2);
    plot_options.ylabel = false;
    plotERSProi(stc_freq_rois, 'RSClh', operations, plot_options)
    subplot(2,3,3);
    plotERSProi(stc_freq_rois, 'OPAlh', operations, plot_options)
    subplot(2,3,4);
    plot_options.xlabel = true;
    plot_options.ylabel = true;
    plotERSProi(stc_freq_rois, 'PPArh', operations, plot_options)
    subplot(2,3,5);
    plot_options.ylabel = false;
    plotERSProi(stc_freq_rois, 'RSCrh', operations, plot_options)
    subplot(2,3,6);
    plotERSProi(stc_freq_rois, 'OPArh', operations, plot_options)
    encoding_tr = contains(stc_freq_rois.trialinfo.Phase,'Encoding');
    suptitle(sprintf('%s - Mean activity for each ROI during encoding - %d trials', subject, sum(encoding_tr)));
    % Custom color map
    c = colorbar;
    c.Position = [0.93 0.33, 0.015, 0.3];
    c.Label.String = 'Log power wrt prestimulus baseline (dB)';
    c.Label.FontSize = 12;
    
    saveCurrentFig(fullfile(study_config.figures_folder,'ERSPs',arch),...
        sprintf('%s_GA-ERSProis_inter_encoding', subject), {'png'}, [1600,800]);
    
    %% Test trials
    % reuse already existing structures
    operations.preStimBaseline = true;
    plot_options.type = 'Intersection';
    plot_options.style = 'Test';
    test_tr = contains(stc_freq_rois.trialinfo.Phase,'Test');
    plot_options.events.time_ms = [0, mean(stc_freq_rois.trialinfo.MainSequenceDurations_ms(test_tr,2))];
    plot_options.events.labels = {{'Intersection','reached'},{'Average','Movement Start'}};
    
    figure;
    subplot(2,3,1);
    plot_options.xlabel = false;
    plot_options.ylabel = true;
    plotERSProi(stc_freq_rois, 'PPAlh', operations, plot_options)
    subplot(2,3,2);
    plot_options.ylabel = false;
    plotERSProi(stc_freq_rois, 'RSClh', operations, plot_options)
    subplot(2,3,3);
    plotERSProi(stc_freq_rois, 'OPAlh', operations, plot_options)
    subplot(2,3,4);
    plot_options.xlabel = true;
    plot_options.ylabel = true;
    plotERSProi(stc_freq_rois, 'PPArh', operations, plot_options)
    subplot(2,3,5);
    plot_options.ylabel = false;
    plotERSProi(stc_freq_rois, 'RSCrh', operations, plot_options)
    subplot(2,3,6);
    plotERSProi(stc_freq_rois, 'OPArh', operations, plot_options)
    suptitle(sprintf('%s - Mean activity for each ROI during test - %d trials', subject, sum(test_tr)));
    % Custom color map
    c = colorbar;
    c.Position = [0.93 0.33, 0.015, 0.3];
    c.Label.String = 'Log power wrt prestimulus baseline (dB)';
    c.Label.FontSize = 12;
    
    saveCurrentFig(fullfile(study_config.figures_folder,'ERSPs',arch),...
        sprintf('%s_GA-ERSProis_inter_test', subject), {'png'}, [1600,800]);
    %}
        
    %% Subject-level stats
options.pairing = 'off';
for r = 1:numel(ROIs)
    plot_params.ROI = ROIs{r};
    ERSPs = stc_freq_rois.(sprintf('%s_pow',ROIs{r}));
    ERSPs_base_surrog = nan(size(ERSPs_all));
    
        fprintf('Computing statistics for %s - %s\n', subject, ROIs{r});
        plot_params.saveFigName_ave = sprintf('%s_%s_%sAve', subject, ROIs{r}, comparison);
        plot_params.saveFigName_diff = sprintf('%s_%s_%sDiff', subject, ROIs{r}, comparison);
        
        if isfile(fullfile(plot_params.saveFigFolder, plot_params.saveFigName_ave))
            continue
        end
        
        %% Prepare ERSPs
        % Baseline per subject and per condition
        for fld = 1:numel(options.fields)
            switch comparison
                case 'condition'
                    field_trials = contains(stc_freq_rois_all.trialinfo.Condition, options.fields{fld});
            end
            
            prestimBaseline = mean(mean(ERSPs_all(field_trials&subj_trials,:,base_inds),3),1);
            ERSPs_all(field_trials&subj_trials,:,:) = ERSPs_all(field_trials&subj_trials,:,:)...
                ./repmat(prestimBaseline,[sum(field_trials&subj_trials),1,length(stc_freq_rois_all.time)]);
            
            % Create baseline surrogate trials
            %ERSPs_base_surrog(field_trials&subj_trials,:,base_inds) = ERSPs_all(field_trials&subj_trials,:,base_inds);
            for frq = 1:length(stc_freq_rois_all.freq)
                base_distrib = ERSPs_all(field_trials&subj_trials,frq,base_inds);
                base_distrib = base_distrib(:);
                % Simply draw some samples from the distribution
                draw = ceil(length(base_distrib)*...
                    rand(sum(field_trials&subj_trials),length(base_inds)));
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
    
    %% Prepare for MSanalysis
    if merge4MSanalysis
        % Operations
        for r = 1:numel(ROIs)
            % Average each trial over time:
            singletrialbase = mean(stc_freq_rois.(sprintf('%s_pow',ROIs{r})),4);
            % Single trial normalization (for pow only)
            stc_freq_rois.(sprintf('%s_pow',ROIs{r})) = stc_freq_rois.(sprintf('%s_pow',ROIs{r}))./...
                repmat(singletrialbase, [1,1,1,length(stc_freq_rois.time)]);
            % dB transformation:
            %stc_freq_rois.(sprintf('%s_pow',ROIs{r})) = 10*log10(stc_freq_rois.(sprintf('%s_pow',ROIs{r})));
        end
        
        if ~exist('stc_freq_rois_all','var')
            stc_freq_rois_all = struct('cfg',stc_freq_rois.cfg,'method','rawtrial_roiAve',...
                'freq',stc_freq_rois.freq,'time',stc_freq_rois.time);
            stc_freq_rois_all.trialinfo = stc_freq_rois.trialinfo;
            for r = 1:numel(ROIs)
                % Average over all dipoles in the ROI
                stc_freq_rois_all.(sprintf('%s_pow',ROIs{r})) = squeeze(mean(stc_freq_rois.(sprintf('%s_pow',ROIs{r})),1));
                % Store the black screen baseline information (roi and subject specific)
                stc_freq_rois_all.blackScreenBase.(subject).(ROIs{r}) = squeeze(mean(mean(mean(stc_freq_rois_b.(sprintf('%s_pow',ROIs{r})),4),1),2));
            end
        else
            stc_freq_rois_all.trialinfo = [stc_freq_rois_all.trialinfo; stc_freq_rois.trialinfo];
            for r = 1:numel(ROIs)
                stc_freq_rois_all.(sprintf('%s_pow',ROIs{r})) = cat(1,stc_freq_rois_all.(sprintf('%s_pow',ROIs{r})),...
                    squeeze(mean(stc_freq_rois.(sprintf('%s_pow',ROIs{r})),1)));
                stc_freq_rois_all.blackScreenBase.(subject).(ROIs{r}) = squeeze(mean(mean(mean(stc_freq_rois_b.(sprintf('%s_pow',ROIs{r})),4),1),2));
            end
        end
    end
    
    clear stc_freq_rois_b stc_freq_rois
end

if merge4MSanalysis
    save(fullfile(N.searchFolder_4arch_rej,'ERSPsources_ROIs'), 'stc_freq_rois_all', '-v7.3');
end