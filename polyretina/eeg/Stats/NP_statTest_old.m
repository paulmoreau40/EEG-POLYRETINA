function [clustered_stats, all_clusters, stats_surrog, pairwiseStats, history] = NP_statTest_old(data_conds, options_stats)
% Non parametrical statistical test to compare multiple (2 or more) EEG epoched
% conditions (typically ERSPs of 1 IC)
% 1. Collect the trials of the experimental conditions in a single set.
% 2. Permute datasets between conditions.
%       The result of this procedure is called a random partition.
% 3. Calculate the test statistic on this random partition.
% 4. Repeat steps 2 and 3 a large number of times to construct a histogram of the test statistics.
% 5. From the test statistic that was actually observed and the histogram in step 4,
%       calculate the proportion of random partitions that resulted in a larger test statistic
%       than the observed one. This proportion is called the Monte Carlo p-value.
% 6. If the p-value is smaller than the critical alpha-level (0.05),
%       then conclude that the data in the two experimental conditions are significantly different.
%
% Based on:
% Maris, E. & Oostenveld, R. Nonparametric statistical testing of EEG- and MEG-data. J. Neurosci. Methods 164, 177?190 (2007).
%
% Inputs:
% data_conds    - Data. Cell of matrices, 1 per condition. The last
%                   dimension of the matrices corresponds to the repetitions
%                   along which you want the data average (typically subjects
%                   in MS analysis.)
% options_stats - options structure containing at least the following fields:
%                   'fields': cell of strings containing the name of each
%                   condition (as inputed in data_conds)
%                   'pairing': do you want to compute paired statistics or
%                   not ('on'|'off')
%                   'N-reps': Number of repetitions for the permutation
%                   analysis.
%                   'removeSmallestClusters': boolean indicating if you
%                   want to exclude some very small clusters right away
%                   (relatively to the size of the dataset). They will
%                   probably never be significant in the end anyway.
%
% Outputs:
%   cluster_stats   - Table grouping a information summary about each cluster.
%   all_clusters    - Span of all non-corrected significant clusters (chan x time)
%                       in the original data (output of clusteredStats).
%   stats_surrog    - Maximum statistical values issued by the random partition.
%   pairwiseStats   - If the number of conditions is 3, this function will
%                       compute all pairwise comparisons between conditions
%                       on the same permutation occurences as for the 3
%                       conditions analysis. This structure summarizes the
%                       results. Can be used for post-hoc statistics.


n_conds = length(options_stats.fields);
%n_subjs = length(data_conds.Subjects);
if strcmp(options_stats.pairing, 'on')
    n_trials = size(data_conds.(options_stats.fields{1}),3);
else
    n_trials = [];
    for c = 1:n_conds
        n_trials = [n_trials, size(data_conds.(options_stats.fields{c}),3)];
    end
end
%options_stats.n_subjs = n_subjs;

% Do the test on original data:
disp('Testing original data')
data4stats = cell(1,n_conds);
for c = 1:n_conds
    if contains(data_conds.BaselineModel, '_Gain_Log')
        data4stats{c} = 10.*log10(data_conds.(options_stats.fields{c}));
    else
        data4stats{c} = data_conds.(options_stats.fields{c});
    end
end

[~, all_statVal, all_clusters, options_stats] = clusteredStats(data4stats, options_stats);

pairwiseStats = struct();
if n_conds == 3
    % Do pairwise comparisons too
    options_stats2 = options_stats;
    for p = 1:3
        inds2keep = setdiff(1:3,p);
        options_stats2.fields = options_stats.fields(inds2keep);
        pairwiseStats.(['Pair' num2str(p)]).Conditions = options_stats.fields(inds2keep);
        switch options_stats.model
            case 'classic'
                [~, all_statVal_pair, all_clusters_pair, ~] = clusteredStats(data4stats(inds2keep), options_stats2);
                pairwiseStats.(['Pair' num2str(p)]).StatValsData = all_statVal_pair;
                pairwiseStats.(['Pair' num2str(p)]).ClustersData = all_clusters_pair;
            case 'mixedEffects'
                % Pairwise comparisons have already been computed
                pairwiseStats.(['Pair' num2str(p)]).StatValsData = all_statVal.(['Pair' num2str(p)]);
                pairwiseStats.(['Pair' num2str(p)]).ClustersData = all_clusters.(['Pair' num2str(p)]);
        end
    end
    if strcmp(options_stats.model, 'mixedEffects')
        % Change the output to respect later code
        all_statVal = all_statVal.Main;
        all_clusters = all_clusters.Main;
    end
end

disp('Computing permutation statistics')
if isempty(all_statVal)
    warning('No cluster found in the original data')
    clustered_stats = [];
    stats_surrog = [];
    history = 1:n_trials*n_conds;
else
    if strcmp(options_stats.pairing, 'on')
        n_tot = n_trials*n_conds;
        maxReps_th = factorial(n_conds)^n_trials;
        maxReps_th = maxReps_th - 1; % Already one permutation represented by the original data
        maxReps = min(options_stats.N_reps, maxReps_th);
        if maxReps<options_stats.N_reps
            disp('Theoretical maximum number of repetitions reached, computing exact p_value')
        end
        
        if n_conds == 3
            stats_surrog = zeros(maxReps,4);
        else
            stats_surrog = zeros(maxReps,1);
        end
        history = 1:n_tot;
        
        progressBar = waitbar(0,'Initializing', 'Name','NP statistical test');
        for rep = 1:maxReps
            if options_stats.reusePerms
                permutation_inds = options_stats.permutations(rep,:);
            else
                permutation_inds = newRandomPartition(history, n_tot, 'on', n_conds, n_trials);
            end
            
            switch options_stats.model
                case 'classic'
                    % Create new data4stats, keep options_stats
                    [data4stats, ~, history] = createPermutationData(data_conds, permutation_inds,...
                        n_conds, n_trials, 'on', history, options_stats);
                    [stats_surrog(rep,1),~,~, ~] = clusteredStats(data4stats, options_stats);
                case 'mixedEffects'
                    % Create new data4stats, create new options_stats
                    [data4stats, options_stats_perm, history] = createPermutationData(data_conds, permutation_inds,...
                        n_conds, n_trials, 'on', history, options_stats);
                    [stats_surrog_temp,~,~,~] = clusteredStats(data4stats, options_stats_perm);
            end
            
            if n_conds == 2 && strcmp(options_stats.model, 'mixedEffects')
                stats_surrog(rep,1) = stats_surrog_temp;
            end
            if n_conds == 3
                % Do pairwise comparisons too
                options_stats2 = options_stats;
                for p = 1:3
                    inds2keep = setdiff(1:3,p);
                    options_stats2.fields = options_stats.fields(inds2keep);
                    switch options_stats.model
                        case 'classic'
                            [stats_surrog(rep,p+1),~,~, ~] =...
                                clusteredStats(data4stats(inds2keep), options_stats2);
                        case 'mixedEffects'
                            stats_surrog(rep,p+1) = stats_surrog_temp.(['Pair' num2str(p)]);
                    end
                end
                if strcmp(options_stats.model, 'mixedEffects')
                    stats_surrog(rep,1) = stats_surrog_temp.Main;
                end
            end
            
            waitbar(rep/maxReps,progressBar,sprintf('Permutation %4d complete', rep))
        end
        delete(progressBar);
    else
        %n_tot = n_subjs*n_conds;
        n_tot = sum(n_trials);
        maxReps_th = 1;
        %for c = n_conds:-1:2
        %    maxReps_th = maxReps_th*nchoosek(n_subjs*c,n_subjs);
        %end
        for c = 1:n_conds-1
            maxReps_th = maxReps_th*nchoosek(sum(n_trials(c:end)),n_trials(c));
        end
        maxReps_th = maxReps_th - 1; %Already one as the original data
        
        maxReps = min(options_stats.N_reps, maxReps_th);
        if maxReps<options_stats.N_reps
            disp('Theoretical maximum number of repetitions reached, computing exact p_value')
        end
        
        if n_conds == 3
            stats_surrog = zeros(maxReps,4);
        else
            stats_surrog = zeros(maxReps,1);
        end
        history = 1:n_tot;
        
        progressBar = waitbar(0,'Initializing', 'Name','NP statistical test');
        for rep = 1:maxReps
            if options_stats.reusePerms
                permutation_inds = options_stats.permutations(rep,:);
            else
                permutation_inds = newRandomPartition(history, n_tot, 'off', n_conds, n_trials);
            end
            
            switch options_stats.model
                case 'classic'
                    % Create new data4stats, keep options_stats
                    [data4stats, ~, history] = createPermutationData(data_conds, permutation_inds,...
                        n_conds, n_trials, 'off', history, options_stats);
                    [stats_surrog(rep,1),~,~, ~] = clusteredStats(data4stats, options_stats);
                case 'mixedEffects'
                    % Create new data4stats, create new options_stats
                    [data4stats, options_stats_perm, history] = createPermutationData(data_conds, permutation_inds,...
                        n_conds, n_trials, 'off', history, options_stats);
                    [stats_surrog_temp,~,~, ~] = clusteredStats(data4stats, options_stats_perm);
            end
            
            
            if n_conds == 2 && strcmp(options_stats.model, 'mixedEffects')
                stats_surrog(rep,1) = stats_surrog_temp;
            end
            if n_conds == 3
                % Do pairwise comparisons too
                options_stats2 = options_stats;
                for p = 1:3
                    inds2keep = setdiff(1:3,p);
                    options_stats2.fields = options_stats.fields(inds2keep);
                    switch options_stats.model
                        case 'classic'
                            [stats_surrog(rep,p+1),~,~, ~] =...
                                clusteredStats(data4stats(inds2keep), options_stats2);
                        case 'mixedEffects'
                            stats_surrog(rep,p+1) = stats_surrog_temp.(['Pair' num2str(p)]);
                    end
                end
                if strcmp(options_stats.model, 'mixedEffects')
                    stats_surrog(rep,1) = stats_surrog_temp.Main;
                end
            end
            
            waitbar(rep/maxReps,progressBar,sprintf('Permutation %4d complete', rep))
        end
        delete(progressBar);
    end
    
    all_pVals = all_statVal; % just to initialize the size
    for cl = 1:length(all_statVal)
        all_pVals(cl) = sum(stats_surrog(:,1)>abs(all_statVal(cl)))/maxReps;
    end
    
    n_clust = length(all_pVals);
    ClusterID = transpose(1:n_clust);
    ClusterStatVal = all_statVal;
    ClusterPval = all_pVals;
    
    % Extracting median time sample and freq of the cluster
    medSamp = zeros(n_clust,1);
    medTime = zeros(n_clust,1);
    medFreq = zeros(n_clust,1);
    for cl = 1:n_clust
        [clust_freqs, clust_samps] = find(all_clusters == cl);
        medSamp(cl) = round(median(clust_samps));
        medTime(cl) = data_conds.Times(round(median(clust_samps)));
        medFreq(cl) = data_conds.Freqs(round(median(clust_freqs)));
    end
    
    clustered_stats = table(ClusterID, ClusterStatVal, ClusterPval, medSamp, medTime, medFreq);
    
    if n_conds ==3
        for p = 1:3
            n_clust = length(pairwiseStats.(['Pair' num2str(p)]).StatValsData);
            all_pVals_pair = zeros(n_clust,1);
            for cl = 1:n_clust
                all_pVals_pair(cl) = sum(stats_surrog(:,p+1) >...
                    abs(pairwiseStats.(['Pair' num2str(p)]).StatValsData(cl)))/maxReps;
            end
            pairwiseStats.(['Pair' num2str(p)]).PValsData = all_pVals_pair;
            
            pairwiseStats.(['Pair' num2str(p)]).MaxStatValReps = stats_surrog(:,p+1);
        end
        stats_surrog = stats_surrog(:,1);
    end
end

    function permutation_inds = newRandomPartition(history, n_tot, pairing, n_conds, n_trials)
        % Randomly create a new partition and check whether it has already
        % been drawn thanks to the history variable
        % when pairing is 'on', n_trials should be an int (number of trials
        % in each condition assumed to be the same)
        % when pairing is 'off', n_trials should be an array of ints
        % (number of trials per condition)
        
        already_drawn = true;
        while already_drawn
            if strcmp(pairing,'on')
                % Create a random partition (depending on the number of conditions)
                permutation_inds = zeros(1,n_tot);
                
                ind_trials = [];
                for cnd = 1:n_conds
                    ind_trials = [ind_trials; n_trials*(cnd-1)+1:n_trials*cnd];
                end
                
                for cnd = n_conds:-1:1
                    draw_conds = randi(cnd, 1, n_trials);
                    for tr = 1:n_trials
                        permutation_inds(n_trials*(cnd-1)+tr) = ind_trials(draw_conds(tr),tr);
                        ind_trials(draw_conds(tr),tr) = NaN;
                    end
                    if cnd>1
                        ind_trials(isnan(ind_trials)) = [];
                        ind_trials = reshape(ind_trials, cnd-1, n_trials);
                    end
                end
            else
                permutation_inds = randperm(n_tot);
            end
            
            % Check if this partition is already in history:
            inHistory_acc = true(size(history,1),1);
            for cnd = 1:n_conds-1
                if strcmp(pairing,'on')
                    first = n_trials*(cnd-1)+1;
                    last = n_trials*cnd;
                else
                    if cnd == 1
                        first = 1;
                    else
                        first = sum(n_trials(1:cnd-1))+1;
                    end
                    last = sum(n_trials(1:cnd));
                end
                cond_inds = sort(permutation_inds(first:last));
                inHistory = ismember(history(:,first:last), cond_inds, 'rows');
                inHistory_acc = inHistory_acc & inHistory;
                if sum(inHistory_acc) == 0
                    already_drawn = false;
                    break
                end
            end
        end
    end

    function [data4stats, options_stats, history] = createPermutationData(data_conds, permutation_inds, n_conds, n_trials, pairing, history, options_stats)
        % Select data from data_conds according to the permutation_inds
        % and update the history of permutations
        % when pairing is 'on', n_trials should be an int (number of trials
        % in each condition assumed to be the same)
        % when pairing is 'off', n_trials should be an array of ints
        % (number of trials per condition)
        
        data4stats = cell(1,n_conds);
        
        newhist_line = [];
        for cnd = 1:n_conds
            if strcmp(pairing,'on')
                first = n_trials*(cnd-1)+1;
                last = n_trials*cnd;
            else
                if cnd == 1
                    first = 1;
                else
                    first = sum(n_trials(1:cnd-1))+1;
                end
                last = sum(n_trials(1:cnd));
            end
            cond_inds = sort(permutation_inds(first:last));
            
            % Create the data for statistics
            data_cond = zeros(size(data_conds.(options_stats.fields{cnd})));
            if strcmp(options_stats.model, 'mixedEffects')
                trInfo_cond = options_stats.(sprintf('trialsInfo_%s',options_stats.fields{cnd}));
            end
            
            if strcmp(pairing, 'on')
                max_trials = n_trials;
            else
                max_trials = n_trials(cnd);
            end
            
            for tr = 1:max_trials
                if strcmp(pairing, 'on')
                    field = ceil(cond_inds(tr)/n_trials);
                    trial = mod(cond_inds(tr),n_trials);
                    if trial == 0
                        trial = n_trials;
                    end
                else
                    field = find(cumsum(n_trials)>=cond_inds(tr),1);
                    if field == 1
                        trial = cond_inds(tr);
                    else
                        trial = cond_inds(tr) - sum(n_trials(1:field-1));
                    end
                end
                
                if contains(data_conds.BaselineModel, '_Gain_Log')
                    data_cond(:,:,tr) = 10.*log10(data_conds.(options_stats.fields{field})(:,:,trial));
                else
                    data_cond(:,:,tr) = data_conds.(options_stats.fields{field})(:,:,trial);
                end
                
                if strcmp(options_stats.model, 'mixedEffects')
                    if strcmp(options_stats.MEterm, 'Baseline')
                        trInfo_cond(tr,:) = options_stats.(sprintf('trialsInfo_%s',options_stats.fields{field}))(trial,:);
                    else
                        cols2update = ~strcmp(trInfo_cond.Properties.VariableNames,options_stats.MEterm);
                        trInfo_cond(tr,cols2update) = options_stats.(sprintf('trialsInfo_%s',options_stats.fields{field}))(trial,cols2update);
                    end
                end
            end
            
            data4stats{cnd} = data_cond;
            if strcmp(options_stats.model, 'mixedEffects')
                options_stats.(sprintf('trialsInfo_%s',options_stats.fields{cnd})) = trInfo_cond;
            end
            
            % Store the successful permutation in history
            newhist_line = [newhist_line, cond_inds];
        end
        history = [history; newhist_line];
    end
end
