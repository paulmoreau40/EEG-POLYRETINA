function [stats_all] =...
    ERSPstats(erspCluster_comp, options, params, n_targets, diff_inds)
% The ERSP input is such as is can be averaged (power or z-score but no dB)
n_conds = length(options.fields);

stats_all = struct();

%% Condition vs Baseline
for cnd = 1:n_conds
    %continue
    fprintf('Computing statistics of %s condition vs baseline\n',options.fields{cnd})
    if contains(erspCluster_comp.BaselineModel, '_Gain_Log')
        %data4stats_pow = cell(1,2);
        %data4stats_pow{1} = erspCluster_comp.(options.fields{cnd});
        %data4stats_pow{2} = erspCluster_comp.(sprintf('%s_base',options.fields{cnd}));
        
        data4stats_dB = cell(1,2);
        data4stats_dB{1} = 10.*log10(erspCluster_comp.(options.fields{cnd}));
        data4stats_dB{2} = 10.*log10(erspCluster_comp.(sprintf('%s_base',options.fields{cnd})));
        
        data4stats = data4stats_dB;
    else
        data4stats = cell(1,2);
        data4stats{1} = erspCluster_comp.(options.fields{cnd});
        data4stats{2} = erspCluster_comp.(sprintf('%s_base',options.fields{cnd}));
    end
    
    % options for baseline comparison
    optionsBase = options;
    optionsBase.fields = {options.fields{cnd},sprintf('%s_base',options.fields{cnd})};
    optionsBase.pairing = 'on';
    if isfield(options, 'permutations')
        optionsBase.permutations = options.permutations.(options.fields{cnd});
    end
    
    if strcmp(optionsBase.model,'mixedEffects')
        optionsBase.pairing = 'off';
        optionsBase.MEterm = 'Baseline';
        optionsBase.formula = erspCluster_comp.baseFormula;
        for f = 1:numel(optionsBase.fields)
            optionsBase.(sprintf('trialsInfo_%s',optionsBase.fields{f})) = erspCluster_comp.(sprintf('trialsInfo_%s',options.fields{cnd}));
        end
    end
    
    masks = struct();
    masks.n_masks = 1;
    masks.nameMask1 = 'Clustered T-stats Masks';
    
    %%%%% Permutation statistics (not clustered) %%%%%%
    %stats_perm_struct = statcond(data4stats, 'paired', options2.pairing,...
    %    'method', 'perm', 'naccu', options.N_reps, 'alpha', 0.05,...
    %    'structoutput', 'on','verbose', 'off');
    %mask_perm = stats_perm_struct.mask;
    
    %%%%% Clustered statistics %%%%%
    [clustered_stats_table, statistical_clusters, stats_surrog, ~, permutations] =...
        NP_statTest(erspCluster_comp, optionsBase);
    
    [alpha_test, optionsBase] = alphaMultTargets(0.05, n_targets, params, optionsBase);
    
    signif_clusters = clustered_stats_table.ClusterID(clustered_stats_table.ClusterPval <= alpha_test);
    clustered_mask = false(size(statistical_clusters));
    for signCl = signif_clusters'
        clustered_mask = clustered_mask | statistical_clusters == signCl;
    end
    masks.mask1 = clustered_mask;
    
    %stats_all.(options.fields{cnd}).sampledStats = stats_perm_struct;
    stats_all.(options.fields{cnd}).clusteredStats.permutations = permutations;
    stats_all.(options.fields{cnd}).clusteredStats.Ttest.results = clustered_stats_table;
    stats_all.(options.fields{cnd}).clusteredStats.Ttest.clusters = statistical_clusters;
    stats_all.(options.fields{cnd}).clusteredStats.Ttest.MaxStatValReps = stats_surrog;
    
    stats_all.(options.fields{cnd}).Masks = masks;
    stats_all.(options.fields{cnd}).options = optionsBase;
end

%% Between conditions
fprintf('Computing statistics between conditions\n')
if contains(erspCluster_comp.BaselineModel, '_Gain_Log')
    data4stats_pow = cell(1,n_conds);
    data4stats_dB = cell(1,n_conds);
    %     figure;
    %     hold on;
    for c = 1:n_conds
        data4stats_pow{c} = erspCluster_comp.(options.fields{c});
        %         subplot(2,n_conds,c)
        %         histogram(data4stats_pow{c}(:))
        %         title([options.fields{c} ' - power'])
        %
        data4stats_dB{c} = 10.*log10(erspCluster_comp.(options.fields{c}));
        %         subplot(2,n_conds,c+n_conds)
        %         histogram(data4stats_dB{c}(:))
        %         title([options.fields{c} ' - dB'])
    end
    data4stats = data4stats_dB;
else
    data4stats = cell(1,n_conds);
    for c = 1:n_conds
        data4stats{c} = erspCluster_comp.(options.fields{c});
    end
end

%%%%% Permutation statistics (not clustered) %%%%%%
%stats_perm_struct = statcond(data4stats, 'paired', options.pairing,...
%    'method', 'perm', 'naccu', options.N_reps, 'alpha', 0.05,...
%    'structoutput', 'on','verbose', 'off');

%mask_perm = stats_perm_struct.mask;
%             signif_samps = find(stats_perm_struct.pval <= 0.05);
%             for signSamp = signif_samps'
%                 % Post-hoc analysis
%                 data4posthoc = zeros(n_subj, n_conds);
%                 for c = 1:n_conds
%                     cond_data = erspCluster_comp.(options.fields{c});
%                     for s = 1:n_subj
%                         subj_cond_data = squeeze(cond_data(:,:,s));
%                         data4posthoc(s,c) = subj_cond_data(signSamp);
%                     end
%                 end
%
%                 [~,~,stats] = anova1(data4posthoc, options_comp.fields, 'on');
%                 multcompare(stats);
%             end

%%%%% Clustered statistics %%%%%
optionsBet = options;
if isfield(options, 'permutations')
    optionsBet.permutations = options.permutations.BETWEEN;
end

if strcmp(optionsBet.model,'mixedEffects')
    optionsBet.formula = erspCluster_comp.mainFormula;
    for f = 1:numel(optionsBet.fields)
        optionsBet.(sprintf('trialsInfo_%s',optionsBet.fields{f})) = erspCluster_comp.(sprintf('trialsInfo_%s',optionsBet.fields{f}));
    end
end

[clustered_stats_table, statistical_clusters, stats_surrog, pairwiseStats, permutations] =...
    NP_statTest(erspCluster_comp, optionsBet);

[alpha_test, optionsBet] = alphaMultTargets(0.05, n_targets, params, optionsBet);

signif_clusters = clustered_stats_table.ClusterID(clustered_stats_table.ClusterPval <= alpha_test);
clustered_mask = false(size(statistical_clusters));
for signCl = signif_clusters'
    clustered_mask = clustered_mask | statistical_clusters == signCl;
end

%%%%% Masks definition %%%%%
masks = struct();
masks.n_masks = 1;
masks.nameMask1 = 'Clustered T-stats Masks';

if n_conds == 2
    masks.mask1 = clustered_mask;
elseif n_conds == 3
    masks_pairs = cell(1,3);
    
    [alpha_pair, optionsBet] = alphaPairwiseComp(alpha_test, 3, params, optionsBet);
    
    for p = 1:3
        % Find the row corresponding to the comparison (as stored
        % in diff_inds)
        r = 1;
        while r <=3
            if ~isempty(intersect(diff_inds(r,:), p))
                r = r+1;
            else
                break
            end
        end
        signif_clusters_pair = find(pairwiseStats.(['Pair', num2str(p)]).PValsData <= alpha_pair);
        clustered_mask_pair = false(size(statistical_clusters));
        for signCl = signif_clusters_pair'
            clustered_mask_pair = clustered_mask_pair | pairwiseStats.(['Pair', num2str(p)]).ClustersData == signCl;
        end
        masks_pairs{r} = clustered_mask_pair & clustered_mask;
    end
    
    masks.mask1 = masks_pairs;
end

%stats_all.BETWEEN.sampledStats = stats_perm_struct;
stats_all.BETWEEN.clusteredStats.permutations = permutations;
if n_conds == 2
    stats_all.BETWEEN.clusteredStats.Ttest.results = clustered_stats_table;
    stats_all.BETWEEN.clusteredStats.Ttest.clusters = statistical_clusters;
    stats_all.BETWEEN.clusteredStats.Ttest.MaxStatValReps = stats_surrog;
elseif n_conds == 3
    stats_all.BETWEEN.clusteredStats.ANOVA.results = clustered_stats_table;
    stats_all.BETWEEN.clusteredStats.ANOVA.clusters = statistical_clusters;
    stats_all.BETWEEN.clusteredStats.ANOVA.MaxStatValReps = stats_surrog;
    stats_all.BETWEEN.clusteredStats.Pairwise = pairwiseStats;
end
stats_all.BETWEEN.Masks = masks;
stats_all.BETWEEN.options = optionsBet;

    function [alpha, options] = alphaMultTargets(alpha_init, n_targets, params, options)
        % If you want to apply an additionnal multiple comparison corection
        % accounting for the number of parrallel comparisons on the same
        % dataset (example: n_targets = number of ROIs you consider)
        if params.apply_clusterCorr
            options.multClustersCorrected = 'yes';
            
            switch params.correctionType
                case 'bonf'
                    alpha = alpha_init/nchoosek(n_targets,2); %Bonferroni
                case 'dunn'
                    alpha = 1-(1-alpha_init)^(1/nchoosek(n_targets,2)); %Dunn-Sidak
                otherwise
                    error('Unknown correction type');
            end
        else
            options.multClustersCorrected = 'no';
            alpha = alpha_init;
        end
        options.alpha_test = alpha;
    end

    function [alpha, options] = alphaPairwiseComp(alpha_init, n_pairs, params, options)
        if params.apply_pairwiseCorr
            options.pairwiseCompsCorrected = 'yes';
            switch params.correctionType
                % factor 2 appears because it is two-tailed
                case 'bonf'
                    alpha = alpha_init/(2*nchoosek(n_pairs,2)); %Bonferroni
                case 'dunn'
                    alpha = (1-(1-alpha_init)^(1/nchoosek(n_pairs,2)))/2; %Dunn-Sidak
                otherwise
                    error('Unknown correction type');
            end
        else
            options.pairwiseCompsCorrected = 'no';
            alpha = alpha_init;
        end
        options.alpha_pairwise = alpha;
    end
end