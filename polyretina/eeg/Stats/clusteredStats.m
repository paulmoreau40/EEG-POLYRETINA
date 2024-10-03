function [maxstat_val, cluster_vals, clusters, options_stats] = clusteredStats(data, options_stats)
% Gives the ERSP clustered statistics for N conditions according to the options.
% 1. For every freq x time sample, the statistical test compare the N
%       conditions (ex: t-test if 2 conditions, ANOVA if 3 conditions).
% 2. Select all samples whose statistical value is larger than some
%       threshold (based on the test distribution).
% 3. Cluster the selected samples in connected sets on the basis of
%       temporal and spectral adjacency. Possibility to exclude clusters that are too small
%       with respect to the size of the dataset (see options).
% 4. Calculate cluster-level statistics by taking the sum of the statistical values within a cluster.
% Follows the procedure detailed in:
% Maris, E. & Oostenveld, R. Nonparametric statistical testing of EEG- and MEG-data. J. Neurosci. Methods 164, 177?190 (2007).
%
% Inputs:
% data       - Data. Should be a 1 dimensional cell whose length corresponds to
%                   the number of conditions. For each matrix of 1
%                   condition, the last dimension will be used for
%                   permutations.
% options    -
%
% Outputs:
%   maxstat_val     - Largest of the cluster-level statistics (absolute value).
%   cluster_vals    - Computed statistical values for all clusters.
%   clusters        - [freq x time] matrix containing the information
%                       about the clusters (each freq x time pair is
%                       associated with a cluster index - in correspondence
%                       with the cluster_vals vector - or 0 if it does not
%                       belong to any cluster)

n_conds = length(options_stats.fields);

if n_conds < 2
    error('Cannot perform this analysis with only one condition')
end

if n_conds == 2
    alpha_crit = 0.975;
else
    alpha_crit = 0.95;
end

%% Perform the test comparing the samples:
switch options_stats.model
    case 'classic'
        stats = statcond(data, 'paired', options_stats.pairing,...
            'method', 'param', 'structoutput', 'on','verbose', 'off');
        values = stats.stat;
        
        % Define threshold
        if n_conds == 2
            if strcmp(options_stats.pairing, 'on')
                thresh = tinv(alpha_crit, stats.df);
            else
                thresh = tinv(alpha_crit, mean(stats.df));
            end
        else
            thresh = finv(alpha_crit, stats.df(1), stats.df(2));
        end
        
        [maxstat_val, cluster_vals, clusters] = computeClusterStatsTimeFreq(values, thresh, n_conds,...
            options_stats.removeSmallestClusters);
        
    case 'mixedEffects'
        % Mixed-effects model
        [data_me, trInfo] = buildMEmodel(data, options_stats);
        
        %         options_stats2 = options_stats;
        %     for p = 1:3
        %         inds2keep = setdiff(1:3,p);
        %         options_stats2.fields = options_stats.fields(inds2keep);
        %         pairwiseStats.(['Pair' num2str(p)]).Conditions = options_stats.fields(inds2keep);
        %         [~, all_statVal_pair, all_clusters_pair, ~] = clusteredStats(data4stats(inds2keep), options_stats2);
        %
        %         pairwiseStats.(['Pair' num2str(p)]).StatValsData = all_statVal_pair;
        %         pairwiseStats.(['Pair' num2str(p)]).ClustersData = all_clusters_pair;
        %     end
        
        MEterm = options_stats.MEterm;
        formula = options_stats.formula;
        fields = options_stats.fields;
        main_values = nan(size(data_me,1),size(data_me,2));
        siz = size(main_values);
        main_thresholds = nan(siz);
        if n_conds == 3
            pair1 = [2,3];
            values_pair1 = nan(siz);
            %thresholds_pair1 = nan(siz);
            pair2 = [1,3];
            values_pair2 = nan(siz);
            %thresholds_pair2 = nan(siz);
            pair3 = [1,2];
            values_pair3 = nan(siz);
            %thresholds_pair3 = nan(siz);
            thresholds_pairs = nan(siz);
        end
        
        ppm = ParforProgressbar(numel(main_values), 'showWorkerProgress', true,...
            'title', 'Computing LME for each data point');        
        % Loop over each data point
        parfor xy = 1:numel(main_values)
            [x,y] = ind2sub(siz,xy);            
            % Create appopriate data structure
            data_xy = squeeze(data_me(x,y,:));
            tbl = [trInfo, table(data_xy, 'VariableNames', {'Data'})];
            
            % Fit lme
            lme = fitlme(tbl, formula, 'DummyVarCoding', 'effects');
            if n_conds == 2
                line = find(contains(lme.CoefficientNames, MEterm)...
                    & ~contains(lme.CoefficientNames,':'));
                main_value_temp = table2array(lme.Coefficients(line,4));
                main_threshold_temp = tinv(alpha_crit, table2array(lme.Coefficients(line,5)));
            else
                stats = anova(lme);
                % store values
                line = find(strcmp(stats.Term, MEterm));
                main_value_temp = stats.FStat(line);
                main_threshold_temp = finv(alpha_crit, stats.DF1(line), stats.DF2(line));
                
                % Pairwise statistics
                cond1 = find(contains(lme.CoefficientNames, MEterm)...
                    & contains(lme.CoefficientNames, fields{1})...
                    & ~contains(lme.CoefficientNames,':'));
                vect1 = zeros(1,lme.NumCoefficients);
                vect1(1) = 1; % Due to the 'effects' var coding
                vect1(cond1) = 1;
                cond2 = find(contains(lme.CoefficientNames, MEterm)...
                    & contains(lme.CoefficientNames, fields{2})...
                    & ~contains(lme.CoefficientNames,':'));
                vect2 = zeros(1,lme.NumCoefficients);
                vect2(1) = 1; % Due to the 'effects' var coding
                vect2(cond2) = 1;
                %                     cond3 = find(contains(lme.CoefficientNames, MEterm)...
                %                         & contains(lme.CoefficientNames, options_stats.fields{3})...
                %                         & ~contains(lme.CoefficientNames,':'));
                vect3 = zeros(1,lme.NumCoefficients);
                vect3(1) = 1; % Due to the 'effects' var coding
                vect3(cond1) = -1;
                vect3(cond2) = -1;
                
                [~,value_pair1_temp, ~] = coefTest(lme,vect3-vect2);
                [~,value_pair2_temp, ~] = coefTest(lme,vect3-vect1);
                [~,value_pair3_temp, ~] = coefTest(lme,vect2-vect1);
                threshold_pairs_temp = tinv(alpha_crit, lme.DFE);
            end
            main_values(xy) = main_value_temp;
            main_thresholds(xy) = main_threshold_temp;
            
            if n_conds == 3
                values_pair1(xy) = value_pair1_temp;
                values_pair2(xy) = value_pair2_temp;
                values_pair3(xy) = value_pair3_temp;
                thresholds_pairs(xy) = threshold_pairs_temp;
            end
            ppm.increment();
        end
        delete(ppm);
        
        %{
                ppm = ParforProgressbar(size(data_me,2), 'showWorkerProgress', true,...
            'title', 'Computing LME for each data point');
        parfor y = 1:size(data_me,2)
            main_values_y = nan(size(data_me,1),1);
            main_thresholds_y = nan(size(data_me,1),1);
            
            if n_conds == 3
                values_pair1_y = nan(size(data_me,1),1);
                %thresholds_pair1_y = nan(size(data_me,1),1);
                values_pair2_y = nan(size(data_me,1),1);
                %thresholds_pair2_y = nan(size(data_me,1),1);
                values_pair3_y = nan(size(data_me,1),1);
                %thresholds_pair3_y = nan(size(data_me,1),1);
                thresholds_pairs_y = nan(size(data_me,1),1);
            end
            
            for x = 1:size(data_me,1)
                % Create appopriate data structure
                data_xy = squeeze(data_me(x,y,:));
                tbl = [trInfo, table(data_xy, 'VariableNames', {'Data'})];
                
                % Fit lme
                lme = fitlme(tbl, formula, 'DummyVarCoding', 'effects');
                if n_conds == 2
                    line = find(contains(lme.CoefficientNames, MEterm)...
                        & ~contains(lme.CoefficientNames,':'));
                    main_values_y(x) = table2array(lme.Coefficients(line,4));
                    main_thresholds_y(x) = tinv(alpha_crit, table2array(lme.Coefficients(line,5)));
                else
                    stats = anova(lme);
                    % store values
                    line = find(strcmp(stats.Term, MEterm));
                    main_values_y(x) = stats.FStat(line);
                    main_thresholds_y(x) = finv(alpha_crit, stats.DF1(line), stats.DF2(line));
                    
                    % Pairwise statistics
                    cond1 = find(contains(lme.CoefficientNames, MEterm)...
                        & contains(lme.CoefficientNames, fields{1})...
                        & ~contains(lme.CoefficientNames,':'));
                    vect1 = zeros(1,lme.NumCoefficients);
                    vect1(1) = 1; % Due to the 'effects' var coding
                    vect1(cond1) = 1;
                    cond2 = find(contains(lme.CoefficientNames, MEterm)...
                        & contains(lme.CoefficientNames, fields{2})...
                        & ~contains(lme.CoefficientNames,':'));
                    vect2 = zeros(1,lme.NumCoefficients);
                    vect2(1) = 1; % Due to the 'effects' var coding
                    vect2(cond2) = 1;
                    %                     cond3 = find(contains(lme.CoefficientNames, MEterm)...
                    %                         & contains(lme.CoefficientNames, options_stats.fields{3})...
                    %                         & ~contains(lme.CoefficientNames,':'));
                    vect3 = zeros(1,lme.NumCoefficients);
                    vect3(1) = 1; % Due to the 'effects' var coding
                    vect3(cond1) = -1;
                    vect3(cond2) = -1;
                    
                    [~,values_pair1_y(x), ~] = coefTest(lme,vect3-vect2);
                    %thresholds_pair1_y(x) = tinv(alpha_crit, df1);
                    [~,values_pair2_y(x), ~] = coefTest(lme,vect3-vect1);
                    %thresholds_pair2_y(x) = tinv(alpha_crit, df2);
                    [~,values_pair3_y(x), ~] = coefTest(lme,vect2-vect1);
                    %thresholds_pair3_y(x) = tinv(alpha_crit, df3);
                    thresholds_pairs_y(x) = tinv(alpha_crit, lme.DFE);
                end
            end
            main_values(:,y) = main_values_y;
            main_thresholds(:,y) = main_thresholds_y;
            
            if n_conds == 3
                values_pair1(:,y) = values_pair1_y;
                %thresholds_pair1(:,y) = thresholds_pair1_y;
                values_pair2(:,y) = values_pair2_y;
                %thresholds_pair2(:,y) = thresholds_pair2_y;
                values_pair3(:,y) = values_pair3_y;
                %thresholds_pair3(:,y) = thresholds_pair3_y;
                thresholds_pairs(:,y) = thresholds_pairs_y;
            end            
            ppm.increment();
        end
        delete(ppm);
        %}
        main_thresh = mean(main_thresholds,'all');        
        if n_conds == 2
            [maxstat_val, cluster_vals, clusters] = computeClusterStatsTimeFreq(main_values, main_thresh, n_conds,...
                options_stats.removeSmallestClusters);
        else
            [maxstat_val.Main, cluster_vals.Main, clusters.Main] = computeClusterStatsTimeFreq(main_values, main_thresh, n_conds,...
                options_stats.removeSmallestClusters);
            thresh_pairs = mean(thresholds_pairs,'all');
            [maxstat_val.Pair1, cluster_vals.Pair1, clusters.Pair1] = computeClusterStatsTimeFreq(values_pair1, thresh_pairs, n_conds,...
                options_stats.removeSmallestClusters);
            [maxstat_val.Pair2, cluster_vals.Pair2, clusters.Pair2] = computeClusterStatsTimeFreq(values_pair2, thresh_pairs, n_conds,...
                options_stats.removeSmallestClusters);
            [maxstat_val.Pair3, cluster_vals.Pair3, clusters.Pair3] = computeClusterStatsTimeFreq(values_pair3, thresh_pairs, n_conds,...
                options_stats.removeSmallestClusters);
        end
end

% % Select stats above threshold
% if n_conds == 2
%     aboveTh = abs(values)>thresh;
% else
%     aboveTh = values>thresh;
% end

% %% Form clusters based on temporal/frequency adjacency
% clusters = zeros(size(values));
% cluster_count = 0;
% eq_clusters = [];
% for i = 1:size(values,1)
%     for j = 1:size(values,2)
%         if aboveTh(i,j)
%             % Check if one of the neighbours already belongs to a cluster
%             if j>1 && clusters(i,j-1)~=0 && i>1 && clusters(i-1,j)~=0
%                 % Conflicting clusters?
%                 if clusters(i,j-1) == clusters(i-1,j)
%                     % No
%                     clusters(i,j) = clusters(i-1,j);
%                 else
%                     % Yes, put them in the equivalent clusters list
%                     eq_clusters = [eq_clusters;[clusters(i-1,j),clusters(i,j-1)]];
%                     clusters(i,j) = clusters(i-1,j);
%                 end
%             elseif j>1 && clusters(i,j-1)~=0
%                 clusters(i,j) = clusters(i,j-1);
%             elseif i>1 && clusters(i-1,j)~=0
%                 clusters(i,j) = clusters(i-1,j);
%             else
%                 % Create a new cluster
%                 cluster_count = cluster_count+1;
%                 clusters(i,j) = cluster_count;
%             end
%         end
%     end
% end
%
% % Merge equivalent clusters for the neighbours of this channel
% if ~isempty(eq_clusters)
%     while ~isempty(intersect(unique(eq_clusters(:,1)), unique(eq_clusters(:,2))))
%         % Sort the equivalent clusters such as the lowest id is always in the
%         % first column
%         eq_clusters = sort(eq_clusters,2);
%
%         clusters2merge = unique(eq_clusters(:));
%
%         for cl = clusters2merge(2:end)' % Force the vector to be a line
%             if sum(eq_clusters(:) == cl) > 1
%                 lowestEq = find(eq_clusters(:,2) == cl,1);
%                 if ~isempty(lowestEq)
%                     duplicates = find(eq_clusters(:) == cl);
%                     IndtoAvoid = size(eq_clusters,1) + lowestEq;
%
%                     IndsToReplace = setdiff(duplicates,IndtoAvoid);
%                     eq_clusters(IndsToReplace) = eq_clusters(lowestEq,1);
%                 end
%             end
%         end
%
%         % Remove lines where both clusters are equal
%         eq_clusters(eq_clusters(:,1)==eq_clusters(:,2),:) = [];
%     end
%
%     for eq = 1:size(eq_clusters,1)
%         clusters(clusters == eq_clusters(eq,2)) = eq_clusters(eq,1);
%     end
% end

% if cluster_count == 0
%     maxstat_val = NaN;
%     cluster_vals = [];
% else
%     cluster_count = length(unique(clusters))-1; % There still is 0 values in the matrix
%
%     clust_inds = sort(unique(clusters));
%     clust_inds = clust_inds(2:end);
%
%     if options_stats.removeSmallestClusters
%         removedClust_count = 0;
%         %% Remove clusters that are too small (less than 0.1% of the samples)
%         for cl = 1:cluster_count
%             if numel(find(clusters == clust_inds(cl))) <= critSize
%                 % Remove the cluster
%                 clusters(clusters == clust_inds(cl)) = 0;
%                 removedClust_count = removedClust_count + 1;
%             else
%                 % Renumber the cluster from 1 to cluster_count
%                 clusters(clusters == clust_inds(cl)) = cl - removedClust_count;
%             end
%         end
%         cluster_count = cluster_count - removedClust_count;
%     end
%
%     %% Give the final clusters an index between 1 and cluster_count
%     for cl = 1:cluster_count
%         clusters(clusters == clust_inds(cl)) = cl;
%     end
%
%     %% Calculate cluster levels stats
%     cluster_vals = zeros(cluster_count,1);
%     for cl = 1:cluster_count
%         cluster_vals(cl) = sum(values(clusters == cl));
%     end
%
%     %% Take the largest
%     [~,c_max] = max(abs(cluster_vals));
%
%     % if length(c_max)>1 % multiple maximums
%     %     %Take the largest cluster
%     %     clusters_size = sum(clusters == c_max);
%     %     [~,c_max] = max(clusters_size);
%     %     if length(c_max)>1
%     %         warning('Still multiple maxima')
%     %     end
%     % end
%
%     maxstat_val = cluster_vals(c_max);
%     maxstat_val = abs(maxstat_val(1)); %In case there are multiple maxima
% end

    function [data_me, trInfo] = buildMEmodel(data, options_stats)
        data_me = [];
        for f = 1:numel(options_stats.fields)
            data_me = cat(3,data_me,data{f});
            trInfo_cond = options_stats.(sprintf('trialsInfo_%s', options_stats.fields{f}));
            if f==1
                if strcmp(options_stats.MEterm,'Baseline')
                    if contains(options_stats.fields{f},'base')
                        trInfo = [table(true(size(trInfo_cond,1),1),'VariableNames',{'Baseline'}), trInfo_cond];
                    else
                        trInfo = [table(false(size(trInfo_cond,1),1),'VariableNames',{'Baseline'}), trInfo_cond];
                    end
                else
                    trInfo = trInfo_cond;
                end
            else
                if strcmp(options_stats.MEterm,'Baseline')
                    if contains(options_stats.fields{f},'base')
                        trInfo_temp = [table(true(size(trInfo_cond,1),1),'VariableNames',{'Baseline'}), trInfo_cond];
                    else
                        trInfo_temp = [table(false(size(trInfo_cond,1),1),'VariableNames',{'Baseline'}), trInfo_cond];
                    end
                    trInfo = [trInfo; trInfo_temp];
                else
                    trInfo = [trInfo; trInfo_cond];
                end
            end
        end
    end

    function [maxstat_val, cluster_vals, clusters] = computeClusterStatsTimeFreq(values, thresh, n_conds, rmSmallClusters)
        % Select stats above threshold
        if n_conds == 2
            aboveTh = abs(values)>thresh;
        else
            aboveTh = values>thresh;
        end
        
        %% Form clusters based on temporal/frequency adjacency
        clusters = zeros(size(aboveTh));
        cluster_count = 0;
        eq_clusters = [];
        for i = 1:size(aboveTh,1)
            for j = 1:size(aboveTh,2)
                if aboveTh(i,j)
                    % Check if one of the neighbours already belongs to a cluster
                    if j>1 && clusters(i,j-1)~=0 && i>1 && clusters(i-1,j)~=0
                        % Conflicting clusters?
                        if clusters(i,j-1) == clusters(i-1,j)
                            % No
                            clusters(i,j) = clusters(i-1,j);
                        else
                            % Yes, put them in the equivalent clusters list
                            eq_clusters = [eq_clusters;[clusters(i-1,j),clusters(i,j-1)]];
                            clusters(i,j) = clusters(i-1,j);
                        end
                    elseif j>1 && clusters(i,j-1)~=0
                        clusters(i,j) = clusters(i,j-1);
                    elseif i>1 && clusters(i-1,j)~=0
                        clusters(i,j) = clusters(i-1,j);
                    else
                        % Create a new cluster
                        cluster_count = cluster_count+1;
                        clusters(i,j) = cluster_count;
                    end
                end
            end
        end
        % Merge equivalent clusters for the neighbours of this channel
        if ~isempty(eq_clusters)
            while ~isempty(intersect(unique(eq_clusters(:,1)), unique(eq_clusters(:,2))))
                % Sort the equivalent clusters such as the lowest id is always in the
                % first column
                eq_clusters = sort(eq_clusters,2);
                
                clusters2merge = unique(eq_clusters(:));
                
                for cl = clusters2merge(2:end)' % Force the vector to be a line
                    if sum(eq_clusters(:) == cl) > 1
                        lowestEq = find(eq_clusters(:,2) == cl,1);
                        if ~isempty(lowestEq)
                            duplicates = find(eq_clusters(:) == cl);
                            IndtoAvoid = size(eq_clusters,1) + lowestEq;
                            
                            IndsToReplace = setdiff(duplicates,IndtoAvoid);
                            eq_clusters(IndsToReplace) = eq_clusters(lowestEq,1);
                        end
                    end
                end
                
                % Remove lines where both clusters are equal
                eq_clusters(eq_clusters(:,1)==eq_clusters(:,2),:) = [];
            end
            
            for eq = 1:size(eq_clusters,1)
                clusters(clusters == eq_clusters(eq,2)) = eq_clusters(eq,1);
            end
        end
        
        if cluster_count == 0
            maxstat_val = NaN;
            cluster_vals = [];
        else
            cluster_count = length(unique(clusters))-1; % There still is 0 values in the matrix
            
            clust_inds = sort(unique(clusters));
            clust_inds = clust_inds(2:end);
            
            if rmSmallClusters
                removedClust_count = 0;
                %% Remove clusters that are too small (less than 0.1% of the samples)
                for cl = 1:cluster_count
                    if numel(find(clusters == clust_inds(cl))) <= critSize
                        % Remove the cluster
                        clusters(clusters == clust_inds(cl)) = 0;
                        removedClust_count = removedClust_count + 1;
                    else
                        % Renumber the cluster from 1 to cluster_count
                        clusters(clusters == clust_inds(cl)) = cl - removedClust_count;
                    end
                end
                cluster_count = cluster_count - removedClust_count;
            end
            
            %% Give the final clusters an index between 1 and cluster_count
            for cl = 1:cluster_count
                clusters(clusters == clust_inds(cl)) = cl;
            end
            
            %% Calculate cluster levels stats
            cluster_vals = zeros(cluster_count,1);
            for cl = 1:cluster_count
                cluster_vals(cl) = sum(values(clusters == cl));
            end
            
            %% Take the largest
            [~,c_max] = max(abs(cluster_vals));
            
            % if length(c_max)>1 % multiple maximums
            %     %Take the largest cluster
            %     clusters_size = sum(clusters == c_max);
            %     [~,c_max] = max(clusters_size);
            %     if length(c_max)>1
            %         warning('Still multiple maxima')
            %     end
            % end
            
            maxstat_val = cluster_vals(c_max);
            maxstat_val = abs(maxstat_val(1)); %In case there are multiple maxima
        end
    end
end
