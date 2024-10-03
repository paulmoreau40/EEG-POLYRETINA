function [data_heatmap_anova_conditionvscondition, data_heatmap_anova_conditionvscondition_cluster_id] = format_for_heatmap_conditionvcondition_anova(clustered_stats_table_conditionvcondition, statistical_clusters_conditionvcondition)

% Get ID of clusters where the values are statistically significant and
% those that are not
valid_cluster_id = clustered_stats_table_conditionvcondition.ClusterID(find(clustered_stats_table_conditionvcondition.ClusterPval < 0.05));
invalid_cluster_id = clustered_stats_table_conditionvcondition.ClusterID(find(clustered_stats_table_conditionvcondition.ClusterPval > 0.05));

% Setting the invalid clusters back to 0
for invalid_cluster = 1:length(invalid_cluster_id)
    statistical_clusters_conditionvcondition(find(statistical_clusters_conditionvcondition == invalid_cluster_id(invalid_cluster))) = 0;
end

% Returning map with cluster ids, so that can be used for intersection
data_heatmap_anova_conditionvscondition_cluster_id = statistical_clusters_conditionvcondition;

% Setting the values of the valud cluster to their statistical values value
for valid_cluster = 1:length(valid_cluster_id)
    statistical_clusters_conditionvcondition(find(statistical_clusters_conditionvcondition == valid_cluster)) = clustered_stats_table_conditionvcondition.ClusterStatVal(clustered_stats_table_conditionvcondition.ClusterID == valid_cluster);
end

data_heatmap_anova_conditionvscondition = statistical_clusters_conditionvcondition;





end
