function [data_heatmap_conditionvbase] = format_for_heatmap_baselinecorrected_dB(clustered_stats_table_conditionvbase,statistical_clusters_conditionvbase, spectrum_trial_dB_1, spectrum_trial_dB_2)
% Format the data so that can be used for the heatmap to represent
% statistically significant clusters


% 1. Computing the difference in dB between both spectrum after taking the
% average across participants

% Converting dB signals into linear so that we can take the average
spectrum_trial_1_averaged_linear = mean(10.^(spectrum_trial_dB_1./10),3);
spectrum_trial_2_averaged_linear = mean(10.^(spectrum_trial_dB_2./10),3);

% Transforming them back in dB : 
spectrum_trial_1_averaged_dB = 10*log10(spectrum_trial_1_averaged_linear);
spectrum_trial_2_averaged_dB = 10*log10(spectrum_trial_2_averaged_linear);
% Taking the difference in dB
spectrum_relative_averaged_dB = spectrum_trial_1_averaged_dB - spectrum_trial_2_averaged_dB;

% 2. Setting all of the non_statistically_significant clusters to 0 and
% keep only valid clusters, whose values are the difference in dB

% Get ID of clusters where the values are statistically significant and those that are not
valid_cluster_id = clustered_stats_table_conditionvbase.ClusterID(find(clustered_stats_table_conditionvbase.ClusterPval < 0.05));
invalid_cluster_id = clustered_stats_table_conditionvbase.ClusterID(find(clustered_stats_table_conditionvbase.ClusterPval >= 0.05));

% Setting the invalid clusters back to 0
for invalid_cluster = 1:length(invalid_cluster_id)
    statistical_clusters_conditionvbase(find(statistical_clusters_conditionvbase == invalid_cluster_id(invalid_cluster))) = 0;
end

% Retrieving all of the indices where in statistical clusters map where null
indices_non_significant = find(statistical_clusters_conditionvbase == 0);

% Setting the dB difference map to 0 where not statistically significant
spectrum_relative_averaged_dB(indices_non_significant) = 0;

data_heatmap_conditionvbase = spectrum_relative_averaged_dB;

end
