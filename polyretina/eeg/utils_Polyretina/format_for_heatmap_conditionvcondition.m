function [data_heatmap_20v45] = format_for_heatmap_conditionvcondition(data_heatmap_anova_conditionvscondition_cluster_id, pairwise_stats_conditionvcondition, alpha_pair, spectrum_trial_20, spectrum_trial_45) %data_heatmap_20v110, data_heatmap_45v110
% Function to create heatmap data to compae each field of view
% 
% INPUT
% data_heatmap_anova_conditionvscondition_cluster_id            - data of heatmap obtained from anova, used for interesection
% pairwise_stats_conditionvcondition                            - pairwise statistics structure used to get comparision condition per condition
% alpha_pair                                                    - alpha  value which is corrected for the pairwise comparison
%
% OUTPUT
% data_heatmap_20v45                                            - heatmap data with clusters comparing the condition 20 vs 45 FoV
% data_heatmap_20v110                                           - heatmap data with clusters comparing the condition 20 vs 110 FoV
% data_heatmap_45v110                                           - heatmap data with clusters comparing the condition 45 vs 110 FoV

% Looping over each pair of the pairwise statistics structure
for pair = 1:3
    
    current_pair = pairwise_stats_conditionvcondition.(['Pair' num2str(pair)]);
    
    % Seeing which condition we are under to fill the correct heatmap
    if strcmp(current_pair.Conditions{1}, 'FoV45') && strcmp(current_pair.Conditions{2}, 'FoV110')
        % 1. Only keep clusters which are significant with corrected alpha
        data_heatmap_45v110 = current_pair;
        data_heatmap_45v110.ClustersData(data_heatmap_45v110.PValsData > alpha_pair) = 0;
        
        % 2. Intersect with anova cluster map: retrieving clustes which are common to both
        indices_anova_clusters_location = find(data_heatmap_anova_conditionvscondition_cluster_id~=0);
        indices_ttest_clusters_location = find(data_heatmap_45v110.ClustersData~=0);
        indices_valid_clusters = intersect(indices_anova_clusters_location,indices_ttest_clusters_location);
        
        % 3. Compute relative dB spectrum
            % Take average over participants
        spectrum_trial_45_averaged = mean(spectrum_trial_45, 3);
        spectrum_trial_110_averaged = mean(spectrum_trial_110, 3);
            % Convert to dB
        spectrum_trial_45_averaged_dB = 10*log10(spectrum_trial_45_averaged);
        spectrum_trial_110_averaged_dB = 10*log10(spectrum_trial_110_averaged);
            % Computing relative spectrum
        spectrum_relative_dB = spectrum_trial_45_averaged_dB - spectrum_trial_110_averaged_dB;
        
            % Setting everything that is not part of the intersection to 0
        spectrum_relative_dB(setdiff(1:end,indices_valid_clusters)) = 0;
        
        data_heatmap_45v110 = spectrum_relative_dB;
        
    elseif strcmp(current_pair.Conditions{1}, 'FoV20') && strcmp(current_pair.Conditions{2}, 'FoV110')
        % 1. Only keep clusters which are significant with corrected alpha
        data_heatmap_20v110 = current_pair;
        data_heatmap_20v110.ClustersData(data_heatmap_20v110.PValsData > alpha_pair) = 0;
        
        % 2. Intersect with anova cluster map: retrieving clustes which are common to both
        indices_anova_clusters_location = find(data_heatmap_anova_conditionvscondition_cluster_id~=0);
        indices_ttest_clusters_location = find(data_heatmap_20v110.ClustersData~=0);
        indices_valid_clusters = intersect(indices_anova_clusters_location,indices_ttest_clusters_location);

        % 3. Compute relative dB spectrum
            % Take average over participants
        spectrum_trial_20_averaged = mean(spectrum_trial_20, 3);
        spectrum_trial_110_averaged = mean(spectrum_trial_110, 3);
            % Convert to dB
        spectrum_trial_20_averaged_dB = 10*log10(spectrum_trial_20_averaged);
        spectrum_trial_110_averaged_dB = 10*log10(spectrum_trial_110_averaged);
            % Computing relative spectrum
        spectrum_relative_dB = spectrum_trial_20_averaged_dB - spectrum_trial_110_averaged_dB;
        
            % Setting everything that is not part of the intersection to 0
        spectrum_relative_dB(setdiff(1:end,indices_valid_clusters)) = 0;
        
        data_heatmap_20v110 = spectrum_relative_dB;

    elseif strcmp(current_pair.Conditions{1}, 'FoV20') && strcmp(current_pair.Conditions{2}, 'FoV45')
        % 1. Only keep clusters which are significant with corrected alpha
        data_heatmap_20v45 = current_pair;
        data_heatmap_20v45.ClustersData(data_heatmap_20v45.PValsData > alpha_pair) = 0;
        
        % 2. Intersect with anova cluster map: retrieving clustes which are common to both
        indices_anova_clusters_location = find(data_heatmap_anova_conditionvscondition_cluster_id~=0);
        indices_ttest_clusters_location = find(data_heatmap_20v45.ClustersData~=0);
        indices_valid_clusters = intersect(indices_anova_clusters_location,indices_ttest_clusters_location);

        % 3. Compute relative dB spectrum
            % Take average over participants
        spectrum_trial_20_averaged = mean(spectrum_trial_20, 3);
        spectrum_trial_45_averaged = mean(spectrum_trial_45, 3);
            % Convert to dB
        spectrum_trial_20_averaged_dB = 10*log10(spectrum_trial_20_averaged);
        spectrum_trial_45_averaged_dB = 10*log10(spectrum_trial_45_averaged);
            % Computing relative spectrum
        spectrum_relative_dB = spectrum_trial_20_averaged_dB - spectrum_trial_45_averaged_dB;
        
            % Setting everything that is not part of the intersection to 0
        spectrum_relative_dB(setdiff(1:end,indices_valid_clusters)) = 0;
        
        data_heatmap_20v45 = spectrum_relative_dB;

    end
end


end


