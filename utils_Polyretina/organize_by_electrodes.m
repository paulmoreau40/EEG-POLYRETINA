function [reorganized_data_heatmap, new_electrode_labels] = organize_by_electrodes(data_heatmap, old_electrode_labels, sort_alphabetically)
% Re-organizing the electrodes by Brain Regions of interest

num_electrodes = length(old_electrode_labels);

reorganized_data_heatmap = [];

% Defining the electrodes of the brain regions concerned
frontal_electrodes = {'Z1', 'Z2', 'Z3', 'Z4', 'Z5',...
                      'R1', 'R2', 'R3', 'R4', 'R5', 'R6',...
                      'RR1', 'RR2', 'RR3','RR4', 'RR5',...
                      'RD1', 'RD2', 'RD3',...
                      'RC1', 'RC2', 'RC3',...
                      'RB1', 'RB2',...
                      'RA1', 'RA2',...
                      'RE1',...
                      'L1', 'L2', 'L3', 'L4', 'L5', 'L6',...
                      'LL1', 'LL2', 'LL3', 'LL4', 'LL5',...
                      'LD1', 'LD2', 'LD3',...
                      'LC1', 'LC2', 'LC3',...
                      'LB1', 'LB2', ...
                      'LA1', 'LA2',...
                      'LE1'};
parietal_electrodes = {'Z6', 'LL6', 'RR6', 'R7', 'L7', 'RR7', 'LL7', 'LA3', 'RA3', 'LA4', 'RA4', 'R8', 'L8'};
central_electrodes = {'Z8', 'LL8', 'RR8', 'LA5', 'RA5', 'L9', 'R9', 'Z9', 'LL9', 'RR9', 'L10', 'R10', 'RR10', 'LL10'};
occipital_electrodes = {'Z10', 'Z11', 'Z12', 'Z13', 'Z14',...
                        'R11', 'R12', 'R13', 'R14',...
                        'L11', 'L12', 'L13', 'L14',...
                        'LL11', 'LL12', 'LL13',...
                        'RR11', 'RR12', 'RR13', ...
                        'RE4', 'LE4'};
temporal_electrodes = {'RB3', 'RB4', 'RB5', 'RB6',...
                       'RC4', 'RC5', 'RC6', 'RC7',...
                       'RD4', 'RD5', 'RD6', 'RD7',...
                       'RE2', 'RE3', 'Rm',...
                       'LB3', 'LB4', 'LB5', 'LB6',...
                       'LC4', 'LC5', 'LC6', 'LC7', ...
                       'LD4', 'LD5', 'LD6', 'LD7',...
                       'LE2', 'LE3', 'Lm'};

% Separating heatmap data according to brain RoI
frontal_heatmap_data = [];
parietal_heatmap_data = [];
occipital_heatmap_data = [];
central_heatmap_data = [];
temporal_heatmap_data = [];

% Defining new structure containing strings of new electrode setup:
new_electrode_labels = {};
electrode_labels_frontal = {};
electrode_labels_parietal = {};
electrode_labels_occipital = {};
electrode_labels_central = {};
electrode_labels_temporal = {};


% Looping over every electrode
for electrode_index = 1:num_electrodes
    
    % Getting name of current electrode
    current_electrode_name = old_electrode_labels{electrode_index};
    
    % Seeing which brain area it corresponds to
    if ismember(current_electrode_name, frontal_electrodes)        
        % Updating frontal heatmap data:
        frontal_heatmap_data = [frontal_heatmap_data; data_heatmap(electrode_index, :)];
        % Updating label order of electrodes:
        if isempty(electrode_labels_frontal)
            electrode_labels_frontal{1} = current_electrode_name;
        else
            electrode_labels_frontal{end+1} = current_electrode_name;
        end
        
    elseif ismember(current_electrode_name,parietal_electrodes)
        % Updating parietal heatmap data:
        parietal_heatmap_data = [parietal_heatmap_data; data_heatmap(electrode_index, :)];
        % Updating label order of electrodes:
        if isempty(electrode_labels_parietal)
            electrode_labels_parietal{1} = current_electrode_name;
        else
            electrode_labels_parietal{end+1} = current_electrode_name;
        end        
    elseif ismember(current_electrode_name,occipital_electrodes)
        % Updating parietal heatmap data:
        occipital_heatmap_data = [occipital_heatmap_data; data_heatmap(electrode_index, :)];
        % Updating label order of electrodes:
        if isempty(electrode_labels_occipital)
            electrode_labels_occipital{1} = current_electrode_name;
        else
            electrode_labels_occipital{end+1} = current_electrode_name;
        end        
        
    elseif ismember(current_electrode_name, central_electrodes)
        % Updating parietal heatmap data:
        central_heatmap_data = [central_heatmap_data; data_heatmap(electrode_index, :)];
        % Updating label order of electrodes:
        if isempty(electrode_labels_central)
            electrode_labels_central{1} = current_electrode_name;
        else
            electrode_labels_central{end+1} = current_electrode_name;
        end        
        
    elseif ismember(current_electrode_name, temporal_electrodes)
        % Updating parietal heatmap data:
        temporal_heatmap_data = [temporal_heatmap_data; data_heatmap(electrode_index, :)];
        % Updating label order of electrodes:
        if isempty(electrode_labels_temporal)
            electrode_labels_temporal{1} = current_electrode_name;
        else
            electrode_labels_temporal{end+1} = current_electrode_name;
        end           
    end
        
end

% If true, then sorts electrodes alphabetically 
if sort_alphabetically
    % Sort frontal electrodes
    [sorted_labels_frontal, sorted_indices_frontal] = sort(electrode_labels_frontal);
    % Sort parietal electrodes
    [sorted_labels_parietal, sorted_indices_parietal] = sort(electrode_labels_parietal);
    % Sort central electrodes
    [sorted_labels_central, sorted_indices_central] = sort(electrode_labels_central);
    % Sort temporal electrodes
    [sorted_labels_temporal, sorted_indices_temporal] = sort(electrode_labels_temporal);
    % Sort occipital electrodes
    [sorted_labels_occipital, sorted_indices_occipital] = sort(electrode_labels_occipital);

    % Concatenating name of electrodes and electrode data by group
    reorganized_data_heatmap = [frontal_heatmap_data(sorted_indices_frontal,:); parietal_heatmap_data(sorted_indices_parietal,:); central_heatmap_data(sorted_indices_central,:); temporal_heatmap_data(sorted_indices_temporal,:); occipital_heatmap_data(sorted_indices_occipital,:)];

    for current_electrode = 1:length(sorted_labels_frontal)
        if isempty(new_electrode_labels)
            new_electrode_labels{1} = sorted_labels_frontal{current_electrode};
        else
            new_electrode_labels{end+1} = sorted_labels_frontal{current_electrode};
        end
    end
    for current_electrode = 1:length(sorted_labels_parietal)
        new_electrode_labels{end+1} = sorted_labels_parietal{current_electrode};
    end
    for current_electrode = 1:length(sorted_labels_central)
        new_electrode_labels{end+1} = sorted_labels_central{current_electrode};
    end
    for current_electrode = 1:length(sorted_labels_temporal)
        new_electrode_labels{end+1} = sorted_labels_temporal{current_electrode};
    end
    for current_electrode = 1:length(sorted_labels_occipital)
        new_electrode_labels{end+1} = sorted_labels_occipital{current_electrode};
    end

else
    % Concatenating name of electrodes and electrode data by group
    reorganized_data_heatmap = [frontal_heatmap_data; parietal_heatmap_data; central_heatmap_data; temporal_heatmap_data; occipital_heatmap_data];

    for current_electrode = 1:length(electrode_labels_frontal)
        if isempty(new_electrode_labels)
            new_electrode_labels{1} = electrode_labels_frontal{current_electrode};
        else
            new_electrode_labels{end+1} = electrode_labels_frontal{current_electrode};
        end
    end
    for current_electrode = 1:length(electrode_labels_parietal)
        new_electrode_labels{end+1} = electrode_labels_parietal{current_electrode};
    end
    for current_electrode = 1:length(electrode_labels_central)
        new_electrode_labels{end+1} = electrode_labels_central{current_electrode};
    end
    for current_electrode = 1:length(electrode_labels_temporal)
        new_electrode_labels{end+1} = electrode_labels_temporal{current_electrode};
    end
    for current_electrode = 1:length(electrode_labels_occipital)
        new_electrode_labels{end+1} = electrode_labels_occipital{current_electrode};
    end
end

end
