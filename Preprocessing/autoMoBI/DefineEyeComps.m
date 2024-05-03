function [eye_comps] = DefineEyeComps(classification_res, thresholds, crit, save_path_figs, filename)
%% Computes the eye components predicted by ICLabel according to the given criterion
% Also plots the distribution of eye components if save_path_figs is not empty
%
% Inputs:
% classification_res        - classification for all ICs, all classes,
%                               given by ICLabel. [n_IC x n_classes] matrix
% thresholds                - thresholds for each class (necessary for the
%                               uniqueness criterion)
% crit                      - criterion for selecting eye components.
%                               possible values: 'all', 'first', 'unique'
% save_path_figs            -  path for saving figures (replace by [] if you do not want to save)
% filename                  -  filename (necessary if you
%                               want to save figures)
%
% Outputs:
% eye_comps                 - the selection of components that are defined
%                               as eyes according to the criterion.

%% Some useful calculations
% 3 should be the column for eyes
eye_above_th = find(classification_res(:,3) > thresholds(3));
[~, most_probable_class] = max(classification_res,[], 2);
eye_most_probable = find(most_probable_class == 3);

classes_count = zeros(size(classification_res,1),1);
for c = 1:length(thresholds)
    classes_count = classes_count + (classification_res(:,c) > thresholds(c));
end

%% Eye components definition
switch crit
    case 'all'
        % set eye component if it is above threshold
        eye_comps = eye_above_th;
    case 'first'
        % set eye component if it is above threshold & most probable class
        eye_comps = intersect(eye_most_probable, eye_above_th);
    case 'unique'
        % set eye component if it is the only class predicted by IClabel
        eye_comps = intersect(find(classes_count==1), eye_above_th);
    otherwise
        error('Criterion for selecting eye components not understood')
end

%% plots
if ~isempty(save_path_figs)
    %{
    figure
    hold on
    yyaxis left
    plot([1,2], [length(eye_above_th), length(eye_most_probable)], 'lineWidth', 2)
    ylabel('eye IC count')
    yyaxis right
    plot([1,2], 100*[mean(classification_res(eye_above_th,3)), mean(classification_res(eye_most_probable,3))],...
        'lineWidth', 2)
    ylabel('mean % certainty of being eye IC')
    xlim([0,3])
    xticks(1:2)
    xticklabels({'Above threshold','Most probable'});
    xlabel('Criterion type')
    %}
    if floor(100*thresholds(3)) == ceil(100*thresholds(3)/5)*5
        edges = ceil(100*thresholds(3)/5)*5:5:100;
        ticks = edges(1:2:end);
    else
        edges = [floor(100*thresholds(3)), ceil(100*thresholds(3)/5)*5:5:100];
        ticks = edges(2:2:end);
    end
    figure
    hold on
    histogram(100*classification_res(eye_above_th,3), edges, 'FaceAlpha', 1)
    histogram(100*classification_res(eye_most_probable,3), edges, 'FaceAlpha', 0.7)
    if strcmp(crit,'unique')
        histogram(100*classification_res(eye_comps,3), edges, 'FaceAlpha', 0.4)
    end
    legend({['Above threshold: ',num2str(length(eye_above_th)),' ICs'],...
        ['Best prediction: ',num2str(length(eye_most_probable)),' ICs'],...
        ['Unique label: ',num2str(length(eye_comps)),' ICs']})
    ylabel('Eye IC count')
    xlabel('% certainty with ICLabel')
    xticks(ticks)
    title({['Eye ICs distribution after AMICA1'],...
        ['Threshold = ' num2str(thresholds(3)*100) '%']})
    saveCurrentFig(save_path_figs, filename, {'png'}, [])
end
end

