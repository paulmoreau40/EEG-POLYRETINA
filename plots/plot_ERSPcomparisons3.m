function plot_ERSPcomparisons3(erspCluster, fields, params_plot, masks)
n_conds = length(fields);

switch erspCluster.BaselineModel
    case 'SingleTrial_Gain_Log'
        ColorBarLegend = 'Log Power wrt Baseline (dB)';
    case 'FullTBNorm_Additive'
        ColorBarLegend = 'Power difference with Baseline (z-score)';
    case 'FullTBNorm_Gain_Log'
        ColorBarLegend = 'Log Power wrt Baseline (dB)';
end

%% ERSPs per condition
% Keep the data in one var
Data = zeros(length(erspCluster.Freqs),length(erspCluster.Times), n_conds);
for f = 1:n_conds
    Data(:,:,f) = mean(erspCluster.(fields{f}),3);
end
if contains(erspCluster.BaselineModel, '_Gain_Log')
    % Convert back to dB:
    Data = 10.*log10(Data);
end

titles = fields;
plot_options = struct();
plot_options.singlePlot = false;
plot_options.plotEvents = true;
plot_options.timewarped = false;
if params_plot.autoClimsCommon
    plot_options.clims = [min(min(Data(:)),0), max(max(Data(:)),0)];
elseif isfield(params_plot, 'minPowPlot_conds')
    plot_options.clims = [params_plot.minPowPlot_conds, params_plot.maxPowPlot_conds];
end
plot_options.cLegend = ColorBarLegend;

%% Custom color map
N_colors_standard = 512;
myCmap = asymColorMapWhiteZero(plot_options.clims, N_colors_standard);
set(0,'DefaultFigureColormap',myCmap)

plot_options.ylabel = true;

if exist('masks', 'var')
    for f = 2:n_conds
        if masks.(fields{f-1}).n_masks ~= masks.(fields{f-1}).n_masks
            error('All conditions should have the same number of masks')
        end
    end
    
    for m = 1:masks.(fields{1}).n_masks
        figure
        for f = 1:n_conds
            subplot(n_conds,1,f)
            if ~iscell(masks.(fields{f}).(['mask',num2str(m)]))
                plot_options.mask = masks.(fields{f}).(['mask',num2str(m)]);
            end
            
            plot_options.title = titles{f};
            if f == n_conds
                plot_options.nameEvents = true;
            else
                plot_options.nameEvents = false;
            end
            ersp_fieldAve = mean(erspCluster.(fields{f}),3);
            if contains(erspCluster.BaselineModel, '_Gain_Log')
                % Convert back to dB:
                ersp_fieldAve = 10.*log10(ersp_fieldAve);
            end
            
            makeERSPplot(ersp_fieldAve, erspCluster.Times,...
                erspCluster.Freqs, params_plot.principalFreqValues, params_plot.events, plot_options);
            
            if ~isfield(plot_options, 'clims')
                c = colorbar;
                c.Position = [0.92 0.33*(n_conds-f)+0.09, 0.015, 0.15];
                c.Label.String = ColorBarLegend;
                c.Label.FontSize = 11;
            end
        end
        
        if isfield(plot_options, 'clims')
            c = colorbar;
            c.Position = [0.92 0.25, 0.02, 0.5];
            c.Label.String = ColorBarLegend;
            c.Label.FontSize = 12;
        end
        
        sgtitle([num2str(params_plot.ROI),...
            ' - Average Masked with ', masks.(fields{1}).(['nameMask', num2str(m)])])
        
        saveCurrentFig(params_plot.saveFigFolder, [params_plot.saveFigName_ave '-Mask', num2str(m)], {'fig'}, []);
    end
else
    figure
    for f = 1:n_conds
        subplot(n_conds,1,f)
        plot_options.title = titles{f};
        if f == n_conds
            plot_options.nameEvents = true;
        else
            plot_options.nameEvents = false;
        end
        
        ersp_fieldAve = mean(erspCluster.(fields{f}),3);
        if contains(erspCluster.BaselineModel, '_Gain_Log')
            % Convert back to dB:
            ersp_fieldAve = 10.*log10(ersp_fieldAve);
        end
        
        makeERSPplot(ersp_fieldAve, erspCluster.Times,...
            erspCluster.Freqs, params_plot.principalFreqValues, params_plot.events, plot_options);
        if ~isfield(plot_options, 'clims')
            c = colorbar;
            c.Position = [0.92 0.33*(n_conds-f)+0.09, 0.015, 0.15];
            c.Label.String = ColorBarLegend;
            c.Label.FontSize = 11;
        end
    end
    
    if isfield(plot_options, 'clims')
        c = colorbar;
        c.Position = [0.92 0.25, 0.02, 0.5];
        c.Label.String = ColorBarLegend;
        c.Label.FontSize = 12;
    end
    
    sgtitle([num2str(params_plot.ROI), ' - Average'])
    
    saveCurrentFig(params_plot.saveFigFolder, [params_plot.saveFigName_ave, '-Unmasked'], {'fig'}, []);
end

%% Differences
% Keep the diff data in one var
n_diffs = size(params_plot.diff_inds,1);
diffData = zeros(length(erspCluster.Freqs),length(erspCluster.Times), n_diffs);
diff_titles = cell(1, n_diffs);
for f = 1:n_diffs
    switch erspCluster.BaselineModel
        case 'SingleTrial_Gain_Log'
            % Difference of averages
            phase1_all = 10.*log10(mean(erspCluster.(fields{params_plot.diff_inds(f,1)}),3));
            phase2_all = 10.*log10(mean(erspCluster.(fields{params_plot.diff_inds(f,2)}),3));
        case 'FullTBNorm_Additive'
            % Difference of averages
            phase1_all = mean(erspCluster.(fields{params_plot.diff_inds(f,1)}),3);
            phase2_all = mean(erspCluster.(fields{params_plot.diff_inds(f,2)}),3);
        case 'FullTBNorm_Gain_Log'
            % Difference of averages
            phase1_all = 10.*log10(mean(erspCluster.(fields{params_plot.diff_inds(f,1)}),3));
            phase2_all = 10.*log10(mean(erspCluster.(fields{params_plot.diff_inds(f,2)}),3));
    end
    
    diffData(:,:,f) = phase1_all - phase2_all;
    diff_titles{f} = [titles{params_plot.diff_inds(f,1)}, ' > ', titles{params_plot.diff_inds(f,2)}];
end

plot_options = struct();
plot_options.singlePlot = false;
plot_options.plotEvents = true;
plot_options.timewarped = false;
if params_plot.autoClimsCommon
    plot_options.clims = [min(min(diffData(:)),0), max(max(diffData(:)),0)];
elseif isfield(params_plot, 'minPowPlot_diff')
    plot_options.clims = [params_plot.minPowPlot_diff, params_plot.maxPowPlot_diff];
end
plot_options.cLegend = ColorBarLegend;

%% Custom color map
N_colors_standard = 512;
myCmap = asymColorMapWhiteZero(plot_options.clims, N_colors_standard);
set(0,'DefaultFigureColormap',myCmap)

plot_options.ylabel = true;

if exist('masks', 'var')
    for m = 1:masks.BETWEEN.n_masks
        if ~iscell(masks.BETWEEN.(['mask',num2str(m)]))
            plot_options.mask = masks.BETWEEN.(['mask',num2str(m)]);
        end
        
        figure
        if n_diffs == 1
            plot_options.singlePlot = true;
            plot_options.XLabel = 'Time (ms)';
            plot_options.title = {[num2str(params_plot.ROI),...
                ' - Differences Masked with ', masks.BETWEEN.(['nameMask', num2str(m)])],diff_titles{f}};
            plot_options.nameEvents = true;
            makeERSPplot(diffData(:,:,f), erspCluster.Times,...
                erspCluster.Freqs, params_plot.principalFreqValues, params_plot.events, plot_options);
            if isfield(plot_options, 'clims')
                c = colorbar;
                c.Position = [0.92 0.25, 0.02, 0.5];
                c.Label.String = ColorBarLegend;
                c.Label.FontSize = 12;
            end
        else
            for f = 1:n_diffs
                subplot(n_diffs,1,f)
                if iscell(masks.BETWEEN.(['mask',num2str(m)]))
                    plot_options.mask = masks.BETWEEN.(['mask',num2str(m)]){f};
                end
                
                plot_options.title = diff_titles{f};
                if f == n_diffs
                    plot_options.nameEvents = true;
                else
                    plot_options.nameEvents = false;
                end
                
                makeERSPplot(diffData(:,:,f), erspCluster.Times,...
                    erspCluster.Freqs, params_plot.principalFreqValues, params_plot.events, plot_options);
                if ~isfield(plot_options, 'clims')
                    c = colorbar;
                    c.Position = [0.92 0.33*(n_diffs-f)+0.09, 0.015, 0.15];
                    c.Label.String = ColorBarLegend;
                    c.Label.FontSize = 11;
                end
            end
            
            if isfield(plot_options, 'clims')
                c = colorbar;
                c.Position = [0.92 0.25, 0.02, 0.5];
                c.Label.String = ColorBarLegend;
                c.Label.FontSize = 12;
            end
            
            sgtitle([num2str(params_plot.ROI),...
                ' - Differences Masked with ', masks.BETWEEN.(['nameMask', num2str(m)])])
        end
        saveCurrentFig(params_plot.saveFigFolder,...
            [params_plot.saveFigName_diff, '-Mask', num2str(m)], {'fig'}, []);
    end
else
    figure
    for f = 1:n_diffs
        subplot(n_diffs,1,f)
        plot_options.title = diff_titles{f};
        makeERSPplot(diffData(:,:,f), erspCluster.Times,...
            erspCluster.Freqs, params_plot.principalFreqValues, params_plot.events, plot_options);
        if ~isfield(plot_options, 'clims')
            c = colorbar;
            c.Position = [0.92 0.33*(n_diffs-f)+0.09, 0.015, 0.15];
            c.Label.String = ColorBarLegend;
            c.Label.FontSize = 11;
        end
    end
    
    if isfield(plot_options, 'clims')
        c = colorbar;
        c.Position = [0.92 0.25, 0.02, 0.5];
        c.Label.String = ColorBarLegend;
        c.Label.FontSize = 12;
    end
    
    suptitle([num2str(params_plot.ROI), ' - Differences Unmasked'])
    
    saveCurrentFig(params_plot.saveFigFolder,...
        [params_plot.saveFigName_diff, '-Unmasked'], {'fig'}, []);
end
set(0,'DefaultFigureColormap',parula)
end