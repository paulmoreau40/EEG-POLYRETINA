function makeERSPplot(ERSP, times, freqs, principalFreqValues, events, options)
% Takes ERSP input in its final form (no data computation in this function)

hold on;
if isfield(options, 'mask')
    if isfield(options, 'clims')
        % Problem when clims and alphData at the same time
        %imagesc(times, freqs, ERSP, options.clims, 'AlphaData', options.mask);
        im = imagesc(times, freqs, ERSP, options.clims);
        im.AlphaData = options.mask;
    else
        imagesc(times, freqs, ERSP, 'AlphaData', options.mask);
    end
    set(gca,'Color',[0.5,0.5,0.5])
else
    if isfield(options, 'clims')
        imagesc(times, freqs, ERSP, options.clims);
    else
        imagesc(times, freqs, ERSP);
    end
end

% Create a 'grid' to make reading frequencies easier
for f = principalFreqValues
    yline(f, '--');
end
ylim([freqs(1), freqs(end)])

if options.singlePlot
    if options.plotEvents
        % Add events information
        for ev = 1:numel(events.time_ms)
            xl = xline(events.time_ms(ev), 'label', events.labels{ev});
            xl.LineWidth = 1.5;
            xl.FontSize = 12;
            xl.LabelVerticalAlignment = 'top';
            xl.LabelHorizontalAlignment = 'center';
        end
    end
    
    if options.timewarped
        xticks(events.time_ms);
        xticklabels([]);
        xlim([-1000, times(end)])
    else
        % No events, normal time data (no timewarping)
        xlim([times(1), times(end)])
        xlabel(options.XLabel, 'FontSize', 12)
    end
    
    % For the tick labels size
    ax = gca;
    ax.FontSize = 12;
    ylabel('Frequency (Hz)', 'FontSize', 12)
    title(options.title, 'FontSize', 14)
    
    %% Custom color map
    c = colorbar;
    if isfield(options, 'clims')
        clims = options.clims;
    else
        clims = [min(ERSP(:)), max(ERSP(:))];
    end
    N_colors_standard = 512;
    myCmap = asymColorMapWhiteZero(clims, N_colors_standard);
    set(gca,'Colormap',myCmap)
    c.Label.String = options.cLegend;
    c.Label.FontSize = 12;
else
    if options.plotEvents
        if isfield(options, 'nameEvents') && options.nameEvents && ~options.timewarped
            % Add event lines with labels
            for ev = 1:numel(events.time_ms)
                xl = xline(events.time_ms(ev), 'label', events.labels{ev});
                xl.LineWidth = 1.5;
                xl.FontSize = 12;
                xl.LabelVerticalAlignment = 'top';
                xl.LabelHorizontalAlignment = 'center';
            end
        else
            % Just add event lines
            for ev = 1:numel(events.time_ms)
                xl = xline(events.time_ms(ev));
                xl.LineWidth = 1.5;
            end
            if options.timewarped
                if isfield(options, 'nameEvents') && options.nameEvents
                    % Add labels at the bottom
                    xticks(events.time_ms);
                    xticklabels(events.labels);
                    xtickangle(15);
                else
                    xticks(events.time_ms);
                    xticklabels([]);
                end
            end
        end
    end
    xlim([times(1), times(end)])
    
    % For the tick labels size
    ax = gca;
    ax.FontSize = 11;
    if options.ylabel
        ylabel('Frequency (Hz)', 'FontSize', 12)
    end
    
    if isfield(options, 'XLabel')
        xlabel(options.XLabel, 'FontSize', 12)
    end
    title(options.title, 'FontSize', 12)
end
end