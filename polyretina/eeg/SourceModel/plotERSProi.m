function plotERSProi(stc_freq_rois, roi, opes, opts)
times_ms = stc_freq_rois.time*1000;
freqs = stc_freq_rois.freq;
field = sprintf('%s_pow',roi);

principalFreqValues = [4,8,12,20,30];
%% options struct
plot_options.title = roi;
%options.mask;
plot_options.clims = opts.clims;
plot_options.singlePlot = opts.singlePlot; % Adds more information on the plot if single
plot_options.timewarped = false;
if opts.xlabel
    plot_options.XLabel = 'Time (ms)';
end
switch opts.type
    case 'Baseline'
        plot_options.plotEvents = false;
    case 'Observation'
        plot_options.plotEvents = true;
    otherwise
        error('Unknown type')
end
% Following options only if singlePlot == false
plot_options.ylabel = opts.ylabel;
plot_options.nameEvents = opts.nameEvents;


%% Data operations
if opes.singleTrialNorm
    % single trial normalization
    singletrialbase = mean(stc_freq_rois.(field),4);
    ERSP = stc_freq_rois.(field)./repmat(singletrialbase, [1,1,1,length(times_ms)]);
else
    ERSP = stc_freq_rois.(field);
end

switch opts.style
    case 'GrandAverage'
        trials = 1:size(stc_freq_rois.trialinfo,1);
    case 'Encoding'
        trials = contains(stc_freq_rois.trialinfo.Phase,'Encoding');
    case 'Test'
        trials = contains(stc_freq_rois.trialinfo.Phase,'Test');
    otherwise
        error('Unknown Style')
end
ERSP = ERSP(:,trials,:,:);

if opes.preStimBaseline
    preStimBase = mean(ERSP(:,:,:,times_ms<0),[2,4]);
    ERSP = ERSP./repmat(preStimBase, [1,size(ERSP,2),1,size(ERSP,4)]);
    plot_options.cLegend = 'Log power wrt prestimulus baseline (dB)';
elseif ~isempty(opes.baseline)
    % External baseline provided by the user
    siz = size(squeeze(opes.baseline));
    if siz == size(squeeze(ERSP(:,1,:,1)))
        ERSP = ERSP./repmat(opes.baseline, [1,size(ERSP,2),1,size(ERSP,4)]);
        plot_options.cLegend = 'Log power wrt external baseline (dB)';
    else
        error('This type of baseline is no supported')
    end
else
    plot_options.cLegend = 'Log power (dB)';
end

% Average over all dipoles
ERSP = squeeze(mean(ERSP,1));
% dB transformation
ERSP = 10*log10(ERSP);
% Average over all trials
ERSP = squeeze(mean(ERSP,1));

%% Actual figure (no figure creation to allow subplotting)
makeERSPplot(ERSP, times_ms, freqs, principalFreqValues, opts.events, plot_options)
end