function plotERProi(stc_time_rois,roi,style,ylims,n_samp_smooth)
times = stc_time_rois.time;
field = sprintf('%s_pow',roi);

% Actual figure (no figure creation to allow subplotting)
hold on
pow_roiAve = squeeze(mean(stc_time_rois.(field),1));
switch style
    case 'GrandAverage'
        pow_roiAve_trMean = squeeze(movmean(mean(pow_roiAve,1),n_samp_smooth));
        pow_roiAve_trSD = squeeze(movmean(std(pow_roiAve,[],1),n_samp_smooth));
        line = plot(times, pow_roiAve_trMean, 'LineWidth',2);
        patch([times,fliplr(times)],[pow_roiAve_trMean-pow_roiAve_trSD,...
            fliplr(pow_roiAve_trMean+pow_roiAve_trSD)],...
            line.Color, 'EdgeColor', 'none', 'FaceAlpha',0.3);
        yline(0,'--k');
        xline(0,'--k', 'Label', 'Intersection');
        
    case 'Phase-wise'
        encoding_trials = contains(stc_time_rois.trialinfo.Phase,'Encoding');
        pow_roiAve_encodingMean = squeeze(movmean(mean(pow_roiAve(encoding_trials,:),1),n_samp_smooth));
        pow_roiAve_encodingSD = squeeze(movmean(std(pow_roiAve(encoding_trials,:),[],1),n_samp_smooth));
        pow_roiAve_testMean = squeeze(movmean(mean(pow_roiAve(~encoding_trials,:),1),n_samp_smooth));
        pow_roiAve_testSD = squeeze(movmean(std(pow_roiAve(~encoding_trials,:),[],1),n_samp_smooth));
        
        line1 = plot(times, pow_roiAve_encodingMean, 'LineWidth',2);
        line2 = plot(times, pow_roiAve_testMean, 'LineWidth',2);
        patch([times,fliplr(times)],[pow_roiAve_encodingMean-pow_roiAve_encodingSD,...
            fliplr(pow_roiAve_encodingMean+pow_roiAve_encodingSD)],...
            line1.Color, 'EdgeColor', 'none', 'FaceAlpha',0.3)
        patch([times,fliplr(times)],[pow_roiAve_testMean-pow_roiAve_testSD,...
            fliplr(pow_roiAve_testMean+pow_roiAve_testSD)],...
            line2.Color, 'EdgeColor', 'none', 'FaceAlpha',0.3)
        
        yline(0,'--k');
        xline(0,'--k', 'Label', 'Intersection');
        
        legend([line1,line2],{sprintf('Encoding (N=%d)',sum(encoding_trials)),...
            sprintf('Test (N=%d)',sum(~encoding_trials))})
        
    case 'Condition-wise'
        down_trials = contains(stc_time_rois.trialinfo.Condition,'DOWN');
        up_trials = contains(stc_time_rois.trialinfo.Condition,'UP');
        mix_trials = contains(stc_time_rois.trialinfo.Condition,'MIX');
        
        pow_roiAve_downMean = squeeze(movmean(mean(pow_roiAve(down_trials,:),1),n_samp_smooth));
        pow_roiAve_downSD = squeeze(movmean(std(pow_roiAve(down_trials,:),[],1),n_samp_smooth));
        pow_roiAve_upMean = squeeze(movmean(mean(pow_roiAve(up_trials,:),1),n_samp_smooth));
        pow_roiAve_upSD = squeeze(movmean(std(pow_roiAve(up_trials,:),[],1),n_samp_smooth));
        pow_roiAve_mixMean = squeeze(movmean(mean(pow_roiAve(mix_trials,:),1),n_samp_smooth));
        pow_roiAve_mixSD = squeeze(movmean(std(pow_roiAve(mix_trials,:),[],1),n_samp_smooth));
        
        line1 = plot(times, pow_roiAve_downMean, 'LineWidth',2);
        line2 = plot(times, pow_roiAve_upMean, 'LineWidth',2);
        line3 = plot(times, pow_roiAve_mixMean, 'LineWidth',2);
        patch([times,fliplr(times)],[pow_roiAve_downMean-pow_roiAve_downSD,...
            fliplr(pow_roiAve_downMean+pow_roiAve_downSD)],...
            line1.Color, 'EdgeColor', 'none', 'FaceAlpha',0.3)
        patch([times,fliplr(times)],[pow_roiAve_upMean-pow_roiAve_upSD,...
            fliplr(pow_roiAve_upMean+pow_roiAve_upSD)],...
            line2.Color, 'EdgeColor', 'none', 'FaceAlpha',0.3)
        patch([times,fliplr(times)],[pow_roiAve_mixMean-pow_roiAve_mixSD,...
            fliplr(pow_roiAve_mixMean+pow_roiAve_mixSD)],...
            line3.Color, 'EdgeColor', 'none', 'FaceAlpha',0.3)
        
        yline(0,'--k');
        xline(0,'--k', 'Label', 'Intersection');
        
        legend([line1,line2,line3],{sprintf('Down (N=%d)',sum(down_trials)),...
            sprintf('Up (N=%d)',sum(up_trials)), sprintf('Mix (N=%d)',sum(mix_trials))})
        
    otherwise
        error('Unknown style')
end

ylim(ylims);
title(roi);
end