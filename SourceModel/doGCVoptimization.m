function opt_lambda = doGCVoptimization(G, data, domain, P, plotParams)

[U,s,~] = csvd(G);
trials = 1:size(data.trialinfo,1);
chans = 1:numel(data.label);
times = find(0<=data.time);

[~,~,~,regparams] = gcv(U,s,ones(numel(chans),1),'Tikh',true);

switch domain
    case 'time'
        dat = permute(data.trial(:,:,times),[2,1,3]);
        lambdas_optim = nan(length(trials), length(times));
        G_optim = nan(length(trials), length(times));
        
        ppm = ParforProgressbar(length(trials), 'showWorkerProgress', true,...
            'title', sprintf('%s - Computing Generalized Cross Validation to optimize lambda', plotParams.subject));
        % One sample at a time
        parfor tr = trials
            %fprintf('Trial %d: ',tr)
            %tic;
            lambdas_optim_tr = nan(1, length(times));
            G_optim_tr = nan(1, length(times));
            for ti = 1:length(times)
                M = squeeze(dat(chans,tr,ti));
                if ~isempty(P)
                    [lambdas_optim_tr(1,ti),G_optim_tr(1,ti),~,~] = gcv(U,s,P*M,'Tikh',true);
                else
                    [lambdas_optim_tr(1,ti),G_optim_tr(1,ti),~,~] = gcv(U,s,M,'Tikh',true);
                end
            end
            lambdas_optim(tr,:) = lambdas_optim_tr;
            G_optim(tr,:) = G_optim_tr;
            %elapsed = toc;
            %fprintf('%.2f minutes\n',elapsed/60)
            ppm.increment();
        end
        delete(ppm);
        
    case 'freq'
        freqs = 1:length(data.freq);
        dat = permute(data.fourierspctrm(:,:,:,times), [2 1 3 4]);
        lambdas_optim = nan(length(trials), length(freqs), length(times));
        G_optim = nan(length(trials), length(freqs), length(times));
        
        ppm = ParforProgressbar(length(trials), 'showWorkerProgress', true,...
            'title', sprintf('%s - Computing Generalized Cross Validation to optimize lambda', plotParams.subject));
        % One sample at a time
        parfor tr = trials
            %fprintf('Trial %d: ',tr)
            %tic;
            lambdas_optim_tr = nan(1, length(freqs), length(times));
            G_optim_tr = nan(1, length(freqs), length(times));
            for fr = freqs
                for ti = 1:length(times)
                    M = squeeze(dat(chans,tr,fr,ti));
                    if ~isempty(P)
                        [lambdas_optim_tr(1,fr,ti),G_optim_tr(1,fr,ti),~,~] = gcv(U,s,P*M,'Tikh',true);
                    else
                        [lambdas_optim_tr(1,fr,ti),G_optim_tr(1,fr,ti),~,~] = gcv(U,s,M,'Tikh',true);
                    end
                end
            end
            lambdas_optim(tr,:,:) = lambdas_optim_tr;
            G_optim(tr,:,:) = G_optim_tr;
            %elapsed = toc;
            %fprintf('%.2f minutes\n',elapsed/60)
            ppm.increment();
        end
        delete(ppm);
    otherwise
        error('Unknown domain')
end

perc_highbound = 100*sum(lambdas_optim > regparams(2),'all')/numel(lambdas_optim);
perc_lowbound = 100*sum(lambdas_optim < regparams(end-1),'all')/numel(lambdas_optim);
fprintf('%.1f%% of optimal lambda values ceiled to high search boundary\n',perc_highbound);
fprintf('%.1f%% of optimal lambda values floored to low search boundary\n',perc_lowbound);
if perc_highbound + perc_lowbound > 5
    warning('More than 5%% of optimal lambda values were restricted by search boundaries')
end

if plotParams.do
    %% Figures
    [edges, opt_lambda] = plotLambdaGdistrib();
    plotOptLambdaGA(opt_lambda, edges);
    if strcmp(domain,'freq')
        plotOptLambdaAvePerFactor(edges);
    end
else
    opt_lambda = mean(lambdas_optim(:));
end


    function [edges, opt_lambda]= plotLambdaGdistrib()
        figure;
        suptitle(sprintf('%s',plotParams.subject));
        subplot(1,2,1)
        hold on;
        edges = logspace(floor(min(log10([lambdas_optim(:);regparams(end)]))),...
            ceil(max(log10([lambdas_optim(:);regparams(1)]))),100);
        histogram(lambdas_optim(:),edges)
        opt_lambda = mean(lambdas_optim(:));
        [~,opt_lambda_ind] = min(abs(lambdas_optim(:)-opt_lambda));
        xline(opt_lambda,'-.r', 'Label', sprintf('Mean = %.2e', opt_lambda),...
            'LabelHorizontalAlignment','center','LabelOrientation', 'horizontal');
        xline(regparams(end),'-.g','Label', {'Lower search bound',sprintf('(%.1f%%)',perc_lowbound)},...
            'LabelVerticalAlignment','middle','LabelOrientation', 'horizontal');
        xline(regparams(1),'-.g','Label', {'Upper search bound',sprintf('(%.1f%%)',perc_highbound)},...
            'LabelVerticalAlignment','middle','LabelHorizontalAlignment','left','LabelOrientation', 'horizontal');
        xlabel('Optimal lambda', 'Fontsize', 12);
        %xlim([0,prctile(lambdas_optim(:),99)])
        set(gca,'XScale','log')
        ylabel('Count', 'Fontsize', 12);
        switch domain
            case 'time'
                title({'Distribution of optimal lambdas',...
                    sprintf('1 optimum per trial x time point (%d in total)',numel(lambdas_optim))},...
                    'Fontsize', 13);
            case 'freq'
                title({'Distribution of optimal lambdas',...
                    sprintf('1 optimum per trial x freq x time point (%d in total)',numel(lambdas_optim))},...
                    'Fontsize', 13);
        end
        subplot(1,2,2)
        hold on;
        edges_G = logspace(floor(min(log10(G_optim(:)))),ceil(max(log10(G_optim(:)))),100);
        histogram(G_optim(:),edges_G)
        opt_G = G_optim(opt_lambda_ind);
        xline(opt_G,'-.r', 'Label', sprintf('%.2e', opt_G),...
            'LabelOrientation', 'horizontal');
        xlabel('G values', 'Fontsize', 12);
        %xlim([0,prctile(G_optim(:),99)])
        set(gca,'XScale','log')
        ylabel('Count', 'Fontsize', 12);
        switch domain
            case 'time'
                title({'Distribution of G values for optimal lambdas',...
                    sprintf('1 optimum per trial x time point (%d in total)',numel(lambdas_optim))},...
                    'Fontsize', 13);
            case 'freq'
                title({'Distribution of G values for optimal lambdas',...
                    sprintf('1 optimum per trial x freq x time point (%d in total)',numel(lambdas_optim))},...
                    'Fontsize', 13);
        end
        saveCurrentFig(plotParams.saveFolder,...
            sprintf('%s_lambda-G_distributions_%sopt', plotParams.subject, domain), {'png'}, [1600,800]);
        
        %     figure;
        %     scatter(lambdas_optim(:), G_optim(:))
        %     xlabel('Optimal lambda');
        %     set(gca,'XScale','log')
        %     ylabel('G values');
        %     set(gca,'YScale','log')
    end

    function plotOptLambdaGA(opt_lambda, edges)
        figure;
        switch domain
            case 'time'
                figure;
                subplot(1,2,1)
                hold on
                plot(data.time(times), squeeze(mean(lambdas_optim,[1])));
                yline(opt_lambda,'-.r', 'Label', sprintf('%.2e', opt_lambda));
                xlabel('Time (s)');
                ylabel('Mean Lambda');
                ylim([10*edges(1),edges(end)/10])
                set(gca,'YScale','log')
                title({'Optimization across times', 'Grand average'});
                subplot(1,2,2)
                hold on
                plot(trials, squeeze(mean(lambdas_optim,[2])));
                yline(opt_lambda,'-.r', 'Label', sprintf('%.2e', opt_lambda));
                xlabel('Trial (count)');
                ylabel('Mean Lambda');
                ylim([10*edges(1),edges(end)/10])
                set(gca,'YScale','log')
                title({'Optimization across trials', 'Grand average'});
                
            case 'freq'
                figure;
                subplot(1,3,1)
                hold on
                plot(data.time(times), squeeze(mean(lambdas_optim,[1,2])));
                yline(opt_lambda,'-.r', 'Label', sprintf('%.2e', opt_lambda));
                xlabel('Time (s)');
                ylabel('Mean Lambda');
                ylim([10*edges(1),edges(end)/10])
                set(gca,'YScale','log')
                title({'Optimization across times', 'Grand average'});
                subplot(1,3,2)
                hold on
                plot(data.freq(freqs), squeeze(mean(lambdas_optim,[1,3])));
                yline(opt_lambda,'-.r', 'Label', sprintf('%.2e', opt_lambda));
                xlabel('Freq (Hz)');
                ylabel('Mean Lambda');
                ylim([10*edges(1),edges(end)/10])
                set(gca,'YScale','log')
                title({'Optimization across frequencies', 'Grand average'});
                subplot(1,3,3)
                hold on
                plot(trials, squeeze(mean(lambdas_optim,[2,3])));
                yline(opt_lambda,'-.r', 'Label', sprintf('%.2e', opt_lambda));
                xlabel('Trial (count)');
                ylabel('Mean Lambda');
                ylim([10*edges(1),edges(end)/10])
                set(gca,'YScale','log')
                title({'Optimization across trials', 'Grand average'});
        end
        
        saveCurrentFig(plotParams.saveFolder,...
            sprintf('%s_meanlambda_perfactor_%sopt', plotParams.subject, domain), {'png'}, [1600,800]);
    end

    function plotOptLambdaAvePerFactor(edges)
        figure;
        subplot(2,3,1)
        plot(data.time(times), squeeze(mean(lambdas_optim,[1])));
        xlabel('Time (s)');
        ylabel('Mean Lambda');
        ylim([10*edges(1),edges(end)/10])
        set(gca,'YScale','log')
        title({'Optimization across times', 'per frequency'});
        subplot(2,3,4)
        plot(data.time(times), squeeze(mean(lambdas_optim,[2])));
        xlabel('Time (s)');
        ylabel('Mean Lambda');
        ylim([10*edges(1),edges(end)/10])
        set(gca,'YScale','log')
        title({'Optimization across times', 'per trial'});
        subplot(2,3,2)
        plot(data.freq(freqs), squeeze(mean(lambdas_optim,[1])));
        xlabel('Freq (Hz)');
        ylabel('Mean Lambda');
        ylim([10*edges(1),edges(end)/10])
        set(gca,'YScale','log')
        title({'Optimization across frequencies', 'per time'});
        subplot(2,3,5)
        plot(data.freq(freqs), squeeze(mean(lambdas_optim,[3])));
        xlabel('Freq (Hz)');
        ylabel('Mean Lambda');
        ylim([10*edges(1),edges(end)/10])
        set(gca,'YScale','log')
        title({'Optimization across frequencies', 'per trial'});
        subplot(2,3,3)
        plot(trials, squeeze(mean(lambdas_optim,[2])));
        xlabel('Trial (count)');
        ylabel('Mean Lambda');
        ylim([10*edges(1),edges(end)/10])
        set(gca,'YScale','log')
        title({'Optimization across trials', 'per time'});
        subplot(2,3,6)
        plot(trials, squeeze(mean(lambdas_optim,[3])));
        xlabel('Trial (count)');
        ylabel('Mean Lambda');
        ylim([10*edges(1),edges(end)/10])
        set(gca,'YScale','log')
        title({'Optimization across trials', 'per frequency'});
        saveCurrentFig(plotParams.saveFolder,...
            sprintf('%s_meanlambda_perfactor_details_%sopt', plotParams.subject, domain), {'png'}, [1600,800]);
        
        %     figure;
        %     subplot(1,3,1)
        %     plot(data_freq.time(times), reshape(lambdas_optim,[size(lambdas_optim,1)*size(lambdas_optim,2), size(lambdas_optim,3)]));
        %     xlabel('Time (s)');
        %     ylabel('Lambda');
        %     set(gca,'YScale','log')
        %     title({'Evolution of lambda optimization','across consecutive time stamps'});
        %     axis tight
        %     subplot(1,3,2)
        %     plot(data_freq.freq(freqs), reshape(permute(lambdas_optim,[1,3,2]),[size(lambdas_optim,1)*size(lambdas_optim,3), size(lambdas_optim,2)]));
        %     xlabel('Freq (Hz)');
        %     ylabel('Lambda');
        %     set(gca,'YScale','log')
        %     title({'Evolution of lambda optimization','across consecutive freq values'});
        %     axis tight
        %     subplot(1,3,3)
        %     plot(trials, reshape(lambdas_optim,[size(lambdas_optim,1), size(lambdas_optim,2)*size(lambdas_optim,3)]));
        %     xlabel('Trial (count)');
        %     ylabel('Lambda');
        %     set(gca,'YScale','log')
        %     title({'Evolution of lambda optimization','across consecutive trials'});
        %     axis tight
    end
end