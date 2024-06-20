function [alpha, options] = alphaPairwiseComp(alpha_init, n_pairs, params, options)
    if params.apply_pairwiseCorr
        options.pairwiseCompsCorrected = 'yes';
        switch params.correctionType
            % factor 2 appears because it is two-tailed
            case 'bonf'
                alpha = alpha_init/(2*nchoosek(n_pairs,2)); %Bonferroni
            case 'dunn'
                alpha = (1-(1-alpha_init)^(1/nchoosek(n_pairs,2)))/2; %Dunn-Sidak
            otherwise
                error('Unknown correction type');
        end
    else
        options.pairwiseCompsCorrected = 'no';
        alpha = alpha_init;
    end
    options.alpha_pairwise = alpha;
end