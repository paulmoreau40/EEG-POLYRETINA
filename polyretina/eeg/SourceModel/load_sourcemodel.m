function [sourcemodel, cov, cov_chol, dipsROI] = load_sourcemodel(names, chanLabels, params)
%% Inputs:
% names         - folder names structure
% chanLabels    - channel labels {EEG.chanlocs.labels}
% params        - struct containing the following fields:
%   chans2remove            - cell list of electrodes names to remove for leadfield (bad trust in localization)
%   fixed_ori               - boolean to specify whether the leadfield has fixed orientations or not
%   patch_space             - boolean to specify whether to use the patch space when loading the fwd structure
%   cholesky                - boolean to specify whether to use the cholesky decomposition or not
%   normalizeDepth          - boolean to specify whether to normalize with depth
%   normalizeDepthParam     - float between 0 and 1
%   weightWithPatchSize     - boolean to specify whether to weight the leadfield with patch size

chanSelect = true(numel(chanLabels),1);
if isempty(params.chans2remove)
    labels = chanLabels;
else
    for ch = 1:numel(params.chans2remove)
        chanSelect(strcmp(chanLabels, params.chans2remove{ch})) = false;
    end
    labels = chanLabels(chanSelect);
end

%trans = fiff_read_coord_trans(fullfile(names.bemFolder, names.transFile));
%bem = mne_read_bem_surfaces(fullfile(names.bemFolder, names.bemFile),1);
src = mne_read_source_spaces(fullfile(names.bemFolder, names.srcFile),1);
fwd = mne_read_forward_solution(fullfile(names.bemFolder, names.fwdFile), params.fixed_ori, params.patch_space);

%% Compute used triangle information
for s = 1:2
    included_fwd = true(src(s).nuse,1);
    for d = 1:length(included_fwd)
        if isempty(find(fwd.src(s).vertno == src(s).vertno(d),1))
            included_fwd(d) = false;
        end
    end
    positions = src(s).rr(logical(src(s).inuse),:);
    
    orig_tris = src(s).use_tris;
    new_tris = zeros(size(orig_tris),'int32');
    
    verts = src(s).vertno;
    for v = 1:length(verts)
        new_tris(orig_tris == verts(v)) = v;
    end
    
    if s == 1
        included = included_fwd;
        used_pos = positions;
        used_tris = new_tris;
    else
        included = [included; included_fwd];
        used_pos = [used_pos; positions];
        used_tris = [used_tris; new_tris + repmat(src(1).nuse,size(new_tris))];
    end
end

%% Configuration of Source Model with precomputed leadfield
sourcemodel = struct();
sourcemodel.unit = 'mm';
sourcemodel.pos = used_pos;
sourcemodel.tri = used_tris;
sourcemodel.label = labels;
sourcemodel.leadfielddimord = '{pos}_chan_ori';
sourcemodel.inside = included;
conversion = find(included);
if params.fixed_ori
    fullLeadfield = fwd.sol.data;
    reducedLeadfield = fullLeadfield(chanSelect,:);
    
    if params.normalizeDepth
        for d=1:size(reducedLeadfield,2)
            tmplf = reducedLeadfield(:,d);
            % normalize the leadfield by sum of squares of the elements of the leadfield matrix to the power "normalizeparam"
            % this is the same as the Frobenius norm if normalizeDepthParam is 0.5
            nrm = sum(tmplf(:).^2)^params.normalizeDepthParam;
            if nrm>0
                tmplf = tmplf./nrm;
            end
            reducedLeadfield(:,d) = tmplf;
        end
    end
    
    if params.weightWithPatchSize
        error('Not implemented yet')
        for d=1:size(reducedLeadfield,2)
            reducedLeadfield(:,d) = reducedLeadfield(:,d) * patchWeight(d);
        end
    end
    
    % Always return both matrices
    load(fullfile(names.bemFolder, names.covFile), 'cov');
    load(fullfile(names.bemFolder, names.covCholFile), 'cov_chol');
    if params.cholesky
        reducedLeadfield = reducedLeadfield*cov_chol;
    end
%     if params.cholesky
%         % Option 1: Return the cholesky decomposition of the modified covariance matrix
%         % To be multiplied with the leadfield before computing the inverse
%         load(fullfile(names.bemFolder, names.covCholFile), 'cov_chol');
%         cov = cov_chol;
%     else
%         % Option 2: Directly return the modified covariance matrix.
%         % Input for the inverse computation by ft_inverse_mne
%         load(fullfile(names.bemFolder, names.covFile), 'cov');
%     end
    
    % Modify the leadfield and covariance matrix to fit with the
    % 'inside' parameter
    modifLeadfield = zeros(numel(labels),length(included));
    cpt = 0;
    for i = 1:length(included)
        if included(i)
            cpt = cpt+1;
            modifLeadfield(:,i) = reducedLeadfield(:,cpt);
        end
    end
    
    sourcemodel.leadfield = mat2cell(modifLeadfield,numel(labels),ones(length(included),1));
else
    error('Not implemented')
end

%% Find dipoles in each ROI
load(fullfile(names.roiFolder, names.roiSumFile), 'roi_summary');
fields = fieldnames(roi_summary);
for f = 1:numel(fields)
    if contains(fields{f},'lh')
        dipsROI.(fields{f}) = conversion(roi_summary.(fields{f}).dips);
    else
        dipsROI.(fields{f}) = conversion(roi_summary.(fields{f}).dips + double(fwd.src(1).nuse));
    end
end
end