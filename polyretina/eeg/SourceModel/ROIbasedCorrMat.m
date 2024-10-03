%clear all;
close all;
configEEGAFF_Ainhoa;
addpath(genpath('C:\Users\Ainhoa\Documents\M2\MATLAB\codeM2\mne-matlab'))
overwrite = true;

hemi = study_config.recon.hemispheres;
rois = study_config.recon.rois;
res = study_config.recon.resolution;
for subject_ind = subject_inds
    %     if ~exist('ALLEEG','var')
    %         launchEEGLAB;
    %     end
    
    %subject_ind = 22;
    subject = study_config.subjects(subject_ind).id;
    disp(['Subject ' subject]);
    study_config.current_subject = subject_ind;
    
    %% Folders and Files names:
    N = makeFolderFileNames(study_config, subject);
    
    if ~isfile(fullfile(N.roiFolder, N.roiSumFile)) || ~isfile(fullfile(N.bemFolder, N.covFile))...
            || ~isfile(fullfile(N.bemFolder, N.covCholFile)) || overwrite
        
        %% MNE outputs: bem, src and fwd
        trans = fiff_read_coord_trans(fullfile(N.bemFolder, N.transFile));
        src = mne_read_source_spaces(fullfile(N.bemFolder, N.srcFile),...
            study_config.recon.fixed_ori, study_config.recon.patch_space);
        fwd = mne_read_forward_solution(fullfile(N.bemFolder, N.fwdFile),...
            study_config.recon.fixed_ori, study_config.recon.patch_space);
        
        %% Find dipoles associated to each roi
        roi_summary = struct();
        for h = 1:numel(hemi)
            src_trans = load(fullfile(N.bemFolder,sprintf('src%s2mri.mat', hemi{h})));
            
            src_t = mne_transform_source_space_to(fwd.src(h), 5, trans);
            %src_t = src(h);
            dipoles = src_t.rr;
            %dipoles = (trans.trans(1:3,1:3) * src_t.rr')';
            %dipoles = (trans.trans(1:3,:) * [src_t.rr ones(src_t.np,1)]')';
            used_inds = src_t.vertno;
            used_dips_mm = dipoles(used_inds,:)*1000 + src_trans.trans(:,4)';
            
            src_plot = src(h); % Already in MRI space
            %dipoles_plot = src_plot.rr;
            dipoles_plot = src_plot.rr + src_trans.trans(:,4)'./1000;
            verts = src_plot.vertno;
            % Compute used triangles with accurate reference to used dipoles
            orig_tris = src_plot.use_tris;
            used_tris = zeros(size(orig_tris),'int32');
            for v = 1:length(verts)
                used_tris(orig_tris == verts(v)) = v;
            end
            
            figure;
            hold on
            p = patch('Faces', used_tris, 'Vertices', dipoles_plot(verts,:)*1000);
            p.FaceColor = [.5 .5 .5];
            p.EdgeColor = 'black';
            view(3); axis image
            xlabel('X (mm)', 'Fontsize', 12);
            ylabel('Y (mm)', 'Fontsize', 12);
            zlabel('Z (mm)', 'Fontsize', 12);
            title(sprintf('%s - Hemisphere %s', subject, hemi{h}), 'FontSize', 16)
            
            for r = 1:numel(rois)
                roi_file = sprintf('%s_%s_res%dmm_%s.csv',...
                    hemi{h}, lower(rois{r}), res, study_config.recon.roi_creation);
                roi_voxels = table2array(readtable(fullfile(N.roiFolder, roi_file)));
                
                % Plot the ROI on the figure to check for misalignement of reference frames
                switch r
                    case 1
                        s1 = scatter3(roi_voxels(:,1), roi_voxels(:,2), roi_voxels(:,3), '+r');
                    case 2
                        s2 = scatter3(roi_voxels(:,1), roi_voxels(:,2), roi_voxels(:,3), '+b');
                    case 3
                        s3 = scatter3(roi_voxels(:,1), roi_voxels(:,2), roi_voxels(:,3), '+g');
                end
                
                fprintf('Finding source dipoles related to %s in %s (%dmm resolution)\n',...
                    rois{r}, hemi{h}, res)
                voxAssign = nan(size(roi_voxels,1),1);
                voxDist = nan(size(roi_voxels,1),1);
                for v = 1:size(roi_voxels,1)
                    all_diffs = sqrt(sum((used_dips_mm +...
                        repmat(study_config.recon.res_vect - roi_voxels(v,:), size(used_dips_mm,1),1)).^2,2));
                    [min_dist, dip_ind] = min(all_diffs);
                    voxAssign(v) = dip_ind;
                    voxDist(v) = min_dist;
                end
                
                if ~isempty(study_config.recon.dist_limit_mm)
                    fprintf('Using %d mm maximal distance\n', study_config.recon.dist_limit_mm)
                    allowed_dips = voxDist < study_config.recon.dist_limit_mm;
                    roi_dips = unique(voxAssign(allowed_dips));
                    meanDist = mean(voxDist(allowed_dips));
                    stdDist = std(voxDist(allowed_dips));
                else
                    roi_dips = unique(voxAssign);
                    meanDist = mean(voxDist);
                    stdDist = std(voxDist);
                end
                
                roi_summary.(sprintf('%s%s',rois{r},hemi{h})) = struct('voxAssign', voxAssign,...
                    'voxDist', voxDist, 'dips', roi_dips, 'meanDist', meanDist, 'stdDist', stdDist);
                
                fprintf('%d dipoles selected\n',length(roi_dips));
                fprintf('Mean projection distance: %.1f +/- %.1f mm\n', meanDist, stdDist);
            end
            
            % Present the dipoles removed when computing the fwd model, if any
            rem_dips = setdiff(verts,used_inds);
            if isempty(rem_dips)
                legend([s1,s2,s3], study_config.recon.rois);
            else                            
                rem = scatter3(dipoles_plot(rem_dips,1)*1000, dipoles_plot(rem_dips,2)*1000, dipoles_plot(rem_dips,3)*1000, 'w', 'filled');
                legend([rem,s1,s2,s3], [{sprintf('%d dips removed', length(rem_dips))},study_config.recon.rois]);
            end
        end
        close all
        
        cov = zeros(fwd.nsource);
        cov_chol = zeros(fwd.nsource);
        count = 1;
        for h = 1:numel(hemi)
            reorg_dips_h = [];
            
            %% Check if some dipoles are included in more than 1 ROI
            disp('Checking for intersections between ROIs...')
            for r1 = 1:(numel(rois)-1)
                r1field = sprintf('%s%s',rois{r1},hemi{h});
                r1_dips = roi_summary.(r1field).dips;
                for r2 = (1+r1):numel(rois)
                    r2field = sprintf('%s%s',rois{r2},hemi{h});
                    r2_dips = roi_summary.(r2field).dips;
                    
                    shared = intersect(r1_dips, r2_dips);
                    removed_r1 = 0; removed_r2 = 0;
                    if ~isempty(shared)
                        fprintf('Dipoles shared between %s and %s\n', r1field, r2field);
                        for sh = 1:length(shared)
                            instances_r1 = roi_summary.(r1field).voxAssign == shared(sh);
                            minDist_r1 = min(roi_summary.(r1field).voxDist(instances_r1));
                            instances_r2 = roi_summary.(r2field).voxAssign == shared(sh);
                            minDist_r2 = min(roi_summary.(r2field).voxDist(instances_r2));
                            if minDist_r1 < minDist_r2
                                r2_dips = setdiff(r2_dips, shared(sh));
                                removed_r2 = removed_r2+1;
                            elseif minDist_r1 > minDist_r2
                                r1_dips = setdiff(r1_dips, shared(sh));
                                removed_r1 = removed_r1+1;
                            else
                                % Distances are equal: Remove from both
                                r1_dips = setdiff(r1_dips, shared(sh));
                                removed_r1 = removed_r1+1;
                                r2_dips = setdiff(r2_dips, shared(sh));
                                removed_r2 = removed_r2+1;
                            end
                        end
                        fprintf('Removed %d dipoles from %s and %d from %s\n', removed_r1, r1field,...
                            removed_r2, r2field);
                        
                        % Last check
                        shared = intersect(r1_dips, r2_dips);
                        if ~isempty(shared)
                            error('Correction failed');
                        end
                        
                        roi_summary.(r1field).dips = r1_dips;
                        roi_summary.(r2field).dips = r2_dips;
                    end
                end
            end
            disp('Done')
            
            %% Compute distances between dipoles in the ROI
            src_trans = load(fullfile(N.bemFolder,sprintf('src%s2mri.mat', hemi{h})));
            src_t = mne_transform_source_space_to(fwd.src(h), 5, trans);
            dipoles = src_t.rr;
            triangles = src_t.use_tris;
            used_inds = src_t.vertno;
            used_dips_mm = dipoles(used_inds,:)*1000 + src_trans.trans(:,4)';
            for r = 1:numel(rois)
                rfield = sprintf('%s%s',rois{r},hemi{h});
                roi_dips = roi_summary.(rfield).dips;
                if ~isempty(roi_dips)
                    fprintf('Finding neighbours for %s...\n', rfield);
                    
                    % First order neighbours
                    neighs_1st = [];
                    neigh_found = false(length(roi_dips),1);
                    for d1 = 1:(length(roi_dips)-1)
                        vert1 = used_inds(roi_dips(d1));
                        search1 = triangles == vert1;
                        
                        for d2 = (1+d1):length(roi_dips)
                            vert2 = used_inds(roi_dips(d2));
                            search2 = triangles == vert2;
                            
                            if any(any(search1,2) & any(search2,2))
                                neighs_1st = [neighs_1st;[roi_dips(d1),roi_dips(d2)]];
                                neigh_found(d1) = true;
                                neigh_found(d2) = true;
                            end
                        end
                    end
                    fprintf('Removing %d dipoles with no first order neighbour\n',sum(~neigh_found))
                    roi_summary.(rfield).dips_with_neigh = roi_dips(neigh_found);
                    roi_summary.(rfield).neighs_1st = neighs_1st;
                    
                    % Second order neighbours
                    neighs_2nd = [];
                    roi_dips = roi_summary.(rfield).dips_with_neigh';
                    for d = roi_dips
                        % Get all first order neighbours
                        neighbors = [neighs_1st(neighs_1st(:,1) == d,2);neighs_1st(neighs_1st(:,2) == d,1)]';
                        
                        % Start from the neighbors to propagate
                        for n = neighbors
                            candidates = [neighs_1st(neighs_1st(:,1) == n,2);neighs_1st(neighs_1st(:,2) == n,1)]';
                            for c = candidates
                                if ~any(neighbors == c) && d<c &&...
                                        (isempty(neighs_2nd) || ~any(neighs_2nd(:,1)==d & neighs_2nd(:,2)==c))
                                    neighs_2nd = [neighs_2nd;[d,c]];
                                end
                            end
                        end
                    end
                    roi_summary.(rfield).neighs_2nd = neighs_2nd;                   
                    
                    %%%%% Wrong: should be the geodesic distance, not the
                    %%%%% euclidian distance !!!
                    %fprintf('Computing inter-ROI distances for %s...', rfield);
                    %distances = zeros(length(roi_dips));
                    %for d1 = 1:(length(roi_dips)-1)
                    %   vert1 = used_dips_mm(roi_dips(d1),:);
                    %   for d2 = (1+d1):length(roi_dips)
                    %       vert2 = used_dips_mm(roi_dips(d2),:);
                    %       distances(d1,d2) = sqrt(sum((vert1-vert2).^2));
                    %       distances(d2,d1) = distances(d1,d2);
                    %   end
                    %end
                    %roi_summary.(rfield).interDist = distances;
                    
                    %% Compute influence matrix
                    % 1. As the inverse of distance
                    % Change diagonal elements to ones for inversion
                    %distances(1:size(distances,1)+1:end) = 1;
                    %influences = 1./ceil(distances+1);
                    % Remove smallest numbers
                    %influences(influences < study_config.recon.infl_th) = 0;
                    %if any(sum(influences,1)==1)
                    %   warning('Some dipoles within the ROI are isolated')
                    %end
                    
                    % 2. With neighbourhood definition
                    influences = eye(length(roi_dips));
                    for n1 = 1:size(neighs_1st,1)
                        d11 = roi_dips == neighs_1st(n1,1);
                        d12 = roi_dips == neighs_1st(n1,2);
                        influences(d11,d12) = study_config.recon.infl_levels(1);
                        influences(d12,d11) = study_config.recon.infl_levels(1);
                    end
                    
                    for n2 = 1:size(neighs_2nd,1)
                        d21 = roi_dips == neighs_2nd(n2,1);
                        d22 = roi_dips == neighs_2nd(n2,2);
                        influences(d21,d22) = study_config.recon.infl_levels(2);
                        influences(d22,d21) = study_config.recon.infl_levels(2);
                    end
                    
                    disp('Done')
                    
                    %% Compute Cholesky decomposition of covariance matrix
                    disp('Computing cholesky decomposition...');
                    n_dips = length(roi_dips);
                    if all(eig(influences)>0)
                        fprintf('%s influence matrix is positive semi-definite\n', rfield)
                        infl_cho = chol(influences)';
                        
                        cov(count:count+n_dips-1,count:count+n_dips-1) = influences;
                        roi_summary.(rfield).interInfl = influences;
                    else
                        fprintf('%s influence matrix is not positive semi-definite\n', rfield)
                        disp('Replacing negative eigenvalues with slightly positive ones')
                        [V,D] = eig(influences);
                        D(D<0) = study_config.recon.slight_pos_val;
                        influences_pos = V*D/V;
                        infl_cho = chol(influences_pos)';
                        
                        cov(count:count+n_dips-1,count:count+n_dips-1) = influences_pos;
                        roi_summary.(rfield).interInfl = influences_pos;
                    end                    
                    
                    cov_chol(count:count+n_dips-1,count:count+n_dips-1) = infl_cho;
                    reorg_dips_h = [reorg_dips_h,roi_dips];
                    count = count + n_dips;
                else
                    fprintf('Skipping %s (no dipoles included)\n', rfield);
                    %roi_summary.(rfield).interDist = [];
                    roi_summary.(rfield).interInfl = [];
                end
            end
            
            remaining_dips = setdiff(1:src_t.nuse, reorg_dips_h);
            reorg_dips_h = [reorg_dips_h, remaining_dips];
            
            if h == 1
                next = count;
                last = src_t.nuse;
                offset = 0;
                count = src_t.nuse + 1;
            else
                next = count;
                last = fwd.nsource;
                offset = fwd.src(1).nuse;
            end
            
            for d = next:last
                cov(d,d) = 1;
                cov_chol(d,d) = 1;
            end
            
            % Reorganize matices
            new_cov = cov;
            new_chol = cov_chol;
            old_order = 1:length(reorg_dips_h);
            same_order = reorg_dips_h == old_order;
            for i = old_order
                if same_order(i)
                    for j = old_order(~same_order)
                        new_cov(reorg_dips_h(i)+offset, reorg_dips_h(j)+offset) = cov(i+offset, j+offset);
                        new_chol(reorg_dips_h(i)+offset, reorg_dips_h(j)+offset) = cov_chol(i+offset, j+offset);
                    end
                else
                    for j = old_order
                        new_cov(reorg_dips_h(i)+offset, reorg_dips_h(j)+offset) = cov(i+offset, j+offset);
                        new_chol(reorg_dips_h(i)+offset, reorg_dips_h(j)+offset) = cov_chol(i+offset, j+offset);
                    end
                end
            end
            cov = new_cov;
            cov_chol = new_chol;
            clear new_cov new_chol
        end
        
        save(fullfile(N.roiFolder, N.roiSumFile), 'roi_summary');
        save(fullfile(N.bemFolder, N.covFile), 'cov', '-v7.3');
        save(fullfile(N.bemFolder, N.covCholFile), 'cov_chol', '-v7.3');
    end
end