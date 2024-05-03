function EEG = importCustomElectrodes(EEG, study_config)
% Import electrodes' locations found with get_chanlocs from 3D scans
fprintf('Importing channel locations from custom file...\n')

subject = study_config.subjects(study_config.current_subject).id;
elecFolder = [study_config.study_folder study_config.electrodes_folder subject filesep];
file2load = [elecFolder, subject, study_config.indiv_channel_locations_filename];
if isfile(file2load)
    chanlocs = readlocs(file2load,'format',{'labels','X','Y','Z'});
    
    for ch = 1:numel(chanlocs)-3
        line = strcmp({EEG.chanlocs.labels},chanlocs(ch).labels);
        EEG.chanlocs(line).X = chanlocs(ch).X;
        EEG.chanlocs(line).Y = chanlocs(ch).Y;
        EEG.chanlocs(line).Z = chanlocs(ch).Z;
    end
    
    % Put fiducials in chaninfo
    EEG.chaninfo.nodatchans = chanlocs(end-2:end);
    
    % Recompute data in all coordinate frames
    EEG.chanlocs = rmfield(EEG.chanlocs, {'sph_phi','sph_radius','sph_theta','theta','radius'});
    EEG = eeg_checkchanlocs(EEG);
        
    %EEG.chanlocs = readlocs(file2load,'format',{'labels','X','Y','Z'});
    
    for ch = 1:numel(EEG.chanlocs)
        EEG.chanlocs(ch).urchan = ch;
%         if sum(contains(study_config.eog_channels, EEG.chanlocs(ch).labels))>0
%             EEG.chanlocs(ch).type = 'EOG';
%         else
%             EEG.chanlocs(ch).type = 'EEG';
%         end
    end
else
    warning('Could not find the individual electrodes'' location file: skipping import')
end

if study_config.moveElecInwards ~= 0
    EEG = moveElecInwardsEEG(EEG, study_config.moveElecInwards);
    EEG = eeg_checkchanlocs(EEG);
end

figure;
topoplot(ones(1,EEG.nbchan), EEG.chanlocs,...
    'style', 'blank', 'electrodes', 'ptslabels', 'emarker', {'.','r',10,1}, 'plotdisk', 'off');
saveCurrentFig([study_config.figures_folder 'ChannelsImport' filesep],...
    [subject, '_channelsLocation'], {'png'}, [1000,700]);
end