folder = 'F:\EEG-AFF-2DR\figures\ERSPs\bemobil\autoMoBI\MultiSubject\';
list = dir(folder);
for i = 1:numel(list)
    if contains(list(i).name,'.fig')
        figName = list(i).name(1:end-4);
        if ~exist(fullfile(folder,[figName, '.png']), 'file')
            openfig(fullfile(folder,list(i).name));
            saveCurrentFig(folder, figName,{'png'},[1000,800]);
        end
    end
end