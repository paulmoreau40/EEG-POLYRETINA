% Needed to start eeglab on linux
path= 'C:\Program Files\MATLAB\R2019a\toolbox\matlab\graph2d';
%rmpath(path)
addpath(genpath(path), '-begin')
% Should return 'C:\Program Files\MATLAB\R2019a\toolbox\matlab\graph2d\axis.p
which axis

% Commands:
eeglab;
pop_editoptions('option_storedisk', 1, 'option_savetwofiles', 1,...
    'option_saveversion6', 0, 'option_single', 0, 'option_memmapdata', 0,...
    'option_eegobject', 0, 'option_computeica', 1, 'option_scaleicarms', 1,...
    'option_rememberfolder', 1, 'option_donotusetoolboxes', 0, 'option_checkversion', 0);
%     pop_editoptions('option_storedisk', 0, 'option_savetwofiles', 1,...
%         'option_saveversion6', 0, 'option_single', 0, 'option_memmapdata', 0,...
%         'option_eegobject', 0, 'option_computeica', 1, 'option_scaleicarms', 1,...
%         'option_rememberfolder', 1, 'option_donotusetoolboxes', 0, 'option_checkversion', 1, 'option_chat', 1);

%path2load_xdf = 'C:\Users\Alexandre\Desktop\code\FourStreetsCity\xdf';
%addpath(genpath(path2load_xdf))

runMoBI = false;
if runMoBI
    runmobilab;
    
    % Correct path to avoid issues with the workspace (wrong strjoin function is used)
    path= 'C:\Program Files\MATLAB\R2019a\toolbox\matlab\strfun';
    %rmpath(path)
    addpath(genpath(path), '-begin')
    % Should return 'C:\Program Files\MATLAB\R2017a\toolbox\matlab\strfun\strjoin.m
    which strjoin
   
    % Remove mobilab cleanline function from path (interferes with the clean line plugin)
    %rmpath(['C:\Users\Alexandre\Desktop\code'...
    %    filesep 'eeglab14_1_0b' filesep 'plugins' filesep 'mobilab'...
    %    filesep 'dependency' filesep 'cleanline']);
end