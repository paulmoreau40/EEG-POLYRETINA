function Answers = ReadAnswers_EEGAFF(config)

filename = fullfile('/Users/ilaria/Documents/MATLAB/code/EEG-AFF-2', 'Distance_estimation.xlsx');

opts = spreadsheetImportOptions("NumVariables", 7);

% Specify sheet and range
opts.Sheet = config.subjects(config.current_subject).id;
opts.DataRange = "A3:G28";

% Specify column names and types
opts.VariableNames = ["Essai", "Block0", "Comments0", "Block1", "Comments1", "Block2", "Comments2"];
opts.SelectedVariableNames = ["Essai", "Block0", "Comments0", "Block1", "Comments1", "Block2", "Comments2"];
opts.VariableTypes = ["double", "double", "string", "double", "string", "double", "string"];
opts = setvaropts(opts, [3, 5, 7], "WhitespaceRule", "preserve");
opts = setvaropts(opts, [3, 5, 7], "EmptyFieldRule", "auto");

% Import the data
Answers = readtable(filename, opts, "UseExcel", false);

% filename = fullfile('/Users/ilaria/Documents/MATLAB/code/EEG-AFF-2', 'Distance_estimation.xlsx');
% Participants = table2struct(readtable(filename));
% 
% Answers = [];
% 
% for p=1:38
%     if p<10
%         tab = strcat("P00",num2str(p));
%         answer = table2struct(readtable(filename,'Sheet',tab));
%         Answers.(strcat("P00",num2str(p))) =  answer;
%         for i=2:length(Answers.(strcat("P00",num2str(p))))
%             Answers.(strcat("P00",num2str(p)))(i).(strcat("P00",num2str(p))) = str2double(Answers.(strcat("P00",num2str(p)))(i).(strcat("P00",num2str(p))));
%             Answers.(strcat("P00",num2str(p)))(i).Var2 = str2double(Answers.(strcat("P00",num2str(p)))(i).Var2);
%             Answers.(strcat("P00",num2str(p)))(i).Var4 = str2double(Answers.(strcat("P00",num2str(p)))(i).Var4);
%             Answers.(strcat("P00",num2str(p)))(i).Var6 = str2double(Answers.(strcat("P00",num2str(p)))(i).Var6);
%         end
%     else
%         tab = strcat("P0",num2str(p));
%         answer = table2struct(readtable(filename,'Sheet',tab));
%         Answers.(strcat("P0",num2str(p))) =  answer;        
%         for i=2:length(Answers.(strcat("P0",num2str(p))))
%             Answers.(strcat("P0",num2str(p)))(i).(strcat("P0",num2str(p))) = str2double(Answers.(strcat("P0",num2str(p)))(i).(strcat("P0",num2str(p))));
%             Answers.(strcat("P0",num2str(p)))(i).Var2 = str2double(Answers.(strcat("P0",num2str(p)))(i).Var2);
%             Answers.(strcat("P0",num2str(p)))(i).Var4 = str2double(Answers.(strcat("P0",num2str(p)))(i).Var4);
%             Answers.(strcat("P0",num2str(p)))(i).Var6 = str2double(Answers.(strcat("P0",num2str(p)))(i).Var6);
%         end
%     end
% end
% 
% %% Save the output structure
% subject = config.subjects(config.current_subject).id;
% dirname = [config.study_folder config.preprocessing_folder];
% fname = [subject '_ReadAnswers.mat'];
% save(fullfile(dirname, fname), 'Answers')
end