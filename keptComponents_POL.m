% One section per pipeline
% Name the struct with the pipeline name
% Use manual for ambiguous cases

%% Bemobil - autoMoBI pipeline
autoMoBI_KC = struct('ID', [], 'Brain', [], 'BrainWithNoise', [], 'NoiseWithBrain',[], 'Doubts',[]);

autoMoBI_KC(1) = struct('ID','P002',...
    'Brain',[2,6,13,16,27,43,49,50,59],...
    'BrainWithNoise',[7,10,11,12,17,35,37,45,48,51,53,55:58,61,63,65,69,70:73,75,78],...
    'NoiseWithBrain',[1,3,4,9,14,24,26,28,44,46,47,60,64,74,77,83,85,86,89,105],...
    'Doubts',[10,17,45,63,75,85,105]);
 
autoMoBI_KC(2) = struct('ID','P006',...
    'Brain',[1,2,3,5,6,9,15,16,23,36,37,54,55,58],...
    'BrainWithNoise',[10,12,22,24,35,43,46,57,65,68,71,75,79,81,82,85],...
    'NoiseWithBrain',[26,28,40,42,44,47,52,66,69,72,73,76,77,80,84,89,97,106,107],...
    'Doubts',[10,52,65,72,82]);

% autoMoBI_KC(2) = struct('ID','CIG03',...
%     'Brain',[18,27,33,43,52,53,59],...
%     'BrainWithNoise',[48,66:67,73],...
%     'NoiseWithBrain',[29,51,65,68:70,72,76:80,84],...
%     'Doubts',[48,66,68,73]);
% % Notes for CIG03:
% % Ambiguous cases resolved with manual
% 
% autoMoBI_KC(3) = struct('ID','CAN04',...
%     'Brain',[11,19,20,28,38,43,50,52,76,77,81,97],...
%     'BrainWithNoise',[59,62,66,96],...
%     'NoiseWithBrain',[21,27,33,61,72,85,89,92,98,101,102,105:106,108,110:113,115],...
%     'Doubts',[59,62,66,105,108]);
% % Notes for CAN04:
% % revoir 59,62,66,73,89,95,101,104,105,108
% % Ambiguous cases resolved with manual
% 
% autoMoBI_KC(4) = struct('ID','AUT05',...
%     'Brain',[16,23,30,34,37,39,40,42,48,51,58,72],...
%     'BrainWithNoise',[19,29,50,73:74],...
%     'NoiseWithBrain',[21,22,31,33,38,52,57,61,65,70,71,76,78,81:82,84,90,91,93,97,100,103,104],...
%     'Doubts',[29,38,70,73,81]);
% % Notes for AUT05:
% % Revoir 29,38,70,73,81
% % Ambiguous cases resolved with manual
% 
% autoMoBI_KC(5) = struct('ID','ANG06',...
%     'Brain',[5,7,9,10,12:14,18,19,21,24,27:28,33,50,79],...
%     'BrainWithNoise',[25,29,34,52,59,99],...
%     'NoiseWithBrain',[15,31,55,63,70,75,83,87,88,102,104,108:109,112,114:115],...
%     'Doubts',[29,34,55,88,99]);
% % Notes for ANG06: not enough IC (73)
% %
% % Ambiguous cases resolved with manual
% 
% autoMoBI_KC(6) = struct('ID','CRE07',...
%     'Brain',[9,13,20,21,22,24,27,30,40,42,63],...
%     'BrainWithNoise',[10,11,23,26,36,38,41,44,46,49,52,54,65,77],...
%     'NoiseWithBrain',[3,4,15,28,31,32,34,47,57,59,68,75,81,84,87,92,94,96,98,99,101:108,110:113,116,120,123,124],...
%     'Doubts',[23,26,28,38,44,52,54,59,81]);
% % Notes for CRE07:
% % Revoir 23,26,28,36,38,44,52,54,59,81,88
% % Ambiguous cases resolved with manual
% 
% autoMoBI_KC(7) = struct('ID','ABE08',...
%     'Brain',[5,7,8,13,16,21,27,28,35,36,45,50,60,61,96],...
%     'BrainWithNoise',[19,37,53,92,94,98],...
%     'NoiseWithBrain',[15,22,30,32,62,67,74,77,80,82,84,85,87,88,90,93,97,99,102:103,105,106,108,111,112],...
%     'Doubts',[77,98,99]);
% % Notes for ABE08:
% % Revoir 30,32,53,77,98,99,102,103
% % Ambiguous cases resolved with manual
% 
% autoMoBI_KC(8) = struct('ID','COL09',...
%     'Brain',[6,7,15:21,30,38,46,58,71,91],...
%     'BrainWithNoise',[10,37,39,56,65,75,77,85,94,96],...
%     'NoiseWithBrain',[8,24,31,42,45,49,68,78,80,82,89:90,98,101,103,105,106,108,111,112,113,119,121],...
%     'Doubts',[37,39,75,90,101,103,105]);



% Notes for COL09:
% Revoir 8,31,32,37,39,75,89,90,101,103,105,106
% Ambiguous cases resolved with manual



%% Bemobil - APP pipeline
app_KC = struct('ID', [], 'Brain', [], 'BrainWithNoise', [], 'NoiseWithBrain',[]);
app_KC(1) = struct('ID','...','Brain',[],...
    'BrainWithNoise',[],...
    'NoiseWithBrain',[]);
% Notes for ...:
%
% Ambiguous cases resolved with manual

%% Bemobil - ASR pipeline
asr_KC = struct('ID', [], 'Brain', [], 'BrainWithNoise', [], 'NoiseWithBrain',[]);
asr_KC(1) = struct('ID','...','Brain',[],...
    'BrainWithNoise',[],...
    'NoiseWithBrain',[]);
% Notes for ...:
%
% Ambiguous cases resolved with manual






