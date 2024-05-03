% import the data of all participants
filename = fullfile('/Users/ilaria/Documents/MATLAB/data/EEG-AFF', 'MergedDatabase.xlsx');
AllParticipants = table2struct(readtable(filename));


for i=1:length(AllParticipants)
    AllParticipants(i).Responses = str2double(AllParticipants(i).Responses);
end

% split the groups into 2 separate Databases
old = AllParticipants(strcmp({AllParticipants.Group},'O'));
young = AllParticipants(strcmp({AllParticipants.Group},'Y'));

% number of old/young participants
n_old = length(unique(split([old.ID],"P")));
n_young = length(unique(split([young.ID],"P")));

% split condition trials
Short_young = young(strcmp({young.Condition},'Short'));
Medium_young = young(strcmp({young.Condition},'Medium'));
Long_young = young(strcmp({young.Condition},'Long'));
Short_old = old(strcmp({old.Condition},'Short'));
Medium_old = old(strcmp({old.Condition},'Medium'));
Long_old = old(strcmp({old.Condition},'Long'));

% 1 plot, 3 sublot pour les conditions
for c=1:3 % 3 condition
    s = subplot(3, 1, c);
    hold on
    if c==1 % short
        answers_old = [Short_old.Responses]';
        answers_young = [Short_young.Responses]';
        % statistical ttest
        [h_short,p_short]=ttest2(answers_old,answers_young);
        %short zone area
        y = [2 2 3 3];
    elseif c==2 % medium
        answers_old = [Medium_old.Responses]';
        answers_young = [Medium_young.Responses]';
        % statistical ttest
        [h_medium,p_medium]=ttest2(answers_old,answers_young);
        %medium zone area
        y = [4.5 4.5 5.5 5.5];
    elseif c==3 % long
        answers_old = [Long_old.Responses]';
        answers_young = [Long_young.Responses]';
        % statistical ttest
        [h_long,p_long]=ttest2(answers_old,answers_young);
        %long zone area
        y = [7 7 8 8];
    end
    
    % create column with responses
    answers = [answers_young;answers_old];
    x = [0.5 3.5 3.5 0.5];
    p = patch(x,y,'green','EdgeColor',[1 1 1]);
    alpha(p,0.09)
    
    Groups = repelem(["Young";"Old"], [length(answers_young);length(answers_old)]);
    
    results = table(Groups,answers);
    GroupOrder = {'Young','Old'};
    results.Groups = categorical(results.Groups,GroupOrder);
    
    set(gcf, 'Position', [100,100,1500,700])
    
    % 3 boxplot for the 3 zones and mean distances
    boxplot(results.answers,results.Groups, 'Labels',{'Young','Old'})
    %title('Condition Comaraison');
    ax = gca;
    ax.FontSize = 12;
    
    if c==1
        title('SHORT ZONE')
    elseif c==2
        title('MEDIUM ZONE')
    elseif c==3
        title('LONG ZONE')
    end
    
    
end

%% save fig in the current folder
subject = 'Group';
saveCurrentFig(['/Users/ilaria/Documents/MATLAB/data/EEG-AFF' filesep],...
    [subject, '_boxplot'], {'png'}, [1000,700]);

% 3 plot (1 per condition), 1 sublot per block
%% statistical analysis
% all distances old vs young
[h_global,p_global]=ttest2([old.Responses]',[young.Responses]');

% leaning ? block 1 vs bock 3 old
b1_old = old([old.Block]==1);
b3_old = old([old.Block]==3);
[ho,po] = ttest2([b1_old.Responses]',[b3_old.Responses]');

% learning ? block 1 vs bock 3 young
b1_young = young([young.Block]==1);
b3_young = young([young.Block]==3);
[hy,py] = ttest2([b1_young.Responses]',[b3_young.Responses]');

% tailed ttest to see if the absolute error is less in B1 then B3 => learning
[ho_larning,po_learning] = ttest2([b1_old.Err_absolue]',[b3_old.Err_absolue]')
[hy_larning,py_learning] = ttest2([b1_young.Err_absolue]',[b1_young.Err_absolue]')
% one tailed
[ho_larning,po_learning] = ttest2([b1_old.Err_absolue]',[b3_old.Err_absolue]','Tail','Right')
[hy_larning,py_learning] = ttest2([b1_young.Err_absolue]',[b1_young.Err_absolue]','Tail','Right')



