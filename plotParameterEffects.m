% Initial state of the agents
summaryData = nan(11,6);
k = 1
for i = 1:4
    %%
colonyNumber = i; % 4 colonies numbered 1 - 4
colony = load(['data/allDataCol' num2str(colonyNumber) '.mat']);
brood = relabelBroodObject(colony.broodPre);
numFrames = size(colony.preNest,1);
numBees = size(colony.preNest,2);
tags = colony.orTagTreat;
estimatedData = load(['data/Essential_Info_Col_' num2str(colonyNumber) '.mat']);

%%
estimatedTransProb = estimatedData.transProb_Post;
    for tagNumber = 1:3 % untreated, sucrose control, low dose, high dose
        %%
        summaryData(k,1) = estimatedTransProb{1,tagNumber+1}(2,1);
        summaryData(k,2) = estimatedTransProb{1,tagNumber+1}(1,2);
        summaryData(k,3) = estimatedTransProb{2,tagNumber+1}(2,1);
        summaryData(k,4) = estimatedTransProb{2,tagNumber+1}(1,2); 
        summaryData(k,5) = tagNumber; %Which treatment group?
        summaryData(k,6) = i; %Which colony?
        k = k+1;
    end
end
%%
summaryData = array2table(summaryData);
summaryData.Properties.VariableNames = {'AtoI_Unbumped', 'ItoA_Unbumped', 'AtoI_Bumped', 'ItoA_Bumped', 'treatment', 'colony'};

%%
subplot(2,2,1);
boxplot(summaryData.AtoI_Unbumped, summaryData.treatment)
title('AtoIUnbumped');

subplot(2,2,2);
boxplot(summaryData.ItoA_Unbumped, summaryData.treatment);
title('ItoAUnbumped');

subplot(2,2,3);
boxplot(summaryData.AtoI_Bumped, summaryData.treatment);
title('AtoIBumped');

subplot(2,2,4);
boxplot(summaryData.ItoA_Bumped, summaryData.treatment);
title('ItoABumped');