%% navigate to folder containing data
cd(uigetdir());
load('allDataCol1.mat');

broodPre = relabelBroodObject(broodPre);
broodPost = relabelBroodObject(broodPost);

preNest = calculateActivityMatrix(preNest);
postNest = calculateActivityMatrix(postNest);
%% agent movement plotting


figure(1);
xlm = [0 0.25];
ylm = [0 0.2];
for i = 2000:3000
    %%
    plotCoordinatesAndBrood(preNest, broodPre,i);
    hold on;
    addBumpLabels(preNest, i, 0.01);
    xlim(xlm);
    ylim(ylm);
    drawnow
    hold off
end


%% plot effects by group, pre and post
figure(1)
grp = orTagTreat == 1;
plotComparisons(preNest(:,grp,:), postNest(:,grp,:), broodPre, broodPost)

figure(2)
grp = orTagTreat == 3;
plotComparisons(preNest(:,grp,:), postNest(:,grp,:), broodPre, broodPost)
