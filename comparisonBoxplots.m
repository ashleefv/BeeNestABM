function comparisonBoxplots(data1, data2, brood, grouplist)
    %Input two "nestData" shaped datasets, brood object (to pass to summary function) and a tagged group, 
    %
    %Generate boxplots 
    %% Generate summary statistics
    [means1 distributions1] = calculateSummaryStatistics(data1, brood, grouplist);
    [means2 distributions2] = calculateSummaryStatistics(data2, brood, grouplist);
    %% boxplots!
    %Activity
    figure(1);
    subplot(2,1,1);
boxplot(distributions1.activity, grouplist);
    title('Activity Level');
subplot(2,1,2);
boxplot(distributions2.activity, grouplist);

    figure(2);
    subplot(2,1,1);
boxplot(distributions1.distanceToNestmates, grouplist);
    title('Distance To Nestmates');
subplot(2,1,2);
boxplot(distributions2.distanceToNestmates, grouplist);

    figure(3);
    subplot(2,1,1);
boxplot(distributions1.distanceToBrood, grouplist);
    title('Distance To Brood');
subplot(2,1,2);
boxplot(distributions2.distanceToBrood, grouplist);
