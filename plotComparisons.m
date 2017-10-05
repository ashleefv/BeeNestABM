function [] = plotComparisons(nestData1, nestData2, brood1, brood2)
    %% inputs:
    % nestData1 and nestData2, m x n x 5 matrices, of form "preNest" (m x n x 5)
    %
    % brood1 and brood2: brood, formatted as output of 'relabelBroodObject',
    % for the first and group respectively
    
    
    %Create grouping variables
    x1 = repelem(0,size(nestData1,2));
    x2 = repelem(1,size(nestData1,2));
    grps = [x1 x2];
    
    
    %% activity (portion of time active)
    subplot(2,2,1);
    act1 = nanmean(nestData1(:,:,5));
    act2 = nanmean(nestData2(:,:,5));
    act = [act1 act2];
    boxplot(act,grps)
    title('Portion of time active');
    
    
    %% mean instantaneous distance from social center
    subplot(2,2,2);
    socDist1 = calculateMeanDistanceFromSocialCenter(nestData1);
    socDist2 = calculateMeanDistanceFromSocialCenter(nestData2);
    
    socDist = [socDist1 socDist2];
    boxplot(socDist,grps);
    title('Mean distance from soc center');
    
    
    %% mean inst. distance from queen
    subplot(2,2,3);
    queenDist1 = calculateMeanDistanceFromQueen(nestData1);
    queenDist2 = calculateMeanDistanceFromQueen(nestData2);
    queenDist = [queenDist1 queenDist2];
    boxplot(queenDist, grps);
    title('Mean distance from queen');
    
    
    %% portion of time on the nest structure
    subplot(2,2,4);
    porTimeOnNest1 = calculatePortionOfTimeOnNest(nestData1, brood1, 0.01);
    porTimeOnNest2 = calculatePortionOfTimeOnNest(nestData2, brood2, 0.01);
    porTimeOnNest = [porTimeOnNest1 porTimeOnNest2];
    boxplot(porTimeOnNest, grps);
    title('Portion of time on nest structure');
    
    %% Add network plots?
    
end