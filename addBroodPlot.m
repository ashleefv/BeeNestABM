function [] = addBroodPlot(brood)
    %% input: brood object, formatted as output of "relabelBroodObject: output
    brInd = brood(:,3) == 1;
    fpInd = brood(:,3) == 2;
    epInd = brood(:,3) == 3;
    
    
    %% Plot
    %Set graphical parameters
    size = 300;
    FaceAlpha = 0.2;
    EdgeAlpha = 0.5;
    
    %Plot
    scatter(brood(brInd,1), brood(brInd,2), size, 'MarkerFaceAlpha',FaceAlpha, 'MarkerEdgeAlpha', EdgeAlpha, 'MarkerFaceColor', [0.9 0.9 0.3], 'MarkerEdgeColor', 'k');
    hold on
    scatter(brood(epInd,1), brood(epInd,2), size, 'MarkerFaceAlpha',0, 'MarkerEdgeAlpha', EdgeAlpha, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k');
    scatter(brood(fpInd,1), brood(fpInd,2), size, 'MarkerFaceAlpha',FaceAlpha, 'MarkerEdgeAlpha', EdgeAlpha, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k');
