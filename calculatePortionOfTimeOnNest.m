function [porTimeOnNest onNest] = calculatePortionOfTimeOnNest(nestData, brood, distThresh)
    %% inputs:
    % nestData: tracked nest data matrix, of form "preNest" (m x n x 5)
    % brood: brood, formatted as output of 'relabelBroodObject'
    %
    %Output:
    %porTimeOnNest: portion of time spent physically located on nest,
    
    minDistToNest = nan(size(nestData(:,:,1)));
    for i = 1:size(nestData,2)
        %%
        xy = permute(nestData(:,i,1:2), [1 3 2]);

        minDistToNest(:,i) =  min(pdist2(xy, brood(:,1:2)),[],2);
    end
    
    onNest = double(minDistToNest < distThresh);
    onNest(isnan(minDistToNest)) = NaN;
    porTimeOnNest = nanmean(onNest);