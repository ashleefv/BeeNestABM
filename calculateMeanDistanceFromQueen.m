function distQueen = calculateMeanDistanceFromQueen(nestData)
    %%
    %Input: nestData
    % m x n x l matrix, where m = nframes, n = nbees, and the first
    %two elements of 'l' are the x and y coordinates, respectively
    % Queen must be first column!
    %
    % Output:
    % distSocCenter = mean instantaneous distance from the nest "social center"
    
    qpx = nestData(:,1,1);
    qpy = nestData(:,1,2);
    %%
    xd = bsxfun(@minus, nestData(:,:,1), qpx);
    yd = bsxfun(@minus, nestData(:,:,2), qpy);
    
    totDists = sqrt(xd.^2 + yd.^2);
    
    distQueen = nanmean(totDists);
