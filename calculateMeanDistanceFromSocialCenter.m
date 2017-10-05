function distSocCenter = calculateMeanDistanceFromSocialCenter(nestData)
    %%
    %Input: nestData
    % m x n x l matrix, where m = nframes, n = nbees, and the first
    %two elements of 'l' are the x and y coordinates, respectively
    %
    % Output:
    % distSocCenter = mean instantaneous distance from the nest "social center"
    
    xm = nanmean(nanmean(nestData(:,:,1)));
    ym = nanmean(nanmean(nestData(:,:,2)));
    
    xd = nestData(:,:,1) - xm;
    yd = nestData(:,:,2) - ym;
    
    totDists = sqrt(xd.^2 + yd.^2);
    
    distSocCenter = nanmean(totDists);
