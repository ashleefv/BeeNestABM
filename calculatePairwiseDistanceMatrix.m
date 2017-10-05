function distMat = calculatePairwiseDistanceMatrix(nestData)
%Calculates pairwise distance for an MxNx2 video data, where M is the
%number of frames in the video, N is the number of bees, and the third
%dimension represents x and y coordinates. Tags is the taglist, indexed in
%the same order as the columns of nestData
%
%nestDat = preNest;
%vis = 0;
distMat = nan(size(nestData,2), size(nestData,2), size(nestData,1));

for i = 1:size(nestData,1)
    %%
    ind = ~isnan(nestData(i,:,1));
    curData = nestData(i,ind,:);
    dists = squareform(pdist([curData(:,:,1)' curData(:,:,2)']));
    distMat(ind,ind,i) = dists;
end