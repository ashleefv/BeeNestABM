function nestData = calculateActivityMatrix(nestData)
    %% Inputs:
        %Takes matrix of the form 'preNest' as input
        % m x n x 3
        %m: number of frames
        %n: number of bees
    %    
    %
    % Outputs
    %    %returns a m x n x 5 matrix, where elements 4 and 5 in the third dimension
    %    %are speed matrix and activity matrix, respectively
    %%
    fps = 2;
    thresh = 10^-3.9;
    preDiffVel = abs(diff(nestData(:,:,1:2)));
    preVels = sqrt(preDiffVel(:,:,1).^2 + preDiffVel(:,:,2).^2);
    preVels = preVels*fps; %Correct for frame rate
    activityMat = preVels > thresh;
    activityMat = double(activityMat);
    activityMat(isnan(preVels)) = NaN;
    speedMat = preVels;
    
    blankRow = nan(1,size(nestData,2));
    speedMat = [blankRow ; speedMat];
    activityMat = [blankRow ; activityMat];
    
    
    %% Append to "nestData"
    nestData(:,:,4) = speedMat;
    nestData(:,:,5) = activityMat;
