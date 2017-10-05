function [] = addBumpLabels(nestData, frame, distThresh)
    %% Inputs
    % nestData: tracked nest data matrix, of form "preNest" (m x n x 5)
    % frame: what frame/timestep number to plot
    %distThresh: distance threshold for a legitimate "interaction"
    dist = calculatePairwiseDistanceMatrix(nestData, frame); %calculatedistances
    dist = triu(dist);
    dist(dist == 0) = NaN;
    
    intMat = dist < distThresh;
    nbees = size(intMat,1);
    for i = 1:nbees
        
        for j = 1:nbees
            
            if intMat(i,j) == 1
                
                plot(nestData(frame, [i j], 1), nestData(frame, [i j], 2), 'r-');
                
                hold on
            end
            
        end
    end
