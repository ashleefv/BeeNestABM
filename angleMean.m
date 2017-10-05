function ang = angleMean(weights,angles)
% https://en.wikipedia.org/wiki/Mean_of_circular_quantities
numBees = size(angles,2);
ang = zeros(1,numBees);
for i = 1:numBees
    % Deal with the case when an angle is set to NaN because there is no
    % object of a certain type within the cutoff radius
    keep_idx = isfinite(angles(:,i));
    % Recalculate the angleMean without being weighted by the term
    % that had no objects within the cutoff radius
    truncatedWeights = weights(keep_idx,i);
    truncatedAngles = angles(keep_idx,i);
    sinAvg = dot(truncatedWeights,sin(truncatedAngles))./sum(truncatedWeights); 
    cosAvg = dot(truncatedWeights,cos(truncatedAngles))./sum(truncatedWeights);
    temp_ang = atan(sinAvg/cosAvg);
    if cosAvg<0
        temp_ang = temp_ang+pi;
    else
        if sinAvg<0
           temp_ang = temp_ang+2*pi; 
        end
    end
    
%     while temp_ang < 0 
%         temp_ang = temp_ang+2*pi;
%     end
%     while temp_ang > 2*pi
%         temp_ang = temp_ang - 2*pi;
%     end
    ang(1,i) = temp_ang;
end
