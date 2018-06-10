function resultant_angle = angleObj_app(currentDistanceToObj,DeltaX_Obj,DeltaY_Obj,cutoffRadius)

% sum all the x distances and y distances to create a resultant vector for
% those pairs of points that lie within the cutoffRadius from the current
% position
numBees = size(currentDistanceToObj,1);
% numObj = size(currentDistanceToObj,2);
DeltaX_total = zeros(1,numBees);
DeltaY_total = zeros(1,numBees);
resultant_angle = zeros(1,numBees);
for i = 1:numBees
    DeltaX_Obj(i,currentDistanceToObj(i,:)> cutoffRadius(i,1))=0;
    DeltaY_Obj(i,currentDistanceToObj(i,:)> cutoffRadius(i,1))=0;
    DeltaX_total(1,i) = sum(DeltaX_Obj(i,:),2);
    DeltaY_total(1,i) = sum(DeltaY_Obj(i,:),2);
    if abs(DeltaY_total(1,i)) > 0 || abs(DeltaX_total(1,i)) > 0
        temp_ang = atan(DeltaY_total(1,i)/DeltaX_total(1,i));
        if DeltaX_total(1,i)<0
            temp_ang = temp_ang+pi;
        else
            if DeltaY_total(1,i)<0
               temp_ang = temp_ang+2*pi; 
            end
        end
        resultant_angle(1,i) = temp_ang;
    else % there is no object within the cutoff radius of the ith bee
        resultant_angle(1,i) = NaN;
    end
end

end