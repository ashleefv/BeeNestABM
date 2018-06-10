function [DeltaX_Obj,DeltaY_Obj] = Delta_Obj_app(position,ObjPosition)
% calculates the final - initial x and y coordinates (Delta X and Delta y)
% for the position relative to some object
numBees = length(position);
numObj = length(ObjPosition);
DeltaX_Obj = zeros(numBees,numObj);
DeltaY_Obj = zeros(numBees,numObj);

for i = 1:numBees % each row
    for j = 1:numObj 
        DeltaX_Obj(i,j) = ObjPosition(j,1)-position(i,1);
        DeltaY_Obj(i,j) = ObjPosition(j,2)-position(i,2);
    end
end
end