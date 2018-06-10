function brood = relabelBroodObject_app(brood)
    %% Inputs
    %brood: n x 3 matrix, where n = # of mapped brood element
    %
    % Output:
    %brood: n x 3 matrix, where 1st and 2nd columns are x and y coordinates,
    %and 3rd column contains labels for element type:
    %   1 = brood (eggs, larvae, and pupae)
    %   2 = full food pots
    %   3 = empty pots + wax cover
    
    broodNum = str2num(char(brood(:,3)));
    brInd = ismember(broodNum, [1 2 3]);
    fpInd = ismember(broodNum, [5 6]);
    epInd = ismember(broodNum, [4 7]);
    brood(brInd,3) = 1;
    brood(fpInd,3) = 2;
    brood(epInd,3) = 3;