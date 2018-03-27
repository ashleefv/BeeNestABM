function Activity_Prob_Dist_avg = compute_Activity_Prob_Dist_allColonies(colony_2_include,Pre_Post_Indctr)
%% 
% Inputs:
%    ---colony_2_include---
%    An array of Colony-IDs
%    Size: m x 1
%          where m is the number of colonies being combined to
%          compute the estimates of speed distributions
%    ---Pre_Post_Indctr---
%    A binary 0/1 variable
%    0: Consider Pre-Exposure Data
%    1: Consider Post-Exposure Data
%
% Outputs:
%    ---Activity_Prob_Dist---
%    A matrix defining the Activity_Prob_Dist_Pre or _Post averaged over
%    the colonies
%
for colonyNumber = colony_2_include
    file2load = strcat('data/Essential_Info_Col_',num2str(colonyNumber),'.mat');
    %
    if (Pre_Post_Indctr == 0)
         load(file2load,'Activity_Prob_Dist_Pre');
         Activity_Prob_Dist(colonyNumber,:,:) = Activity_Prob_Dist_Pre;
    elseif  (Pre_Post_Indctr == 1)
        load(file2load,'Activity_Prob_Dist_Post');
        Activity_Prob_Dist(colonyNumber,:,:) = Activity_Prob_Dist_Post;
    end
end
%
Activity_Prob_Dist_avg(:,:) = mean(Activity_Prob_Dist,1);