function transProb_Post = calculateBumpTransitionProb(postNest, distMatPost, orTagTreat)
    %
% Counts the total number of time-frames and bees in the dataset
% Number of time-frames can be different between pre- and post- exposure
% datasets

BeeBodyThreshold = 0.01;

Total_Frames = size(distMatPost,3);
Number_of_Bees = size(distMatPost,1);
%
% Computes the activity levels in the individual bees
postNest_MDFD = calculateActivityMatrix(postNest);
actvty_Storage_Post = (postNest_MDFD(:,:,5))';
%
% The first row of the dataset corresponds to the queen
Activity_QUEEN_Post = actvty_Storage_Post(1,:);
%
% Stores acitivity of the other members of the group
act_Store_post_others = actvty_Storage_Post(2:Number_of_Bees,:);
%
% Extracts the exposure categories of the bees into a single vector
% Categories: 
% (0)  -- Not touched at all during the experiment
% (1)  -- Control Group: treated with Sucrose
% (2)  -- Semi-Target Group: treated with a small amount of neonicotinoid
% (3)  -- Target Group: treated with maximum amount of neonicotinoid
Exposure_Condn_others = orTagTreat(2:Number_of_Bees,:);
%
% Extracts the acitivity levels in individual categories
Activity_OTHER_unRM = act_Store_post_others(Exposure_Condn_others==0,:);
Activity_OTHER_Cntl = act_Store_post_others(Exposure_Condn_others==1,:); 
Activity_OTHER_mEXP = act_Store_post_others(Exposure_Condn_others==2,:); 
Activity_OTHER_fEXP = act_Store_post_others(Exposure_Condn_others==3,:); 
% Storing acitivity values of individual categories in separate vectors
Activity_QUEEN_Post = Activity_QUEEN_Post(:);
Activity_OTHER_unRM = Activity_OTHER_unRM(:);
Activity_OTHER_Cntl = Activity_OTHER_Cntl(:); 
Activity_OTHER_mEXP = Activity_OTHER_mEXP(:); 
Activity_OTHER_fEXP = Activity_OTHER_fEXP(:); 
% Cleans up the N-a-N entries
Activity_QUEEN_Post(isnan(Activity_QUEEN_Post)==1) = [];
Activity_OTHER_unRM(isnan(Activity_OTHER_unRM)==1) = [];
Activity_OTHER_Cntl(isnan(Activity_OTHER_Cntl)==1) = [];
Activity_OTHER_mEXP(isnan(Activity_OTHER_mEXP)==1) = [];
Activity_OTHER_fEXP(isnan(Activity_OTHER_fEXP)==1) = [];
%
% Computes the Probability of being active in individual categories
Activity_Prob_Dist_Post = NaN(2,5);
Activity_Prob_Dist_Post(:,1) = ([sum(Activity_QUEEN_Post==0) sum(Activity_QUEEN_Post==1)]')/length(Activity_QUEEN_Post);
Activity_Prob_Dist_Post(:,2) = ([sum(Activity_OTHER_unRM==0) sum(Activity_OTHER_unRM==1)]')/length(Activity_OTHER_unRM);
Activity_Prob_Dist_Post(:,3) = ([sum(Activity_OTHER_Cntl==0) sum(Activity_OTHER_Cntl==1)]')/length(Activity_OTHER_Cntl);
Activity_Prob_Dist_Post(:,4) = ([sum(Activity_OTHER_mEXP==0) sum(Activity_OTHER_mEXP==1)]')/length(Activity_OTHER_mEXP);
Activity_Prob_Dist_Post(:,5) = ([sum(Activity_OTHER_fEXP==0) sum(Activity_OTHER_fEXP==1)]')/length(Activity_OTHER_fEXP);
%
% Stores exposure categories inside a matrix
Exposure_Condn_Post = orTagTreat*ones(1,Total_Frames-1);
%
% Computes whether a pair of bees are bumping into each other
bumped_Storage = zeros(Number_of_Bees,Total_Frames);
for time_idx = 1:Total_Frames
    %
    current_Distance = distMatPost(:,:,time_idx);
    bumped_Storage(:,time_idx) = bump(BeeBodyThreshold,current_Distance);
    %
end
%
% Stores the bump-indicators inside a vector
bumped_Vector = bumped_Storage(:,1:Total_Frames-1);
bumped_Vector = bumped_Vector(:);
% Stores the pre-bump activity level (active or inactive) inside a vector
act_pr_Vector = actvty_Storage_Post(:,1:Total_Frames-1);
act_pr_Vector = act_pr_Vector(:);
% Stores the post-bump activity level (active or inactive) inside a vector
act_fu_Vector = actvty_Storage_Post(:,2:Total_Frames);
act_fu_Vector = act_fu_Vector(:);
% Stores the category of  an individual inside a vector
pstcde_Vector = Exposure_Condn_Post(:);
%
% From the vectors of pre- and post-bump activity levels computes a
% variable to identify the transition (one of the following four types:
% 0->0, 0->1, 1->0, 1->1
Activity_Transition_Indicator = 2*act_pr_Vector + act_fu_Vector;
% Cleans up the N-a-N entries 
Activity_Transition_Indicator(isnan(2*act_pr_Vector + act_fu_Vector)==1) = [];
bumped_Vector(isnan(2*act_pr_Vector + act_fu_Vector)==1) = [];
pstcde_Vector(isnan(2*act_pr_Vector + act_fu_Vector)==1) = [];
% Computes the percentage of data lost due to this clean up
Data_Thrw_Post = 100*(1 - length(Activity_Transition_Indicator)/length(act_fu_Vector));
%
% Computes the Transition Probabilities for both pre- and post-bump
% conditions
transProb_Post = cell(2,4);
condn_bmp_Post = NaN(2,4);
condn_ExP_Post = NaN(2,4);
%
for bump_indi_val = 0:1
    for expsre_indi_val = 0:3
        %
        condn_bmp_Post(bump_indi_val+1,expsre_indi_val+1) = bump_indi_val;
        condn_ExP_Post(bump_indi_val+1,expsre_indi_val+1) = expsre_indi_val;
        %
        temp_Activity_Transition = Activity_Transition_Indicator(4*bumped_Vector + pstcde_Vector == 4*bump_indi_val + expsre_indi_val);
        %
        trans_prob_TEMP(1,:) = [sum(temp_Activity_Transition==0) sum(temp_Activity_Transition==1)];
        trans_prob_TEMP(1,:) = trans_prob_TEMP(1,:)/sum(trans_prob_TEMP(1,:));
        %
        trans_prob_TEMP(2,:) = [sum(temp_Activity_Transition==2) sum(temp_Activity_Transition==3)];
        trans_prob_TEMP(2,:) = trans_prob_TEMP(2,:)/sum(trans_prob_TEMP(2,:));
        %
        transProb_Post{bump_indi_val+1,expsre_indi_val+1} = trans_prob_TEMP;
        %
    end
    %
end
%
%
%Rearrange
transProb