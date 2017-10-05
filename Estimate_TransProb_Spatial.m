clear; clc; close all
%
% Change the COLONY ID to select which colony data you want to process
Colony_ID = 4;
% Bees are considered to be bumping onto each other when their distance is
% below the "BeeBodyThreshold"
BeeBodyThreshold = 0.01;
%
DistToDetect_Nest_OnOff = 1*BeeBodyThreshold;
% Loads the relevant data file
% The data file should be inside "..\data" folder
% This relative location should be maintained for proper functioning
file2load = strcat('data\allDataCol',num2str(Colony_ID),'.mat');
load(file2load);
%
%%
% Counts the total number of time-frames and bees in the dataset
% Number of time-frames can be different between pre- and post- exposure
% datasets
Total_Frames = size(distMatPre,3);
Number_of_Bees = size(distMatPre,1);
%
% Computes the activity levels in the individual bees
preNest_MDFD = calculateActivityMatrix(preNest);
actvty_Storage_Pre = (preNest_MDFD(:,:,5))';
%
%
% Computes whether a pair of bees are bumping into each other
bumped_Storage = zeros(Number_of_Bees,Total_Frames);
for time_idx = 1:Total_Frames
    %
    current_Distance = distMatPre(:,:,time_idx);
    bumped_Storage(:,time_idx) = bump(BeeBodyThreshold,current_Distance);
    %
end
%
% Stores the bump-indicators inside a vector
bumped_Vector = bumped_Storage(:,1:Total_Frames-1);
bumped_Vector = bumped_Vector(:);
% Stores the pre-bump activity level (active or inactive) inside a vector
act_pr_Vector = actvty_Storage_Pre(:,1:Total_Frames-1);
act_pr_Vector = act_pr_Vector(:);
% Stores the post-bump activity level (active or inactive) inside a vector
act_fu_Vector = actvty_Storage_Pre(:,2:Total_Frames);
act_fu_Vector = act_fu_Vector(:);
% From the vectors of pre- and post-bump activity levels computes a
% variable to identify the transition (one of the following four types:
% 0->0, 0->1, 1->0, 1->1
Activity_Transition_Indicator_Pre = 2*act_pr_Vector + act_fu_Vector;
%
%
% Copies the Activity Transition Indicators into two variables
% One for those close to the nest structures
% And another one for those away from the nest structures
Activity_Transition_OnNest_Pre = Activity_Transition_Indicator_Pre;
Activity_Transition_OffNest_Pre = Activity_Transition_Indicator_Pre;
%
bumped_Vector_OnNest = bumped_Vector;
bumped_Vector_OffNest = bumped_Vector;
%
[WeDontCareHere,ON_or_OFF_Nest] = calculatePortionOfTimeOnNest(preNest_MDFD, broodPre, DistToDetect_Nest_OnOff);
% The entries of the vector "ON_or_OFF_Nest" are either 0 or 1. 
% 0 -->> Away from the Nest 
% 1 -->> Close to the Nest
ON_or_OFF_Nest = ON_or_OFF_Nest';
ON_or_OFF_Nest = ON_or_OFF_Nest(:,1:Total_Frames-1);
ON_or_OFF_Nest = ON_or_OFF_Nest(:);
ON_or_OFF_Nest(isnan(ON_or_OFF_Nest)) = 0;
%
% Masks the irrelevant entries
Activity_Transition_OnNest_Pre(ON_or_OFF_Nest==0) = NaN;
Activity_Transition_OffNest_Pre(ON_or_OFF_Nest==1) = NaN;
%
% Removes the N-a-N entries 
nan_Indi_1 = isnan(Activity_Transition_OnNest_Pre);
nan_Indi_2 = isnan(Activity_Transition_OffNest_Pre);
%
Activity_Transition_OnNest_Pre(nan_Indi_1) = [];
Activity_Transition_OffNest_Pre(nan_Indi_2) = [];
%
bumped_Vector_OnNest(nan_Indi_1) = [];
bumped_Vector_OffNest(nan_Indi_2) = [];
%
% Computes the Transition Probabilities for both pre- and post-bump
% conditions
transProb_Pre_OnNest = cell(2,1);
condn_bmp_Pre_OnNest = NaN(2,1);
%
transProb_Pre_OffNest = cell(2,1);
condn_bmp_Pre_OffNest = NaN(2,1);
%
for bump_indi_val = 0:1
    %
    condn_bmp_Pre_OnNest(bump_indi_val+1,1) = bump_indi_val;
    condn_bmp_Pre_OffNest(bump_indi_val+1,1) = bump_indi_val;
    %
    % Selects the rows corresponding to appropriate values (0:pre-bump or
    % 1:post-bump) of the bump indicator variable
    temp_Activity_Transition_OnNest = Activity_Transition_OnNest_Pre(bumped_Vector_OnNest == bump_indi_val);
    temp_Activity_Transition_OffNest = Activity_Transition_OffNest_Pre(bumped_Vector_OffNest == bump_indi_val);
    %
    % The first row corresponds to the transitions from "Inactive State"
    transProb_Pre_OnNest{bump_indi_val+1}(1,:) = [sum(temp_Activity_Transition_OnNest==0) sum(temp_Activity_Transition_OnNest==1)];
    transProb_Pre_OnNest{bump_indi_val+1}(1,:) = transProb_Pre_OnNest{bump_indi_val+1}(1,:)/sum(transProb_Pre_OnNest{bump_indi_val+1}(1,:));
    transProb_Pre_OffNest{bump_indi_val+1}(1,:) = [sum(temp_Activity_Transition_OffNest==0) sum(temp_Activity_Transition_OffNest==1)];
    transProb_Pre_OffNest{bump_indi_val+1}(1,:) = transProb_Pre_OffNest{bump_indi_val+1}(1,:)/sum(transProb_Pre_OffNest{bump_indi_val+1}(1,:));
    %
    % The first row corresponds to the transitions from "Active State"
    transProb_Pre_OnNest{bump_indi_val+1}(2,:) = [sum(temp_Activity_Transition_OnNest==2) sum(temp_Activity_Transition_OnNest==3)];
    transProb_Pre_OnNest{bump_indi_val+1}(2,:) = transProb_Pre_OnNest{bump_indi_val+1}(2,:)/sum(transProb_Pre_OnNest{bump_indi_val+1}(2,:));
    transProb_Pre_OffNest{bump_indi_val+1}(2,:) = [sum(temp_Activity_Transition_OffNest==2) sum(temp_Activity_Transition_OffNest==3)];
    transProb_Pre_OffNest{bump_indi_val+1}(2,:) = transProb_Pre_OffNest{bump_indi_val+1}(2,:)/sum(transProb_Pre_OffNest{bump_indi_val+1}(2,:));
    %
end
%
%
%%
%
% Counts the total number of time-frames and bees in the dataset
% Number of time-frames can be different between pre- and post- exposure
% datasets
Total_Frames = size(distMatPost,3);
Number_of_Bees = size(distMatPost,1);
%
% Computes the activity levels in the individual bees
postNest_MDFD = calculateActivityMatrix(postNest);
actvty_Storage_Post = (postNest_MDFD(:,:,5))';
%
% % Stores exposure categories inside a matrix
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
Activity_Transition_Indicator_Post = 2*act_pr_Vector + act_fu_Vector;
%
%
% Copies the Activity Transition Indicators into two variables
% One for those close to the nest structures
% And another one for those away from the nest structures
Activity_Transition_OnNest_Post = Activity_Transition_Indicator_Post;
Activity_Transition_OffNest_Post = Activity_Transition_Indicator_Post;
%
bumped_Vector_OnNest_Post = bumped_Vector;
bumped_Vector_OffNest_Post = bumped_Vector;
%
pstcde_Vector_OnNest_Post = pstcde_Vector;
pstcde_Vector_OffNest_Post = pstcde_Vector;
%
[WeDontCareHere,ON_or_OFF_Nest_Post] = calculatePortionOfTimeOnNest(postNest_MDFD, broodPost, DistToDetect_Nest_OnOff);
% The entries of the vector "ON_or_OFF_Nest" are either 0 or 1. 
% 0 -->> Away from the Nest 
% 1 -->> Close to the Nest
ON_or_OFF_Nest = ON_or_OFF_Nest';
ON_or_OFF_Nest = ON_or_OFF_Nest(:,1:Total_Frames-1);
ON_or_OFF_Nest = ON_or_OFF_Nest(:);
ON_or_OFF_Nest(isnan(ON_or_OFF_Nest)) = 0;
%
% Masks the irrelevant entries
Activity_Transition_OnNest_Post(ON_or_OFF_Nest==0) = NaN;
Activity_Transition_OffNest_Post(ON_or_OFF_Nest==1) = NaN;
%
% Removes the N-a-N entries 
nan_Indi_1 = isnan(Activity_Transition_OnNest_Post);
nan_Indi_2 = isnan(Activity_Transition_OffNest_Post);
%
Activity_Transition_OnNest_Post(nan_Indi_1) = [];
Activity_Transition_OffNest_Post(nan_Indi_2) = [];
%
bumped_Vector_OnNest_Post(nan_Indi_1) = [];
bumped_Vector_OffNest_Post(nan_Indi_2) = [];
%
pstcde_Vector_OnNest_Post(nan_Indi_1) = [];
pstcde_Vector_OffNest_Post(nan_Indi_2) = [];
%
%
% Computes the Transition Probabilities for both pre- and post-bump
% conditions
transProb_Post_OnNest = cell(2,4);
condn_bmp_Post_OnNest = NaN(2,4);
condn_ExP_Post_OnNest = NaN(2,4);
%
transProb_Post_OffNest = cell(2,4);
condn_bmp_Post_OffNest = NaN(2,4);
condn_ExP_Post_OffNest = NaN(2,4);
%
%
for bump_indi_val = 0:1
    for expsre_indi_val = 0:3
        %
        condn_bmp_Post_OnNest(bump_indi_val+1,expsre_indi_val+1) = bump_indi_val;
        condn_bmp_Post_OffNest(bump_indi_val+1,expsre_indi_val+1) = bump_indi_val;
        %
        condn_ExP_Post_OnNest(bump_indi_val+1,expsre_indi_val+1) = expsre_indi_val;
        condn_ExP_Post_OffNest(bump_indi_val+1,expsre_indi_val+1) = expsre_indi_val;
        %
        temp_Activity_Transition_OnNest = Activity_Transition_OnNest_Post(4*bumped_Vector_OnNest_Post + pstcde_Vector_OnNest_Post == 4*bump_indi_val + expsre_indi_val);
        temp_Activity_Transition_OffNest = Activity_Transition_OffNest_Post(4*bumped_Vector_OffNest_Post + pstcde_Vector_OffNest_Post == 4*bump_indi_val + expsre_indi_val);
        %
        trans_prob_TEMP_OnNest(1,:) = [sum(temp_Activity_Transition_OnNest==0) sum(temp_Activity_Transition_OnNest==1)];
        trans_prob_TEMP_OnNest(1,:) = trans_prob_TEMP_OnNest(1,:)/sum(trans_prob_TEMP_OnNest(1,:));
        %
        trans_prob_TEMP_OffNest(1,:) = [sum(temp_Activity_Transition_OffNest==0) sum(temp_Activity_Transition_OffNest==1)];
        trans_prob_TEMP_OffNest(1,:) = trans_prob_TEMP_OffNest(1,:)/sum(trans_prob_TEMP_OffNest(1,:));
        %
        trans_prob_TEMP_OnNest(2,:) = [sum(temp_Activity_Transition_OnNest==2) sum(temp_Activity_Transition_OnNest==3)];
        trans_prob_TEMP_OnNest(2,:) = trans_prob_TEMP_OnNest(2,:)/sum(trans_prob_TEMP_OnNest(2,:));
        %
        trans_prob_TEMP_OffNest(2,:) = [sum(temp_Activity_Transition_OffNest==2) sum(temp_Activity_Transition_OffNest==3)];
        trans_prob_TEMP_OffNest(2,:) = trans_prob_TEMP_OffNest(2,:)/sum(trans_prob_TEMP_OffNest(2,:));
        %
        transProb_Post_OnNest{bump_indi_val+1,expsre_indi_val+1} = trans_prob_TEMP_OnNest;
        transProb_Post_OffNest{bump_indi_val+1,expsre_indi_val+1} = trans_prob_TEMP_OffNest;
        %
    end
    %
end
%
%
%%
%
% saves the result inside a file
% Desitination file name: TransProb_Spatial_Info_*.mat
% Destinatio folder (relative location): ..\data
file_2_save = strcat('data\TransProb_Spatial_Info_',num2str(Colony_ID),'.mat');
save(file_2_save, 'transProb_Pre_OnNest', 'condn_bmp_Pre_OnNest', ...
    'transProb_Pre_OffNest', 'condn_bmp_Pre_OffNest', ...
    'transProb_Post_OnNest', 'condn_bmp_Post_OnNest', ...
    'condn_ExP_Post_OnNest', 'transProb_Post_OffNest', ...
    'condn_bmp_Post_OffNest', 'condn_ExP_Post_OffNest');
%