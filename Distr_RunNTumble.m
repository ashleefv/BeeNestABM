clear; close all; clc
%
load('E:\Dropbox\pesticideProject\data\allDataCol1.mat');

%% Alt
% cd(uigetdir());
% 
% load('sampleDataPre.mat');
%%
Total_Frames = size(postNest,1);
Number_of_Bees = size(postNest,2);
%
Trajectory_Storage = cell(Number_of_Bees,1);
Intr_Frame_Distr = NaN(Total_Frames-1, Number_of_Bees);
Turn_Angle_Distr = NaN(Total_Frames-2, Number_of_Bees);
%
Mat_Num = [0 -1; 1 0];
Mat_Den = eye(2);
%
for bee_idx = 1:Number_of_Bees
    %
    traj_Current = [postNest(:,bee_idx,1), postNest(:,bee_idx,2)];
    Trajectory_Storage{bee_idx} = traj_Current;
    %
    plot(traj_Current(:,1), postNest(:,2));
    %
    rel_displcmnt = diff(traj_Current);
    %
    distance_betn_datapt = sqrt(sum(rel_displcmnt.^2,2));
    %
    Num_atan2_turn_Angle = diag(rel_displcmnt(1:Total_Frames-2,:)*Mat_Num*(rel_displcmnt(2:Total_Frames-1,:))');
    Den_atan2_turn_Angle = diag(rel_displcmnt(1:Total_Frames-2,:)*Mat_Den*(rel_displcmnt(2:Total_Frames-1,:))');
    trnAngle_betn_datapt = -atan2d(Num_atan2_turn_Angle,Den_atan2_turn_Angle);
    %
    Intr_Frame_Distr(:,bee_idx) = distance_betn_datapt;
    Turn_Angle_Distr(:,bee_idx) = trnAngle_betn_datapt;
    %
end
%
figure(1)
bin_trn_angl = -175:10:175;
n_count = hist(Turn_Angle_Distr(:),bin_trn_angl);
n_count_normalized = n_count/max(n_count);
polarplot(bin_trn_angl*(pi/180),n_count_normalized,'ro-')
%
figure(2)
hist(Intr_Frame_Distr)