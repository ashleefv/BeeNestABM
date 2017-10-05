%Add orientation data and clean a bit
%manually select correct 'allData' file
clear all
[filename pathname] = uigetfile();
load([pathname filename]);

cd(pathname);

%
preFile = dir('*nest*Pre*tracked.mat');
load(preFile(1).name);
preTracks = trackingData;
postFile = dir('*nest*Post*tracked.mat');
load(postFile(1).name);
postTracks = trackingData;

figure(1)
plot(preTracks(:,:,1), preTracks(:,:,2));
figure(2);
plot(postTracks(:,:,1), postTracks(:,:,2));
%%

x1 = preTracks(:,:,1);
y1 = preTracks(:,:,2);

x2 = preTracks(:,:,3);
y2 = preTracks(:,:,4);

xd = x2 - x1;
yd = y2 - y1;

[preTheta rho] = cart2pol(xd,yd);
%preTheta = fixShortNanGaps(preTheta,10);
preNest(:,:,1) = fixShortNanGaps(preNest(:,:,1),10);
preNest(:,:,2) = fixShortNanGaps(preNest(:,:,2),10);

preNest(:,:,3) = preTheta;

%%
% Repeat for post-data
x1 = postTracks(:,:,1);
y1 = postTracks(:,:,2);

x2 = postTracks(:,:,3);
y2 = postTracks(:,:,4);

xd = x2 - x1;
yd = y2 - y1;

[postTheta rho] = cart2pol(xd,yd);
postTheta = fixShortNanGaps(postTheta,10);
postNest(:,:,1) = fixShortNanGaps(postNest(:,:,1),10);
postNest(:,:,2) = fixShortNanGaps(postNest(:,:,2),10);
% Append to appropriate scaled data
postNest(:,:,3) = postTheta;
%% Rename variables
nestTracks = preNest;
brood = broodPre;
%cd(uigetdir());
save('allData.mat', 'preNest', 'postNest', 'preFeederData', 'preFeederTimes', 'postFeederData', 'postFeederTimes','tags', 'mass', 'orTagTreat', 'distMatPre', 'distMatPost', 'broodPre', 'broodPost', 'nestBackImagePost', 'nestBackImagePre');