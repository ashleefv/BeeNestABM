% simulation function
% Output:
%   one simulated set of nestSimulationData for frames corresponding to 
%   the experimental data set.  
function simulation
% function nestSimulationData = simulation
clear; close all; clc

%% Step 1: Defining the initial state & parameters
tic
% Initial state of the agents
colonyNumber = 1; % 4 colonies numbered 1 - 4
colony = load(['data\allDataCol' num2str(colonyNumber) '.mat']);
brood = relabelBroodObject(colony.broodPre);
numFrames = size(colony.preNest,1);
numBees = size(colony.preNest,2);
tags = colony.orTagTreat;
% number of time points
totalTimePoints = 250;
vis = 1;

nestSimulationData = zeros(totalTimePoints,numBees,6); 
estimatedData = load(['data\Essential_Info_Col_' num2str(colonyNumber) '.mat']);

exposure_state = 'pre'; % 'pre' or 'post'



% The initial state for X & Y coordinates is sampled one time from a normal 
% distribution fitted to the observed data for the colony specified by colonyNumber
% All repeated simulations use the same initial coordinates unless
% Estimate_TransProb is run again
nestSimulationData(1,:,1) = estimatedData.X_pos_initial; %initialize as random finite value
nestSimulationData(1,:,2) = estimatedData.Y_pos_initial; %initialize as random finite value
% The orientation of the bee bodies is not used, but preserved for keeping
% the state variable structure the same as that for the experimental data
nestSimulationData(1,:,3) = zeros(size(estimatedData.Y_pos_initial)); % 

% The initial state for velocity is randomly sampled from the Weibul distribution
nestSimulationData(1,:,4) = zeros(size(estimatedData.Y_pos_initial));% velocity
velocityPDF = estimatedData.PDF_Dist_Speed_Log;
 % Velocity sampling from the velocityPDF distribution
initialVelocity = 10.^random(velocityPDF,numBees,1); % m/s
nestSimulationData(1,:,4) = initialVelocity;

% The initial state for activity is based on the proportion of the entire
% experiment duration that bees in a cohort for the specified colony 
% were active or inactive
initialActivity = ones(size(estimatedData.Y_pos_initial));
% Initial activity for cohorts (columns: 1 queen, 2 untreated, 3 control
% sucrose, 4 low dose, 5 high dose)
% first row is inactive probability 
% second row is active probability
if strcmp(exposure_state, 'pre')
   estimatedActiveProb = estimatedData.Activity_Prob_Dist_Pre(2,:);
elseif strcmp(exposure_state, 'post')
   estimatedActiveProb = estimatedData.Activity_Prob_Dist_Post(2,:);
end
randomActive = rand(size(initialActivity));
activeProbBees = zeros(size(initialActivity));
activeProbBees = assignCohortParameters(estimatedActiveProb,activeProbBees,tags);
initialActivity(activeProbBees<randomActive) = 0;
nestSimulationData(1,:,5) = initialActivity;% activity
% Frequency of time resolution 
frequency = 2; % Hz
% Time resolution per step
dt = 1/frequency; % time step, 1/(2 Hz) = 0.5 s
T = size(nestSimulationData,1)-1; % number of time steps
final_time = dt*T; % s

% Static environment for the given colony
%   brood: n x 3 matrix, where 1st and 2nd columns are x and y coordinates,
%   and 3rd column contains labels for element type:
%   1 = brood (eggs, larvae, and pupae)
%   2 = full food pots
%   3 = empty pots + wax cover
broodPosition = brood(brood(:,3)==1,1:2);
fullFoodPosition = brood(brood(:,3)==2,1:2); 
emptyFoodPosition = brood(brood(:,3)==3,1:2); 

initialAngle = 2*pi*rand(size(estimatedData.Y_pos_initial));
nestSimulationData(1,:,6) = initialAngle;% angle 


% Transition probability matrix estimated from colony 1 data
% from Estimate_TransProb.m and the corresponding output
% Essential_Info_Col_*.mat
% transProb_Pre is a 2x1 cell: first entry is matrix without bump and
% second entry is matrix with bump
% A = active (activity = 1)
% I = inactive (activity = 0)
% [II IA % Probability that I stays I and that I transitions to A
%  AI AA] % Probability that A transitions to I and that A stays A
% each bee gets its own AtoI and ItoA _Bumped and _Unbumped
    AtoI_Unbumped = zeros(size(initialActivity));
    ItoA_Unbumped = zeros(size(initialActivity));
    AtoI_Bumped = zeros(size(initialActivity));
    ItoA_Bumped = zeros(size(initialActivity));
if strcmp(exposure_state, 'pre')
    % all bees get the same bump values as there is only one cohort
    estimatedTransProb = estimatedData.transProb_Pre;
    AtoI_Unbumped(:) = estimatedTransProb{1,1}(2,1);
    ItoA_Unbumped(:) = estimatedTransProb{1,1}(1,2);
    AtoI_Bumped(:) = estimatedTransProb{2,1}(2,1);
    ItoA_Bumped(:) = estimatedTransProb{2,1}(1,2);
elseif strcmp(exposure_state, 'post')
    % there are 4 cohorts of bump values as the queen is part of the
    % untreated group
    estimatedTransProb = estimatedData.transProb_Post;
    for tagNumber = 1:4 % untreated, sucrose control, low dose, high dose
        AtoI_Unbumped(tags == tagNumber-1) = estimatedTransProb{1,tagNumber}(2,1);
        ItoA_Unbumped(tags == tagNumber-1) = estimatedTransProb{1,tagNumber}(1,2);
        AtoI_Bumped(tags == tagNumber-1) = estimatedTransProb{2,tagNumber}(2,1);
        ItoA_Bumped(tags == tagNumber-1) = estimatedTransProb{2,tagNumber}(1,2);    
    end
    % the queeen
        AtoI_Unbumped(1) = estimatedTransProb{1,1}(2,1);
        ItoA_Unbumped(1) = estimatedTransProb{1,1}(1,2);
        AtoI_Bumped(1) = estimatedTransProb{2,1}(2,1);
        ItoA_Bumped(1) = estimatedTransProb{2,1}(1,2);        
end
toc
'end of initialization'
% 
%% Step 2: Covering each time step
% Time loop
tic
for timestep=1:totalTimePoints 
	% What happens at each time step?
    % Update agents
    %% Step 3: Covering each agent
    % Agent loop
    
        %% calculate whether each bee is on the nest structure or not for current timestep
        [scrap onNestCurrentFrame] = calculatePortionOfTimeOnNest(nestSimulationData(timestep,:,:), brood, 0.01);
        onNestCurrentFrame = logical(onNestCurrentFrame); %Convert from double to logical
        
        %Write to memory
        onNest(timestep,:) = onNestCurrentFrame;
        
    nestSimulationData(timestep+1,:,:) = rules(dt,nestSimulationData(timestep,:,:),...
        broodPosition,emptyFoodPosition,fullFoodPosition,...
        AtoI_Unbumped,ItoA_Unbumped,AtoI_Bumped,ItoA_Bumped,velocityPDF,exposure_state,tags);
end
toc
'end of time loop'

%update nestSimulationData with onNest
    nestSimulationData = nestSimulationData(1:totalTimePoints,:,:);
    nestSimulationData(:,:,7) = onNest;
        cols = nan(numBees,3);
    
    for  zz = 1:numBees
        if tags(zz) == 0
            cols(zz,:) = [0 1 0];
        elseif tags(zz) == 1
            cols(zz,:) = [0 1 0];
        elseif tags(zz) == 2
            cols(zz,:) = [0 0 1];
        elseif tags(zz) == 3
            cols(zz,:) = [1 0 0];
        end
    end
    
%% Step 5: Final processing
% Outputs and final processing
% avg_x = sum(Agents.xstate)/length(Agents.xstate);
% avg_y = sum(Agents.ystate)/length(Agents.ystate);
% 
tic
if vis == 1
    figure(1)
    % % for video
    % myVideo = VideoWriter('preAttracted.avi','uncompressed avi');
    % open(myVideo);
    for timestep = 1:totalTimePoints 
        %%
        plotCoordinatesAndBrood(nestSimulationData, brood, timestep,cols);
    % % for video
    %     frame = getframe(gcf);
    %     writeVideo(myVideo,frame);
        drawnow
    end
    toc
    'end of visualization'
    % % for video
    % close(myVideo);
end
[means distributions]=calculateSummaryStatistics(nestSimulationData,brood)

end