
%% Step 1: Defining the initial state & parameters
% Initial state of the agents
colonyNumber = 1; % 4 colonies numbered 1 - 4
colony = load(['allDataCol1_app.mat']);
brood = relabelBroodObject_app(colony.broodPre);
numFrames = size(colony.preNest,1);
numBees = size(colony.preNest,2);
tags = colony.orTagTreat;
nestSimulationData = zeros(numFrames,numBees,7); 
estimatedData = load(['Essential_Info_Col_1_app.mat']);

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
   estimatedActiveProb = estimatedData.Activity_Prob_Dist_Pre(2,:);
randomActive = rand(size(initialActivity));
activeProbBees = zeros(size(initialActivity));
activeProbBees(tags == 0) = estimatedActiveProb(2); % untreated
activeProbBees(1) = estimatedActiveProb(1); % queen
activeProbBees(tags == 1) = estimatedActiveProb(3); % sucrose control
activeProbBees(tags == 2) = estimatedActiveProb(4); % low dose
activeProbBees(tags == 3) = estimatedActiveProb(5); % high dose
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
    % all bees get the same bump values as there is only one cohort
    estimatedTransProb = estimatedData.transProb_Pre;
    AtoI_Unbumped(:) = estimatedTransProb{1,1}(2,1);
    ItoA_Unbumped(:) = estimatedTransProb{1,1}(1,2);
    AtoI_Bumped(:) = estimatedTransProb{2,1}(2,1);
    ItoA_Bumped(:) = estimatedTransProb{2,1}(1,2);

handles.totalTimePoints = totalTimePoints;    
handles.brood = brood;

% Update handles structure
guidata(hObject, handles);    