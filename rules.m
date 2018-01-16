function output = rules(dt,currentState,broodPosition,...
    emptyFoodPosition,fullFoodPosition,...
    AtoI_Unbumped,ItoA_Unbumped,AtoI_Bumped,ItoA_Bumped,velocityPDF,exposure_state,tags,varargin)


%% Rules

%% Step 4: Updating agent i at time t

%% Parse currentState matrix into constituent variables
position(:,1:2) = currentState(1,:,1:2); % locations in m
currentVelocity(:,1) = currentState(1,:,4);
currentActivity(:,1) = currentState(1,:,5);
currentAngle(:,1) = currentState(1,:,6);

numBees = length(position);

% Nest chamber size
nestMaxX = 0.25; % m
nestMaxY = 0.2; % m
sensoryRadius = 10000; %How many body lengths can the bee pay attention to stimuli within?>
updatedState = currentState; % initialize

%% Parameters
% Each bee gets its own value, though most will be the same for all or at
% least for all within a treatment cohort

% COHORT parameters are vectors of [queen untreated sucroseControl lowDose highDose]
% for each parameters values that will later be assigned via the tags on
% the bees
if strcmp(exposure_state,'pre')
    COHORTvelocityPerturbationAlwaysAccept = [0.2, 0.2, 0.2, 0.2, 0.2]; 
    % case A & B: no attraction to nestmates, all nest structures lumped,
    % optimization using all summary statistics, pre-exposure
    COHORTenvironmentalStimuliWeight = [0.06675, 0.06675, 0.06675, 0.06675, 0.06675];
    % case C & D: no attraction to nestmates, all nest structures lumped,
    % optimization using all summary statistics EXCEPT interactionRate, pre-exposure 
%     COHORTenvironmentalStimuliWeight = [0.09056, 0.09056, 0.09056, 0.09056, 0.09056];
    COHORTlambdaBrood = [1, 1, 1, 1, 1];
    COHORTlambdaFullFood = [0, 0, 0, 0, 0];
    COHORTlambdaEmptyFood = [0, 0, 0, 0, 0];
    COHORTperturbationAngle = [pi/4, pi/4, pi/4, pi/4, pi/4]; 
    COHORTvelocityPerturbationMightAcceptProb = [0.1, 0.1, 0.1, 0.1, 0.1]; 
    COHORTBeeBodyThreshold = [0.01, 0.01, 0.01, 0.01, 0.01]; % m
    COHORTcutoffRadius = sensoryRadius*COHORTBeeBodyThreshold;
elseif strcmp(exposure_state,'post')
    COHORTvelocityPerturbationAlwaysAccept = [0.2, 0.2, 0.2, 0.2, 0.2]; 
    % case A: estimated coef 4 and 5 distinctly using all summary statistics
%     COHORTenvironmentalStimuliWeight = [0.06675, 0.06675, 0.06675, 0.05014, 0.03489];
    % case B: estimated coef 5 distinctly and set coef 4 = to
    % control/untreated, using all summary statistics
    COHORTenvironmentalStimuliWeight = [0.06675, 0.06675, 0.06675, 0.06675, 0.0317];
    % case C: estimated coef 4 and 5 distinctly, using all summary statistics EXCEPT interactionRate
%     COHORTenvironmentalStimuliWeight = [0.09056, 0.09056, 0.09056, 0.05014, 0.03489]; %NOT OPTIMIZED
    % case D: estimated coef 5 distinctly and set coef 4 = to
    % control/untreated, using all summary statistics EXCEPT interactionRate
%     COHORTenvironmentalStimuliWeight = [0.09056, 0.09056, 0.09056, 0.09056, 0.077875];
    % case E: same as B except that the objective function only considers
    % high dose treated bees instead of the entire population
%     COHORTenvironmentalStimuliWeight = [0.06675, 0.06675, 0.06675, 0.06675, 0.0519];
    COHORTlambdaBrood = [1, 1, 1, 1, 1];
    COHORTlambdaFullFood = [0, 0, 0, 0, 0];
    COHORTlambdaEmptyFood = [0, 0, 0, 0, 0];
    COHORTperturbationAngle = [pi/4, pi/4, pi/4, pi/4, pi/4]; % pi/4 = 45 degrees
    COHORTvelocityPerturbationMightAcceptProb = [0.1, 0.1, 0.1, 0.1, 0.1]; 
    COHORTBeeBodyThreshold = [0.01, 0.01, 0.01, 0.01, 0.01]; % m
    COHORTcutoffRadius = sensoryRadius*COHORTBeeBodyThreshold;
end

if nargin > 12
    coefficients = varargin{1};
    paramcase = varargin{2};
    if paramcase == 0
%         paramcase = 0; pre-exposure Cases A, B, C, & D       
        COHORTenvironmentalStimuliWeight = coefficients(1).*ones(1,5);
    elseif paramcase == 1
%        paramcase = 1;
        COHORTenvironmentalStimuliWeight = coefficients(1).*ones(1,5);
        COHORTlambdaBrood = coefficients(2).*ones(1,5);
        COHORTlambdaFullFood = coefficients(3).*zeros(1,5);
        COHORTlambdaEmptyFood = coefficients(4).*zeros(1,5);
    elseif paramcase == 2
%         paramcase = 2; post-exposure estimate
        COHORTenvironmentalStimuliWeight(4:5) = [coefficients(1) coefficients(2)];     

    elseif paramcase == 3
        % case B & D: post-exposure estimate
        % coef 5 distinctly and set coef 4 = to control/untreated
        COHORTenvironmentalStimuliWeight(5) = coefficients(1); 
    elseif paramcase == 4
        COHORTenvironmentalStimuliWeight = coefficients(1).*ones(1,5); % untreated
        COHORTenvironmentalStimuliWeight(5) = coefficients(2); % treated
    end
    if nargin > 14 % editable size of the domain
        nestSize = varargin{3};
        nestMaxX = nestSize(1); % cm
        nestMaxY = nestSize(2); % cm
    end
end
% Initialize all parameters to zero for each bee
velocityPerturbationAlwaysAccept = zeros(numBees,1); 
lambdaBrood = zeros(numBees,1);
lambdaFullFood = zeros(numBees,1);
lambdaEmptyFood = zeros(numBees,1);
lambdaBees = zeros(numBees,1);
environmentalStimuliWeight = zeros(numBees,1);
perturbationAngle = zeros(numBees,1); 
velocityPerturbationMightAcceptProb = zeros(numBees,1); 
BeeBodyThreshold = zeros(numBees,1); 
cutoffRadius = zeros(numBees,1);

% Set parameter values

% Fraction within +/- currentVelocity in which new velocity samples are 
% always acceptable
velocityPerturbationAlwaysAccept = assignCohortParameters(COHORTvelocityPerturbationAlwaysAccept,velocityPerturbationAlwaysAccept,tags);

% Lambda weights for strength of attraction to nestmates (other bees),
% brood, empty food pots, and full food pots
lambdaBrood = assignCohortParameters(COHORTlambdaBrood,lambdaBrood,tags);
lambdaFullFood = assignCohortParameters(COHORTlambdaFullFood,lambdaFullFood,tags);
lambdaEmptyFood = assignCohortParameters(COHORTlambdaEmptyFood,lambdaEmptyFood,tags);
lambdaBees(:) = 0; %1-(lambdaBrood+lambdaFullFood+lambdaEmptyFood);

environmentalStimuliWeight = assignCohortParameters(COHORTenvironmentalStimuliWeight,environmentalStimuliWeight,tags);

perturbationAngle = assignCohortParameters(COHORTperturbationAngle,perturbationAngle,tags);

% This fraction of samples outside of the velocityPerturbationAlwaysAccept 
% window will be allowed to deviate to a proposedVelocity from the sampled 
% velocityPDF. Those outside this fraction and also outside the 
% velocityPerturbationAlwaysAccept window will remain at their currentVelocity
velocityPerturbationMightAcceptProb = assignCohortParameters(COHORTvelocityPerturbationMightAcceptProb,velocityPerturbationMightAcceptProb,tags);

% Spatial resolution that constitutes the distance at which bees are
% considered to overlap or bump into each other
BeeBodyThreshold = assignCohortParameters(COHORTBeeBodyThreshold,BeeBodyThreshold,tags);

% Cutoff radius for neighborhood sensing influence
cutoffRadius = assignCohortParameters(COHORTcutoffRadius,cutoffRadius,tags);

%% Update pairwise distances, Delta X & Y, and angles of resultant vectors for the current time step
% Input:
%   current position vectors
% Output:
%   updated pairwise distance calculations

currentDistanceToBees = pdist2(position,position); % m
currentDistanceToBrood = pdist2(position,broodPosition);% m
currentDistanceToFullFood = pdist2(position,fullFoodPosition);% m
currentDistanceToEmptyFood = pdist2(position,emptyFoodPosition);% m

% Input:
%   current position vectors
% Output:
%   updated difference in x and y coordinates to objects relative to bee positions
[DeltaX_Bees,DeltaY_Bees] = Delta_Obj(position,position);% m
[DeltaX_Brood,DeltaY_Brood] = Delta_Obj(position,broodPosition);% m
[DeltaX_FullFood,DeltaY_FullFood] = Delta_Obj(position,fullFoodPosition);% m
[DeltaX_EmptyFood,DeltaY_EmptyFood] = Delta_Obj(position,emptyFoodPosition);% m

% Input:
%   current position vectors, DeltaX and DeltaY vectors
% Output:
%   updated angles between bees and the angle of the resultant vector of 
%   the distance from all attractors of a given category (other nestmates,
%   brood, full full pots, and empty food pots)
angleBees = angleObj(currentDistanceToBees,DeltaX_Bees,DeltaY_Bees,cutoffRadius);% radians
angleBrood = angleObj(currentDistanceToBrood,DeltaX_Brood,DeltaY_Brood,cutoffRadius);% radians
angleFullFood = angleObj(currentDistanceToFullFood,DeltaX_FullFood,DeltaY_FullFood,cutoffRadius);% radians
angleEmptyFood = angleObj(currentDistanceToEmptyFood,DeltaX_EmptyFood,DeltaY_EmptyFood,cutoffRadius);% radians

%% Check Bump
% If the bee is in the nest chamber & less than BeeBodyThreshold from
% another bee, those two bees are close enough to be in contact. Excludes 
% self-contact for bees in the nest chamber and excludes bees outside the
% chamber

% Input: 
%   default transition probability vector
%   current position of each bee and the threshold distance between a pair
%   of bees
% Output:
%   updated transition probability vector
bumpedStorage = bump(BeeBodyThreshold,currentDistanceToBees);

%% Define Transition Probability Vector 
% Probabilities for state transitions from active to inactive and inactive
% to active either with or without being bumped
transitionVector = zeros(size(bumpedStorage));
transitionVector(bumpedStorage == 0 & currentActivity == 1) = AtoI_Unbumped(bumpedStorage == 0 & currentActivity == 1);
transitionVector(bumpedStorage == 0 & currentActivity == 0) = ItoA_Unbumped(bumpedStorage == 0 & currentActivity == 0);
transitionVector(bumpedStorage == 1 & currentActivity == 1) = AtoI_Bumped(bumpedStorage == 1 & currentActivity == 1);
transitionVector(bumpedStorage == 1 & currentActivity == 0) = ItoA_Bumped(bumpedStorage == 1 & currentActivity == 0);


%% Switch States
% Use the output transition probability vector from the check bump rule to
% determine whether the active agents will stay active or become inactive
% and whether the inactive agents will stay inactive or become active

% Input:
%   Current activitivity indicator
% Output: 
%   Updated activitiy indicator at the current time step

% Check that the agent is in the domain. NaN cannot switch states
randomTransition = rand(size(currentActivity)); % sample from uniform distribution [0,1]
updatedActivity = currentActivity; % initialize
% if the agent is in the domain and the random transition probability is 
% less than the threshold in the transition vector, then the activity state 
% switch happens
% logical index for which bees switch
switch_idx = (isfinite(position(:,1))) & (randomTransition < transitionVector); 
updatedActivity(switch_idx) = 1-currentActivity(switch_idx); % switch activity state index
 
%% Reenter the Nest Chamber
% This effect is neglected in this version of the model. Only bees that
% stay in the nest are considered in this model.
% % A bee may have left the nest chamber or might have been temporily
% % visually blocked from view of the camera for a window of time. Both of 
% % these scenarios result in NaN for the bee's position for that period. We
% % need to consider rules for such bees' departures and returns, or we need 
% % to cull those bees from the data set and only consider bees tracked 
% % explicitly during the entire collection time frame or those that can be 
% % reconstructured after brief orientation-based departures from the visible field

%% Move
% If the bee is active and in the nest, then compute the movement rules

% Input:
%   Updated activity state
%   Current position 
%   Current velocityPDF 
%   currentDistanceToBrood
%   currentDistanceToFullFood
%   currentDistanceToEmptyFood
%   currentDistanceToBees
%   weight parameters: lambda values and environmentalStimuliWeight
% Output:
%   Updated position

% logical index for which bees move
move_idx = (updatedActivity == 1) & (isfinite(position(:,1)));

% initialize such that if the sampled velocity not accepted, then the value 
% will default to the currentVelocity
updatedVelocity = currentVelocity;

% Velocity sampling from the velocityPDF distribution
proposedVelocity = 10.^random(velocityPDF,numBees,1); % m/s

% 100% acceptance within +/- velocityPerturbationAlwaysAccept of the
% currentVelocity
% idx is the index for the bees that satisfy the criteria
velocityPerturbationAlwaysAcceptAccept_idx = ...
    (proposedVelocity <= (1+velocityPerturbationAlwaysAccept).*currentVelocity) ...
    & (proposedVelocity >= (1-velocityPerturbationAlwaysAccept).*currentVelocity);

updatedVelocity(velocityPerturbationAlwaysAcceptAccept_idx) = proposedVelocity(velocityPerturbationAlwaysAcceptAccept_idx);

% velocityPerturbationMightAcceptProb acceptance rate outside of the window around
% the currentVelocity
% consider all that were not already accepted
velocityPerturbationMightAccept_idx = 1-velocityPerturbationAlwaysAcceptAccept_idx; 
velocityRandom = rand(size(proposedVelocity));
% set the velocityRandom for the already accepted where they will always
% fail the probability check below so that they are not manipulated again
velocityRandom(~velocityPerturbationMightAccept_idx) = 1; 
% idx is the index for the bees that satisfy the criteria
velocityPerturbationMightAccept_idx  = (velocityRandom < velocityPerturbationMightAcceptProb);
updatedVelocity(velocityPerturbationMightAccept_idx) = proposedVelocity(velocityPerturbationMightAccept_idx);

% find the bees that are moving that switched their activity this
% iteration 
switchToActive_idx = switch_idx & move_idx; 
% theese bees must now accept their sampled velocity rather than
% staying at zero with a high probability
updatedVelocity(switchToActive_idx) = proposedVelocity(switchToActive_idx);

% for bees that are not moving (inactive), set velocity to 0
updatedVelocity(~move_idx) = 0;

stepsize = updatedVelocity'.*dt;

% Random walk angle
%randomWalkAngle = currentAngle - pi/4+ (currentAngle + pi/4 -
%(currentAngle - pi/4) )*rand(1,numBees); % random perturbation within 
% +/- 45 degrees from currentAngle
randomWalkAngle = currentAngle' - perturbationAngle'+ 2*perturbationAngle'.*rand(1,numBees); %
% simplified math of the eqn above with generalized +/- perturbationAngle
% from currentAngle
% randomWalkAngle = 2*pi*rand(1,numBees); %random directions
% Net angle from environmental attractors
environWeights = [lambdaBees'; lambdaBrood'; lambdaFullFood'; lambdaEmptyFood'];
environAngles = [angleBees; angleBrood; angleFullFood; angleEmptyFood];
netEnvironAngle = angleMean(environWeights,environAngles);

% Net angle from random walk and environmental attractors
weights = [environmentalStimuliWeight'; 1-environmentalStimuliWeight'];
angles = [netEnvironAngle; randomWalkAngle];
updatedAngle = angleMean(weights,angles);

poX = move_idx'.*stepsize.*cumsum([zeros(1,numBees); cos(updatedAngle)]); % before and after for each bee
poY = move_idx'.*stepsize.*cumsum([zeros(1,numBees); sin(updatedAngle)]);
% poX(ang == NaN) = 0;
% poY(ang == NaN) = 0;
moveDistance = [poX', poY']; % column order: before & after x then before & after y starting from origin
% beforeAfterPosition tells where the coordinates of the bees before and
% after move
beforeAfterPosition = [position(:,1) position(:,1) position(:, 2) position(:,2)]+moveDistance;
% % Plotting for debugging: see where the bees go after a single time step
% plot(beforeAfterPosition(:,1:2)',beforeAfterPosition(:,3:4)','x-')
% hold on
% plot(position(:,1)',position(:,2)','o')
% hold off
% legend('1','2','3')
% figure(2)
% plot(poX, poY,'x-')
% legend('1','2','3')

updatedPosition = [position(:,1)+moveDistance(:,2) position(:,2)+moveDistance(:,4)];
%% Check if hit the wall after movement
% Determine if the bee hit the wall. If it hits the wall, truncate
% it final position to be the wall in that direction.
updatedPosition(updatedPosition(:,1)>nestMaxX, 1) = nestMaxX;
updatedPosition(updatedPosition(:,1)<0, 1) = 0;
updatedPosition(updatedPosition(:,2)>nestMaxY,2) = nestMaxY;
updatedPosition(updatedPosition(:,2)<0,2) = 0;

%% Update state
% Store updatedPosition into updatedState
updatedState(1,:,1:2) = updatedPosition; % PLACEHOLDER for cellarray? of updatedState and desired output for comparison, visualization
updatedState(1,:,4) = updatedVelocity;
updatedState(1,:,5) = updatedActivity;
updatedState(1,:,6) = updatedAngle';

%% Calculate any local variables that need to be passed as output

%% Output
output = updatedState;
