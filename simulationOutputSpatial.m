function nestSimulationData = simulationOutput(colony, estimatedData, exposure_state, totalTimePoints, motionParams, vis, varargin)
    %Function to take colony data and return simulated data
    %
    %Inputs:
    %   colony: same object as in "simulation.m", but now loaded externally to
    %       avoid repetition
    %
    %   estimatedData: bump transition probabilities calculated directly for
    %       various data sets
    %
    %   exposure_state: what kind of simulation to run? Options are 'pre', 'post',
    %       'postSpont', and 'postInt'. "Int" and "Spont" stand for driving only
    %       bumped and spontaneous parameter changes, respectively.
    %
    %   totalTimePoints: number of time points
    %
    %   vis: Boolean if true plot visualization (slower)
    %
    %   motionParams: specified motion parameter estimated:
    %       Rows(n = 4) are AtoI_bumped, ItoA_bumped, AtoI_unbumped, and
    %       ItoA_unbumped
    %       Columns (n = 4) are for 4 groups (untreated, control, 0.1 ng, and
    %       1.0 ng)
    %       Sheets (n = 2) are for on nest and off nest-specific parameters
    %
    %Output:
    %   nestSimulationData
    %   Same as out of "calculateActivityMatrix.m", except 6th sheet is a
    %   boolean of whether the bees are on the nest or not
    %
    if nargin > 6
        coefficients = varargin{1};
        param_estim = 1;
    else
        param_estim = 0;
    end
    %% Step 1: Defining the initial state & parameters
    if vis == 1
        tic
    end
    % Initial state of the agents
    brood = relabelBroodObject(colony.broodPre);
    numFrames = size(colony.preNest,1);
    numBees = size(colony.preNest,2);
    tags = colony.orTagTreat;
    % number of time points
    %totalTimePoints = 250;
    
    nestSimulationData = zeros(totalTimePoints,numBees,6);
    onNest = nan(totalTimePoints,numBees);
    %exposure_state = 'post'; % 'pre' or 'post' or
    
    
    
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
    
    % ON NEST PARAMETERS
    AtoI_Unbumped_OnNest = zeros(size(initialActivity));
    ItoA_Unbumped_OnNest = zeros(size(initialActivity));
    AtoI_Bumped_OnNest = zeros(size(initialActivity));
    ItoA_Bumped_OnNest = zeros(size(initialActivity));
    
%     if strcmp(exposure_state, 'pre')
%         % all bees get the same bump values as there is only one cohort
%         estimatedTransProb = estimatedData.transProb_Pre;
%         AtoI_Unbumped_OnNest(:) = motionParams(1,1);
%         ItoA_Unbumped_OnNest(:) = motionParams(2,1);
%         AtoI_Bumped_OnNest(:) = motionParams(3,1);
%         ItoA_Bumped_OnNest(:) = motionParams(4,1);
%         %         % the queeen
%         %         AtoI_Unbumped_OnNest(1) = estimatedTransProb{1,1}(2,1);
%         %         ItoA_Unbumped_OnNest(1) = estimatedTransProb{1,1}(1,2);
%         %         AtoI_Bumped_OnNest(1) = estimatedTransProb{2,1}(2,1);
%         %         ItoA_Bumped_OnNest(1) = estimatedTransProb{2,1}(1,2);
        
        
        
        %FULL MODEL
   % elseif strcmp(exposure_state, 'post')
        % there are 4 cohorts of bump values as the queen is part of the
        % untreated group
        estimatedTransProb = estimatedData.transProb_Post;
        for tagNumber = 1:4 % untreated, sucrose control, low dose, high dose
            AtoI_Unbumped_OnNest(tags == tagNumber-1) = motionParams(1,tagNumber);
            ItoA_Unbumped_OnNest(tags == tagNumber-1) = motionParams(2,tagNumber);
            AtoI_Bumped_OnNest(tags == tagNumber-1) = motionParams(3,tagNumber);
            ItoA_Bumped_OnNest(tags == tagNumber-1) = motionParams(4,tagNumber);
        end
        %         % the queeen
        %         AtoI_Unbumped_OnNest(1) = estimatedTransProb{1,1}(2,1);
        %         ItoA_Unbumped_OnNest(1) = estimatedTransProb{1,1}(1,2);
        %         AtoI_Bumped_OnNest(1) = estimatedTransProb{2,1}(2,1);
        %         ItoA_Bumped_OnNest(1) = estimatedTransProb{2,1}(1,2);
        
  %  end
    
    
    % OFF NEST PARAMETERS
    AtoI_Unbumped_OffNest = zeros(size(initialActivity));
    ItoA_Unbumped_OffNest = zeros(size(initialActivity));
    AtoI_Bumped_OffNest = zeros(size(initialActivity));
    ItoA_Bumped_OffNest = zeros(size(initialActivity));
    
%     if strcmp(exposure_state, 'pre')
%         % all bees get the same bump values as there is only one cohort
%         estimatedTransProb = estimatedData.transProb_Pre;
%         AtoI_Unbumped_OffNest(:) = motionParams(1,1,1); %Take data from first "sheets" of motionParams
%         ItoA_Unbumped_OffNest(:) = motionParams(2,1,1);
%         AtoI_Bumped_OffNest(:) = motionParams(3,1,1);
%         ItoA_Bumped_OffNest(:) = motionParams(4,1,1);
%         %         % the queeen
%         %         AtoI_Unbumped_OffNest(1) = estimatedTransProb{1,1}(2,1);
%         %         ItoA_Unbumped_OffNest(1) = estimatedTransProb{1,1}(1,2);
%         %         AtoI_Bumped_OffNest(1) = estimatedTransProb{2,1}(2,1);
%         %         ItoA_Bumped_OffNest(1) = estimatedTransProb{2,1}(1,2);
%         
        
        
        %FULL MODEL
%    elseif strcmp(exposure_state, 'post')
        % there are 4 cohorts of bump values as the queen is part of the
        % untreated group
        estimatedTransProb = estimatedData.transProb_Post;
        for tagNumber = 1:4 % untreated, sucrose control, low dose, high dose
            AtoI_Unbumped_OffNest(tags == tagNumber-1) = motionParams(1,tagNumber,2);
            ItoA_Unbumped_OffNest(tags == tagNumber-1) = motionParams(2,tagNumber,2);
            AtoI_Bumped_OffNest(tags == tagNumber-1) = motionParams(3,tagNumber,2);
            ItoA_Bumped_OffNest(tags == tagNumber-1) = motionParams(4,tagNumber,2);
        end
        %         % the queeen
        %         AtoI_Unbumped_OffNest(1) = estimatedTransProb{1,1}(2,1);
        %         ItoA_Unbumped_OffNest(1) = estimatedTransProb{1,1}(1,2);
        %         AtoI_Bumped_OffNest(1) = estimatedTransProb{2,1}(2,1);
        %         ItoA_Bumped_OffNest(1) = estimatedTransProb{2,1}(1,2);
        
%    end
    
    if vis == 1
        toc
        'end of initialization'
    end
    %
    %% Step 2: Covering each time step
    % Time loop
    if vis == 1
        tic
    end
    for timestep=1:totalTimePoints
        % What happens at each time step?
        % Update agents
        
        %% calculate whether each bee is on the nest structure or not for current timestep
        [scrap onNestCurrentFrame] = calculatePortionOfTimeOnNest(nestSimulationData(timestep,:,:), brood, 0.01);
        onNestCurrentFrame = logical(onNestCurrentFrame); %Convert from double to logical
        
        %Write to memory
        onNest(timestep,:) = onNestCurrentFrame;
        %% Assign transition parameters based on spatial location
        
        %Initalize empty transition probability vectors
        AtoI_Unbumped = zeros(size(initialActivity));
        ItoA_Unbumped = zeros(size(initialActivity));
        AtoI_Bumped = zeros(size(initialActivity));
        ItoA_Bumped = zeros(size(initialActivity));
        
        %Assign values to transition prob vectors based on current spatial state
        %for each bee
        AtoI_Unbumped(onNestCurrentFrame) = AtoI_Unbumped_OnNest(onNestCurrentFrame); %For bees that are on the nest, assigns values from "OnNest" object to appropriate slots
        AtoI_Unbumped(~onNestCurrentFrame) = AtoI_Unbumped_OffNest(~onNestCurrentFrame); %For bee off nest, assign values from "OffNest" object, ...
        
        ItoA_Unbumped(onNestCurrentFrame) = ItoA_Unbumped_OnNest(onNestCurrentFrame);
        ItoA_Unbumped(~onNestCurrentFrame) = ItoA_Unbumped_OffNest(~onNestCurrentFrame);
        
        AtoI_Bumped(onNestCurrentFrame) = AtoI_Bumped_OnNest(onNestCurrentFrame);
        AtoI_Bumped(~onNestCurrentFrame) = AtoI_Bumped_OffNest(~onNestCurrentFrame);
        
        ItoA_Bumped(onNestCurrentFrame) = ItoA_Bumped_OnNest(onNestCurrentFrame);
        ItoA_Bumped(~onNestCurrentFrame) = ItoA_Bumped_OffNest(~onNestCurrentFrame);
        %% Step 3: Covering each agent
        % Agent loop
        if param_estim == 0
            nestSimulationData(timestep+1,:,:) = rules(dt,nestSimulationData(timestep,:,:),...
                broodPosition,emptyFoodPosition,fullFoodPosition,...
                AtoI_Unbumped,ItoA_Unbumped,AtoI_Bumped,ItoA_Bumped,velocityPDF,exposure_state,tags);
        else
            nestSimulationData(timestep+1,:,:) = rules(dt,nestSimulationData(timestep,:,:),...
                broodPosition,emptyFoodPosition,fullFoodPosition,...
                AtoI_Unbumped,ItoA_Unbumped,AtoI_Bumped,ItoA_Bumped,velocityPDF,exposure_state,tags,coefficients);
        end
    end
    if vis == 1
        toc
        'end of time loop'
    end
    %% Step 5: Final processing
    % Outputs and final processing
    % avg_x = sum(Agents.xstate)/length(Agents.xstate);
    % avg_y = sum(Agents.ystate)/length(Agents.ystate);
    %
    
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
    if vis == 1
        tic
        % % for video
        % myVideo = VideoWriter('preAttracted.avi','uncompressed avi');
        % open(myVideo);
        figure(1)
        for timestep = 1:totalTimePoints
            %%
            
            
            plotCoordinatesAndBrood(nestSimulationData, brood, timestep, cols);
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
end