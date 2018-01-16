% simulation function
% Output:
%   one simulated set of nestSimulationData for frames corresponding to 
%   the experimental data set.  
function simulation_attraction(hObject, eventdata, handles)
attraction = handles.attraction;
totalTimePoints = handles.totalTimePoints;

bee_setup_script

%% Step 2: Covering each time step
% Time loop
for timestep=1:totalTimePoints 
	% What happens at each time step?
    % Update agents
    %% Step 3: Covering each agent
    % Agent loop
    nestSimulationData(timestep+1,:,:) = rules(dt,nestSimulationData(timestep,:,:),...
        broodPosition,emptyFoodPosition,fullFoodPosition,...
        AtoI_Unbumped,ItoA_Unbumped,AtoI_Bumped,ItoA_Bumped,velocityPDF,exposure_state,tags,[attraction 0.1 0.1 0.1]);
end
% handles.nestSimulationData = nestSimulationData;
% % Update handles structure
% guidata(hObject, handles); 
% 
% %% Plotting
% totalTimePoints = handles.totalTimePoints;
% brood = handles.brood;
% nestSimulationData = handles.nestSimulationData;
for timestep = 1:totalTimePoints 
    %%
    plotCoordinatesAndBrood(nestSimulationData, brood, timestep);
    drawnow
end