# BeeNestABM
Bee Nest Agent-Based Model: Agent-based model of spatiotemporal distribution of bumblebees in nests

[![DOI](https://zenodo.org/badge/105827383.svg)](https://zenodo.org/badge/latestdoi/105827383)

![Initial position of bees (blue dots), queen (large red circle) and nest structures (other large circles)](thumbnail.PNG)

## Overview
This software features the MATLAB source code for an interactive computational model that can be used to study the localized responses 
of bumblebees to sublethal exposures to a prevalent class of pesticides called neonicotinoids. The code involves an agent-based stochastic model for the interactions between 
and movements of individual bees within a nest and the nest-related disruptions that occur due to pesticide exposure. The dynamic states of the bees are stored in a matrix, the default data structure 
of MATLAB. Agent-based modeling allows for building understanding of the colony scale impacts of multiple interacting factors that affect numerous individuals in close proximity and how those change upon pesticide exposure. The scienfic significance is that the model solved in the software focuses 
on the effects of pesticides that occur in bumblebees over short time scales (hours to days) inside a single nest. These temporal and spatial scales are appropriate for modeling the effects of 
neonicitinoid pesticides that account for colony size and interactions between exposed and unexposed individuals. The short and local scales allow the model to explicitly consider neighbor interactions and individual bee interactions with structures inside an isolated
environment without confounding external factors. Novel research results using this software for scientific applications have been obtained (Crall et al. 2018a, 2018b).

![A) The BeeNestABM model tracks bumblebee activity and motility using empirically estimated probabilities for transitions between active (mobile) and inactive states. B) The location of a bee in relation to the structures such as brood and food pots within the nest influence the transition probabilities and the orientation of bee movement through a combination of random walk and attraction toward the nest structures. C) The transition probabilities contain a component that considers whether the transition is occuring spontaneously or due to social modulation upon collison with a neighboring bee.](BeeNestABM.png)
To facilitate model reuse and reproducible computational research, we have followed the Overview, Design concepts, and Details (ODD) protocol (Grimm et al. 2006, 2010) for standardizing communication of agent-based models. The purpose of this agent-based model is to track the movements of individual bees within a nest chamber on relatively short time scales (a few seconds to less than one day) considering interactions with nestmates and nest structures such as food or brood pots. The entities are the bumblebees. 

## Entities, state variables, and scales
The entities are characterized by the following state variables that are updated at each time step: the x- and y-coordinates, the velocity, the activity state, and the directional heading angle. The states are stored in the nestSimulationData 3-dimensional MATLAB matrix. The simulation environment also includes nest stuctures--full and empty food pots and brood pots--dispersed throughout the nest. The model considers collectives of bees that are on or off of the nest based on their locations relative to the nest structures. The default length scales are the nestMaxX and nestMaxY settings set to 25 cm and 20 cm, respectively for the dimensions of the nest. The time step is 1/frequency, where the default frequency is 2 Hz, and a simuation runs for  for totalTimePoints, which has typical values between 600 and 7200 for 5-60 min of simulated time. 

## Process overview and scheduling

  1. Checking whether or not bees are close enough for social interactions.
  2. Transitioning between active (moving) and inactive (stationary) states.
  3. Moving with some velocity and at some angle heading (0-2 pi in 2D) from the current position. 
    - velocity: Bees that were previously inactive are set at a velocity selected from the empirically determined distribution of bee velocities. Bees that were already moving have two options for velocities: if the resampled velocity from the distribution is +/- 20% of the current velocity, it is updated, but if the sample velocity is outside this range, then velocities have a 10% probability of switching to the sample value and 90 % probability of staying at the current velocity. 
    - angle: First, an angle perturbed within +/- pi/4 of the current angle heading is calculated as the random walk angle. Next, the pairwise distances between a bee and the collection of other nestmates and the nest structures (itemized by brood, empty food pots, and full food pots) are calculated along with the angles of the resultant total distance vectors for each category. The attractions to the nestmates and the nest structures are specified as attraction weights called lambda (the current values have all three types of nest structures lumped together and no attraction toward nestmates). These are collectively referred to as the  environmental stimuli. The net angle for the environmental stimuli is the mean in polar coordinates of the angles between each bee and the nestmates and structures in the nest weighted by their attraction lambda values. The net environmental stimuli angle and the random walk angle are combined by a weighted mean in poloar coordinates with the weight of the environmental stimuli set to different weights for different cohorts (default values stored in rules.m) and the random walk deviations weighted at 1-environmental stimuli weight. 
  4. Truncating a bee's movement if it hits the wall.
  5. Updating all the states for the next iterations.

## Design concepts

### Basic principles
A key concept of the model is to include rules for behaviors inside the nest that can be modulated upon sublethal exposure to pesticides.

## Authors
James D. Crall, Department of Organismic and Evolutionary Biology, Harvard University, jcrall@oeb.harvard.edu

Biswadip Dey, Department of Mechanical & Aerospace Engineering, Princeton University, biswadip@princeton.edu

Ashlee N. Ford Versypt, School of Chemical Engineering, Oklahoma State University, ashleefv@okstate.edu

## Main files
   
* inSilicoExp_Working.m.
    Script to adjust input parameters for spontaneous and social transitions between 
    active/inactive states, run the simulation via simulationOutputSpatial.m, and 
    analyze the statistics via calculateSummaryStatistics.m.

OR

* optimizerShellSpatial.m
    Function to estimate the unknown parameters (attraction coefficients) using 
    empirically derived values for the transition probabilities and 
    simulationOutputSpatial.m comparing the model output to the data via all 
    metrics from calculateSummaryStatistics.m.

OR

* localSensitivity.m
  Runs the full model with the optimized parameter cases when sensitivity = 0 and 
  analyzes the statistics via calculateSummaryStatistics.m.. numberTrials indicates how many repeated simulations to run. For   numberTrials > 1, means and standard deviations are calculated for each output metric.
  
OR

* BeeAttractionApp/BeeAttraction.m
  MATLAB graphical user interface for educational purposes. It uses an older version of the ABM model where the bees are attracted environmental stimuli weight "attraction" is variable and is comprised of 10% each from the three types of nest structures and 70% from attraction to nestmates. This is no longer used in the production runs for scientific purposes; however, it does serve as a nice example of an agent-based model for educational outreach use.

## Dependent files
* simulationOutputSpatial.m.
   The function runs the model defined by rules.m for a user-defined number of time steps with the default parameters. This file is intended to be called by other scripts to repeat the simulation for determining statistics or to test the impacts of setting different parameter sets for parameter estimation, sensitivity analysis, or in silico experiments. Set vis = 1 to generate plots and videos of the output. Set vis = 0 to run without figures.
   
* rules.m.
NOT intended to be called directly. This function contains all of the rules of the agent-based model that are called at each time step. These are listed under Process overview and scheduling section above.
  
* calculateSummaryStatistics.m.
    Given a nestData object and the coordinates of the nest structures (brood and 
    food pots), the function returns the means over all   time and all bees and the 
    time averaged values for each of the bees for the following calculated metrics: 
    activity (moving or not), distance to other bees, distance to nest structures, 
    distance to the social center of the nest, and portion of time spent on the nest.
    
## References
* Crall, J.A., Dey, B., Ford Versypt, A. N., Colony size and pesticide exposure: using an agent-based model to explore social buffering of neonicotinoid exposure in bumblebees, Submitted to Frontiers in Ecology and Evolution, 2018.
* Crall, J.A., Switzer, C., Oppenheimer, R., Combes, S., Pierce, N., De Bivort, B., Ford Versypt, A. N., Dey, B., Brown, A., Eyster, M., Guerin, C. Chronic neonicotinoid exposure disrupts bumblebee nest behavior, social networks, and thermoregulation, Submitted to Science, 2018.
* Grimm, V., Berger, U., Bastiansen, F., Eliassen, S., Ginot, V., Giske, J., Goss-Custard, J., Grand, T., Heinz, S., Huse, G., Huth, A., Jepsen, J.U., Jorgensen, C., Mooij, W.M., Muller, B., Pe'er, G., Piou, C., Railsback, S.F., Robbins, A.M., Robbins, M.M., Rossmanith, E., Ruger, N., Strand, E., Souissi, S., Stillman, R.A., Vabo, R., Visser, U., DeAngelis, D.L., 2006. A standard protocol for describing individual-based and agent-based models, Ecol. Model. 198, 115-126. 
* Grimm, V., Berger, U., DeAngelis, D.L., Polhill, J.G., Giske, J., Railsback, S.F., 2010. The ODD protocol: A review and first update, Ecol. Model. 221, 2760-2768. 
