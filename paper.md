---
title: 'BeeNestABM: An open-source agent-based model of spatiotemporal distribution of bumblebees in nests'
tags:
  - toxicology
  - biology
  - mathematical model
  - MATLAB
authors:
 - name: Ashlee N. Ford-Versypt
   orcid: 0000-0001-9059-5703
   affiliation: 1
 - name: James D. Crall
   affiliation: 2
 - name: Biswadip Dey
   affiliation: 3
   
affiliations:
 - name: School of Chemical Engineering, Oklahoma State University
   index: 1
 - name: Department of Organismic and Evolutionary Biology, Harvard University
   index: 2 
 - name: Department of Mechanical & Aerospace Engineering, Princeton University
   index: 3  
date: 18 April 2018
bibliography: paper.bib
---

# Summary

This software features the MATLAB source code for an interactive computational model that can be used to study the localized responses 
of bumblebees to sublethal exposures to a prevalent class of pesticides called neonicotinoids. The code involves an agent-based stochastic model for the interactions between 
and movements of individual bees within a nest and the nest-related disruptions that occur due to pesticide exposure. The dynamic states of the bees are stored in a matrix, the default data structure 
of MATLAB. Agent-based modeling allows for building understanding of the colony scale impacts of multiple interacting factors that affect numerous individuals in close proximity and how those change upon pesticide exposure. The scienfic significance is that the model solved in the software focuses 
on the effects of pesticides that occur in bumblebees over short time scales (hours to days) inside a single nest, which includes a much 
smaller spatial region with finer resolution than that considered in the state-of-the-art
[@Thorbek2016]. These temporal and spatial scales are appropriate for modeling the effects of 
neonicotinoid pesticides that account for colony size and interactions between exposed and unexposed individuals. The short and local scales allow the model to explicitly consider neighbor interactions and individual bee interactions with structures inside an isolated
environment without confounding external factors. 

![A) The BeeNestABM model tracks bumblebee activity and motility using empirically estimated probabilities for transitions between active (mobile) and inactive states. B) The location of a bee in relation to the structures such as brood and food pots within the nest influence the transition probabilities and the orientation of bee movement through a combination of random walk and attraction toward the nest structures. C) The transition probabilities contain a component that considers whether the transition is occuring spontaneously or due to social modulation upon collison with a neighboring bee.](BeeNestABMtransitions.png)

# References
