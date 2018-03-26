---
title: 'BeeNestABM: An open-source agent-based model of spatiotemporal distribution of bumblebees in nests'
tags:
  - toxicology
  - biology
  - mathematical model
  - MATLAB
authors:
 - name: Ashlee N. Ford Versypt
   orcid: 0000-0001-9059-5703
   affiliation: 1,2
 - name: James R. Crall
   affiliation: 3
 - name: Biswadip Dey
   affilitation: 4
   
affiliations:
 - name: School of Chemical Engineering, Oklahoma State University
   index: 1
 - name: Interdisciplinary Toxicology Program, Oklahoma State University
   index: 2  
 - name: Department of Organismic and Evolutionary Biology, Harvard University
   index: 3 
 - name: Department of Mechanical & Aerospace Engineering, Princeton University
   index: 4  
date: 26 March 2018
bibliography: paper.bib
---

# Summary

This software features the MATLAB source code for an interactive computational model that can be used to study the localized responses 
of bumblebees to sublethal exposures to a prevalent class of pesticides called neonicotinoids. The code involves an agent-based stochastic model for the interactions 
and movements between individual bees within a nest. The dynamic states of the bees are stored in a matrix, the default data structure 
of MATLAB. Agent-based modeling allows for understanding the impacts of multiple interacting factors on individuals at the colony 
scale where numerous individuals are in close proximity. The scienfic significance is that the model solved in the software focuses 
on the effects of pesticides that instantiate in bumblebees over shorter time scales (hours to days) inside a single nest, which includes a much 
smaller spatial region with finer resolution and far fewer individual bees than those considered in the state-of-the-art
@Thorbek2016. These temporal and spatial scales allow the model to consider the mechanisms underlying the negative effects of 
neonicitinoid pesticides due to colony size and interactions between exposed and unexposed individuals. The short and local scales 
allow the model to explicitly consider neighbor interactions and individual bee interactions with nest structures in an enclosed 
environment isolated from confounding external factors. The model developed here focuses on nest-related disruptions due to the 
pesticide exposure. Novel research results using this software for scientific applications have been obtained [@Crall2018a and @Crall2018b].

![GUI screenshot](thumbnail.png)

# References
