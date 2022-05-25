# Modeling polygenic antibiotic resistance evolution in biofilms

## Introduction
This code was used to perform simulations and to generate results and figures published in article:  Modeling polygenic antibiotic resistance evolution in biofilms. 
It contains main file with classed defining populations (planktonic or biofilm) and experiments functions: biofilm_function.py. This file needs to be loaded by jupyter notebooks in order to perform the simulations. 


## File description
Jupyter notebooks contain the code used to parametrize the simulations to simulate various experiments: 

### IllustrationOfEffects.ipynb

This Jupyter notebook contains code to simulate single treatment at given concentration, starting with mixed population (50:50 wt to mutant). Results of these simulations are shown in Figure 5. 
Note that a change of the name and concentration is required to simulate various scenarios (different panels of figure 5). 

### EvolutinaryExperiment_passages.ipynb
Code in this Jupyter notebook is used to simulate a passage experiment. For various scenarios (plankton, biofilm, or individual mechanisms), benefit and cost associated with the biofilm lifestyle needs to be defined. Results of these simulations are shown in Figure 6.

The same main function is used to simulate two populations - biofilm and plankton, that are connected via dispersal. Various rates of dispersal and attachement can be specified by supplying AdRe parameter to the main function. 
Results of these simulations are shown in Figure 8.

Code in this Jupyter notebook is used 


### Repeated treatment.ipynb
Code in this Jupyter notebook is used to simulate a 7 day treatment, allways with the same concentration, for various concentrations. For various scenarios (plankton, biofilm, or individual mechanisms), benefit and cost associated with the biofilm lifestyle needs to be defined. Results of these simulations are shown in Figure 7.

### Zip files

Zip files contain either figures (Figures.zip) or simulation results used to generate these figures. 

