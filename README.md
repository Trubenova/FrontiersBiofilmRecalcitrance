# Modeling polygenic antibiotic resistance evolution in biofilms

## Introduction
This code was used to perform simulations and to generate results and figures published in the article:  Modeling polygenic antibiotic resistance evolution in biofilms. 
It contains the main file with classed defining populations (planktonic or biofilm) and functions defining treatments: biofilm_function.py. This file needs to be loaded by Jupyter notebooks in order to perform the simulations. 


## File description
Jupyter notebooks contain the code used to parametrize the simulations to simulate various experiments: 

### IllustrationOfEffects.ipynb

This Jupyter notebook contains code to simulate a single treatment at a given concentration, starting with a mixed population (50:50 wt to mutant). The results of these simulations are shown in Figure 5. 
Note that a change of the name and concentration is required to simulate various scenarios (different panels of figure 5). 

### EvolutinaryExperiment_passages.ipynb
Code in this Jupyter notebook is used to simulate a passage experiment. For various scenarios (plankton, biofilm, or individual mechanisms), the benefit and the cost associated with the biofilm lifestyle need to be defined. The results of these simulations are shown in Figure 6.

The same main function is used to simulate two populations - biofilm and plankton, that are connected via dispersal. Various rates of dispersal and attachment can be specified by supplying AdRe parameter to the main function. 
The results of these simulations are shown in Figure 8.

### Repeated treatment.ipynb
Code in this Jupyter notebook is used to simulate a 7-day treatment, always with the same concentration, for various concentrations. For various scenarios (plankton, biofilm, or individual mechanisms), the benefit and the cost associated with the biofilm lifestyle need to be defined. The results of these simulations are shown in Figure 7.

### FinalFiguresFrontiers 
This is the original file with the code used to create all the figures in the manuscript, based on the created data files (in zip). Relevant and cleaned parts of this notebook were copied into the notebooks described above. 

### Zip files

Zip files contain either figures (Figures.zip) or simulation results used to generate these figures. 

