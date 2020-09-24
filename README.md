# The-spatial-stochastic-model-of-the-circadian-clock
NetLogo code to simulate the effect of the cytoplasmic congestion in the mammalian circadian clock. Please see the info panel of the code for further description of the model.

WHAT IS IT?

This model was extended from the previous mathematical model of the mammalian circadian clock (Kim et al., 2012, MSB). To investigate the effect of the crowdedness of cytoplasmic contents on the circadina clock, this model was developed.

The model simulation suggests that the spatial regulation of PER protein, the repressor of the circadian clock, leads to its sharp switch-like phosphorylation and thus sharp nuclear translocation. This nonlinear nuclear entry plays the critial role in generating the robust circadian rhythms. Please see Beesley et al., for details description of the model.

HOW TO USE IT

Step 1. Please set the initial number of agents.
Atot: Total number of activators.
M: The initial number of Per mRNA.
HYPOPER: The initial number of hypophosphorylated PER.
HYPERPER: The initial number of hyperphosphorylated PER.

O: Number of cytopalsmic obstacles. Please choose the value of 0 to regulate the crowdedness of the cytoplasmic contents. 
Examples 
- The model with 150 cytoplasmic obstacles is the normal cell. 
- The model with 275 cytoplasmic obstacles is the overcrowded cell. 
- The model with more than 300 cytoplasmoc obstacles is the extremely overcrowded cell (i.e. adipocyte). Thus, the model do not simulate the rhythmic PER expression.

Step 2. Please set the parameters, which are described below. The current parameter setting makes the model to simulate the rhythmic PER expression of the normal cell.

pa1: Reaction probability for Per mRNA production for each time step (i.e. tick). 

pa2: Reaction probability for PER protein translation for each time step. 

pd1: Reaction probability for Per mRNA degradation for each time step. 

pd2: Reaction probability for hypophos. PER degradtion for each time step. 

pd3: Reaction probability for hyperphos. PER degradtion for each time step.

Kd: Dissociation constant between hyperphos. PER and activator. 

D: Movement step size of PermRNA and PER protein for each time step. 

Dobs: Movement step size of cytoplasmic obstacles for each time step. 

padvec: Probability that the PER protein is advected to the peri-nucleus by the cytoplasmic flux. 
pim: Probability that hyperphos.PER in the cytoplasm is imported to the nucleus for each time step.

Step 3. Load the reation probabilities for hyperphosphorylation and dephosphrylation by pushing “load” buttom.

Step 4. Set the initial condition of the model by pushing “setup” buttom.

Step 5. Run the simulation by pushing “go” buttom.
