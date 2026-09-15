We are going to write a library of python helper functions for working through microbial physiological thermodynamic calculations. Please review the AGENTS.md file found in this directory. 

Two different backend libraries, equilibrator API or pygcc can be used. And are set when the library is first loaded. Please review the following:
 
https://gitlab.com/equilibrator/equilibrator-api https://equilibrator.readthedocs.io/en/latest/

https://bitbucket.org/Tutolo-RTG/pygcc/ 
https://pygcc.readthedocs.io/ 

Please research python programs for balancing chemical reactions that will work well with the above packages, that can calculate formal redox states of atoms, and that handle half reactions and electrons well. Use open source libraries when that is good design practice. 

Research libraries to implement the following chemical figures and plots and interactive plots in a jupyter notebook. Use pandas, matplotlib, seaborn when possible, but don’t hesitate to use a good 

The starting point will be a unbalanced reaction, and optionally, a pH condition and with specific concentrations of other reaction components (rather than 1 M). 

For each reaction, we will first make a figure that shows the two half reactions, one above the other aligned at the arrows. There should be a connecting arrow indicating the transfer of however number of electrons would be transferred. For the redox pair the electron donor should be in the oxidative direction and for the electron acceptor the reaction should be written in the reductive direction. The equations should be balanced for electrons, charge, and atoms. The default normalization is to an electron pair, but this can be overridden to normalize to electron donor, acceptor, H2, or any other component. The oxidative and reductive equations should be labeled on the far left. The formal redox state of the main atom should be above (top half reaction) and below (second half reaction) the major reactants and products (the average if there are more than one atom), aligned properly. On the far right report two values for each half reaction: $$E^{\circ \prime}$$ and $$E^\circ$$, ie at pH=7 and at [H+]=1M. 

Next, we will have a function to plot a redox tower with the two half reactions indicated. Default should be at $$E^{\circ \prime}$$. The title should report the overall $$\Delta G^{\circ \prime}$$ of the reaction. Have both a static form saved as a SVG and PNG and an interactive form, with slides for factors that would change the voltage and gibbs free energy (pH, temperature, concentrations, and gas solubility, etc).  

Lastly, have a function that generates an interactive diagram that plots the gibbs free energy in Y and pH, temperature, H2, or other constituent concentrations on the X axes (selected by a pulldown, like gapminder). There should be a button to save the figure as an svg. 

