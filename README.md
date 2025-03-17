# HBT-Calculations

This is the repository that I had been using for calculating and fitting the HBT Correlation Functions to find the three HBT radii and plot them. This is part of my Masters project. 

I have used a hybrid model involving hydrodynamic evolution with hadronic transport, to simulate the collision dynamics of QCD conserved charges. The main components of this hybrid model include: GICG (Glauber Initial Condition Generator), 3+1-D hydrodynamics code MUSIC, ISS (Isochronous Spectra Sampler) for sampling particles from freeze-out using Cooper-Frye formula, and URQMD for late stage hadronic interactions. Using this model, I eventually get a particle data set containing event-by-event information of each particle's multiplicities, as well as their 4-position and 4-momentum values. 

In the code HBT_Correlations.cpp, I extracted data from this particle data set, calculated the 3-D correlation function of identical particle pairs, and stored them in separate output files. 

In the code HBT_corr_func_fit.cpp, I utilized these output files to fit into the correlation function expression and obtain the three HBT radii as functions of the transverse momentum k_T. I then plotted the three HBT radii with k_T. 
