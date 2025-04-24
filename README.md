# Mathematical-Modelling-of-Infectious-Disease
This repository is designed for my MBiol project, which models an epidemic in hypothetical populations with dynamic vulnerable compartments. The project contains 4 modelling exercises with different model structures and/or assumptions to understand the role of time-limited NPIs on disease dynamics:
- Exercise 1: SIR model
- Exercise 2: SIRS model with no immune protection against death
- Exercise 3: SIRS model with perfect immune protection against death
- Exercise 4: SIRS model with perfect immune protection against death + a non-linear relationship between intervention coverage & efficacy in the non-vulnerable sub-population. 

The "Dissertation_code" folder contains the following R (or R markdown) scripts for reproducing data and figures in the dissertation:

"population_structures": This code is for calculating the population-specific parameter alpha (rate of transition to vulnerable state) and long-term equilibrium values of S and Sv compartments in a disease-free model. To calculate alpha, input population-specific parameters mu & mu_v. The resulting parameters are then used as inputs in a disease-free model to find the long-term equilibrium values of S and Sv for each population typology, which will later serve as starting conditions in the disease models. 

"exc1-3_full_code": This R markdown file contains code reproducing figures in exercises 1-3. The code runs with parameters for Population 1 (Age-only framework). To simulate the other 5 population typologies, please replace some of the parameters and starting conditions (those to be replaced are shown in the script). In specific, this script contains:
- Code for simulating SIR dynamics. This generates heat maps for exercise 1.
- Code for simulating SIRS dynamics. This generates heat maps for exercise 2 and 3. For exercise 2, simply run the code as how it is. For exercise 3, change theta in parms1 from 0 to 1 to test different assumptions on the level of immune protection against death.

"ex1-3_time_series_phase_plot": This generates time series plot for peak Iv & cumulative disease deaths, as well as phase plots for exercise 1-3. To test different intervention scenarios, please change the intervention-related parameters to reflect the intervention-of-interest.

"ex2_calculate_Rt": This calculates real-time effective reproduction number Rt for any intervention scenario (exercise 2). 

"ex4_nonlinearity": This generates heat maps for exercise 4, which assumes a non-linear relationship between intervention coverage and efficacy in the non-vulnerable sub-population. Multiple heat maps will be generated to reflect different intervention start days and durations. 

"ex4_time_series_phase_plot": This generates time series plots and phase plots for any intervention-of-interest in exercise 4.

"ex4_comparison_across_6_populations": This generates peak Iv and cumulative disease deaths heat maps for all 6 population typologies, under assumptions of exercise 4.

Exploratory R codes written at early stages of the project are not included in this repository because most of them involve incorrect parameter values, unrealistic model settings (e.g. permanent interventions) or uninteresting results. 





