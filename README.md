# reactor-sim
Repository for Final Year Design Project python simulation

**Demo:** [![Streamlit App](https://static.streamlit.io/badges/streamlit_badge_black_white.svg)](https://share.streamlit.io/reoneo97/reactor-sim/app.py)
## Task 
Design a reactor for esterification of palmitic acid with isopropanol to produce isopropyl palmitate (IPP)

## How to Run - Online Application:
[![Streamlit App](https://static.streamlit.io/badges/streamlit_badge_black_white.svg)](https://share.streamlit.io/reoneo97/reactor-sim/app.py) for an older version of the build. 

Alternatively clone the repo, install requirements and use `streamlit run app.py` 

## Repo Overview
- opt.py: Optimization code using Bayesian Optimization
- app.py: Entrypoint for Streamlit Application
- main.py: Script for manual simulation of code

- app/: Contains code for streamlit dashboards
- simulator/: Root module for simulation code 
   - reactor/: Sub-module for reactor class files
   - const.py: Location of all thermodynamic/kinetic and costing constants
   - energy.py: Python Classes for energy blaance
   - reaction.py: Python Classes for mass balance
   - metrics.py: Functions for calculating reactor/auxillary unit costsQ   
## Reactors 
Reactors should be held in isothermal conditions. Esterification reaction is an equilibrium reaction and if the conditions for the reactor are adiabatic, the backward reaction will be favoured with time (Le Chatlier's Principle). There are 3 reaction configurations for this project:

1. CSTR using p-Toluene Sulfonic Acid catalyst
2. PFR using p-Toluene Sulfonic Acid
3. Packed bed Reactor using ZnA/SG 
## References
