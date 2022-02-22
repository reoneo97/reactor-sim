# reactor-sim
Repository for Final Year Design Project python simulation

**Demo:** [![Streamlit App](https://static.streamlit.io/badges/streamlit_badge_black_white.svg)](https://share.streamlit.io/reoneo97/reactor-sim/app.py)
## Task 
Design a reactor for esterification of palmitic acid with isopropanol to produce isopropyl palmitate (IPP)


## TODO:
1. Metrics 
To add in design metrics (other than conversion) such as costs, purity, heat recovered etc. This will allow for optimization to be conducted
2. Heating Duty calculations
   - Reaction Heat
3. Non-isothermal consideration
   - Fluid Properties - Cp values for the various compounds
   - Interface Mass Transfer
   - Properties for heating vessel - Heat transfer coefficient
   - Solve for temperature profile
4.  Momentum Balance
   - Fluid friction factors
   - Drag force based on material
   - Calculate Pressure Drop - Can become a metric
5. Mass Balance
   - No idea why this will be an issue but can look at textbook to see 
6. Side Reactions
   - Kinetics for Linoeic Acid
## How to Run:
[![Streamlit App](https://static.streamlit.io/badges/streamlit_badge_black_white.svg)](https://share.streamlit.io/reoneo97/reactor-sim/app.py) for an older version of the build. 

Alternatively clone the repo, install requirements and use `streamlit run app.py` 

### Specifications:
- Target IPP Purity: 99% 
- Plant Capacity: 600 kilo Tonne
- Feedstock: 4% Weight Palmitic Acid
- Palmitic Acid Weight = 24 kilo tonne
# Reaction Information

$$
\ce{C15H31COOH  + C3H7OH<--> C19H38O2 + H2O}
$$

### Esterification

Mechanism is well-known and is a nucleophilic acyl substitution of palmitic acid
Reaction is 1st order with respect to both Palmitic Acid, Isopropanol and the catalyst H



## Thermodynamics
Reaction is a reversible one. 
Endothermic reaction -> Favoured by high temperatures 
Equilibrium pushed forward by
- Letting isopropanol be in excess
- $\ce{H2O}$ to be removed constantly 

## Kinetic/Thermodynamic Reaction Constants
*Obtained from Chin (2006)*


##
Reactor Design Equations
![](pics/2022-01-20-17-52-09.png)
![](pics/2022-01-20-17-52-23.png)


## Reactor 
Reactors should be held in isothermal conditions. Esterification reaction is an equilibrium reaction and if the conditions for the reactor are adiabatic, the backward reaction will be favoured with time (Le Chatlier's Principle)

1. CSTR using p-Toluene Sulfonic Acid catalyst
2. PFR using p-Toluene Sulfonic Acid
3. Packed bed Reactor using ZnA/SG 
## References

https://www.chemanalyst.com/Pricing-data/isopropyl-alcohol-31