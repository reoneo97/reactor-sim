
# Introduction

## Task 
Design a reactor for esterification of palmitic acid with isopropanol to produce isopropyl palmitate (IPP)

### Specifications:
Target IPP Purity: 99%
Plant Capacity: 600 kilo Tonne
Feedstock: 4% Weight Palmitic Acid

Palmitic Acid Weight = 24 kilo tonne
# Reaction Information

$$
\ce{C15H31COOH  + C3H7OH<--> C19H38O2 + H2O}
$$

### Esterification

Mechanism is well-known and is a nucleophilic acyl substitution of palmitic acid
Reaction is 1st order with respect to both Palmitic Acid, Isopropanol and the catalyst H

## Rate Equations. 

## Simulation Methodology

First simulation was conducted using ideal assumptions. The purpose of these simulations are to derive a reliable baseline for the reactor design and create good estimates of the conversions and whether it is possible to meet the design specifications that have been set. 

For the ideal case:
- Reactor is completely isothermal - Uniform temperature profile both radially and axially
- No radial concentration profile
- Perfect Mixing for CSTR


Simulation of reactors were done differently for each individual reactor. In the first phase, simulation was conducted on ideal reactors. This was done by setting a system of ordinary differential equations together with the mass balance equations of the different reactors being simulated. 

In the ideal case, the first key assumption
## Simulation of Ideal Reactors

## Simulation of 


## Assumptions of Simulation

1. No diffusive mass transfer element. All material changes are due to either convection or reaction
   - Convection only happens in the z-direction
   - Only sources 
2. No Backmixing
3. No Radial Material Flow - The only flow in the radial direction is heat flow
4. Constant Heat of Reaction
5. C_p as a linear function of temperature 
6. Work from Shear Stress does not affect energy balance
Concentration Gradient
- There will be a concentration variation in the z-direction. Makes sense because this is due to difference in residence time
- There is also a concentration gradient in the r-direction. This concentration gradient arises from the temperature gradient. Lower temperature - low rate, 
7. Negligible Mixing Enthalpy
8. Constant Heat of Reaction 

## ODE solution
Model the reactor in cylindrical coordinates $r, \theta, z$
- We assume that $\theta$ direction is not important and everything is uniform w.r.t $\theta$ 

Then we discretize the reactor into 100 intervals. r/R and l/L. Perform mass balance and energy balance on each cylindrical rectangle using Euler's method. 


## Optimization

Idea is to optimize the reactor individually, 
- Potential for recycle stream will increase the performance and reduction in costs even more but for now not very necessary 
Metrics
- Capital Costs
  - Based on L, R, material
  - Pressure not an important variable
    - Set it at a certain value so the reaction will not be in a VLE state
- Operating Cost
  - M - IPA Cost
    - Use Counter Current heating
    - To simplify the calculation, we can first fix the temperature at the outermost radial slice. Given that radial slice, the heat profile should be the same. 
    - 
Parameters
## Energy Balance Equations

## Report Bookmarks
 Page 598 - 603: Non-isothermal design with radial variation


**Laminar Flow:** $U_z = 2U_0\bigg[ 1- \Big(\frac{r}{R} \Big)^2\bigg]$



## Mass balance:
$$
V\frac{dC_A}{dt} = v_0 (C_{A,i-1} - C_{A,i}) - r_AV   
$$

## Energy Balance: 


$$
\dot{m}C_p\frac{dT}{dt} = 
$$
Note: Everything is done in the mol basis

## Conduction


$$
\begin{aligned}
   q_z &= -k_e \frac{\partial T}{\partial z}  \\
   &\approx -k_e(\frac{T_{z+\Delta z} - T_z}{\Delta z}) \\
   q_r &= -k_e \frac{\partial T}{\partial r}\\  
   &\approx -k_e(\frac{T_{r+\Delta r} - T_r}{\Delta z} ) \\ 

\end{aligned}
$$


## Assumption Validity
- Error analysis portion that looks at whether some of the assumptions being amde are valid or not. 


## Optimization Metrics

IPA Cost: 1045 per metric tonne in Dec 2021 -https://www.chemanalyst.com/Pricing-data/isopropyl-alcohol-31

IPP Profit: 7500 per MT -  https://www.alibaba.com/product-detail/Isopropyl-Palmitate-Factory-Supply-99-Isopropyl_1600274796462.html?spm=a2700.7724857.normal_offer.d_title.48cf30cewSSBa2&s=p
  - For the given purity 


Azeotropic Distillation of IPA - https://www.cheric.org/PDF/KJChE/KC23/KC23-1-0001.pdf 
  - Allow for recycling of IPA

## Capital Costing
Towler Textbook - 1104 / 1323
![](2022-02-28-10-50-01.png)


## Heat Transfer Costing
Towler - 913 / 1323

Low Pressure Steam occurs at 6 bar - Don't want large pressure differences due to safety concerns 
- Conservative Estimate of Heat Transfer Coefficient to be 4000 W m2C. 
  - Can possibly do the actual Nusselt number calculations for a better estimate but no point for now 

Dimpled Jacket - However, their one weakness is fatigue failure if subjected to rapid thermal cycling, which typically happens if the same jacket is quickly switched from a heating mode to a cooling mode or vice versa.

4500 W/M2C  