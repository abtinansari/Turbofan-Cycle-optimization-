# Turbofan_Cycle_optimization
A simple package to calculate and optimize the Specific Thrust of unMixed flow dual spool Turbofan Engines based on the method described in:
"Mattingly, J. D., Boyer, K. M., & von Ohain, H. (2006). Elements of propulsion: gas turbines and rockets"


The repository contains 6 codes: 
1- Constant_Specific_Heat.m:
This function uses the HP compressor pressure ratio as an input and computes the SFC (working fluid is modeled as Constant Specific Heat(CSH))
2- Modified_Specific_Heat_Constraint.m: 
uses HP compressor pressure ratio, By-pass ratio, and Fan pressure ratio as inputs to compute the Specific Thrust for an UnMixed flow dual spool Turbofan Engine (working fluid is  modeled as Modified Specific Heat) 
3- Tamb.m:
computes the Standard Tempreture in the atmosphere as a function of height
4- enthalpy.m:
computes enthalpy of the combustible mixture as a function of Tempreture & fuel-air ratio
5- dedmov.m:
Differential Evolution (DE) package for a stochastic minimization of continuous space functions that may be non-differentiable,nonlinear and multimodal    
6- Differential_Evolution.m:
This function initializes variables required for DE optimization

                
