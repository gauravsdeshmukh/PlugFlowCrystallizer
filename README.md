# PlugFlowCrystallizer
Mathematical modelling of a plug flow crystallizer using the FlowPy fluid flow solver (refer to my FlowPy repository for the code).

The crystalization population balance only includes the term for crystal growth; no nucleation, aggregation or breakage. 

The excel file contains solubility data of a lot of salts and has been pulled from Wikipedia. The solubility data for every compound 
has been regressed and is available in the form of polynomial equations of temperature. Note however that the data does not include the
molecular weights, densities, heats of crystallization and crystallization kinetic parameters of the crystals, those have to be entered 
manually.

The program works over as another layer on the FlowPy fluid flow code as follows:

1. FlowPy calculates the velocity and pressure distribution. 
2. Using the obtained velocity profile, the temperatures and concentration at each point are calculated from the enthalpy balance and species balance equations.
3. The solubility is calculated at every point (since temperature varies at every point).
4. From the concentration and solubility (C and C*), the driving force and hence the crystal growth rate G is calculated.
5. The population balance equation is solved using the method of moments. The zeroth moment is the number of crystals while the first moment is the length of crystals. These two ODEs are then solved by inputting the appropirate growth rates and boundary conditions. 
6. The final profiles for pressure, temperature, concentration and crystal length are plotted.

Note that this program is still in a nascent stage and that it does not yet completely support physically real parameter inputs (like extremely low diffusivity and thermal diffusivity). Further adjustments in CFL selection and time step calculation would lead to guaranteed convergence. 
