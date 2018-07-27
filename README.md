# PlugFlowCrystallizer
Mathematical modelling of a plug flow crystallizer using the FlowPy fluid flow solver (refer to my FlowPy repository for the code).

The crystalization population balance only includes the term for crystal growth; no nucleation, aggregation or breakage. 

The excel file contains solubility data of a lot of salts and has been pulled from Wikipedia. The solubility data for every compound 
has been regressed and is available in the form of polynomial equations of temperature. Note however that the data does not include the
molecular weights, densities, heats of crystallization and crystallization kinetic parameters of the crystals, those have to be entered 
manually.

