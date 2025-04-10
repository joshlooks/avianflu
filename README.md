# avianflu
Code and data associated with the paper: Introducing a framework for within-host dynamics and mutations modelling of H5N1 influenza infection in humans.

Primary authors (contributed equally): Joshua Looker, Daniel Higgins, Robert Sunnucks.
Other authors (contributed equally): Jonathan Carruthers, Thomas Finnie, Matt J. Keeling, Edward M. Hill

## Running the models
The `Models` directory contains a selection of different models investigated in this project. The final 'version' can be found using `Models.diffusion_advection`.

The branching process model code and run scripts can be found in the `bpm_final.ipynb` and `bpm_final_shift.ipynb` Jupyter notebooks. They call functions from the `Mortality.py` and `simfromcsv.py` files and makes use of `log_adaptive_added_10.csv`, `log_adaptive_shift_10.csv` (final parameter posteriors). Sensitivity analysis to a lower case fatality rate are found in the same Jupyter notebooks.

The `Simulation` directory contains wrappers to the `scipy.odeint` library designed to make running the model easier. 

`article_plotting_final.ipynb`,`article_plotting_shift.ipynb`, `model_sim.py` and `article_plotting_mort_final.ipynb` contain more plotting code used to create figures in the publication (all figures can be found in the `Plots` directory).

## ABC code
ABC code used to fit the ODE model can be found in the `ABC` directory and used methods from Minter and Retkute, 2019 [https://doi.org/10.1016/j.epidem.2019.100368]. The resultant parameter posterior can be found in the `Posterior` directory.

## Results
Final and intermediary results (for the lifespans distribution, effective replication number, etc.) can be found in the `Results` directory as csv files.