# avianflu
Code and data associated with the paper: Introducing a framework for intra-host dynamics and mutations modelling of H5N1 influenza infection in humans.

Authors (contributed equally): Joshua Looker, Daniel Higgins, Robert Sunnucks.

## Running the models
The `Models` directory contains a selection of different models investigated in this project. The final 'version' can be found using `Models.diffusion_advection`.

The branching process model code and run scripts can be found in the `bpm.ipynb` Jupyter notebook. This calls functions from the `Mortality.py` and `simfromcsv.py` files and makes use of `data.csv` (a version of the final parameter posterior).

The `Simulation` directory contains wrappers to the `scipy.odeint` library designed to make running the model easier. 

`stoch_model.R` contains a stochastic tau-leaping implementation of the model that was not used in the final results but indicates the level of model uncertainty in our results.

## ABC code
This can be found in the `ABC` directory and used methods from Minter and Retkute, 2019 [https://doi.org/10.1016/j.epidem.2019.100368]. The resultant parameter posterior can be found in the `Posterior` directory.

## Results
Final and intermediary results (for the lifespans distribution, effective replication number, etc.) can be found in the `Results` directory as csv files. The plots used in the paper are included as pdfs in the `Plots`directory.