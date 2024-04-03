from Simulation import simulate
from Models import *
import numpy as np

ts = np.linspace(0, 100, 1)
inits = np.zeros(10)
simulate.solve_model(dual_target_patch.model, ts, inits)