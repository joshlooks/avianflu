from Simulation import simulate
from Models import *
import numpy as np
from matplotlib import pyplot as plt

ts = np.linspace(0, 1, 10)
inits = np.array([1e9, 0, 0, 0, 40000, 1e9, 0, 0, 0, 0, 0, 0])
# ["lambda_u", "lambda_l", "b_u", "b_l", "g", "c", "d", "p", "kappa", "gamma", "k", "w", "delta", "f", "r", "t_1", "t_2"]
p = [0, 0.015, 1e-7, 1e-7, 4.0, 2.0, 10.0, 1, 0.045, 0.13, 10, 10, 1, 1e-6, 0.27, 1, 1]
results = simulate.solve_model(single_target_patch.model, inits, ts, (p,))

#Eclipsed, Infected and Dead cells in URT
plt.plot(results[:,1])
plt.plot(results[:,2])
plt.plot(results[:,3])
plt.show()

plt.plot(results[:,4])
plt.show()