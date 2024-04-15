from Simulation import simulate
from Models import diffusion_advection
import numpy as np
from matplotlib import pyplot as plt

ts = np.linspace(0, 1, 10)
inits = np.array([1e9, 0, 0, 0, 40000, 1e9, 0, 0, 0, 0, 0, 0])
# Parameters b_u, b_l, g, c, d, p_u, p_l, gamma, k, f, r, D, a
# States are U1, E1, I1, V1, U2, E2, I2, V2, X
p = [0, 0, 4, 2, 5.2, 0, 0, 0.15, 20, 2.8e-6, 0.27, 0, 0]
results = simulate.solve_model(diffusion_advection.model, inits, ts, (p,))

#Eclipsed, Infected and Dead cells in URT
plt.plot(results[:,1])
plt.plot(results[:,2])
plt.plot(results[:,3])
plt.show()

plt.plot(results[:,4])
plt.show()