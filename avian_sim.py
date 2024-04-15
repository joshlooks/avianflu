from Simulation import simulate
from Models import diffusion_advection
import numpy as np
from matplotlib import pyplot as plt

ts = np.linspace(0, 30, 1000)
inits = np.array([1e9, 0, 0, 40000, 1e9, 0, 0, 0, 0])
# Parameters b_u, b_l, g, c, d, p_u, p_l, gamma, k, f, r, D, a
# States are U1, E1, I1, V1, U2, E2, I2, V2, X
p = [1.9e-9, 1.9e-7, 4, 2, 5.2, 1, 1, 0.15, 20, 2.8e-6, 0.27, 1, 0.01]
results = simulate.solve_model(diffusion_advection.model, inits, ts, (p,))

#Eclipsed, Infected in URT
plt.plot(ts,results[:,1],label="EU")
plt.plot(ts,results[:,2],label='IU')
plt.legend()
plt.show()

plt.plot(ts, results[:,0],label="SU")
plt.legend()
plt.show()
plt.plot(ts, results[:,4],label="SL")
plt.legend()
plt.show()


#VU and VL
plt.plot(ts,results[:,3],label="VU")
plt.legend()
plt.show()
plt.plot(ts,results[:,7],label="VL")
plt.legend()
plt.show()