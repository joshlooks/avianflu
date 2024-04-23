import numpy as np
from matplotlib import pyplot as plt
from Models import diffusion_advection
from Simulation import simulate
from scipy.optimize import least_squares

def sol_minimiser(theta):
    Cell_Count_U = 4e8
    Cell_Count_L = 6.25e9
    inits = np.array([Cell_Count_U, 0, 0, 1.3e3, Cell_Count_L, 0, 0, 0, 0])
    p = [theta[0], theta[1], 4, 2, 5.2, theta[2], theta[3], 0, 20, 0.56*2.8e-6 / 7, 0.27 / 7, theta[4], theta[5]]
    ts = np.linspace(0,10,1001)
    x = 100 * np.array([4,4,5,5,5,6,6,6,6,6,6,7,7,7,7,8])
    results = simulate.solve_model(diffusion_advection.model, inits, ts, (p,))
    #Lresults = results[x,3] + results[x,7]
    Lresults = results[x,3]
    y = 10**np.array([7.2,5.40,5.1,7.05,7.7,7.68,7.85,7.73,7.01,6.073,4.45,5.894,5.5356,4.325,3.959,6.548])
    return y - Lresults

Cell_Count_U = 4e8
Cell_Count_L = 6.25e9
test = least_squares(sol_minimiser, [1.5e-8,1.5e-6,5e7/(Cell_Count_U),5e7/(Cell_Count_L),0.01,0.1],bounds=(0,np.inf))

Cell_Count_U = 4e8
Cell_Count_L = 6.25e9
ts = np.linspace(0, 10, 1000)

inits = np.array([Cell_Count_U, 0, 0, 1.3e3, Cell_Count_L, 0, 0, 0, 0])
# Parameters b_u, b_l, g, c, d, p_u, p_l, gamma, k, f, r, D, a
# States are U1, E1, I1, V1, U2, E2, I2, V2, X
theta = test.x
p = [theta[0], theta[1], 4, 2, 5.2, theta[2], theta[3], 0, 20, 0.56*2.8e-6 / 7, 0.27 / 7, theta[4], theta[5]]
results = simulate.solve_model(diffusion_advection.model, inits, ts, (p,))

x = np.array([4,4,5,5,5,6,6,6,6,6,6,7,7,7,7,8])
y = 10**np.array([7.2,5.40,5.1,7.05,7.7,7.68,7.85,7.73,7.01,6.073,4.45,5.894,5.5356,4.325,3.959,6.548])
plt.semilogy(ts, results[:,3], label="VU", color="#001380")
plt.semilogy(ts, results[:,7], label="VL", color="#800000")
plt.semilogy(ts, results[:,3] + results[:,7], label="Total")
plt.semilogy(x, y, "o", color='k', label="H5N1 Data")
plt.xlabel("Days post infection", fontsize=12)
plt.ylabel(r"Viral Titre (TCID$_{50}$/ml)", fontsize=12)
plt.legend(fontsize=12)
plt.ylim([10**0,10**8])
plt.savefig("Plots/lss_low_start.pdf")
plt.show()

test2 = least_squares(sol_minimiser, [1.5e-8,1.5e-6,5e7/(Cell_Count_U),5e7/(Cell_Count_L),1,1],bounds=(0,np.inf))

Cell_Count_U = 4e8
Cell_Count_L = 6.25e9
ts = np.linspace(0, 10, 1000)

inits = np.array([Cell_Count_U, 0, 0, 1.3e3, Cell_Count_L, 0, 0, 0, 0])
# Parameters b_u, b_l, g, c, d, p_u, p_l, gamma, k, f, r, D, a
# States are U1, E1, I1, V1, U2, E2, I2, V2, X
theta = test2.x
p = [theta[0], theta[1], 4, 2, 5.2, theta[2], theta[3], 0, 20, 0.56*2.8e-6 / 7, 0.27 / 7, theta[4], theta[5]]
results2 = simulate.solve_model(diffusion_advection.model, inits, ts, (p,))

x = np.array([4,4,5,5,5,6,6,6,6,6,6,7,7,7,7,8])
y = 10**np.array([7.2,5.40,5.1,7.05,7.7,7.68,7.85,7.73,7.01,6.073,4.45,5.894,5.5356,4.325,3.959,6.548])
plt.semilogy(ts, results2[:,3], label="VU", color="#001380")
plt.semilogy(ts, results2[:,7], label="VL", color="#800000")
plt.semilogy(ts, results2[:,3] + results2[:,7], label="Total")
plt.semilogy(x, y, "o", color='k', label="H5N1 Data")
plt.xlabel("Days post infection", fontsize=12)
plt.ylabel(r"Viral Titre (TCID$_{50}$/ml)", fontsize=12)
plt.legend(fontsize=12)
plt.ylim([10**0,10**8])
plt.savefig("Plots/lss_high_start.pdf")
plt.show()

print(f"Poster values: {0.5*np.sum(np.power(sol_minimiser([1.5e-8,1.5e-6,5e7/(Cell_Count_U),5e7/(Cell_Count_L),0.01,0.1]),2))}")
print(f"Low D,a seed: {test.cost}")
print(f"[b_u,b_l,p_u,p_l,D,a]: {test.x}")
print(f"1, 1 seed: {test2.cost}")
print(f"[b_u,b_l,p_u,p_l,D,a]: {test2.x}")
print(f"Starting cells: {[Cell_Count_U, Cell_Count_L]}")
print(f"Dead cells lowest LSS: {[Cell_Count_U-results[-1,0],Cell_Count_L-results[-1,7]]}")
print(f"Dead cells high seed: {[Cell_Count_U-results2[-1,0],Cell_Count_L-results2[-1,7]]}")