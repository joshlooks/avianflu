from Simulation import simulate
from Models import diffusion_advection
import numpy as np
from matplotlib import pyplot as plt

Cell_Count_U = 4e8
Cell_Count_L = 6.25e9
ts = np.linspace(0, 10, 1000)

inits = np.array([Cell_Count_U, 0, 0, 1.3e3, Cell_Count_L, 0, 0, 0, 0])
# Parameters b_u, b_l, g, c, d, p_u, p_l, gamma, k, f, r, D, a
# States are U1, E1, I1, V1, U2, E2, I2, V2, X
p = [1.5e-8, 1.5e-6, 4, 2, 5.2, 5e7/(Cell_Count_U), 5e7/(Cell_Count_L), 0, 20, 0.56*2.8e-6 / 7, 0.27 / 7, 0.01, 0.1]
results = simulate.solve_model(diffusion_advection.model, inits, ts, (p,))

experiment = False

if (experiment):
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
else:
    #VU and VL
    x = np.array([4,4,5,5,5,6,6,6,6,6,6,7,7,7,7,8])
    y = 10**np.array([7.2,5.40,5.1,7.05,7.7,7.68,7.85,7.73,7.01,6.073,4.45,5.894,5.5356,4.325,3.959,6.548])
    plt.semilogy(ts, results[:,3], label="VU", color="#001380")
    plt.semilogy(ts, results[:,7], label="VL", color="#800000")
    plt.semilogy(x, y, "o", color='k', label="H5N1 Data")
    plt.xlabel("Days post infection", fontsize=12)
    plt.ylabel(r"Viral Titre (TCID$_{50}$/ml)", fontsize=12)
    plt.legend(fontsize=12)
    plt.ylim([10**0,10**8])
    plt.savefig("Plots/VUvVL.pdf",dpi=1200,bbox_inches='tight')
    plt.close()

    #Dead and infected cells
    plt.semilogy(ts, results[:,2], label="IU", color="#001380")
    plt.semilogy(ts, results[:,6], label="IL", color="#800000")
    plt.semilogy(ts, Cell_Count_U - results[:,0], label="DU", color="#6680FF")
    plt.semilogy(ts, Cell_Count_L - results[:,4], label="DL", color="#FF8080")
    plt.xlabel("Days post infection", fontsize=12)
    plt.ylabel("Number of cells", fontsize=12)
    plt.legend(fontsize=12)
    plt.savefig("Plots/ID.pdf",dpi=1200,bbox_inches='tight')