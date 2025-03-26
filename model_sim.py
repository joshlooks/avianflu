from Simulation import simulate
from Models import diffusion_advection
import numpy as np
import Mortality
import matplotlib.pyplot as plt
from scipy import integrate
from pathlib import Path

delta_t = 0.001

path = Path(__file__).parent.absolute() / 'log_adaptive_added_10.csv'

df = np.genfromtxt(path,delimiter=",")[1:]

numrows = np.shape(df)[0]
numcols = np.shape(df)[1]

TimeList = []
VList = []
VUList = []
VLList = []
TUList = []
TLList = []
EUList = []
ELList = []
IUList = []
ILList = []


#need to check rows vs columns
for i in range(0, numrows):
    params =  [df[i,0], df[i,1], 4, 2, 5.2, df[i,2], df[i,3], df[i,4], 20, 0.56*2.8e-7/7, 0.27/7, df[i,5], df[i,6]]

    Cell_Count_U = 4e8
    Cell_Count_L = 6.25e9

    ts = np.linspace(0, 30, int(30 / delta_t) + 1)

    inits = np.array([Cell_Count_U, 0, 0, 1.3e3, Cell_Count_L, 0, 0, 0, 0])
    # Parameters b_u, b_l, g, c, d, p_u, p_l, gamma, k, f, r, D, a
    # States are U1, E1, I1, V1, U2, E2, I2, V2, X
    results = simulate.solve_model(diffusion_advection.model, inits, ts, (params,))

    #New upper Viral particles, N_U, d/dt (N_U) = pU I_U
    #New lower Viral particles, N_L, d/dt (N_L) = pL I_L
    #Dead upper Viral particles, D_U, d/dt (D_U) = c V_U
    #Dead lower Viral particles, D_U, d/dt (D_L) = (c + kX) V_L

    #Need to select correct parameters
    pU = params[5]
    pL = params[6]
    c = params[3]
    betaU = params[0]
    betaL = params[1]
    gamma = params[7]
    k = params[8]
    I_U = results[:,2]
    I_L = results[:,6]
    T_U = results[:,0]
    T_L = results[:,4]
    V_U = results[:,3]
    V_L = results[:,7]
    X = results[:,-1]


    VUList.append(results[:,3])
    VLList.append(results[:,7])
    VList.append(results[:,3] + results[:,7])

    TUList.append(results[:,0])
    EUList.append(results[:,1])
    IUList.append(results[:,2])
    TLList.append(results[:,4])
    ELList.append(results[:,5])
    ILList.append(results[:,6])

    TimeList.append(ts)

#contains all the lifespans from all the simulations
LifespansOutput = Mortality.ViralLifespanDist(VList, TimeList, 13/18, 10**4)
Lifespans = LifespansOutput[0]
Death_List = LifespansOutput[1]


T_upper = []
T_lower = []
T_upper_top = []
T_lower_top = []
T_upper_bot = []
T_lower_bot = []

V_upper = []
V_upper_top = []
V_upper_bot = []
V_lower = []
V_lower_top = []
V_lower_bot = []
Overall_Time = []


i = 0
reached_end = False

while i < len(ts) and reached_end == False:
    myTU = []
    myTL = []
    myVU = []
    myVL = []

    for j in range(0,numrows):
        if ts[i] <= Lifespans[j]:
            myTU.append(TUList[j][i])
            myTL.append(TLList[j][i])
            myVU.append(VUList[j][i])
            myVL.append(VLList[j][i])
    
    

    if len(myTU) > 0:
        Overall_Time.append(ts[i])

        T_upper.append(np.median(np.array(myTU)))
        T_upper_top.append(np.percentile(np.array(myTU),97.5))
        T_upper_bot.append(np.percentile(np.array(myTU),2.5))

        T_lower.append(np.median(np.array(myTL)))
        T_lower_top.append(np.percentile(np.array(myTL),97.5))
        T_lower_bot.append(np.percentile(np.array(myTL),2.5))

        V_upper.append(np.median(np.array(myVU)))
        V_upper_top.append(np.percentile(np.array(myVU),97.5))
        V_upper_bot.append(np.percentile(np.array(myVU),2.5))

        V_lower.append(np.median(np.array(myVL)))
        V_lower_top.append(np.percentile(np.array(myVL),97.5))
        V_lower_bot.append(np.percentile(np.array(myVL),2.5))
    else:
        reached_end = True

    i += 1


#Set up plot
fig1, ax1 = plt.subplots()
ax1.plot(Overall_Time, T_upper, label='Median Upper Dynamics')
ax1.fill(Overall_Time + list(reversed(Overall_Time)), T_upper_bot + list(reversed(T_upper_top)),label='95% Pointwise PI of Upper Dynamics',alpha=0.2)
ax1.plot(Overall_Time, T_lower, label='Median Lower Dynamics')
ax1.fill(Overall_Time + list(reversed(Overall_Time)), T_lower_bot + list(reversed(T_lower_top)),label='95% Pointwise PI of Lower Dynamics',alpha=0.2)
ax1.set_ylabel("Uninfected Target Cells",fontsize=12)
ax1.set_xlabel("Days post infection",fontsize=12)
plt.yscale('log')
plt.savefig('Plots/TargetCellsUpperLower.pdf', format = 'pdf', dpi=2400,bbox_inches='tight')
plt.show()



#Set up plot
fig2, ax2 = plt.subplots()
params = {'mathtext.default': 'regular' }          
plt.rcParams.update(params)
ax2.plot(Overall_Time, V_upper, label='Median Upper Dynamics')
ax2.fill(Overall_Time + list(reversed(Overall_Time)), V_upper_bot + list(reversed(V_upper_top)),label='95% Pointwise PI of Upper Dynamics',alpha=0.2)
ax2.plot(Overall_Time, V_lower, label='Median Lower Dynamics')
ax2.fill(Overall_Time + list(reversed(Overall_Time)), V_lower_bot + list(reversed(V_lower_top)),label='95% Pointwise PI of Lower Dynamics',alpha=0.2)
ax2.set_ylabel('Viral Titre (TCID$_{50}$/ml)',fontsize=12)
ax2.set_xlabel("Days post infection",fontsize=12)
ax2.legend()
plt.yscale('log')
plt.savefig('Plots/VirionsUpperLower.pdf', format = 'pdf', dpi=2400,bbox_inches='tight')
plt.show()