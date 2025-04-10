from Simulation import simulate
from Models import diffusion_advection
import numpy as np
import Mortality
import matplotlib.pyplot as plt
from scipy import integrate
from pathlib import Path



def getBPMparams(gen,file='log_adaptive_added_10.csv', m_value=0.53):

    delta_t = 0.001

    path = Path(__file__).parent.absolute() / file

    df = np.genfromtxt(path,delimiter=",")[1:]

    numrows = np.shape(df)[0]
    numcols = np.shape(df)[1]


    TList = []
    VList = []
    VUList = []
    VLList = []
    rep_U_List = []
    kill_U_List = []
    rep_L_List = []
    kill_L_List = []

    R_List = []


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


        dNU = pU * I_U
        dNL = pL * I_L
        dDU = c * V_U + gamma * betaU * T_U * V_U
        dDL = (c + k * X) * V_L + gamma * betaL * T_L * V_L

        NU_List = []
        NL_List = []
        DU_List = []
        DL_List = []
        #for i in range(0,len(ts)):
        #    NU_List.append(integrate.trapezoid(dNU[0:i+1], x=ts[0:i+1]))
        #    NL_List.append(integrate.trapezoid(dNL[0:i+1], x=ts[0:i+1]))
        #    DU_List.append(integrate.trapezoid(dDU[0:i+1], x=ts[0:i+1]))
        #    DL_List.append(integrate.trapezoid(dDL[0:i+1], x=ts[0:i+1]))
        #NU = np.array(NU_List)
        #NL = np.array(NL_List)
        #DU = np.array(DU_List)
        #DL = np.array(DL_List)

        NU = integrate.cumulative_trapezoid(dNU, ts, initial=0)
        NL = integrate.cumulative_trapezoid(dNL, ts, initial=0)
        DU = integrate.cumulative_trapezoid(dDU, ts, initial=0)
        DL = integrate.cumulative_trapezoid(dDL, ts, initial=0)


        
        #viral growth rate in upper and lower
        rep_U_List.append(1+ (NU[1:] - NU[0:-1]) / V_U[0:-1])
        rep_L_List.append(np.concatenate((np.ones(1),(1+ (NL[2:] - NL[1:-1]) / V_L[1:-1]))))

        
        #viral death rate in upper and lower
        kill_U_List.append((DU[1:] - DU[0:-1]) / V_U[0:-1])
        kill_L_List.append(np.concatenate((np.zeros(1),((DL[2:] - DL[1:-1]) / V_L[1:-1]))))


        VUList.append(results[:,3])
        VLList.append(results[:,7])
        VList.append(results[:,3] + results[:,7])
        TList.append(ts)

    #contains all the lifespans from all the simulations
    LifespansOutput = Mortality.ViralLifespanDist(VList, TList, m_value, 10**4)
    Lifespans = LifespansOutput[0]
    Death_List = LifespansOutput[1]

    timesteps_per_gen = int(gen / delta_t)


    for i in range(0, numrows):
        reps = []
        R = (np.array(rep_U_List[i]) * (1 - np.array(kill_U_List[i])) * np.array(VUList[i][0:-1]) + np.array(rep_L_List[i]) * (1 - np.array(kill_L_List[i])) * np.array(VLList[i][0:-1])) / np.array(VList[i][0:-1])
        count = 0
        while count < (len(R) - timesteps_per_gen):
            rval = 1
            for j in range(count, count + timesteps_per_gen):
                rval *= R[j]
            reps.append(rval)
            count += timesteps_per_gen
        R_List.append(reps)
        
    #Every element of R_list is itself a list, containing the r values needed for the branching process


    #List of lifespans, list of r lists
    return [Lifespans, R_List, Death_List]


