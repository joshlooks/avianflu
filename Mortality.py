import numpy as np
import matplotlib.pyplot as plt
def logint(v,t):
    #we return the integral of the logarithm
    return np.trapz(np.log(v),t)



def findM(V_lists, T_lists, mort_rate):
    #create an empty list to store integrals
    Int_Lists = []
    list_size = len(T_lists)

    #We store all of the entire integral
    for i in range(0,list_size):
        Int_Lists.append(logint(V_lists[i],T_lists[i]))
    
    #We sort the list
    Int_Lists.sort()

    #We find the index location to seperate between alive and dead
    boundaryloc = (list_size-1) * mort_rate
    boundaryind = int(boundaryloc)

    #We use linear interpolation to find the M value
    alpha = boundaryloc - boundaryind
    M_val = Int_Lists[boundaryind] + alpha * (Int_Lists[boundaryind+1] - Int_Lists[boundaryind])
    return M_val



def logintlist(v,t):
    #we create a list of the integral of the logarithm
    #each element of the list is a slightly longer integral
    integlist = []
    for i in range(1,len(t)+1):
        integlist.append(logint(v[0:i], t[0:i]))
    return integlist

def lifespan(v,t,M):
    #If the integral never exceeds M, the host never dies, so the whole time range is the viral lifespan
    if M > logint(v,t):
        return t[-1]
    else:
        #otherwise, we work out when the host dies, and return that
        mylist = logintlist(v,t)
        ind = np.where(mylist >= M)[0][0]
        alpha = (M - mylist[ind-1]) / (mylist[ind] - mylist[ind-1])
        timeofdeath = t[ind-1] + alpha * (t[ind] - t[ind-1])
        return timeofdeath

#A function to cut off simulations upon going below a certain value
def cutoff(v,t,endval):
    #if it ends off above that value, we return the whole thing
    if v[-1] > endval:
        return [v,t]
    else:
        #otherwise we calculate the index where it goes below the value
        #this is slightly more complicated as we must account for possibly starting below the cutoff
        exceedindex = np.where(v>endval)[0][0]
        endindices = np.where(v<endval)[0]
        thisindex = np.where(endindices > exceedindex)[0][0]
        endindex = endindices[thisindex]

        alpha = (endval - v[endindex - 1]) / (v[endindex] - v[endindex - 1])

        t_final = t[endindex - 1] + alpha * (t[endindex] - t[endindex - 1])

        truncv = np.zeros(endindex + 1)
        truncv[0:endindex] = v[0:endindex] + 0
        truncv[-1] = endval

        trunct = np.zeros(endindex + 1)
        trunct[0:endindex] = t[0:endindex] + 0
        trunct[-1] = t_final

        #we return the truncated lists
        return [truncv, trunct]


def ViralLifespanDist(V_List, T_List, mortality_rate,end_viral_load):
    #This function will take in simulation results, as well as some other parameters
    #And will return the distribution of viral lifespans

    #cut off simulations when they reach the end viral load value
    numsims = len(V_List)
    newV = []
    newT = []
    for i in range(0,numsims):
        Vi, Ti = cutoff(V_List[i], T_List[i], end_viral_load)
        newV.append(Vi)
        newT.append(Ti)

    #calculate a value of M
    M_val = findM(newV,newT,mortality_rate)

    #get lifespans for each simulation
    lifespans = []
    for i in range(0,numsims):
        lifespans.append(lifespan(newV[i],newT[i],M_val))

    #return the lifespans
    return lifespans



# #Testing our function

# #Make some generic T and V
# T = np.linspace(0,50,1000)
# V = T * np.exp(1-T) + 1/(T+2)**10

# Ts = []
# Vs = []
# #Make copies, with varying magnitudes
# N = 10000
# for i in range(0,N):
#     Ts.append(T)
#     randomval = np.random.normal(loc=8.5, scale=1)
#     Vs.append(10**randomval * V)


# #Plot these
# for i in range(0,len(Vs)):
#     plt.semilogy(Ts[i],Vs[i])
# plt.show()

# #Plot these with cutoff
# for i in range(0,len(Vs)):
#     V,T = cutoff(Vs[i],Ts[i],10**4)
#     plt.semilogy(T,V)
# plt.semilogy(np.linspace(0,30,10), 10**4 * np.ones(10))
# plt.show()

# #Calculate lifespans
# mydist = ViralLifespanDist(Vs,Ts,0.5,10**4)

# plt.hist(mydist, bins = int(np.sqrt(N)))
# plt.show()

