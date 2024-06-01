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
    boundaryloc = (list_size-1) * (1-mort_rate)
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
        return [t[-1],False]
    else:
        #otherwise, we work out when the host dies, and return that
        mylist = logintlist(v,t)
        ind = np.where(mylist >= M)[0][0]
        alpha = (M - mylist[ind-1]) / (mylist[ind] - mylist[ind-1])
        timeofdeath = t[ind-1] + alpha * (t[ind] - t[ind-1])
        return [timeofdeath,True]

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
    print(M_val)

    #get lifespans for each simulation
    death_list = []
    lifespans = []
    for i in range(0,numsims):
        lifespanoutput = lifespan(newV[i],newT[i],M_val)
        lifespans.append(lifespanoutput[0])
        death_list.append(lifespanoutput[1])

    #return the lifespans
    return [lifespans,death_list]