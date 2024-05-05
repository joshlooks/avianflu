import numpy as np
from scipy.stats import loguniform, norm, truncnorm, multivariate_normal
from numpy.random import choice
import ABC.tvmn
from copy import copy

def abc_truncated(dist_fun, N_part = 1000, N_gen = 7, tol = [1e17,1e17,1e17,1e17,1.5e16,1.4e16,1.35e16], std = [1e-9,1e-9,1e-2,1e-2,1e-5,5e-3,5e-3]):
    # NB: HYPERPARAMETERS:
    # Distance function, tolerance and std for kernel all effect results
    # N_part and N_gen just improve posterior resolution/accuracy
    wold = np.zeros(N_part)
    wnew = np.zeros(N_part)
    resold = np.zeros((7, N_part))
    resnew = np.zeros((7, N_part))
    lb = np.array([1e-9,1e-9,1e-4,1e-4,1e-6,1e-4,1e-4])
    ub = np.array([1e-5,1e-5,1,1,1e-2,1,1])
    
    # For each generation
    for g in range(N_gen):
        i = 0
        print(f"Beginning iteration {g}")
        # Until we have required number of particles in this generation
        while i < N_part:
            # If first generation sample from prior
            if g == 0:
                b_u = loguniform.rvs(lb[0],ub[0])
                b_l = loguniform.rvs(lb[1],ub[1])
                p_u = loguniform.rvs(lb[2],ub[2])
                p_l = loguniform.rvs(lb[3],ub[3])
                gamma = loguniform.rvs(lb[4],ub[4])
                D = loguniform.rvs(lb[5],ub[5])
                a = loguniform.rvs(lb[6],ub[6])
                theta = np.array([b_u, b_l, p_u, p_l, gamma, D, a])
            # Else sample from previous generation
            else:
                # Sample based on particle weights
                p = choice(N_part,size=1,p=wold)[0]
                # Perturb via independent normal distributions
                # a, b = (a_trunc - loc) / scale, (b_trunc - loc) / scale
                theta = np.array([truncnorm(loc=resold[k,p],scale=std[k],a=(lb[k]-resold[k,p])/std[k],b=(ub[k]-resold[k,p])/std[k]).rvs(1)[0] for k in range(len(theta))])
            # Check if theta within prior range
            if ((lb < theta) & (theta < ub)).all():
                # Calculate distance metric
                dist = dist_fun(theta)
                # If particle in acceptable range
                if dist <= tol[g]:
                    # Store particle and update weights
                    resnew[:,i] = copy(theta)
                    w1 = np.prod([loguniform.pdf(theta[k],a=lb[k],b=ub[k]) for k in range(len(theta))])
                    if g == 0:
                        w2 = 1
                    else:
                        w2 = 0
                        for k in range(N_part):
                            w2 += np.prod([norm(resold[k,p],std[k]).pdf(theta[k]) for k in range(len(theta))])
                    wnew[i] = w1/w2
                    i += 1
                    if (i%100 == 0):
                        print(f"Particle {i} done")
        resold = copy(resnew)
        wold = wnew/np.sum(wnew)
    return resnew   

def abc(dist_fun, N_part = 1000, N_gen = 7, tol = [1e17,1e17,1e17,1e17,1.5e16,1.4e16,1.35e16], std = np.array([1e-9,1e-9,1e-2,5e-3,5e-6,1e-3,1e-3])):
    # wide std: std = np.array([1e-9,1e-9,1e-2,1e-2,1e-5,5e-3,5e-3])
    # NB: HYPERPARAMETERS:
    # Distance function, tolerance and std for kernel all effect results
    # N_part and N_gen just improve posterior resolution/accuracy
    wold = np.zeros(N_part)
    wnew = np.zeros(N_part)
    resold = np.zeros((7, N_part))
    resnew = np.zeros((7, N_part))
    lb = np.array([1e-9,1e-9,1e-4,1e-4,1e-6,1e-4,1e-4])
    ub = np.array([1e-7,1e-7,1,1,5e-5,5e-2,5e-2])
    
    # For each generation
    for g in range(N_gen):
        i = 0
        print(f"Beginning generation {g}")
        # Until we have required number of particles in this generation
        while i < N_part:
            # If first generation sample from prior
            if g == 0:
                b_u = loguniform.rvs(lb[0],ub[0])
                b_l = loguniform.rvs(lb[1],ub[1])
                p_u = loguniform.rvs(lb[2],ub[2])
                p_l = loguniform.rvs(lb[3],ub[3])
                gamma = loguniform.rvs(lb[4],ub[4])
                D = loguniform.rvs(lb[5],ub[5])
                a = loguniform.rvs(lb[6],ub[6])
                theta = np.array([b_u, b_l, p_u, p_l, gamma, D, a])
            # Else sample from previous generation
            else:
                # Sample based on particle weights
                p = choice(N_part,size=1,p=wold)[0]
                # Perturb via independent normal distributions
                theta = np.array([norm(resold[k,p],std[k]).rvs(1)[0] for k in range(len(theta))])
                # theta = np.array([truncnorm(loc=resold[k,p],scale=std[k],a=(lb[k]-resold[k,p])/std[k],b=(ub[k]-resold[k,p])/std[k]).rvs(1)[0] for k in range(len(theta))])
                # theta = ABC.tvmn.TruncatedMVN(mu=resold[:,p],cov=np.diag(std**2),lb=lb,ub=ub).sample(1)
                # theta = multivariate_normal(mean = resold[:,p], cov = cov, allow_singular=True).rvs()
            # Check if theta within prior range
            if ((lb < theta) & (theta < ub)).all():
                # Calculate distance metric
                dist = dist_fun(theta)
                # If particle in acceptable range
                if dist <= tol[g]:
                    # Store particle and update weights
                    resnew[:,i] = copy(theta)
                    w1 = np.prod([loguniform.pdf(theta[k],a=lb[k],b=ub[k]) for k in range(len(theta))])
                    if g == 0:
                        w2 = 1
                    else:
                        w2 = 0
                        for k in range(N_part):
                            w2 += wold[k]*np.prod([norm(resold[j,k],std[j]).pdf(theta[j]) for j in range(len(theta))])
                            # w2 += wold[k]*multivariate_normal(mean = resold[:,k], cov = cov, allow_singular = True).pdf(theta)
                    wnew[i] = w1/w2
                    i += 1
                    if (i%100 == 0):
                        print(f"Particle {i} done")
        resold = copy(resnew)
        wold = wnew/np.sum(wnew)
    return resnew  