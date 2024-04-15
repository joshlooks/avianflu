# Model code to run 
import numpy as np
from collections import namedtuple

params = namedtuple('params', ["b_u", "b_l", "g", "c", "d", "p_u", "p_l", "gamma", "k", "f", "r", "D", "a"])
def model(y, t, p):
    # Parameters b_u, b_l, g, c, d, p_u, p_l, gamma, k, f, r, D, a
    # States are U1, E1, I1, V1, U2, E2, I2, V2, X
    p = params(*p)
    dydt = np.zeros(len(y)) + 0.0
    dydt[0] = - p.b_u*y[0]*y[3]
    dydt[1] = p.b_u*y[0]*y[3] - p.g*y[1]
    dydt[2] = p.g*y[1] - p.d*y[2]
    dydt[3] = p.p_u*y[2] - p.c*y[3] - p.gamma*p.b_u*y[0]*y[3] - p.D*(y[3]-y[7]) + p.a*y[7]
    dydt[4] = - p.b_l*y[4]*y[7]
    dydt[5] = p.b_l*y[4]*y[7] - p.g*y[6]
    dydt[6] = p.g*y[6] - p.d*y[7]
    dydt[7] = p.p_l*y[6] - p.c*y[7] - p.gamma*p.b_l*y[5]*y[9] - p.k*y[9]*y[11] + p.D*(y[3]-y[7]) - p.a*y[7]
    dydt[8] = p.f*y[7] + p.r*y[8]
    return dydt