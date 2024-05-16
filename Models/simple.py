# Model code to run 
import numpy as np
from collections import namedtuple

params = namedtuple('params', ["b", "g", "c", "d", "p"])
def model(y, t, p):
    # Parameters b, g, c, d, p
    # States are U1, E1, I1, V1
    p = params(*p)
    dydt = np.zeros(len(y)) + 0.0
    #dSu/dt = beta_u*Su*Vu
    dydt[0] = - p.b*y[0]*y[3]
    #dEu/dt = beta_u*Su*Vu - g*Eu
    dydt[1] = p.b*y[0]*y[3] - p.g*y[1]
    #dIu/dt = g*Eu - d*Iu
    dydt[2] = p.g*y[1] - p.d*y[2]
    #dVu/dt = p_u*Iu - c*Vu
    dydt[3] = p.p*y[2] - p.c*y[3]
    #dDu/dt = p.c*y[3]
    dydt[4] = p.d*y[2]
    return dydt