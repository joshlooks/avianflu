# Model code to run 
import numpy as np
from collections import namedtuple

params = namedtuple('params', ["lambda_u", "lambda_l", "b_u", "b_l", "g", "c", "d", "p", "kappa", "gamma", "k", "w", "delta", "f", "r", "t_1", "t_2"])
def model(y, t, p):
    # Parameters lambda_u, lambda_l, b_u, b_l, g, c, d, p, kappa, gamma, k, w, delta, f, r, t_1, t_2
    # States are U1, E1, I1, D1, V1, U2, E2, I2, D2, V2, F, X
    p = params(*p)
    dydt = np.zeros(len(y)) + 0.0
    dydt[0] = p.lambda_u*y[3] - p.b_u*y[0]*y[4]
    dydt[1] = p.b_u*y[0]*y[4] - p.g*y[1]
    dydt[2] = p.g*y[1] - p.d*y[2]
    dydt[3] = p.d*y[2] - p.lambda_u*y[3]
    dydt[4] = p.p*y[2] - p.c*y[4] - p.g*p.b_u*y[0]*y[4] - p.t_1*y[4] + p.t_2*y[9]
    dydt[5] = p.lambda_l*y[8] - p.b_l*y[5]*y[9]
    dydt[6] = p.b_l*y[5]*y[9] - p.g*y[6]
    dydt[7] = p.g*y[6] - p.d*y[7]
    dydt[8] = p.d*y[7] - p.lambda_l*y[8]
    dydt[9] = p.p*y[7]/(1+p.kappa*y[10]) - p.c*y[9] - p.gamma*p.b_l*y[5]*y[9] - p.k*y[9]*y[11] - p.t_2*y[9] + p.t_2*y[4]
    dydt[10] = p.w*y[9] - p.delta*y[10]
    dydt[11] = p.f*y[9] - p.r*y[11]
    return dydt