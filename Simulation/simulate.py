import numpy as np
from scipy import integrate

def solve_model(model, ts, inits, params, int_tool = "odeint"):
    if int_tool == "odeint":
        return integrate.odeint(model, inits, ts, args=params)
    else:
        pass
