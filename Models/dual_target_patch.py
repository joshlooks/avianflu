# Model code to run 
import numpy as np

# Model ODEs
def model(y, t, params):
    dydt = np.zeros(len(y))
    return dydt