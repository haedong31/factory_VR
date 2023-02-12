import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

def learning_curve(a,b,x):
    return a*np.power(x,-b)

