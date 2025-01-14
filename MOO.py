# I need a fucking class i cannot do this shit anymore

# okay so we start with imports 
import numpy as np
import matplotlib.pyplot as plt
import PyCrystalField as cef
import scipy
from scipy.optimize import fsolve
import pandas as pd
import lmfit

plt.ion()

class MOO(): 
    def __init__(ion, Bparams, JJ, q):
        # init stuff
        # wrap CF levels
        self.g = cef.LandeGFactor(ion)
        self.ionObj = cef.CFLevels.Bdict(ion, Bparams)
        self.JJ = JJ # 3d vector JJx, JJy, JJz
        self.q = q
        return

    def internal_mag(mag,h, temperature): 
        newh = [h[i] + self.q*self.JJ[i]*mag[i]/

    def mag_mft(temps, field): 

        return magnetization
    


