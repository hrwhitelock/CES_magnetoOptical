# I need a fucking class i cannot do this shit anymore

# okay so we start with imports 
import numpy as np
import matplotlib.pyplot as plt
import PyCrystalField as cef
import scipy
from scipy.optimize import minimize
from scipy.optimize import leastsq
import pandas as pd
import lmfit

plt.ion()

class CES_obj(): 
    def __init__(ion, Bparams):
        # init stuff
        # wrap CF levels
        self.g = cef.LandeGFactor(ion)
        self.ErObj = cef.CFLevels.Bdict(ion, Bparams)
        return

    def PCFmagnetization(temps, field): 
        fieldA = [[b,0,0]for b in field]
        fieldB = [[0,b,0]for b in field]
        fieldC = [[0,0,b]for b in field]
        tempMagA = tempMagB = tempMagC = []
        for temperature in temps: 
            magA = np.array([self.ErObj.magnetization(ion, temperature, f) for f in fieldA]).T[0]*-1
            magB = np.array([self.ErObj.magnetization(ion, temperature, f) for f in fieldA]).T[1]*-1
            magC = np.array([self.ErObj.magnetization(ion, temperature, f) for f in fieldC]).T[2]*-1
            tempMagA.append(magA) # what the actual fuck is this naming holy shit
            tempMagB.append(magB)
            tempMagC.append(magC)

        magnetization = [tempMagA, tempMagB, tempMagC]
        return magnetization

    def MolecularFieldTheory(H, Hth, Mth, lamb): 
            '''Use mean field theory to correct magnetization
        for an exchange strength lamb. H is the number of fields to calculate,
        Hth and Mth are the theoretical single-ion magnetization curves to correct.'''
        n = 10
        newM = np.interp(H, Hth, Mth)
        for i in range(n):
            newH = H + 6*lamb*newM/muB/(gJ)**2
            newM = np.interp(newH,Hth,Mth)
        return newM
    
    def PCFsusceptibility(): 
        #

    def MFTsusceptiblity(): 

