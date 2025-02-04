import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import PyCrystalField as cef
import scipy
from scipy.optimize import minimize
from scipy.optimize import leastsq
from scipy.optimize import fsolve
import pandas as pd
import lmfit
from matplotlib import font_manager
import matplotlib.colors as  mcolors
from matplotlib.collections import LineCollection
font_manager.fontManager.addfont('/Users/hopeless/Library/Fonts/cmunrm.ttf')
plt.rcParams['mathtext.fontset'] = 'stix'
plt.rcParams['font.family'] = 'CMU Serif'
plt.rcParams['axes.unicode_minus'] = False
plt.rcParams.update({'font.size': 16})

plt.ion()

# define some constants
temperature = 2 # in K

muB = 5.7883818012e-2  # meV/T
mu0 = np.pi*4e-7       # T*m/A
kB  = 8.617e-2         # [meV/K]
meVToCm =meVTocCmInv= 8.066 
ion = 'Er3+'

kBT = kB*temperature
gJ = cef.LandeGFactor('Er3+')
q = 6




def lorentzian( wave, amp, cen, wid ):
    return np.array([amp * wid**2 / ( wid**2 + ( x - cen )**2) for x in wave])

def newH(ionObj, H, J, temperature): 
    # okay so let's start by defining the MF ham
    # only doing c direction rn because that is faster
    q = 6
    ion = 'Er3+' # hard coded for now dont' @ me
    n = 4 # iterations
    # first step, call the pcf mag because that's a good place to start
    mag = 0
    for i in range(n): 
        mag = fsolve(cmag, mag, args = (J, H, temperature, ionObj))[0] # fsolve spits out an array - this time its one val
    newh = H+q*J*mag/muB/1.2/1.2
    return newh

def newHAB(ionObj, H, J, temperature): 
    # okay so let's start by defining the MF ham
    # only doing c direction rn because that is faster
    q = 6
    ion = 'Er3+' # hard coded for now dont' @ me
    n = 4 # iterations
    # first step, call the pcf mag because that's a good place to start
    mag = 0
    for i in range(n): 
        mag = fsolve(bmag, mag, args = (J, H, temperature, ionObj))[0] # fsolve spits out an array - this time its one val
    newh = H+q*J*mag/muB/1.2/1.2
    return newh

def diagonalizeC(ionObj, ion, Jz, h, temperature): 
    # first calc effective h
    h = newH(ionObj, h, Jz, temperature)
    JdotB = muB*(h*cef.Operator.Jz(ionObj.J))*cef.LandeGFactor(ion)
    H_cef = np.sum([a*b for a,b in zip(ionObj.O, ionObj.B)], axis=0)
    ionObj.diagonalize(H_cef + JdotB.O) # this is just H = Hcef + Hmag
    return ionObj.eigenvalues 

def diagonalizeAB(ionObj, ion, J, h, temperature): 
    h = newHAB(ionObj, h, J, temperature)
    JdotB = muB*(h*cef.Operator.Jy(ionObj.J))*cef.LandeGFactor(ion)
    H = np.sum([a*b for a,b in zip(ionObj.O, ionObj.B)], axis=0)
    ionObj.diagonalize(H + JdotB.O) # this is just H = Hcef + Hmag
    evals = ionObj.eigenvalues
    return evals


def zeemanSplitLinesC(field, B20, B40, B43, B60, B63, B66, Jz, temperature):     
    # assuming only H||B rn
    # assuming that x is an array
    amp = []#[1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1,]#[0, .15, .15, .2, 0.15,0.15,0.15,0.15,0.07,0.07, .1,.1,.1,.1,.1]
    Bparams =  {'B20': B20, 'B40':B40,'B43': B43, 'B60': B60, 'B63':B63,'B66':B66}
    ionObj = cef.CFLevels.Bdict(ion,Bparams)
    amp = []
    dE = []
    for b in field: 
        evals = diagonalizeC(ionObj, ion, Jz, b, temperature)
        dE_temp =[eval for eval in evals] # this is the spitting if everything is in the GS -> not necessarily true for finite temp
        # first, lets calculate the partition function - the microstates of the system don't change, jsut the number of abs line
        # we can make between them, so we can make this now
        # Z = sum_i (e^(-beta*E_i))
        Z = [np.exp(-Ei/kBT) for Ei in dE_temp]
        Z = sum(Z)
        # now that we have the partition fn, we can calculate probabilities
        p = [1/Z*np.exp(-Ei/kBT) for Ei in dE_temp]
        # temp_amp = [1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1,]
        # okay, so now we want to add lines between non ground states 
        # we want the amplitude to be the probability -> main lines are already determined
        numLines = len(dE_temp)
        for i in range(1,numLines): 
            # skip GS - already have those dE
            for j in range(i+1, numLines):
                temp = dE_temp[j]-dE_temp[i]
                dE_temp.append(temp)
                p.append(p[i])
        dE.append(dE_temp)
        amp.append(p)
    return amp, dE # dE returned in meV

def zeemanSplitC_raman(field, wavenum, B20, B40, B43, B60, B63, B66, Jz, temperature):    
    # assuming that x is an array
    # amp = [amp1, amp2, amp3, amp4, amp5, amp6, amp7, amp8, amp9, amp10]#, amp11, amp12, amp13, amp14, amp15, amp16]
    amp = [1,1,1,1,1,1,0.3, 0.3, 0.3, 0.3]
    dEphonon = 49.3
    phononAmp = 0.499
    phononSig = 0.95
    fun = []
    Bparams =  {'B20': B20, 'B40':B40,'B43': B43, 'B60': B60, 'B63':B63,'B66':B66}
    ionObj = cef.CFLevels.Bdict(ion,Bparams)
    h=0
    evals = diagonalizeC(ionObj, ion, Jz, h, temperature)
    dE =[eval for eval in evals] # this is the spitting if everything is in the GS -> not necessarily true for finite temp
    # dE = dE[0:10] # only want to look at the bottome few lines - lets see if this works??
    amp = [1 for e in dE]
    if dE[2]*meVTocCmInv > 60: 
        fun = np.ones((len(field), len(wavenum)))*dE[2]*meVTocCmInv # penalty for high energy B, this is only used for global fit fun!!
    else: 
        for b in field: 
            evals = diagonalizeC(ionObj, ion, Jz, b, temperature)
            dE =[eval for eval in evals] # this is the spitting if everything is in the GS -> not necessarily true for finite temp
            dE = dE[0:10] # only want to look at the bottome few lines - lets see if this works??
            tempAmp = amp # amp array for thermal lines
            Z = [np.exp(-Ei/kBT) for Ei in dE]
            Z = sum(Z)
            p = [1/Z*np.exp(-Ei/kBT) for Ei in dE] # prob of population in each line
            numLines = len(dE)
            for i in range(1,numLines): # make thermal lines. starting at 1 so we skip gs
                for j in range(i+1, numLines):
                    temp = dE[j]-dE[i]
                    dE.append(temp)
                    tempAmp.append(p[i]*amp[j])
            wid = 1
            centers = dE
            tempfun = lorentzian(wavenum, phononAmp, dEphonon, phononSig)
            # tempfun = lorentzian(wavenum, 0, dEphonon, phononSig)
            for i in range(len(centers)):
                a = tempAmp[i]
                tempfun += lorentzian(wavenum, a, centers[i]*meVToCm, wid)
            fun.append(tempfun)
    return fun

def zeemanSplitC_IR(field, wavenum, B20, B40, B43, B60, B63, B66, Jz, temperature):    
    # assuming that x is an array
    # amp = [amp1, amp2, amp3, amp4, amp5, amp6, amp7, amp8, amp9, amp10]#, amp11, amp12, amp13, amp14, amp15, amp16]
    amp = [0.1,0.3,0.3,0.15,0.2,0.2,0.287, 0.2, 0.135, 0.097]
    dEphonon = 49.3
    phononAmp = 0.499
    phononSig = 0.95
    fun = []
    Bparams =  {'B20': B20, 'B40':B40,'B43': B43, 'B60': B60, 'B63':B63,'B66':B66}
    ionObj = cef.CFLevels.Bdict(ion,Bparams)
    for b in field: 
        evals = diagonalizeC(ionObj, ion, Jz, b, temperature)
        dE =[eval for eval in evals] # this is the spitting if everything is in the GS -> not necessarily true for finite temp
        dE = dE[0:10] # only want to look at the bottome few lines - lets see if this works??
        tempAmp = amp
        Z = [np.exp(-Ei/kBT) for Ei in dE]
        Z = sum(Z)
        p = [1/Z*np.exp(-Ei/kBT) for Ei in dE]
        numLines = len(dE)
        for i in range(1,numLines): 
            for j in range(i+1, numLines):
                temp = dE[j]-dE[i]
                dE.append(temp)
                tempAmp.append(p[i]*amp[j])
        wid = 1
        centers = dE
        tempfun = lorentzian(wavenum, 0, dEphonon, phononSig)
        for i in range(len(centers)):
            a = tempAmp[i]
            tempfun += lorentzian(wavenum, a, centers[i]*meVToCm, wid)
        fun.append(tempfun)
    return fun

def zeemanSplitAB(field, wavenum, B20, B40, B43, B60, B63, B66, Jperp, temperature):    
    # assuming that x is an array
    # amp = [amp1, amp2, amp3, amp4, amp5, amp6, amp7, amp8, amp9, amp10]#, amp11, amp12, amp13, amp14, amp15, amp16]
    amp = [0.1,0.3,0.3,0.15,0.2,0.2,0.287, 0.2, 0.135, 0.097]
    dEphonon = 49.3
    phononAmp = 0.499
    phononSig = 0.95
    fun = []
    Bparams =  {'B20': B20, 'B40':B40,'B43': B43, 'B60': B60, 'B63':B63,'B66':B66}
    ionObj = cef.CFLevels.Bdict(ion,Bparams)
    for b in field: 
        evals = diagonalizeAB(ionObj, ion, Jperp, b, temperature)
        dE =[eval for eval in evals] # this is the spitting if everything is in the GS -> not necessarily true for finite temp
        dE = dE[0:10] # only want to look at the bottome few lines - lets see if this works??
        tempAmp = amp
        Z = [np.exp(-Ei/kBT) for Ei in dE]
        Z = sum(Z)
        p = [1/Z*np.exp(-Ei/kBT) for Ei in dE]
        numLines = len(dE)
        for i in range(1,numLines): 
            for j in range(i+1, numLines):
                temp = dE[j]-dE[i]
                dE.append(temp)
                tempAmp.append(p[i]*amp[j])
        wid = 1
        centers = dE
        tempfun = lorentzian(wavenum, 0, dEphonon, phononSig)
        for i in range(len(centers)):
            a = tempAmp[i]
            tempfun += lorentzian(wavenum, a, centers[i]*meVToCm, wid)
        fun.append(tempfun)
    return fun

def zeemanSplitLinesAB(field, B20, B40, B43, B60, B63, B66, Jperp, temperature):     
    # assuming only H||B rn
    # assuming that x is an array
    amp = []#[1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1,]#[0, .15, .15, .2, 0.15,0.15,0.15,0.15,0.07,0.07, .1,.1,.1,.1,.1]
    Bparams =  {'B20': B20, 'B40':B40,'B43': B43, 'B60': B60, 'B63':B63,'B66':B66}
    ionObj = cef.CFLevels.Bdict(ion,Bparams)
    amp = []
    dE = []
    for b in field: 
        evals = diagonalizeAB(ionObj, ion, Jperp, b, temperature)
        dE_temp =[eval for eval in evals] # this is the spitting if everything is in the GS -> not necessarily true for finite temp
        # first, lets calculate the partition function - the microstates of the system don't change, jsut the number of abs line
        # we can make between them, so we can make this now
        # Z = sum_i (e^(-beta*E_i))
        Z = [np.exp(-Ei/kBT) for Ei in dE_temp]
        Z = sum(Z)
        # now that we have the partition fn, we can calculate probabilities
        p = [1/Z*np.exp(-Ei/kBT) for Ei in dE_temp]
        # temp_amp = [1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1,]
        # okay, so now we want to add lines between non ground states 
        # we want the amplitude to be the probability -> main lines are already determined
        numLines = len(dE_temp)
        for i in range(1,numLines): 
            # skip GS - already have those dE
            for j in range(i+1, numLines):
                temp = dE_temp[j]-dE_temp[i]
                dE_temp.append(temp)
                p.append(p[i])
        dE.append(dE_temp)
        amp.append(p)
    return amp, dE