# helper functions for fitting
import numpy as np
import matplotlib.pyplot as plt
import PyCrystalField as cef
import scipy
from scipy.optimize import minimize
from scipy.optimize import leastsq
import pandas as pd
import lmfit

ion = 'Er3+'
gz = 1.2

temperature = 15 # in K

muB = 5.7883818012e-2  # meV/T
mu0 = np.pi*4e-7       # T*m/A
kB  = 8.617e-2         # [meV/K];
meVToCm =meVTocCmInv= 8.066 
ion = 'Er3+'
kBT =kB*temperature


def gaussian(x, amp,cen , wid):
    return amp * np.exp(-(x-cen)**2 / wid)

def lorentzian( x, amp, cen, wid ):
    return amp * wid**2 / ( wid**2 + ( x - cen )**2)

def diagonalize(ionObj, ion, Field): 
    JdotB = muB*(Field*cef.Operator.Jx(ionObj.J))*cef.LandeGFactor(ion)
    H = np.sum([a*b for a,b in zip(ionObj.O, ionObj.B)], axis=0)
    ionObj.diagonalize(H + JdotB.O) # this is just H = Hcef + Hmag
    return ionObj.eigenvalues 

def zeemanSplitC(field, wavenum, B20, B40, B43, B60, B63, B66, amp1, amp2, amp3, amp4, amp5, amp6, amp7, amp8,width):# amp9, amp10, width): #, amp11, amp12, amp13, amp14, amp15, amp16):     
    # assuming that x is an array
    amp = [amp1, amp2, amp3, amp4, amp5, amp6, amp7, amp8]#, amp9, amp10, amp11, amp12, amp13, amp14, amp15, amp16]
    # amp = np.ones(10)*0.5
    # amp = amp.tolist()
    Bparams =  {'B20': B20, 'B40':B40,'B43': B43, 'B60': B60, 'B63':B63,'B66':B66}
    ionObj = cef.CFLevels.Bdict(ion,Bparams)
    dEphonon = 49.3
    phononAmp = 0.9
    phononSig = 1
    try: # for the case we pass in a single B val
        evals = diagonalizeC(ionObj, ion, field)
        dE =[eval for eval in evals] # this is the spitting if everything is in the GS -> not necessarily true for finite temp
        # first, lets calculate the partition function - the microstates of the system don't change, jsut the number of abs line
        # we can make between them, so we can make this now
        # Z = sum_i (e^(-beta*E_i))
        Z = [np.exp(-Ei/kBT) for Ei in dE]
        Z = sum(Z)
        # now that we have the partition fn, we can calculate probabilities
        p = [(1/Z)*np.exp(-Ei/kBT) for Ei in dE]
        # okay, so now we want to add lines between non ground states 
        # we want the amplitude to be the probability -> main lines are already determined
        numLines = len(dE)
        for i in range(1,numLines): 
            # skip GS - already have those dE
            for j in range(i+1, numLines):
                temp = dE[j]-dE[i]
                dE.append(temp)
                amp.append(p[i]*amp[j])
        # so now we want to create a multi gaussian with centers at the eigenvalues
        wid = width # pulled this out my ass, .5cm-1 resolution roughly
        a = amp[0]
        centers =dE 
        fun = lorentzian(wavenum, a, centers[0]*meVToCm, wid)
        for i in range(len(centers[1:])):
            a = amp[i]
            fun += lorentzian(wavenum, a, centers[i]*meVToCm, wid)
        for i in range(len(dEphonon)):
            fun+= lorentzian(wavenum, phononAmp[i], dEphonon[i], wid)
    except AttributeError : # for when we pass in an array 
        # pyCrystalField has trouble with the operators because its not *just* a function
        # so this error handling lets us have a function that works for single values and array inputs
        fun = []
        for b in field: 
            evals = diagonalizeC(ionObj, ion, b)
            dE =[eval for eval in evals] # this is the spitting if everything is in the GS -> not necessarily true for finite temp
            dE = dE[0:10] # only want to look at the bottome few lines - lets see if this works??

            # first, lets calculate the partition function - the microstates of the system don't change, jsut the number of abs line
            # we can make between them, so we can make this now
            # Z = sum_i (e^(-beta*E_i))
            Z = [np.exp(-Ei/kBT) for Ei in dE]
            Z = sum(Z)
            # now that we have the partition fn, we can calculate probabilities
            p = [1/Z*np.exp(-Ei/kBT) for Ei in dE]
            # okay, so now we want to add lines between non ground states 
            # we want the amplitude to be the probability -> main lines are already determined
            numLines = len(dE)
            for i in range(1,numLines): 
                # skip GS - already have those dE
                for j in range(i+1, numLines):
                    temp = dE[j]-dE[i]
                    dE.append(temp)
                    amp.append(p[i]*amp[j])
            # so now we want to create a multi gaussian with centers at the eigenvalues
            wid = width # pulled this out my ass, .5cm-1 resolution roughly
            centers = dE
            tempfun = lorentzian(wavenum, phononAmp, dEphonon, phononSig)
            for i in range(len(centers)):
                a = amp[i]
                tempfun += lorentzian(wavenum, a, centers[i]*meVToCm, wid)
            fun.append(tempfun)
    return fun


def diagonalizeC(ionObj, ion, Field): 
    JdotB = muB*(Field*cef.Operator.Jz(ionObj.J))*cef.LandeGFactor(ion)
    H = np.sum([a*b for a,b in zip(ionObj.O, ionObj.B)], axis=0)
    ionObj.diagonalize(H + JdotB.O) # this is just H = Hcef + Hmag
    return ionObj.eigenvalues 

def zeemanSplitLinesC(field, B20, B40, B43, B60, B63, B66):     
    # assuming only H||B rn
    # assuming that x is an array
    amp = []#[1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1,]#[0, .15, .15, .2, 0.15,0.15,0.15,0.15,0.07,0.07, .1,.1,.1,.1,.1]
    Bparams =  {'B20': B20, 'B40':B40,'B43': B43, 'B60': B60, 'B63':B63,'B66':B66}
    ionObj = cef.CFLevels.Bdict(ion,Bparams)
    kBT = 10*kB
    try: # for the case we pass in a single B val
        evals = diagonalizeC(ionObj, ion, field)
        dE =[eval for eval in evals] # this is the spitting if everything is in the GS -> not necessarily true for finite temp
        # first, lets calculate the partition function - the microstates of the system don't change, jsut the number of abs line
        # we can make between them, so we can make this now
        # Z = sum_i (e^(-beta*E_i))
        Z = [np.exp(-Ei/kBT) for Ei in dE]
        Z = sum(Z)
        # now that we have the partition fn, we can calculate probabilities
        p = [1/Z*np.exp(-Ei/kBT) for Ei in dE]
        # okay, so now we want to add lines between non ground states 
        # we want the amplitude to be the probability -> main lines are already determined
        numLines = len(dE)
        for i in range(1,numLines): 
            # skip GS - already have those dE
            for j in range(i+1, numLines):
                temp = dE[j]-dE[i]
                dE.append(temp)
        amp.append(p)
        # so now we want to create a multi gaussian with centers at the eigenvalues
    except AttributeError : # for when we pass in an array 
        # pyCrystalField has trouble with the operators because its not *just* a function
        # so this error handling lets us have a function that works for single values and array inputs
        amp = []
        dE = []
        for b in field: 
            evals = diagonalizeC(ionObj, ion, b)
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
                    # temp_amp.append(p[i])
            dE.append(dE_temp)
            amp.append(p)
    return amp, dE

def dummyModelC(B20, B40, B43, B60, B63, B66): 
    Bparams = {'B20': B20, 'B40':B40,'B43': B43, 'B60': B60, 'B63':B63,'B66':B66}
    ionCef = cef.CFLevels.Bdict(ion,Bparams)
    calcA, calcB, calcC = calculateZeemanSpec(ionCef, ion, BfieldC)
    return calcC.T