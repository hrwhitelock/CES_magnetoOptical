# cleaned up code to go with MOS paper

# start w/ imports: 
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
temperature = 10# in K

muB = 5.7883818012e-2  # meV/T
mu0 = np.pi*4e-7       # T*m/A
kB  = 8.617e-2         # [meV/K]
meVToCm =meVTocCmInv= 8.066 
ion = 'Er3+'

kBT = kB*temperature
gJ = cef.LandeGFactor('Er3+')
q = 6

# define functions
# this is split in AB/C 
# AB is just B axis, as our "AB plane" data is actually oriented along B, and thers not much
# anisotropy between A and B


def bmag(mag, J,  h, temperature, ionObj):
    newh = h +q*J*mag/muB/1.2/1.2
    mag = -1*ionObj.magnetization(ion, temperature, [0,newh, 0]).T[1]-mag
    return mag

def cmag(mag, J,  h, temperature, ionObj):
    newh = h +q*J*mag/muB/1.2/1.2
    mag = -1*ionObj.magnetization(ion, temperature, [0,0, newh]).T[2]-mag
    return mag


def newHc(ionObj, H, J, temperature): 
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
    h = newHc(ionObj, h, Jz, temperature)
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

def lorentzian( wave, amp, cen, wid ):
    return np.array([amp * wid**2 / ( wid**2 + ( x - cen )**2) for x in wave])


def zeemanSplitC_raman(field, wavenum, B20, B40, B43, B60, B63, B66, Jz, temperature):    
    # assuming that x is an array
    # amp = [amp1, amp2, amp3, amp4, amp5, amp6, amp7, amp8, amp9, amp10]#, amp11, amp12, amp13, amp14, amp15, amp16]
    amp = [0.1,0.3,0.3,0.15,0.3,0.3,0.3, 0.3, 0.3, 0.3]
    dEphonon = 49.3
    phononAmp = 0.499
    phononSig = 0.95
    fun = []
    Bparams =  {'B20': B20, 'B40':B40,'B43': B43, 'B60': B60, 'B63':B63,'B66':B66}
    ionObj = cef.CFLevels.Bdict(ion,Bparams)
    h=0
    evals = diagonalizeC(ionObj, ion, Jz, h, temperature)
    dE =[eval for eval in evals] # this is the spitting if everything is in the GS -> not necessarily true for finite temp
    dE = dE[0:10] # only want to look at the bottome few lines - lets see if this works??
    if dE[-1]*meVTocCmInv > 50: 
        fun = np.ones((len(field), len(wavenum)))*100 # penalty for high energy B
    else: 
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
            tempfun = lorentzian(wavenum, phononAmp, dEphonon, phononSig)
            # tempfun = lorentzian(wavenum, 0, dEphonon, phononSig)
            for i in range(len(centers)):
                a = tempAmp[i]
                tempfun += lorentzian(wavenum, a, centers[i]*meVToCm, wid)
            fun.append(tempfun)
    return fun

def zeemanSplitAB(field, wavenum, B20, B40, B43, B60, B63, B66, Jperp):    
    amp = [0.1,0.3,0.3,0.15,0.2,0.2,0.287, 0.2, 0.135, 0.097]
    dEphonon = 0
    phononAmp = 0
    phononSig = 0
    fun = []
    Bparams =  {'B20': B20, 'B40':B40,'B43': B43, 'B60': B60, 'B63':B63,'B66':B66}
    ionObj = cef.CFLevels.Bdict(ion,Bparams)
    temperature = 10
    for b in field: 
        evals = diagonalizeAB(ionObj, ion, Jperp, b, temperature)
        dE =[eval for eval in evals]
        dE = dE[0:10] # we only want to look at the bottom 10 lines, as our data is cut off
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

## now we can start fitting
## first we import the C-axis raman data
## I'm choosing to use the raman data because it's cleaner than the IR data
## this also works on the IR data 

# replace with your file path 
fname = '/Users/hopeless/Desktop/LeeLab/data/CsErSe2_data/CsErSe2_raman.csv'
rawData = pd.read_csv(fname, index_col=0, skiprows=0, header=1, delimiter=',')
ramanData = rawData

# put the global minimum at 0 
ramanData = np.log(ramanData)
ramanData = ramanData-ramanData.min(axis=None)
# normalize - global max 
# here we don't want to normalize spectrum by spectrum, because the decay of lines
# carries information 
ramanData = ramanData/ramanData.max(axis=None)

ramanField = [float(b) for b in ramanData.columns.values]
ramanWavenums = [float(i) for i in ramanData.index.values]

fitData = ramanData

dropidx =[]
# here I'm cutting out all the data above 100cm^-1 
# above there CES has a 'dead zone' where we can't see anything

for idx in fitData.index: 
    if idx>100: 
        dropidx.append(idx)
fitData = fitData.drop(labels = dropidx, axis = 0)

field = [float(b) for b in fitData.columns.values]
wavenum = [float(i) for i in fitData.index.values]

# now we pick starting params 
# we're lucky enough to have a starting guess: 
# model 1 from https://doi.org/10.1103/PhysRevB.101.144432

B20 = -3.559e-2
B40 = -3.849e-4
B43 = -1.393e-2
B60 =  3.154e-6
B63 = -4.695e-6
B66 =  3.3815e-5
field = [float(b) for b in fitData.columns.values]
wavenum = [float(i) for i in fitData.index.values]
# now do fit
model = lmfit.Model(zeemanSplitC_raman, independent_vars=['field', 'wavenum'])
params = model.make_params()

params['B20'].set(value= B20, min = -.06, max = 0.06)
params['B40'].set(value= B40, min = -.06, max = 0.06)
params['B43'].set(value= B43, min = -.06, max = 0.06)
params['B60'].set(value= B60, min = -.06, max = 0.06)
params['B63'].set(value= B63, min = -.06, max = 0.06)
params['B66'].set(value= B66, min = -.06, max = 0.06)
params['Jz'].set(value = 0)

z = np.array(fitData.to_numpy()) # gotta do it twice with tuples :((((
z = z.T

# 
resultC = model.fit(z, field=field, wavenum=wavenum, params =params)

print(resultC.fit_report())

###########################################################################################
###########################################################################################
###########################################################################################
# let's do the same fit for the b-axis
# first import b data
# replace with your file location 
bfname = '/Users/hopeless/Desktop/LeeLab/data/CsErSe2_data/Maglab_IRMeasurement_June2022/ProcessingCode/Load1_TrimData/P2_CsEr_100-FIR_RAWAVG.dat'

# clean up AB first
dataB = pd.read_csv(bfname, index_col=0, skiprows=0, header=1)
dataB = dataB.dropna(axis = 0)
dataB = dataB.dropna(axis=1)
dataB = dataB.drop(labels = '-1.1', axis=1)


normSpec = dataB['0.001']/max(dataB['0.001'])*-1
avgSpec = normSpec
for column in dataB.columns: 
    dataB[column] = max(dataB[column]) -dataB[column]
    dataB[column] = dataB[column]/(max(dataB[column])) -normSpec
    avgSpec = avgSpec + dataB[column]

for column in dataB.columns: 
    dataB[column] = dataB[column]-avgSpec/len(dataB.columns)
    dataB[column] = dataB[column]-(sum(dataB[column])/len(dataB[column]))

dataB = dataB.drop(labels='0.001', axis=1) # drop this column because we used it as bg


dataB = dataB/dataB.max(axis=None)

for index in dataB.index:
    if index<20: 
        dataB = dataB.drop(labels=index , axis = 0)
    if index>110:
        dataB = dataB.drop(labels=index , axis = 0)

for col in dataB.columns: 
    for idx in dataB.index: 
        if dataB.loc[idx,col] <0: 
            dataB.loc[idx,col] = 0


# now do the same thing, fit only Jperp
# params from full fit
B20 = -0.03721092 # +/- 3.3632e-04 (0.90%) (init = -0.03559)
B40 = -3.8796e-04 # +/- 8.4839e-07 (0.22%) (init = -0.0003849)
B43 = -0.01406804 # +/- 2.5391e-05 (0.18%) (init = -0.01393)
B60 =  3.1865e-06 # +/- 5.3881e-09 (0.17%) (init = 3.154e-06)
B63 = -3.5930e-06 # +/- 2.4015e-07 (6.68%) (init = -4.695e-06)
B66 =  3.4913e-05 # +/- 3.1308e-07 (0.90%) (init = 3.3815e-05)
Jperp =  0# +/- 1.7886e-05 (0.68%) (init = -0.003)
field = [float(b) for b in dataB.columns.values]
wavenum = [float(i) for i in dataB.index.values]
# now do fit
model = lmfit.Model(zeemanSplitAB, independent_vars=['field', 'wavenum'])
params = model.make_params()

params['B20'].set(value= B20, vary = False)# min = -.06, max = 0.06)
params['B40'].set(value= B40, vary = False)# min = -.06, max = 0.06)
params['B43'].set(value= B43, vary = False)# min = -.06, max = 0.06)
params['B60'].set(value= B60, vary = False)# min = -.06, max = 0.06)
params['B63'].set(value= B63, vary = False)# min = -.06, max = 0.06)
params['B66'].set(value= B66, vary = False)# min = -.06, max = 0.06)
params['Jperp'].set(value = 0)

z = np.array(dataB.to_numpy()) # gotta do it twice with tuples :((((
z = z.T

result = model.fit(z, field=field, wavenum=wavenum, params =params)

print(result.fit_report())
