# okay, so I am now happy with the fit, I am going to generate data and figures

# first, imports
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
font_manager.fontManager.addfont('/Users/hopeless/Library/Fonts/cmunrm.ttf')
plt.rcParams['mathtext.fontset'] = 'stix'
plt.rcParams['font.family'] = 'CMU Serif'
plt.rcParams['axes.unicode_minus'] = False
plt.rcParams.update({'font.size': 16})

plt.ion()

# define constants
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

# first we start with our magnetic calculations. These are much quicker

# let's define the custom mag functions -> later this week i'd like to extend the pcf package
# in the future, I'd like to define these in 3d
# and have the results save as attributes of the object

def bmag(mag, J,  h, temperature, ionObj):
    newh = h +q*J*mag/muB/1.2/1.2
    mag = -1*ionObj.magnetization(ion, temperature, [0,newh, 0]).T[1]-mag
    return mag

def cmag(mag, J,  h, temperature, ionObj):
    newh = h +q*J*mag/muB/1.2/1.2
    mag = -1*ionObj.magnetization(ion, temperature, [0,0, newh]).T[2]-mag
    return mag

def MFTmagC(ionObj, H, J, temperature): 
    # okay so let's start by defining the MF ham
    # only doing c direction rn because that is faster
    q = 6
    ion = 'Er3+' # hard coded for now dont' @ me
    n = 10 # iterations
    # first step, call the pcf mag because that's a good place to start
    magArr = []
    for h in H: 
        mag = 0
        for i in range(n): 
            mag = fsolve(cmag, mag, args = (J, h, temperature, ionObj))[0] # fsolve spits out an array - this time its one val
        magArr.append(mag) # fso
    return magArr

def MFTmagB(ionObj, H, J, temperature): 
    # okay so let's start by defining the MF ham
    # only doing c direction rn because that is faster
    q = 6
    ion = 'Er3+' # hard coded for now dont' @ me
    n = 10 # iterations
    # first step, call the pcf mag because that's a good place to start
    magArr = []
    for h in H: 
        mag = 0
        for i in range(n): 
            mag = fsolve(bmag, mag, args = (J, h, temperature, ionObj))[0] # fsolve spits out an array - this time its one val
        magArr.append(mag) # fso
    return magArr

# now, we define everything for an ionObj generate -> define
# B params, Jz, Jperp in a terminal
# init cf levels obj, call it ionObj, then run this code for that obj to keep things seperate

# # to init allen's terminal: 
# B20 = -3.559e-2
# B40 = -3.849e-4
# B43 = -1.393e-2
# B60 = 3.154e-6
# B63 = -4.695e-6
# B66 = 3.3815e-5
# Jperp = -0.2e-3
# Jz = -2.4e-3
# g = cef.LandeGFactor(ion)
# AllenBparams =  {'B20': B20, 'B40':B40,'B43': B43, 'B60': B60, 'B63':B63,'B66':B66}
# ionObj = cef.CFLevels.Bdict(ion,AllenBparams)

# # to init my terminal
B20 =   -0.03721092
B40 =   -0.00038796
B43 =   -0.01406804
B60 =    3.1865e-06
B63 =   -3.593e-06
B66 =    3.4913e-05
Jperp = -.53070e-03 # +/- 2.6332e-06 (0.50%) (init = 0)
Jz =    -2.63253e-03

# make my er obj
g = cef.LandeGFactor(ion)
myBparams =  {'B20': B20, 'B40':B40,'B43': B43, 'B60': B60, 'B63':B63,'B66':B66}
ionObj = cef.CFLevels.Bdict(ion,myBparams)


# calculate temp dependent magnetization in AB/C

temps = [0.025, 0.045,0.1, 0.171, .25, .35, .45, .543, .827, 1.008, 2, 6, 20]
labels = [str(T)+ 'K' for T in temps]
n = len(temps)
colors = plt.cm.jet(np.linspace(0,1,n))
tempMagB = []
tempMagC = []
H = np.concatenate((np.linspace(0,1,100), np.linspace(1.01,10, 100)))



for temperature in temps: 
    tempMagB.append(MFTmagB(ionObj, H, Jperp, temperature))
    tempMagC.append(MFTmagC(ionObj, H, Jz, temperature))

# now calculate dm/dh
dmdhB =[]
dmdhC = []
for magB, magC in zip(tempMagB, tempMagC): 
    dmdhB.append(np.gradient(magB, H))
    dmdhC.append(np.gradient(magC, H))

# import mpms, scm1 data

Na = 6.02214076e23 
SCF = 1/(1.07828221e24/Na)
# import susceptibility
RawMTdata = np.genfromtxt('/Users/hopeless/Desktop/LeeLab/data/CsErSe2_data/CsErSe2_MTall.dat', 
                       delimiter='\t', unpack=True, skip_header=1)
## Take some averages because it is TOO many data points
CESMTdata = []
for i in range(len(RawMTdata)):
    CESMTdata.append(np.mean(RawMTdata[i].reshape(-1,5), axis=1))

### Import magnetization

CESMHdata = np.genfromtxt('/Users/hopeless/Desktop/LeeLab/data/CsErSe2_data/CsErSe2_MHall.dat', 
                       delimiter='\t', unpack=True, skip_header=1)

## import mpms data
fname20K = '/Users/hopeless/Desktop/LeeLab/data/CsErSe2_data/MPMS/CsErSe2/CsErSe2_HParC_MvsH_20K.txt'
fname6K = '/Users/hopeless/Desktop/LeeLab/data/CsErSe2_data/MPMS/CsErSe2/CsErSe2_HParC_MvsH_6K.txt'
fname2K = '/Users/hopeless/Desktop/LeeLab/data/CsErSe2_data/MPMS/CsErSe2/CsErSe2_HParC_MvsH_2K.txt'

Mdata20K = np.genfromtxt(fname20K, delimiter=',',  unpack=True, skip_header=1)
Mdata6K = np.genfromtxt(fname6K, delimiter=',',  unpack=True, skip_header=1)
Mdata2K = np.genfromtxt(fname2K, delimiter=',',  unpack=True, skip_header=1)

fnameAB20K = '/Users/hopeless/Desktop/LeeLab/data/CsErSe2_data/MPMS/CsErSe2/CsErSe2_HParAB_MvsH_20K.txt'
fnameAB6K = '/Users/hopeless/Desktop/LeeLab/data/CsErSe2_data/MPMS/CsErSe2/CsErSe2_HParAB_MvsH_6K.txt'
fnameAB2K = '/Users/hopeless/Desktop/LeeLab/data/CsErSe2_data/MPMS/CsErSe2/CsErSe2_HParAB_MvsH_2K.txt'

MdataAB20K = np.genfromtxt(fnameAB20K, delimiter=',',  unpack=True, skip_header=1)
MdataAB6K = np.genfromtxt(fnameAB6K, delimiter=',',  unpack=True, skip_header=1)
MdataAB2K = np.genfromtxt(fnameAB2K, delimiter=',',  unpack=True, skip_header=1)

# import dmdh data
fnames = ['/Users/hopeless/Desktop/LeeLab/data/CsErSe2_data/25mKdn022.txt', # '/Users/hopeless/Desktop/LeeLab/data/CsErSe2_data/45mKdn035.txt',
            '/Users/hopeless/Desktop/LeeLab/data/CsErSe2_data/100mKUp021.txt',
            '/Users/hopeless/Desktop/LeeLab/data/CsErSe2_data/174mKDn020.txt', 
            '/Users/hopeless/Desktop/LeeLab/data/CsErSe2_data/249mKUp019.txt', 
            '/Users/hopeless/Desktop/LeeLab/data/CsErSe2_data/345mKDn018.txt', 
            '/Users/hopeless/Desktop/LeeLab/data/CsErSe2_data/422mKUp017.txt', 
            '/Users/hopeless/Desktop/LeeLab/data/CsErSe2_data/543mKDn016.txt',
            '/Users/hopeless/Desktop/LeeLab/data/CsErSe2_data/658mKUp015.txt', 
            '/Users/hopeless/Desktop/LeeLab/data/CsErSe2_data/827mKDn013-14.txt', 
            '/Users/hopeless/Desktop/LeeLab/data/CsErSe2_data/1008mKdn033.txt']
dmdhLabels = ['25mK', '100mK', '174mK', '249mK', '345mK', '422mK', '543mK', '658mK', '827mK', '1008mK']
dmdhField = []
dmdhData =[]

for fname in fnames: 
    temp = np.genfromtxt(fname, delimiter=',',  unpack=True, skip_header=1)
    dmdhField.append(temp[0])
    dmdhData.append(temp[1])


# fname should be defined in terminal
# data is saved in file -> this is redundant, but very easy to deal with
# save data
# fname = 'mag_properties_calculated_whos_parmas.h5'
with h5py.File(calc_fname, 'w') as hdf:
    hdf.create_dataset('tempMagC', data=tempMagC)
    hdf.create_dataset('tempMagB', data=tempMagB)
    hdf.create_dataset('dmdhC', data=dmdhC)
    hdf.create_dataset('dmdhB', data=dmdhB)
    hdf.create_dataset('CESMHdata', data = CESMHdata)
    hdf.create_dataset('magnetization2K', data=Mdata2K)
    hdf.create_dataset('magnetization6K', data=Mdata6K)
    hdf.create_dataset('magnetization20K', data=Mdata20K)
    hdf.create_dataset('dmdhField', data=np.array(dmdhField, dtype=object), dtype=h5py.vlen_dtype(float))
    hdf.create_dataset('dmdhData', data=np.array(dmdhData, dtype=object), dtype=h5py.vlen_dtype(float))
    hdf.create_dataset('dmdhLabels', data = dmdhLabels)
    hdf.create_dataset('temps', data = temps)
    hdf.attrs['B20'] = B20
    hdf.attrs['B40'] = B40
    hdf.attrs['B43'] = B43
    hdf.attrs['B60'] = B60
    hdf.attrs['B63'] = B63
    hdf.attrs['B66'] = B66
    hdf.attrs['Jz'] = Jz
    hdf.attrs['Jperp'] = Jperp
    hdf.attrs['Na'] = Na
    hdf.attrs['SCF'] = SCF

with h5py.File(calc_fname, 'a') as hdf:
    hdf.create_dataset('H', data = H)
    hdf.create_dataset('magnetizationAB2K', data=MdataAB2K)
    hdf.create_dataset('magnetizationAB6K', data=MdataAB6K)
    hdf.create_dataset('magnetizationAb20K', data=MdataAB20K)

# now we do susceptibility

# first define functions
def susceptibilityB(ionObj, fieldVal, temps):
    chi = []
    for temp in temps: 
        # f = np.arange(fieldVal-.5*fieldVal, .5*fieldVal+fieldVal, .05*fieldVal) 
        f = np.linspace(fieldVal-.1, fieldVal+.1,3)
        m = MFTmagB(ionObj, f, Jz, [temp])
        x = np.gradient(m, f) 
        chi.append(x[1])
    return chi

def susceptibilityC(ionObj, fieldVal, temps):
    chi = []
    for temp in temps: 
        # f = np.arange(fieldVal-.5*fieldVal, .5*fieldVal+fieldVal, .05*fieldVal) 
        f = np.linspace(fieldVal-0.1, fieldVal+0.1,3)
        m = MFTmagC(ionObj, f, Jz, [temp])
        x = np.gradient(m, f) 
        chi.append(x[1]) 
    return chi

temps = np.arange(1,300, 1)
fieldVals = [0, 0.1, 1, 3, -6]
susB = []
susC = []
for fieldVal in fieldVals: 
    susB.append(susceptibilityB(ionObj, fieldVal, temps))
    susC.append(susceptibilityC(ionObj, fieldVal, temps))

# load chi(T) data
#AB
fnames = ['/Users/hopeless/Desktop/LeeLab/data/CsErSe2_data/scm1_chi_T_AB_-6T.txt', 
            '/Users/hopeless/Desktop/LeeLab/data/CsErSe2_data/scm1_chi_T_AB_01T.txt', 
            '/Users/hopeless/Desktop/LeeLab/data/CsErSe2_data/scm1_chi_T_AB_1T.txt']
blabels = ['-6T', '0.1', '1T']
data_temps_AB = []
data_sus_AB =[]


for fname in fnames: 
    temp = np.genfromtxt(fname, delimiter=',',  unpack=True, skip_header=1)
    data_temps_AB.append(temp[0])
    data_sus_AB.append(temp[1])

### grab C axis susceptibility
fnames = ['/Users/hopeless/Desktop/LeeLab/data/CsErSe2_data/scm1_chi_T_01T.txt', 
            '/Users/hopeless/Desktop/LeeLab/data/CsErSe2_data/scm1_chi_T_1T.txt', 
            '/Users/hopeless/Desktop/LeeLab/data/CsErSe2_data/scm1_chi_T_3T.txt']
clabels = ['0.1T', '1T', '3T']
data_temps_C = []
data_sus_C =[]


for fname in fnames: 
    temp = np.genfromtxt(fname, delimiter=',',  unpack=True, skip_header=1)
    data_temps_C.append(temp[0])
    data_sus_C.append(temp[1])

# now we calculate low temp sus on the c axis

lowTempFields = [0.358, .116, .37, 5.4, .64]
lowTempSusC =[]
lowTempSusB = []
lowTemps = np.linspace(0,1,100)
for field in lowTempFields: 
    lowTempSusB.append(susceptibilityB(ionObj, fieldVal, lowTemps))
    lowTempSusC.append(susceptibilityC(ionObj, fieldVal, lowTemps))

# load low temp susceptibility data

fnames = ['/Users/hopeless/Desktop/LeeLab/data/CsErSe2_data/00358T0warm24.txt', 
            '/Users/hopeless/Desktop/LeeLab/data/CsErSe2_data/00358T0cool25.txt', 
            '/Users/hopeless/Desktop/LeeLab/data/CsErSe2_data/0116T0warm26.txt',
            '/Users/hopeless/Desktop/LeeLab/data/CsErSe2_data/0116T0cool27.txt', 
            '/Users/hopeless/Desktop/LeeLab/data/CsErSe2_data/037T0warm28.txt', 
            '/Users/hopeless/Desktop/LeeLab/data/CsErSe2_data/037T0cool29.txt',
            '/Users/hopeless/Desktop/LeeLab/data/CsErSe2_data/54T0warm32.txt', 
            '/Users/hopeless/Desktop/LeeLab/data/CsErSe2_data/064T0warm30.txt', 
            '/Users/hopeless/Desktop/LeeLab/data/CsErSe2_data/064T0cool31.txt']
labels = ['.0358T', '.0358T', '.116T', '.116T','.37T', '.37T','5.4T', '0.64T', '0.64T']
data_lowTemps = []
data_lowTempSusC =[]


for fname in fnames: 
    temp = np.genfromtxt(fname, delimiter=',',  unpack=True, skip_header=1)
    data_lowTemps.append(temp[0])
    data_lowTempSusC.append(temp[1])


Na = 6.02214076e23 
SCF = 1/(1.07828221e24/Na)
# fname should be defined in terminal
# data is saved in file -> this is redundant, but very easy to deal with
# save data
# fname = 'susceptibility_calculated_whos_parmas.h5'
with h5py.File(sus_calc_fname, 'w') as hdf:
    hdf.create_dataset('temps', data = temps)
    # hdf.create_dataset('lowTemps', data = lowTemps)
    hdf.create_dataset('susB', data = susB)
    hdf.create_dataset('susC', data = susC)
    hdf.create_dataset('data_temps_C', data = np.array(data_temps_C, dtype=object), dtype=h5py.vlen_dtype(float))
    hdf.create_dataset('data_temps_AB', data = np.array(data_temps_AB, dtype=object), dtype=h5py.vlen_dtype(float))
    hdf.create_dataset('data_sus_C', data = np.array(data_sus_C, dtype=object), dtype=h5py.vlen_dtype(float))
    hdf.create_dataset('data_sus_AB', data = np.array(data_sus_AB, dtype=object), dtype=h5py.vlen_dtype(float))
    hdf.create_dataset('clabels', data = clabels)
    hdf.create_dataset('blabels', data = blabels)
    hdf.create_dataset('fieldVals', data = fieldVals)
    # hdf.create_dataset('lowTempField', data = lowTempFields)
    # hdf.create_dataset('lowTempSusB', data = lowTempSusB)
    # hdf.create_dataset('lowTempSusC', data = lowTempSusC)
    # hdf.create_dataset('lowTempLabels', data = labels)
    # hdf.create_dataset('data_lowTemps', data = np.array(data_lowTemps, dtype=object), dtype=h5py.vlen_dtype(float))
    # hdf.create_dataset('data_lowTempSus', data = np.array(data_lowTempSusC, dtype=object), dtype=h5py.vlen_dtype(float))
    hdf.attrs['B20'] = B20
    hdf.attrs['B40'] = B40
    hdf.attrs['B43'] = B43
    hdf.attrs['B60'] = B60
    hdf.attrs['B63'] = B63
    hdf.attrs['B66'] = B66
    hdf.attrs['Jz'] = Jz
    hdf.attrs['Jperp'] = Jperp
    hdf.attrs['Na'] = Na
    hdf.attrs['SCF'] = SCF

#########################################################################################
# now we scan through some Js
# just for a checkeroni

# grab 2k data
fname2K = '/Users/hopeless/Desktop/LeeLab/data/CsErSe2_data/MPMS/CsErSe2/CsErSe2_HParC_MvsH_2K.txt'
Mdata2K = np.genfromtxt(fname2K, delimiter=',',  unpack=True, skip_header=1)
temperature = 1.8
# J vals
Jvals = np.linspace(0, 2, 5)*Jz
# tempMagC = []
plt.figure()
for j in Jvals: 
    temp = MFTmagC(ionObj, H, j, temperature)
    tempMagC.append(temp)
    plt.plot(H, temp, label = 'JJz =' + str(j))


plt.figure()
for j, temp in zip(Jvals, tempMagC): 
    plt.plot(H, temp, label = 'JJz =' + str(j))

plt.plot(Mdata2K[0], Mdata2K[1]/1.37, 'bo')
plt.grid(True)
plt.xlim(0,8)

with h5py.File('jjz_test.h5', 'w') as hdf:
    hdf.create_dataset('H', data=H)
    hdf.create_dataset('Jvals', data=Jvals)
    hdf.create_dataset('tempMagC', data=np.array(tempMagC))
    hdf.create_dataset('Mdata2K_x', data=Mdata2K[0])
    hdf.create_dataset('Mdata2K_y', data=Mdata2K[1])
    hdf.attrs['temperature'] = temperature
#########################################################################################
#########################################################################################
#########################################################################################
#########################################################################################

# okay, so now we do spectroscopy calculations
# assume ion obj are still initialized

# first, let me define some custom functions
# again, this needs to be built into a class later

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


def zeemanSplitLinesC(field, B20, B40, B43, B60, B63, B66, Jz):     
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
    return amp, dE

def zeemanSplitC_raman(field, wavenum, B20, B40, B43, B60, B63, B66, Jz):    
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
        wid = 0.6
        centers = dE
        tempfun = lorentzian(wavenum, phononAmp, dEphonon, phononSig)
        # tempfun = lorentzian(wavenum, 0, dEphonon, phononSig)
        for i in range(len(centers)):
            a = tempAmp[i]
            tempfun += lorentzian(wavenum, a, centers[i]*meVToCm, wid)
        fun.append(tempfun)
    return fun

def zeemanSplitC_IR(field, wavenum, B20, B40, B43, B60, B63, B66, Jz):    
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
        wid = 0.6
        centers = dE
        tempfun = lorentzian(wavenum, 0, dEphonon, phononSig)
        for i in range(len(centers)):
            a = tempAmp[i]
            tempfun += lorentzian(wavenum, a, centers[i]*meVToCm, wid)
        fun.append(tempfun)
    return fun

def zeemanSplitAB(field, wavenum, B20, B40, B43, B60, B63, B66, Jperp):    
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
        wid = 0.6
        centers = dE
        tempfun = lorentzian(wavenum, 0, dEphonon, phononSig)
        for i in range(len(centers)):
            a = tempAmp[i]
            tempfun += lorentzian(wavenum, a, centers[i]*meVToCm, wid)
        fun.append(tempfun)
    return fun

def zeemanSplitLinesAB(field, B20, B40, B43, B60, B63, B66, Jperp):     
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


# first, load real data

# raman data
fname = '/Users/hopeless/Desktop/LeeLab/data/CsErSe2_data/CsErSe2_raman.csv'
rawData = pd.read_csv(fname, index_col=0, skiprows=0, header=1, delimiter=',')
ramanData = rawData

# take the log - make sure to do this FIRST
ramanData = np.log(ramanData)
ramanData = ramanData-ramanData.min(axis=None)
# normalize 
ramanData = ramanData/ramanData.max(axis=None)

ramanField = [float(b) for b in ramanData.columns.values]
ramanWavenums = [float(i) for i in ramanData.index.values]


# next c-axis IR
cfname = '/Users/hopeless/Desktop/LeeLab/data/CsErSe2_data/Maglab_IRMeasurement_June2022/ProcessingCode/Load2_TrimData/P3_CsEr_100_RAWAVG.dat'

# clean up first
dataC = pd.read_csv(cfname, index_col=0, skiprows=0, header=1)
dataC = dataC.dropna(axis = 0)
dataC = dataC.dropna(axis=1)
dataC = dataC.drop(labels = '-1.1', axis=1)
rawData = dataC


normSpec = dataC['0.001']/max(dataC['0.001'])*-1
avgSpec = normSpec
for column in dataC.columns: 
    dataC[column] = max(dataC[column]) -dataC[column]
    dataC[column] = dataC[column]/(max(dataC[column])) -normSpec
    avgSpec = avgSpec + dataC[column]

for column in dataC.columns: 
    dataC[column] = dataC[column]-avgSpec/len(dataC.columns)
    dataC[column] = dataC[column]-(sum(dataC[column])/len(dataC[column]))

dataC = dataC.drop(labels='0.001', axis=1) # drop this column because we used it as bg


dataC = dataC/dataC.max(axis=None)

IRCfield = [float(b) for b in dataC.columns.values]
IRCwavenums = [float(i) for i in dataC.index.values]


# finally, b-axis ir
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

IRBfield = [float(b) for b in dataB.columns.values]
IRBwavenums = [float(i) for i in dataB.index.values]

# now, generate lines for each, simulated spectrum for each
calc_field = np.arange(0,18, 0.02)
calc_wavenums = np.arange(0,120, 0.1)
temperature = 5
kBT = kB*temperature

ampC,arrC = zeemanSplitLinesC(calc_field, B20, B40, B43, B60, B63, B66, Jz)
arrC = np.array(arrC)
arrC = arrC*meVToCm
arrC = arrC.T

ampC = np.array(ampC)
ampC = ampC.T

ampB,arrB = zeemanSplitLinesAB(calc_field, B20, B40, B43, B60, B63, B66, Jperp)
arrB = np.array(arrB)
arrB = arrB*meVToCm
arrB = arrB.T

ampB = np.array(ampB)
ampB = ampB.T

# now generate lines with no mean field

ampC_nomft, arrC_nomft = zeemanSplitLinesC(calc_field, B20, B40, B43, B60, B63, B66, 0)
arrC_nomft = np.array(arrC_nomft)
arrC_nomft = arrC_nomft*meVToCm
arrC_nomft = arrC_nomft.T

ampC_nomft = np.array(ampC_nomft)
ampC_nomft = ampC_nomft.T

ampB_nomft, arrB_nomft = zeemanSplitLinesAB(calc_field, B20, B40, B43, B60, B63, B66, 0)
arrB_nomft = np.array(arrB_nomft)
arrB_nomft = arrB_nomft*meVToCm
arrB_nomft = arrB_nomft.T

ampB_nomft = np.array(ampB_nomft)
ampB_nomft = ampB_nomft.T

# now generate simulated IR specs

simulated_spec_Raman_C = zeemanSplitC_raman(calc_field, calc_wavenums, B20, B40, B43, B60, B63, B66, Jz)
simulated_spec_Raman_C  = np.array(simulated_spec_Raman_C)

simulated_spec_IR_B = zeemanSplitAB(calc_field, calc_wavenums, B20, B40, B43, B60, B63, B66, Jperp)
simulated_spec_IR_B  = np.array(simulated_spec_IR_B)

simulated_spec_IR_C = zeemanSplitC_IR(calc_field, calc_wavenums, B20, B40, B43, B60, B63, B66, Jz)
simulated_spec_IR_C  = np.array(simulated_spec_IR_C)

# spectroscopy_calc_fname
with h5py.File(spec_calc_fname, 'w') as hdf:
    hdf.create_dataset('linesC', data = arrC)
    hdf.create_dataset('linesB', data = arrB)
    hdf.create_dataset('ampC', data = ampC)
    hdf.create_dataset('ampB', data = ampB)
    hdf.create_dataset('linesC_nomft', data = arrC_nomft)
    hdf.create_dataset('linesB_nomft', data = arrB_nomft)
    hdf.create_dataset('ampC_nomft', data = ampC_nomft)
    hdf.create_dataset('ampB_nomft', data = ampB_nomft)
    hdf.create_dataset('calc_field', data = calc_field)
    hdf.create_dataset('calc_wavenums', data= calc_wavenums)
    hdf.create_dataset('simulated_raman', data = simulated_spec_Raman_C)
    hdf.create_dataset('simulated_IR_B', data = simulated_spec_IR_B)
    hdf.create_dataset('simulated_IR_C', data = simulated_spec_IR_C)
    hdf.create_dataset('ramanData', data=ramanData.values)
    hdf.create_dataset('IR_dataB', data=dataB.values)
    hdf.create_dataset('IR_dataC', data=dataC.values)
    hdf.create_dataset('IR_B_wavenums', data = IRBwavenums)
    hdf.create_dataset('IR_C_wavenums', data = IRCwavenums)
    hdf.create_dataset('raman_wavenums', data = ramanWavenums)
    hdf.create_dataset('IR_B_field', data = IRBfield)
    hdf.create_dataset('IR_C_field', data = IRCfield)
    hdf.create_dataset('raman_field', data = ramanField)
    hdf.attrs['B20'] = B20
    hdf.attrs['B40'] = B40
    hdf.attrs['B43'] = B43
    hdf.attrs['B60'] = B60
    hdf.attrs['B63'] = B63
    hdf.attrs['B66'] = B66
    hdf.attrs['Jz'] = Jz
    hdf.attrs['Jperp'] = Jperp

#########################################################################################
#########################################################################################
#########################################################################################
#########################################################################################
#########################################################################################
# refitting starting with params from model 2
# this is based on a hunch
fitData =ramanData


dropidx =[]
# here I'm cutting out all the phonons I can - the phonon spec is boring and idgaf
for idx in fitData.index: 
    if idx>100: 
        dropidx.append(idx)
fitData = fitData.drop(labels = dropidx, axis = 0)


B20 = 3.114e-2
B40 = -4.718e-4
B43 = 1.259e-2
B60 =  9.324e-7
B63 = 4.715e-5
B66 =  2.011e-5
field = [float(b) for b in fitData.columns.values]
wavenum = [float(i) for i in fitData.index.values]
# now do fit
model = lmfit.Model(zeemanSplitC, independent_vars=['field', 'wavenum'])
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

result = model.fit(z, field=field, wavenum=wavenum, params =params, method =  'ampgo')

print(result.fit_report())

######################################



ampC, arrC = zeemanSplitLinesC(np.linspace(0,15,100), B20, B40, B43, B60, B63, B66, Jz)
arrC = np.array(arrC)
arrC = arrC*meVToCm
arrC = arrC.T

ampC = np.array(ampC)
ampC = ampC.T

# import matplotlib.colors as mcolors
cmap = mcolors.ListedColormap(['cyan'])

fig, ax = plt.subplots()
plt.contourf(ramanField, ramanWavenums, ramanData,100, cmap = 'gnuplot')
plt.xlim(0,14)
plt.ylim(0,120)
plt.clim(-0.3, 1)
plt.xlabel('Field (T)')
plt.ylabel('Energy (cm$^{-1}$)')
plt.colorbar()
plt.title('CsErSe2 H||C with overlayed  calclines\n B20: '+ str(B20)+' B40: '+str(B40)+' B43: ' +str(B43)+ '\n B60: ' +str(B60) + ' B63: ' + str(B63)+ ' B66: ' + str(B66))
for i in range(len(arrC)):
    if i<16: 
        plt.plot(np.linspace(0,15,100), arrC[i], 'c', alpha=1, linewidth= .7)
    if i>=16:  
        alphas = ampC[i]
        # Create a LineCollection
        points = np.array([np.linspace(0,15,100), arrC[i]]).T.reshape(-1, 1, 2)
        segments = np.concatenate([points[:-1], points[1:]], axis=1)
        lc = LineCollection(segments, cmap=cmap, alpha=alphas)
        ax.add_collection(lc)
        ax.autoscale()
plt.show()