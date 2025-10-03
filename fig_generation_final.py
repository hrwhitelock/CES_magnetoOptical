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
# ion = 'Er3+'

kBT = kB*temperature
gJ = cef.LandeGFactor('Er3+')
q = 6

# first we start with our magnetic calculations. These are much quicker

# let's define the custom mag functions -> later this week i'd like to extend the pcf package
# in the future, I'd like to define these in 3d
# and have the results save as attributes of the object

def bmag(mag, J,  h, temperature, ionObj):
    newh = h +q*J*mag/muB/1.2/1.2
    mag = -1*ionObj.magnetization(ion, temperature, [float(newh),0, 0]).T[1]-mag
    return mag

def cmag(mag, J,  h, temperature, ionObj):
    newh = h +q*J*mag/muB/1.2/1.2
    mag = -1*ionObj.magnetization(ion, temperature, [0.0,0.0, float(newh)]).T[2]-mag
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

# to init allen's terminal: 
B20 = -3.559e-2
B40 = -3.849e-4
B43 = -1.393e-2
B60 = 3.154e-6
B63 = -4.695e-6
B66 = 3.3815e-5
Jperp = -0.2e-3
Jz = -2.4e-3
g = cef.LandeGFactor(ion)
AllenBparams =  {'B20': B20, 'B40':B40,'B43': B43, 'B60': B60, 'B63':B63,'B66':B66}
ionObj = cef.CFLevels.Bdict(ion,AllenBparams)

# # to init my terminal
B20 =   -0.03721092
B40 =   -0.00038796
B43 =   -0.01406804
B60 =    3.1865e-06
B63 =   -3.593e-06
B66 =    3.4913e-05
Jperp = -.53070e-03 
Jz =    -2.63253e-03

# make my er obj
g = cef.LandeGFactor(ion)
Bparams =  {'B20': B20, 'B40':B40,'B43': B43, 'B60': B60, 'B63':B63,'B66':B66}
ionObj = cef.CFLevels.Bdict(ion,Bparams)


# calculate temp dependent magnetization in AB/C

temps = [0.025, 0.045,0.1, 0.171, .25, .35, .45, .543, .827, 1.008, 2,3,4,5, 6,7,8,9,10,15, 20]
labels = [str(T)+ 'K' for T in temps]
n = len(temps)
colors = plt.cm.jet(np.linspace(0,1,n))
tempMagB = []
tempMagC = []
H = np.concatenate((np.linspace(0,2,100), np.linspace(2.01,18, 150)))



for temperature in temps: 
    tempMagB.append(MFTmagB(ionObj, H, Jperp, temperature))
    tempMagC.append(MFTmagC(ionObj, H, Jz, temperature))


# real quick plot the magnetization
plt.figure()
for i in range(len(temps)):
    plt.plot(H, tempMagC[i], label = str(temps[i]))

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

calc_fname = 'mag_calc_mft_2025May12.h5'

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
        f = np.linspace(fieldVal-.01, fieldVal+.01,3)
        m = MFTmagB(ionObj, f, Jperp, [temp])
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

temps = np.linspace(.1,30,300) #np.concatenate((np.arange(0.025, 25, .025),np.arange(1,300, 1)))
fieldVals = [0.0, 0.1, 0.2]
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

# # now we calculate low temp sus on the c axis

# lowTempFields = [0.358, .116, .37, .64,1, 2, 3, 5.4]
# lowTempSusC =[]
# lowTempSusB = []
# lowTemps = np.linspace(0,1,100)
# for field in lowTempFields: 
#     lowTempSusB.append(susceptibilityB(ionObj, field, lowTemps))
#     lowTempSusC.append(susceptibilityC(ionObj, field, lowTemps))

# load low temp susceptibility data

# fnames = ['/Users/hopeless/Desktop/LeeLab/data/CsErSe2_data/00358T0warm24.txt', 
#             '/Users/hopeless/Desktop/LeeLab/data/CsErSe2_data/00358T0cool25.txt', 
#             '/Users/hopeless/Desktop/LeeLab/data/CsErSe2_data/0116T0warm26.txt',
#             '/Users/hopeless/Desktop/LeeLab/data/CsErSe2_data/0116T0cool27.txt', 
#             '/Users/hopeless/Desktop/LeeLab/data/CsErSe2_data/037T0warm28.txt', 
#             '/Users/hopeless/Desktop/LeeLab/data/CsErSe2_data/037T0cool29.txt',
#             '/Users/hopeless/Desktop/LeeLab/data/CsErSe2_data/064T0warm30.txt', 
#             '/Users/hopeless/Desktop/LeeLab/data/CsErSe2_data/064T0cool31.txt',
#             '/Users/hopeless/Desktop/LeeLab/data/CsErSe2_data/54T0warm32.txt']
# labels = ['.0358T', '.0358T', '.116T', '.116T','.37T', '.37T', '0.64T', '0.64T', '5.4T']
# data_lowTemps = []
# data_lowTempSusC =[]


# for fname in fnames: 
#     temp = np.genfromtxt(fname, delimiter=',',  unpack=True, skip_header=1)
#     data_lowTemps.append(temp[0])
#     data_lowTempSusC.append(temp[1])


Na = 6.02214076e23 
SCF = 1/(1.07828221e24/Na)
# fname should be defined in terminal
# data is saved in file -> this is redundant, but very easy to deal with
# save data
# fname = 'susceptibility_calculated_whos_parmas.h5'
sus_calc_fname = 'sus_mft_0p48ueVJJz_2025May11.h5'

with h5py.File(sus_calc_fname, 'a') as hdf:
    hdf.create_dataset('temps', data = temps)
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
    # hdf.create_dataset('lowTemps', data = lowTemps)
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
    # hdf.attrs['Na'] = Na
    # hdf.attrs['SCF'] = SCF


#########################################################################################
# now we scan through some Js
# just for a checkeroni

# grab 2k data
fname2K = '/Users/hopeless/Desktop/LeeLab/data/CsErSe2_data/MPMS/CsErSe2/CsErSe2_HParC_MvsH_2K.txt'
Mdata2K = np.genfromtxt(fname2K, delimiter=',',  unpack=True, skip_header=1)
temperature = 2.1
# J vals
Jvals = np.linspace(0.75, 1, 1.25, 1.5)*Jz
tempMagC = []

for j in Jvals: 
    temp = MFTmagC(ionObj, H, j, temperature)
    tempMagC.append(temp)



plt.figure()
for j, temp in zip(Jvals, tempMagC): 
    plt.plot(H, temp, label = 'JJz =' + str(j))

plt.plot(Mdata2K[0], Mdata2K[1]/1.42, 'bo')
plt.grid(True)
plt.xlim(0,8)
plt.ylim(0,8)

with h5py.File('jjz_test_dec20_model.h5', 'w') as hdf:
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
    h = newH(ionObj, h, Jz, temperature)
    JdotB = muB*(h*cef.Operator.Jz(ionObj.J))*cef.LandeGFactor(ion)
    H_cef = np.sum([a*b for a,b in zip(ionObj.O, ionObj.B)], axis=0)
    ionObj.diagonalize(H_cef + JdotB.O) # this is just H = Hcef + Hmag
    return ionObj

def diagonalizeAB(ionObj, ion, J, h, temperature): 
    h = newHAB(ionObj, h, J, temperature)
    JdotB = muB*(h*cef.Operator.Jy(ionObj.J))*cef.LandeGFactor(ion)
    H = np.sum([a*b for a,b in zip(ionObj.O, ionObj.B)], axis=0)
    ionObj.diagonalize(H + JdotB.O) # this is just H = Hcef + Hmag
    # evals = ionObj.eigenvalues
    return ionObj

def transition(ionObj,ii, jj, temperature):
    beta = 1/(8.61733e-2*temperature)  # Boltzmann constant is in meV/K
    Z = sum([np.exp(-beta*en) for en in ionObj.eigenvalues])
    # compute population factor
    ket1 = ionObj.eigenvectors.real[ii]
    ket2 = ionObj.eigenvectors.real[jj]
    pn = np.exp(-beta *ionObj.eigenvalues[ii])/Z
    ax = np.dot(ket1, np.dot(ionObj.opttran.Jx, ket2))**2
    ay = np.dot(ket1, np.dot(ionObj.opttran.Jy, ket2))**2
    az = np.dot(ket1, np.dot(ionObj.opttran.Jz, ket2))**2
    return pn*ax, pn*ay, pn*az

def zeemanSplitLinesC(field, B20, B40, B43, B60, B63, B66, Jz):     
    Bparams =  {'B20': B20, 'B40':B40,'B43': B43, 'B60': B60, 'B63':B63,'B66':B66}
    ionObj = cef.CFLevels.Bdict(ion,Bparams)
    amp = []
    dE = []
    for b in field: 
        evals = diagonalizeC(ionObj, ion, Jz, b, temperature)
        dE_temp =[eval for eval in evals] # this is the spitting if everything is in the GS -> not necessarily true for finite temp
        p = [transition(ionObj, 0,i,temperature)[0]+transition(ionObj, 0,i,temperature)[1] for i in range(len(dE_temp))]
        # p = [ionObj.transitionIntensity(0,i, temperature) for i in range(len(dE_temp))]
        # p[0] = 0
        numLines = len(dE_temp)
        for i in range(1,numLines): 
            for j in range(i+1, numLines):
                temp = dE_temp[j]-dE_temp[i]
                dE_temp.append(temp)
                trans = transition(ionObj, i,j,temperature)
                p.append(trans[0]+trans[1])
        dE.append(dE_temp)
        amp.append(p)
    return amp, dE

def zeemanSplitC(field, wavenum, B20, B40, B43, B60, B63, B66, Jz,temperature):    
    dEphonon = 49.3
    phononAmp = 0.499
    phononSig = 0.95
    fun = []
    amp = []
    energies = []
    Bparams =  {'B20': B20, 'B40':B40,'B43': B43, 'B60': B60, 'B63':B63,'B66':B66}
    ionObj = cef.CFLevels.Bdict(ion,Bparams)
    for b in field: 
        ionObj = diagonalizeC(ionObj, ion, Jz, b, temperature)
        dE =[eval for eval in ionObj.eigenvalues] # this is the spitting if everything is in the GS -> not necessarily true for finite temp
        # tempAmp = [0 if i == 0 else transition(ionObj, 0,i,temperature)[0]+transition(ionObj, 0,i,temperature)[1] for i in range(len(dE_temp))]
        tempAmp = [ionObj.transitionIntensity(0,i, temperature) for i in range(len(dE_temp))]
        numLines = len(dE)
        for i in range(1,numLines): 
            for j in range(i+1, numLines):
                temp = dE[j]-dE[i]
                dE.append(temp)
                tempAmp.append(ionObj.transitionIntensity(i,j,temperature))
                # tempAmp.append(transition(ionObj, i,j,temperature)[0]+transition(ionObj, i,j,temperature)[1])
        wid = 1
        centers = dE
        tempfun = lorentzian(wavenum, 0, dEphonon, phononSig)
        for i in range(len(centers)):
            a = tempAmp[i]
            tempfun += lorentzian(wavenum, a, centers[i]*meVToCm, wid)
        fun.append(tempfun)
        amp.append(tempAmp)
        energies.append(centers)
    return fun, amp, energies


def zeemanSplitAB(field, wavenum, B20, B40, B43, B60, B63, B66, Jperp, temperature):    
    # assuming that x is an array
    # amp = [amp1, amp2, amp3, amp4, amp5, amp6, amp7, amp8, amp9, amp10]#, amp11, amp12, amp13, amp14, amp15, amp16]
    dEphonon = 49.3
    phononAmp = 0.499
    phononSig = 0.95
    fun = []
    amp = []
    energies = []
    Bparams =  {'B20': B20, 'B40':B40,'B43': B43, 'B60': B60, 'B63':B63,'B66':B66}
    ionObj = cef.CFLevels.Bdict(ion,Bparams)
    for b in field: 
        ionObj = diagonalizeAB(ionObj, ion, Jperp, b, temperature)
        dE =[eval for eval in ionObj.eigenvalues] # this is the spitting if everything is in the GS -> not necessarily true for finite temp
        # dE = dE[0:10] # only want to look at the bottome few lines - lets see if this works??
        # tempAmp = [0 if i == 0 else transition(ionObj, 0,i,temperature)[0]+transition(ionObj, 0,i,temperature)[2] for i in range(len(dE_temp))]
        tempAmp = [ionObj.transitionIntensity(0,i, temperature) for i in range(len(dE_temp))]
        numLines = len(dE)
        for i in range(1,numLines): 
            for j in range(i+1, numLines):
                temp = dE[j]-dE[i]
                dE.append(temp)
                # tempAmp.append(ionObj.transitionIntensity(i,j,temperature))
                tempAmp.append(transition(ionObj, i,j,temperature)[0]+transition(ionObj, i,j,temperature)[2])
        wid = 1
        centers = dE
        tempfun = lorentzian(wavenum, 0, dEphonon, phononSig)
        for i in range(len(centers)):
            a = tempAmp[i]
            tempfun += lorentzian(wavenum, a, centers[i]*meVToCm, wid)
        fun.append(tempfun)
        amp.append(tempAmp)
        energies.append(centers)
    return fun, amp, energies



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
        dE_temp =[eval for eval in evals] 
        # p = [0 if i == 0 else transition(ionObj, 0,i,temperature)[0]+transition(ionObj, 0,i,temperature)[2] for i in range(len(dE_temp))]
        p = [ionObj.transitionIntensity(0,i, temperature) for i in range(len(dE_temp))]
        numLines = len(dE_temp)
        for i in range(1,numLines): 
            # skip GS - already have those dE
            for j in range(i+1, numLines):
                temp = dE_temp[j]-dE_temp[i]
                dE_temp.append(temp)
                # p.append(transition(ionObj, i,j,temperature)[0]+transition(ionObj, i,j,temperature)[2])
                p.append(ionObj.transitionIntensity(i,j,temperature))
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
rawData = dataC.copy()


normSpec = dataC['0.001'].copy()/max(dataC['0.001'])*-1
avgSpec = normSpec.copy()
plt.figure()
for column in dataC.columns: 
    dataC[column] = max(dataC[column]) -dataC[column]
    dataC[column] = dataC[column]/(max(dataC[column])) 
    dataC[column] = dataC[column]-normSpec
    plt.plot(IRCwavenums, dataC[column])
    avgSpec = avgSpec + dataC[column]
plt.xlabel('E[cm^-1]')
plt.ylabel('Absorption')

plt.figure()
for column in dataC.columns: 
    dataC[column] = dataC[column]-avgSpec/len(dataC.columns)
    dataC[column] = dataC[column]-(sum(dataC[column])/len(dataC[column]))
    plt.plot(IRCwavenums, dataC[column])

dataC = dataC.drop(labels='0.001', axis=1) # drop this column because we used it as bg


dataC = dataC/dataC.max(axis=None)

IRCfield = [float(b) for b in dataC.columns.values]
IRCwavenums = [float(i) for i in dataC.index.values]

## for figure making
cfname = '/Users/hopeless/Desktop/LeeLab/data/CsErSe2_data/Maglab_IRMeasurement_June2022/ProcessingCode/Load2_TrimData/P3_CsEr_100_RAWAVG.dat'

# clean up first
dataC = pd.read_csv(cfname, index_col=0, skiprows=0, header=1)
dataC = dataC.dropna(axis = 0)
dataC = dataC.dropna(axis=1)
dataC = dataC.drop(labels = '-1.1', axis=1)
rawData = dataC.copy()

plt.figure()
for column in dataC.columns: 
    dataC[column] = max(dataC[column]) -dataC[column]
    dataC[column] = dataC[column]/dataC[column][151.6]
    # dataC[column] = dataC[column]/(max(dataC[column])) -normSpec
    # avgSpec = avgSpec + dataC[column]
    plt.plot(dataC[column])

normSpec = dataC['0.001'].copy()
plt.figure()
for column in dataC.columns: 
    dataC[column] = dataC[column] -normSpec
    # avgSpec = avgSpec + dataC[column]
    plt.plot(dataC[column])
# for column in dataC.columns: 
#     dataC[column] = dataC[column]-avgSpec/len(dataC.columns)
#     dataC[column] = dataC[column]-(sum(dataC[column])/len(dataC[column]))
avgSpec = np.zeros_like(normSpec)
for col in dataC.columns: 
    avgSpec += dataC[col]

avgSpec = avgSpec/len(dataC.columns)
plt.figure()
for col in dataC.columns: 
    dataC[col] = dataC[col] - avgSpec
    plt.plot(dataC[col])

plt.figure()
for col in dataC.columns: 
    dataC[col] = dataC[col] - min(dataC[col])
    plt.plot(dataC[col])

dataC = dataC.drop(labels='0.001', axis=1) # drop this column because we used it as bg


dataC = dataC/dataC.max(axis=None)

IRCfield = [float(b) for b in dataC.columns.values]
IRCwavenums = [float(i) for i in dataC.index.values]

##lets try an average subtraction
cfname = '/Users/hopeless/Desktop/LeeLab/data/CsErSe2_data/Maglab_IRMeasurement_June2022/ProcessingCode/Load2_TrimData/P3_CsEr_100_RAWAVG.dat'

# clean up first
dataC = pd.read_csv(cfname, index_col=0, skiprows=0, header=1)
dataC = dataC.dropna(axis = 0)
dataC = dataC.dropna(axis=1)
dataC = dataC.drop(labels = '-1.1', axis=1)
rawData = dataC.copy()

# plt.figure()
# normSpec = dataC['-1.1'].copy()
# normSpec = max(normSpec)-normSpec
avgSpec = np.zeros_like(dataC['0.001'])
for column in dataC.columns: 
    dataC[column] = max(dataC[column]) -dataC[column]
    avgSpec = avgSpec + dataC[column]
    # plt.plot(dataC[column])

avgSpec = np.zeros_like(normSpec)
for col in dataC.columns: 
    avgSpec += dataC[col]

avgSpec = avgSpec/len(dataC.columns)
plt.figure()
for col in dataC.columns: 
    dataC[col] = dataC[col]/avgSpec
    plt.plot(dataC[col])



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
calc_field = np.arange(0,20, .2)
calc_wavenums = np.arange(0,300, 1)
temperature = 2
kBT = kB*temperature
Jperp = 0

simulated_spec_C, ampC, arrC = zeemanSplitC(calc_field, calc_wavenums, B20, B40, B43, B60, B63, B66, Jz, temperature)
simulated_spec_C  = np.array(simulated_spec_C)
ampC = np.array(ampC).T/np.array(ampC).max(axis=None)
arrC = np.array(arrC).T

# # calc_wavenums = np.linspace(0,400, 400)
simulated_spec_IR_B, ampB, arrB = zeemanSplitAB(calc_field, calc_wavenums, B20, B40, B43, B60, B63, B66, Jperp, temperature)
simulated_spec_IR_B  = np.array(simulated_spec_IR_B)
ampB = np.array(ampB).T
arrB = np.array(arrB).T

# spectroscopy_calc_fname
spec_calc_fname = 'spectroscopy_fgr_all_2025Aug05.h5'
with h5py.File(spec_calc_fname, 'w') as hdf:
    hdf.create_dataset('linesC', data = arrC)
    hdf.create_dataset('linesB', data = arrB)
    hdf.create_dataset('ampC', data = ampC)
    hdf.create_dataset('ampB', data = ampB)
    # hdf.create_dataset('linesC_nomft', data = arrC_nomft)
    # hdf.create_dataset('linesB_nomft', data = arrB_nomft)
    # hdf.create_dataset('ampC_nomft', data = ampC_nomft)
    # hdf.create_dataset('ampB_nomft', data = ampB_nomft)
    hdf.create_dataset('calc_field', data = calc_field)
    hdf.create_dataset('calc_wavenums', data= calc_wavenums)
    hdf.create_dataset('simulated_C', data = simulated_spec_C)
    hdf.create_dataset('simulated_IR_B', data = simulated_spec_IR_B)
    # hdf.create_dataset('simulated_IR_C', data = simulated_spec_IR_C)
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
    hdf.attrs['notes'] = 'recalculated with fgr Jul 30 2025 '

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

field = [float(b) for b in fitData.columns.values]
wavenum = [float(i) for i in fitData.index.values]
# now do fit
temperature = 7
model = lmfit.Model(zeemanSplitC_raman, independent_vars=['field', 'wavenum'])
params = model.make_params()



params['B20'].set(value= 1e-6, min = -.06, max = 0.06)
params['B40'].set(value= 1e-6, min = -.06, max = 0.06)
params['B43'].set(value= 1e-6, min = -.06, max = 0.06)
params['B60'].set(value= 1e-6, min = -.06, max = 0.06)
params['B63'].set(value= 1e-6, min = -.06, max = 0.06)
params['B66'].set(value= 1e-6, min = -.06, max = 0.06)
params['Jz'].set(value = 0, vary=False)
params['temperature'].set(value = 7, vary = False)

z = np.array(fitData.to_numpy()) # gotta do it twice with tuples :((((
z = z.T

result = model.fit(z, field=field, wavenum=wavenum, params =params, method='basin_hopping')

print(result.fit_report())

######################################
# plot the fuck around fits
B20 =  0.02876578 # 2.5225e-04 (0.88%) (init = -0.03001018)
B40 =  5.4242e-05 # 9.5667e-07 (1.76%) (init = -0.00011757)
B43 = -0.00119712 # 2.7421e-05 (2.29%) (init = -0.00384833)
B60 =  8.2943e-07 # 7.2025e-09 (0.87%) (init = 3.9e-07)
B63 = -1.0528e-06 # 1.6958e-07 (16.11%) (init = -1.08e-06)
B66 =  2.3953e-05 # 7.9272e-08 (0.33%) (init = 3.49e-06)
Jz =  -0.02026981 # 2.9081e-04 (1.43%) (init = 0)




calc_field = np.linspace(0,15,100)
calc_wave = np.linspace(0,10,10)
sim, ampC, arrC = zeemanSplitC(calc_field, calc_wave, B20, B40, B43, B60, B63, B66, Jz, temperature)
arrC = np.array(arrC)
# arrC = arrC*meVToCm
arrC = arrC.T
arrC = arrC[1:-1]

ampC = np.array(ampC)
ampC = ampC.T/ampC.max(axis=None)
ampC = ampC[1:-1]
ampC = ampC/ampC.max(axis = None)

# import matplotlib.colors as mcolors
cmap = mcolors.ListedColormap(['cyan'])
field = [float(b) for b in ramanData.columns.values]
wavenums = [float(i) for i in ramanData.index.values]
fig, ax = plt.subplots()
plt.contourf(field, wavenums, ramanData,100, cmap = 'Oranges')
plt.xlim(0,14)
plt.ylim(0,120)
plt.colorbar()
# plt.clim(-0.3, 1)
plt.xlabel('Field (T)')
plt.ylabel('Energy (cm$^{-1}$)')

plt.title('CsErSe2 H||C with overlayed  calclines\n B20: '+ str(B20)+' B40: '+str(B40)+' B43: ' +str(B43)+ '\n B60: ' +str(B60) + ' B63: ' + str(B63)+ ' B66: ' + str(B66)+'\n temp = '+str(temperature))
for i in range(1, len(arrC)):
    # if i<16: 
    # plt.plot(calc_field, arrC[i]*meVToCm, 'c', alpha=1, linewidth= .7)
# if i>=16:  
    alphas = ampC[i]
    # Create a LineCollection
    points = np.array([calc_field, arrC[i]*meVToCm]).T.reshape(-1, 1, 2)
    segments = np.concatenate([points[:-1], points[1:]], axis=1)
    lc = LineCollection(segments, cmap=cmap, alpha=alphas)
    ax.add_collection(lc)
    ax.autoscale()
plt.show()

cmap = mcolors.ListedColormap(['cyan'])
field = IRBfield#[float(b) for b in ramanData.columns.values]
wavenums = IRBwavenums#[float(i) for i in ramanData.index.values]
fig, ax = plt.subplots()
plt.contourf(field, wavenums, dataB,100, cmap = 'Oranges')
plt.xlim(0,14)
plt.ylim(0,120)
plt.colorbar()
# plt.clim(-0.3, 1)
plt.xlabel('Field (T)')
plt.ylabel('Energy (cm$^{-1}$)')

plt.title('CsErSe2 H||C with overlayed  calclines\n B20: '+ str(B20)+' B40: '+str(B40)+' B43: ' +str(B43)+ '\n B60: ' +str(B60) + ' B63: ' + str(B63)+ ' B66: ' + str(B66)+'\n temp = '+str(temperature))
for i in range(len(arrB)):
    if i<16: 
        plt.plot(calc_field, arrB[i], 'c', alpha=1, linewidth= .7)
    if i>=16:  
        alphas = ampB[i]
        # Create a LineCollection
        points = np.array([calc_field, arrB[i]]).T.reshape(-1, 1, 2)
        segments = np.concatenate([points[:-1], points[1:]], axis=1)
        lc = LineCollection(segments, cmap=cmap, alpha=alphas)
        ax.add_collection(lc)
        ax.autoscale()
plt.show()


# temp test
Tarr = [0.01, 0.025, 0.1, 0.25, 0.4, 0.5, 0.75, 1, 2]
arrC = []
arrC_nomft = []
for temperature in Tarr: 
    calc_field = np.linspace(0,15,1000)
    ampC_1, arrC_1 = zeemanSplitLinesC(np.linspace(0,15,1000), B20, B40, B43, B60, B63, B66, Jz, temperature)
    arrC_1 = np.array(arrC_1)
    arrC_1 = arrC_1*meVToCm
    arrC_1 = arrC_1.T
    arrC.append(arrC_1.tolist())

    ampC_no_mft, arrC_no_mft = zeemanSplitLinesC(np.linspace(0,15,1000), B20, B40, B43, B60, B63, B66, 0, temperature)
    arrC_no_mft = np.array(arrC_no_mft)
    arrC_no_mft = arrC_no_mft*meVToCm
    arrC_no_mft = arrC_no_mft.T
    arrC_nomft.append(arrC_no_mft.tolist())




plt.figure()
for j in range(len(Tarr)): 
    arrC_1 = arrC[j]
    arrC_1_no_mft = arrC_nomft[j]
    plt.subplot(3, 3, j+1) 
    plt.title(str(Tarr[j])+'K')
    for i in range(len(arrC_1)):
        if i<16: 
            if i ==0: 
                plt.plot(calc_field, arrC_1[i], 'b', label = str(Tarr[j]))
                plt.plot(calc_field, arrC_1_no_mft[i], 'b', label = str(Tarr[j]) +'no mft')
            plt.plot(calc_field, arrC_1[i], 'b-')
            plt.plot(calc_field, arrC_1_no_mft[i], 'b--')

plt.show()


with h5py.File('zeeman_c_temperature.h5', 'w') as f:
    f.create_dataset('calc_field', data=calc_field)
    f.create_dataset('Tarr', data=Tarr)
    f.create_dataset('arrC', data=np.array(arrC))
    f.create_dataset('arrC_nomft', data=np.array(arrC_nomft))
# plt.legend()
########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
# ookay, let's make some line plots for fig 1
# panel 1, hline at evals
# panel 2, ab plane
# panel 3, c axis
# need new fn that returns no norm
# Define calculation functions
def splitLinesC(field, B20, B40, B43, B60, B63, B66, Jz, temperature):
    Bparams =  {'B20': B20, 'B40':B40,'B43': B43, 'B60': B60, 'B63':B63,'B66':B66}
    ionObj = cef.CFLevels.Bdict(ion,Bparams)
    kBT = temperature*kB
    dE = []
    for b in field: 
        h = newH(ionObj, b, Jz, temperature)
        JdotB = muB*(h*cef.Operator.Jz(ionObj.J))*cef.LandeGFactor(ion)
        H = np.sum([a*b for a,b in zip(ionObj.O, ionObj.B)], axis=0)
        ionObj.diagonalize(H + JdotB.O) # this is just H = Hcef + Hmag
        evals = ionObj.eigenvaluesNoNorm
        dE_temp =[eval for eval in evals] # this is the spitting if everything is in the GS -> not necessarily true for finite temp
        dE.append(dE_temp)
    return dE

def splitLinesAB(field, B20, B40, B43, B60, B63, B66, Jperp, temperature):
    Bparams =  {'B20': B20, 'B40':B40,'B43': B43, 'B60': B60, 'B63':B63,'B66':B66}
    ionObj = cef.CFLevels.Bdict(ion,Bparams)
    kBT = temperature*kB
    dE = []
    for b in field: 
        h = newHAB(ionObj, b, Jperp, temperature)
        JdotB = muB*(h*cef.Operator.Jy(ionObj.J))*cef.LandeGFactor(ion)
        H = np.sum([a*b for a,b in zip(ionObj.O, ionObj.B)], axis=0)
        ionObj.diagonalize(H + JdotB.O) # this is just H = Hcef + Hmag
        evals = ionObj.eigenvaluesNoNorm
        dE_temp =[eval for eval in evals] # this is the spitting if everything is in the GS -> not necessarily true for finite temp
        dE.append(dE_temp)
    return dE

# Simulate and save for each temperature
temperatures = [10]
field = np.linspace(0, 100, 1000)

for T in temperatures:
    ZFevals = splitLinesC([0], B20, B40, B43, B60, B63, B66, 0, T)
    ABevals = np.array(splitLinesAB(field, B20, B40, B43, B60, B63, B66, Jperp, T)).T
    ABevals_nomft = np.array(splitLinesAB(field, B20, B40, B43, B60, B63, B66, 0, T)).T
    Cevals = np.array(splitLinesC(field, B20, B40, B43, B60, B63, B66, Jz, T)).T
    Cevals_nomft = np.array(splitLinesC(field, B20, B40, B43, B60, B63, B66, 0, T)).T

    filename = f'EvsH_mft_T{T}K_1015May12.h5'
    with h5py.File(filename, 'w') as hdf:
        hdf.create_dataset('ZFevals', data=ZFevals)
        hdf.create_dataset('ABevals', data=ABevals)
        hdf.create_dataset('Cevals', data=Cevals)
        hdf.create_dataset('ABevals_nomft', data=ABevals_nomft)
        hdf.create_dataset('Cevals_nomft', data=Cevals_nomft)
        hdf.create_dataset('field', data=field)
        hdf.attrs['B20'] = B20
        hdf.attrs['B40'] = B40
        hdf.attrs['B43'] = B43
        hdf.attrs['B60'] = B60
        hdf.attrs['B63'] = B63
        hdf.attrs['B66'] = B66
        hdf.attrs['Jz'] = Jz
        hdf.attrs['Jperp'] = Jperp
        hdf.attrs['temperature'] = T

#######################################
#######################################
#######################################
#######################################
#######################################
#######################################
#######################################
#######################################
# take a sec and calculate a couple temps w J
# just around 2K, to see how that looks
fname2K = '/Users/hopeless/Desktop/LeeLab/data/CsErSe2_data/MPMS/CsErSe2/CsErSe2_HParC_MvsH_2K.txt'

Mdata2K = np.genfromtxt(fname2K, delimiter=',',  unpack=True, skip_header=1)
temps = [1.8,1.9,2.0,2.1,2.1]

H = np.linspace(0, 8, 30)  # Field values
J_values = [-2.11e-3, -2.05e-3, -1.9e-3, -2.53e-3, -2.63e-3]  # J values

# Initialize HDF5 file
with h5py.File('jjz_temperature_test_dec_19_model.h5', 'w') as hdf:
    hdf.create_dataset('H', data=H)
    hdf.create_dataset('Mdata2K', data=Mdata2K)
    hdf.create_dataset('temps', data=temps)
    
    for J in J_values:
        tempMag = []
        for temp in temps:
            mag = MFTmagC(ionObj, H, J, temp)  # Replace `None` with the actual ion object if needed
            tempMag.append(mag)
        tempMag = np.array(tempMag)
        group = hdf.create_group(f'J_{J}')
        group.create_dataset('mag', data=tempMag)
        group.attrs['J'] = J


#########################################################################################
#########################################################################################
#########################################################################################
#########################################################################################

# now I want to check the eigenvectors before and after the sero energy crossing - is it
# a crossing? or is it an avoided crossing

# ass B, J are already loaded into a terminal

# write a function to spit out eigenvectors (in latex form) at a given field

def printLaTexEigenvectors(ionObj, field, precision = 4):
    '''prints eigenvectors and eigenvalues in the output that Latex can read'''
    h = newH(ionObj, field, Jz, temperature)
    JdotB = muB*(h*cef.Operator.Jz(ionObj.J))*cef.LandeGFactor(ion)
    H_cef = np.sum([a*b for a,b in zip(ionObj.O, ionObj.B)], axis=0)
    ionObj.diagonalize(H_cef + JdotB.O) # this is just H = Hcef + Hmag
    evals = ionObj.eigenvalues 
    
    print('\\begin{table*}\n\\caption{Eigenvectors and Eigenvalues...}')
    print('\\begin{ruledtabular}')
    numev = len(ionObj.eigenvalues)
    print('\\begin{tabular}{c|'+'c'*numev+'}')
    if numev % 2 == 1:
        print('E (meV) &'+' & '.join(['$|'+str(int(kk))+'\\rangle$' for kk in 
                    np.arange(-(numev-1)/2,numev/2)])
            +' \\tabularnewline\n \\hline ')
    else:
        print('E (meV) &'+
            ' & '.join(['$| -\\frac{'+str(abs(kk))+'}{2}\\rangle$' if kk <0
                        else '$| \\frac{'+str(abs(kk))+'}{2}\\rangle$'
                        for kk in np.arange(-(numev-1),numev,2)])
            +' \\tabularnewline\n \\hline ')
    sortinds = ionObj.eigenvalues.argsort()
    sortEVal= np.around(ionObj.eigenvalues[sortinds],3)
    sortEVec= np.around(ionObj.eigenvectors[sortinds],precision)
    for i in range(len(sortinds)):
        print(format(ionObj._Re(sortEVal[i]), '.3f'),'&', 
            ' & '.join([str(eevv) for eevv in ionObj._Re(sortEVec[i])]), '\\tabularnewline')
    print('\\end{tabular}\\end{ruledtabular}')
    print('\\label{flo:Eigenvectors}\n\\end{table*}')



######################################################################
######################################################################
######################################################################
######################################################################
######################################################################
######################################################################
######################################################################
######################################################################
######################################################################
######################################################################

# now lets calc the magnetotropic coeff
Bparams =  {'B20': B20, 'B40':B40,'B43': B43, 'B60': B60, 'B63':B63,'B66':B66}
cser = MOO(ion, Bparams, Jperp, Jz, q)

k_temp_arr =[]
tau_temp_arr =[]
t_arr = [2,7,15, 35, 45,55,70, 90]
phi = 0
theta = np.linspace(-np.pi/64, np.pi/64, 5)
Hval = np.linspace(0.1,20,80)
dth = 0.001
for temperature in t_arr:
    k_arr = []
    tau_arr =[]
    for h in Hval: 
        k,t = cser.magnetoTropic(h, theta, phi, dth, temperature)
        k_arr.append(k)
        tau_arr.append(t)
    k_temp_arr.append(k_arr)
    tau_temp_arr.append(tau_arr)


n = len(Hval)
colors = plt.cm.turbo(np.linspace(0,1,n))
k_fig, k_ax = plt.subplots()
# tau_fig, tau_ax = plt.subplots()
deg = np.linspace(0, 180, len(theta))
k_arr = k_temp_arr[0]
for i, h in enumerate(Hval): 
    k_ax.plot(deg, k_arr[i], label = str(h)+'T', color=colors[i])
    # tau_ax.plot(deg, tau_arr[i],label = str(h)+'T', color=colors[i])
k_ax.set_title(" Magnetotro")
k_ax.set_xlabel("angle [deg]")
k_ax.set_ylabel("k")


# Save to HDF5
filename = 'magnetotropic_phi_45_20jun2025.h5'

with h5py.File(filename, 'w') as f:
    f.create_dataset('t_arr', data=np.array(t_arr))
    f.create_dataset('Hval', data=np.array(Hval))
    f.create_dataset('theta', data=np.array(theta))
    
    # save k_temp_arr and tau_temp_arr as 3D arrays: (temperature index, H index, theta index)
    k_temp = np.array([np.array(k_arr) for k_arr in k_temp_arr]) # shape (len(t_arr), len(Hval), len(theta))
    tau_temp = np.array([np.array(tau_arr) for tau_arr in tau_temp_arr])
    
    f.create_dataset('k_temp', data=k_temp)
    f.create_dataset('tau_temp', data=tau_temp)

print(f'Saved HDF5 to {filename}')




## for peak picking

import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import find_peaks

def interactive_peak_selector(data, wavenums, field, height=0.1, distance=3, tol=1.0):
    """
    Interactively select peaks from a 2D false-color map.
    Uses 1D peak detection (find_peaks) on each spectrum at each field.
    Returns picked (field, wavenumber) coordinates.
    """
    peak_coords = []
    plt.ioff()
    # Loop over each field (column of data)
    for j, B in enumerate(field):
        spectrum = data[:, j]
        peaks, _ = find_peaks(spectrum, height=height, distance=distance)

        # Convert peak indices to coordinates
        for i in peaks:
            peak_coords.append((B, wavenums[i]))

    peak_coords = np.array(peak_coords)

    # === Initialize plot ===
    fig, ax = plt.subplots(figsize=(8, 5))
    c = ax.contourf(field, wavenums, data, 100, cmap='viridis')
    plt.colorbar(c, ax=ax, label='Intensity')

    ax.set_xlabel('Field (T)')
    ax.set_ylabel('Wavenumber (cm$^{-1}$)')
    ax.set_title('Click red dots to select peaks (close window when done)')

    # Plot red dots for auto peaks
    ax.plot(peak_coords[:, 0], peak_coords[:, 1], 'ro', markersize=4, label='Detected peaks')

    picked = []

    def onclick(event):
        if event.inaxes != ax:
            return
        click_x, click_y = event.xdata, event.ydata

        # Find closest peak within tolerance
        dists = np.sqrt((peak_coords[:, 0] - click_x)**2 + (peak_coords[:, 1] - click_y)**2)
        nearest = np.argmin(dists)
        if dists[nearest] < tol:
            peak = tuple(peak_coords[nearest])
            if peak not in picked:
                picked.append(peak)
                ax.plot(*peak, 'wx', markersize=6)  # mark accepted peak
                fig.canvas.draw()

    cid = fig.canvas.mpl_connect('button_press_event', onclick)
    plt.tight_layout()
    plt.show()
    fig.canvas.mpl_disconnect(cid)

    # === Summary ===
    print(f"\nYou picked {len(picked)} peaks:")
    for f, w in picked:
        print(f"  Field = {f:.2f} T, Wavenumber = {w:.2f} cm⁻¹")
    plt.ion()
    return picked


# cut data to only show CEF lines
data = np.array(ramanData)
wavenums = np.array(ramanWavenums)
field = np.array(ramanField)
cutoff = 120 # want to cut above 300, because there's a line from the 5/2 manifold
# Find indices where wavenumber is <= 300
mask = wavenums <= cutoff 

# Apply mask to both wavenums and data
wavenums = wavenums[mask]
data = data[mask, :]  # Keep all fields, just cut rows

# start with main lines
main = interactive_peak_selector(data, wavenums, field, height=0.1, distance=3, tol=1.0)

## now get IR peaks 

data = dataB.to_numpy()
wavenums = np.array(IRBwavenums)
field = np.array(IRBfield)

peak_coords = 
np.savetxt('peak_coords_raman_high_energy.csv',np.array(peak_coords), delimiter = ',')


###########################################################################
# let's make some simulated spectra with fermi's golden rule
# def zeemanSplitC(field, wavenum, B20, B40, B43, B60, B63, B66, Jz,temperature):    
#     # assuming that x is an array
#     # amp = [amp1, amp2, amp3, amp4, amp5, amp6, amp7, amp8, amp9, amp10]#, amp11, amp12, amp13, amp14, amp15, amp16]
#     amp = []
#     dEphonon = 49.3
#     phononAmp = .5
#     phononSig = 0.95
#     fun = []
#     Bparams =  {'B20': B20, 'B40':B40,'B43': B43, 'B60': B60, 'B63':B63,'B66':B66}
#     ionObj = cef.CFLevels.Bdict(ion,Bparams)
#     # dE = dE[0:10] # only want to look at the bottome few lines - lets see if this works??
#     for b in field: 
#         evals = diagonalizeC(ionObj, ion, Jz, b, temperature)
#         dE =[eval for eval in evals] # this is the spitting if everything is in the GS -> not necessarily true for finite temp
#         # dE = dE[0:10] # only want to look at the bottome few lines - lets see if this works??
#         tempAmp = [transition(ionObj, 0,i,temperature)[0]+transition(ionObj, 0,i,temperature)[1] for i in range(len(dE))]
#         tempAmp[0] = 0
#         Z = [np.exp(-Ei/kBT) for Ei in dE]
#         Z = sum(Z)
#         p = [1/Z*np.exp(-Ei/kBT) for Ei in dE]
#         numLines = len(dE)
#         for i in range(1,numLines): 
#             for j in range(i+1, numLines):
#                 temp = dE[j]-dE[i]
#                 dE.append(temp)
#                 tempAmp.append(transition(ionObj, i,j,temperature)[0]+transition(ionObj, i,j,temperature)[1])
#         wid = 1
#         centers = dE
#         tempfun = lorentzian(wavenum, 0, dEphonon, phononSig)
#         # tempfun = lorentzian(wavenum, 0, dEphonon, phononSig)
#         for i in range(len(centers)):
#             a = tempAmp[i]
#             tempfun += lorentzian(wavenum, a, centers[i]*meVToCm, wid)
#         fun.append(tempfun)
#     return fun

sim_field = np.linspace(0,18,180)
sim_wave = np.linspace(0,120,120*5)
sim_data = zeemanSplitC(sim_field, sim_wave, B20, B40, B43, B60, B63, B66, Jz, 10)

plt.figure()
plt.contourf(sim_field, sim_wave, np.array(sim_data).T/np.array(sim_data).max(axis=None), 100)
# plt.contourf(ramanField, ramanWavenums, np.array(sim_data).T/10-np.array(ramanData), 100)
# plt.figure()
# plt.contourf(ramanField, ramanWavenums, ramanData, 100)


import h5py

# Save to HDF5
with h5py.File("fgr_HparC_simulation.h5", "w") as f:
    f.create_dataset("ramanField", data=sim_field)
    f.create_dataset("ramanWavenums", data=sim_wave)
    f.create_dataset("sim_data", data=np.array(sim_data))
    # f.create_dataset("ramanData", data=np.array(ramanData))


# track intensities
H = np.linspace(0,18,200)
inten = []
for h in H: 
    evals = diagonalizeC(ionObj, ion, Jz, b, temperature)
    inten.append(ionObj.transitionIntensity(2,1, temperature))

plt.plot(H, inten)


# def zeemanSplitLinesC(field, B20, B40, B43, B60, B63, B66, Jz):     
#     Bparams =  {'B20': B20, 'B40':B40,'B43': B43, 'B60': B60, 'B63':B63,'B66':B66}
#     ionObj = cef.CFLevels.Bdict(ion,Bparams)
#     amp = []
#     dE = []
#     temperature = 4
#     for b in field: 
#         evals = diagonalizeC(ionObj, ion, Jz, b, temperature)
#         dE_temp =[eval for eval in evals] # this is the spitting if everything is in the GS -> not necessarily true for finite temp
#         tempAmp = [ionObj.transitionIntensity(0,i, temperature) for i in range(len(dE_temp))]
#         # tempAmp = [transition(ionObj, 0,i,temperature)[0]+transition(ionObj, 0,i,temperature)[1] for i in range(len(dE_temp))]
#         numLines = len(dE_temp)
#         for i in range(1,numLines): 
#             # skip GS - already have those dE
#             for j in range(i+1, numLines):
#                 temp = dE_temp[j]-dE_temp[i]
#                 dE_temp.append(temp)
#                 # tempAmp.append()
#                 tempAmp.append(ionObj.transitionIntensity(i,j, temperature))
#         dE.append(dE_temp)
#         amp.append(tempAmp)
#     return amp, dE


calc_field = np.linspace(0,2, 100)
amp, dE = zeemanSplitLinesAB(calc_field, B20, B40, B43, B60, B63, B66, Jz)



arrC = np.array(dE)
arrC = arrC*meVToCm
arrC = arrC.T

ampC = np.array(amp)
ampC = ampC.T
ampC = ampC/ampC.max(axis=None)


cmap = mcolors.ListedColormap(['cyan'])
field = [float(b) for b in ramanData.columns.values]
wavenums = [float(i) for i in ramanData.index.values]
fig, ax = plt.subplots()
plt.contourf(field, wavenums, ramanData,100, cmap = 'Oranges')
plt.xlim(0,14)
plt.ylim(0,120)
plt.colorbar()
# plt.clim(-0.3, 1)
plt.xlabel('Field (T)')
plt.ylabel('Energy (cm$^{-1}$)')

# plt.title('CsErSe2 H||C with overlayed  calclines\n B20: '+ str(B20)+' B40: '+str(B40)+' B43: ' +str(B43)+ '\n B60: ' +str(B60) + ' B63: ' + str(B63)+ ' B66: ' + str(B66)+'\n temp = '+str(temperature))
for i in range(len(arrC)):
    alphas = ampC[i]
    # Create a LineCollection
    points = np.array([calc_field, arrC[i]]).T.reshape(-1, 1, 2)
    segments = np.concatenate([points[:-1], points[1:]], axis=1)
    lc = LineCollection(segments, cmap=cmap, alpha=alphas)
    ax.add_collection(lc)
    ax.autoscale()
plt.show()

plt.figure()
for i in range(1,16): 
    plt.plot(calc_field, ampC[i], label = str(i) +'to ground')

plt.xlabel('H[T]')
plt.ylabel('transition probability')

plt.figure()
for i in range(16,31): 
    plt.plot(calc_field, ampC[i], label = str(i) +'to first')

plt.legend()
plt.xlabel('H[T]')
plt.ylabel('transition probability')


binaryAmp = (ampC > 0).astype(int)

np.savetxt('binary_amp.csv', binaryAmp)
np.savetxt('amplitude.csv', ampC)
np.savetxt('energies.csv', arrC)
np.savetxt('calc_field.csv', calc_field)


cmap = mcolors.ListedColormap(['cyan'])
field = [float(b) for b in ramanData.columns.values]
wavenums = [float(i) for i in ramanData.index.values]
fig, ax = plt.subplots()
plt.contourf(field, wavenums, ramanData,100, cmap = 'Oranges')
plt.xlim(0,14)
plt.ylim(0,120)
plt.colorbar()
# plt.clim(-0.3, 1)
plt.xlabel('Field (T)')
plt.ylabel('Energy (cm$^{-1}$)')

# plt.title('CsErSe2 H||C with overlayed  calclines\n B20: '+ str(B20)+' B40: '+str(B40)+' B43: ' +str(B43)+ '\n B60: ' +str(B60) + ' B63: ' + str(B63)+ ' B66: ' + str(B66)+'\n temp = '+str(temperature))
for i in range(len(arrC)):
    alphas =  binaryAmp[i]
    # Create a LineCollection
    points = np.array([calc_field, arrC[i]]).T.reshape(-1, 1, 2)
    segments = np.concatenate([points[:-1], points[1:]], axis=1)
    lc = LineCollection(segments, cmap=cmap, alpha=alphas)
    ax.add_collection(lc)
    ax.autoscale()
plt.show()




## learned that it was actually a-axis
## gotta redo the fit
## time to cry just a little bit
# check magnetization rq
# JJperp = 0
q = 6
ion = 'Er3+'
Bparams =  {'B20': B20, 'B40':B40,'B43': B43, 'B60': B60, 'B63':B63,'B66':B66}
kys = MOO(ion, Bparams, Jperp, Jz, q)

temps = [ 2, 6, 10]
result = []
for i, T in enumerate(temps):
    res = kys.calc_over_field(T, np.linspace(0,7, 50))
    result.append(res)

# Create figures and axes
# magC_fig, magC_ax = plt.subplots()
magAB_fig, magAB_ax = plt.subplots()
# linesC_fig, linesC_ax = plt.subplots()
# linesAB_fig, linesAB_ax = plt.subplots()
# dmdhC_fig, dmdhC_ax = plt.subplots()
# dmdhAB_fig, dmdhAB_ax = plt.subplots()
n = len(temps)
colors = plt.cm.jet(np.linspace(0,1,n))

for i, T in enumerate(temps): 
    dataA = result[i]['x']
    dataB = result[i]['y']
    # dataC = result[i]['z']
    
    # Magnetization
    # magC_ax.plot(dataC['H'], dataC['m'].T[2], label=f'T={T}', color = colors[i])
    magAB_ax.plot(dataA['H'], dataA['m'].T[0], label=f'T={T}', color = colors[i])
    magAB_ax.plot(dataB['H'], -1*dataB['m'].T[1], label = f'T{T}', color = colors[i], linestyle = '--')
    # dM/dH
    # dmdhC_ax.plot(dataC['H'], np.gradient(dataC['m'].T[2], dataC['H']), label=f'T={T}', color = colors[i])
    # dmdhAB_ax.plot(dataA['H'], np.gradient(dataA['m'].T[0], dataA['H']), label=f'T={T}', color = colors[i])

# Eigenvalues - plot at low temp
dataA = result[0]['x']
dataC = result[0]['z']
for j in range(dataC['evals'].shape[1]): 
    linesC_ax.plot(dataC['H'], dataC['evals'][:, j], label=f'T={T}, lvl={j}', color='black')
    linesAB_ax.plot(dataA['H'], dataA['evals'][:, j], label=f'T={T}, lvl={j}', color='black')

# Optional: add legends
magC_ax.set_title("H||c")
magC_ax.set_xlabel("H [T]")
magC_ax.set_ylabel("M (μB)")

magAB_ax.set_title("H||ab")
magAB_ax.set_xlabel("H[T]")
magAB_ax.set_ylabel("M (μB)")

dmdhC_ax.set_title("H||c")
dmdhC_ax.set_xlabel("H[T]")
dmdhC_ax.set_ylabel("dM/dH (μB/T)")

dmdhAB_ax.set_title("H||ab")
dmdhAB_ax.set_xlabel("H[T]")
dmdhAB_ax.set_ylabel("dM/dH (μB/T)")

linesC_ax.set_title("H||c")
linesC_ax.set_xlabel("H[T]")
linesC_ax.set_ylabel("Energy [meV]")

linesAB_ax.set_title("H||ab")
linesAB_ax.set_xlabel("H[T]")
linesAB_ax.set_ylabel("Energy [meV]")

plt.show()

## get peaks code
peak_coords = []
for j, B in enumerate(field):
    spectrum = data[:, j]
    peaks, _ = find_peaks(spectrum, height=height, distance=distance)
    for i in peaks:
        peak_coords.append((B, wavenums[i]))

peak_coords = np.array(peak_coords)
rejected = []

# === Plot ===
fig, ax = plt.subplots(figsize=(8, 5))
c = ax.contourf(field, wavenums, data, 100, cmap='viridis')
plt.colorbar(c, ax=ax, label='Intensity')

ax.plot(peak_coords[:, 0], peak_coords[:, 1], 'ro', markersize=4, label='Auto Peaks')
ax.set_xlabel('Field (T)')
ax.set_ylabel('Wavenumber (cm$^{-1}$)')
ax.set_title('Click red dots to REJECT peaks (close window when done)')


fname = 'peak_coords_highRaman.csv'
np.savetxt(fname, peak_coords)


## can i actually write down the matrix elements??
import numpy as np
from sympy.physics.wigner import clebsch_gordan

def dipole_operator_q(J, q, Jp=None, reduced_element=1.0):
    """
    Construct D_q dipole operator in the |J, m> basis.
    If Jp != J, it acts between different J multiplets.
    """
    if Jp is None:
        Jp = J  # default: act within the same multiplet

    dim_in = int(2*J+1)
    dim_out = int(2*Jp+1)
    Dq = np.zeros((dim_out, dim_in), dtype=complex)

    m_vals = np.arange(-J, J+1, 1)
    mp_vals = np.arange(-Jp, Jp+1, 1)

    for i, m in enumerate(m_vals):
        for j, mp in enumerate(mp_vals):
            cg = clebsch_gordan(J, m, 1, q, Jp, mp)
            if cg != 0:
                Dq[j, i] = reduced_element * cg

    return Dq

# Define J = 15/2 basis
J = 15/2
mJ = np.arange(-J, J+1, 1)

# Construct Stevens operators O_k^q as matrices
O20 = stevens_operator(J, k=2, q=0)
O40 = stevens_operator(J, k=4, q=0)
# ... other O_k^q ...

# Build crystal field Hamiltonian
Hcf = B20*O20 + B40*O40 + B43*O43 + B60*O60 + B63*O63 + B66*O66

# Diagonalize to get CF eigenstates
eigvals, eigvecs = np.linalg.eigh(Hcf)

# Pick i = ground, f = excited CF state
psi_i = eigvecs[:, 0]
psi_f = eigvecs[:, k]  # some excited state

# Construct dipole operator D_q (rank-1 spherical tensor)
Dq = dipole_operator_q(J, q=0 or ±1)

# Define intermediate state manifold (e.g., J=11/2 basis)
# and repeat similarly to get |n⟩ and energy E_n

# Compute matrix elements
M_in = np.dot(psi_n.conj().T, Dq @ psi_i)
M_fout = np.dot(psi_f.conj().T, Dq @ psi_n)

# Assemble Kramers-Heisenberg term
A_fi = np.sum(M_fout * M_in / (E_i + hw_in - E_n + 1j*Gamma))

from sympy.physics.wigner import clebsch_gordan
import numpy as np

def dipole_J15to11(q, reduced_element=1.0):
    J1 = 15/2
    J2 = 11/2
    m1 = np.arange(-J1, J1+1, 1)
    m2 = np.arange(-J2, J2+1, 1)
    D = np.zeros((len(m2), len(m1)), dtype=complex)

    for i, m in enumerate(m1):
        for j, mp in enumerate(m2):
            cg = clebsch_gordan(J1, m, 1, q, J2, mp)
            if cg != 0:
                D[j, i] = reduced_element * cg
    return D

# Dipole matrices from J=15/2 to J=11/2
D0 = dipole_J15to11(0)
Dp = dipole_J15to11(1)
Dm = dipole_J15to11(-1)

# Raman effective operators in J=15/2 space
R00 = D0.conj().T @ D0
Rpp = Dp.conj().T @ Dp
Rmm = Dm.conj().T @ Dm
Rpm = Dp.conj().T @ Dm
Rmp = Dm.conj().T @ Dp
# etc...

# Then for specific CF eigenstates:
I = np.abs(psi_f.conj().T @ R00 @ psi_i)**2

import numpy as np
from sympy.physics.wigner import clebsch_gordan
from sympy import Rational
import matplotlib.pyplot as plt

# --- Dipole construction ---
def dipole_J15to11(q, reduced_element=1.0):
    J1 = Rational(15, 2)
    J2 = Rational(11, 2)
    m1 = [Rational(i, 2) for i in range(-15, 16, 2)]  # 16 mJ values for J=15/2
    m2 = [Rational(i, 2) for i in range(-11, 12, 2)]  # 12 mJ values for J=11/2

    D = np.zeros((len(m2), len(m1)), dtype=complex)
    for i, m in enumerate(m1):
        for j, mp in enumerate(m2):
            cg = clebsch_gordan(J1, m, 1, q, J2, mp).evalf()
            print(cg)
            if cg != 0:
                D[j, i] = reduced_element * cg
    return D

def dipole_J15to11(q, reduced_element=1.0):
    Ji = Rational(15, 2)  # initial J
    Jf = Rational(11, 2)  # final J
    m_i = [Rational(i, 2) for i in range(-15, 16, 2)]
    m_f = [Rational(i, 2) for i in range(-11, 12, 2)]

    D = np.zeros((len(m_f), len(m_i)), dtype=complex)

    for i, m in enumerate(m_i):
        for j, mp in enumerate(m_f):
            cg = clebsch_gordan(Ji, m, 1, q, Jf, mp).evalf()
            if cg != 0:
                D[j, i] = reduced_element * cg
    return D

# Build dipole matrices
D0 = dipole_J15to11(0)
Dp = dipole_J15to11(1)
Dm = dipole_J15to11(-1)

# Effective Raman operators in J=15/2 manifold
R00 = D0.conj().T @ D0
Rpp = Dp.conj().T @ Dp
Rmm = Dm.conj().T @ Dm

# --- Load your actual eigenstates ---
# Assume ionObj.eigenvectors is a 16x16 complex array with eigenvectors as columns
# Replace this line with your actual object reference:
CF_eigenstates = ionObj.eigenvectors.real   # shape (16, 16)

# Pick initial and final states (indices into eigenvectors)
initial_index = 0   # ground state (column 0)
initial_state = CF_eigenstates[0, :]

# Compute Raman intensities to all final CF states using R00 (q=0 component)
intensities = []
labels = []

for f in range(CF_eigenstates.shape[1]):
    final_state = CF_eigenstates[f, :]
    amp = np.vdot(final_state, R00 @ initial_state)  # ⟨ψ_f|R|ψ_i⟩
    I = np.abs(amp)**2
    intensities.append(I)
    labels.append(f"State {f}")

# --- Plot ---
plt.figure(figsize=(8,4))
plt.bar(labels, intensities)
plt.ylabel("Relative Raman Intensity")
plt.title("Raman Intensities via J=11/2 Intermediate (q=0)")
plt.xticks(rotation=90)
plt.tight_layout()
plt.show()


import numpy as np
from sympy.physics.wigner import clebsch_gordan
from sympy import Rational
import matplotlib.pyplot as plt

# --- Dipole construction for Ji=15/2 -> Jf=13/2 ---
def dipole_J15to13(q, reduced_element=1.0):
    Ji = Rational(15, 2)
    Jf = Rational(13, 2)  # using 13/2 as effective J
    m_i = [Rational(i, 2) for i in range(-15, 16, 2)]
    m_f = [Rational(i, 2) for i in range(-13, 14, 2)]

    D = np.zeros((len(m_f), len(m_i)), dtype=complex)

    for i, m in enumerate(m_i):
        for j, mp in enumerate(m_f):
            cg = clebsch_gordan(Ji, m, 1, q, Jf, mp).evalf()
            print(cg)
            if cg != 0:
                D[j, i] = reduced_element * cg
    return D

# Build dipole matrices with mixing
alpha = 0.1  # mixing factor, adjust as needed
D0 = alpha * dipole_J15to13(0)
Dp = alpha * dipole_J15to13(1)
Dm = alpha * dipole_J15to13(-1)

# Effective Raman operator in J=15/2 space
R00 = D0.conj().T @ D0

# Load your CF eigenstates
CF_eigenstates = ionObj.eigenvectors   # 16x16 numpy array

# Pick ground state as initial
initial_index = 0
initial_state = CF_eigenstates[:, initial_index]

# Compute intensities for all final CF states
intensities = []
labels = []

for f in range(CF_eigenstates.shape[1]):
    final_state = CF_eigenstates[:, f]
    amp = np.vdot(final_state, R00 @ initial_state)
    I = np.abs(amp)**2
    intensities.append(I)
    labels.append(f"State {f}")

# Plot
plt.figure(figsize=(8,4))
plt.bar(labels, intensities)
plt.ylabel("Relative Raman Intensity")
plt.title(f"Raman Intensities with J-mixing factor α={alpha}")
plt.xticks(rotation=90)
plt.tight_layout()
plt.show()

## do over field
alpha = 0.1
N_ground = 16
N_intermediate = 12
fields = np.linspace(0, 10, 1000)

# kB = 0.695  # cm^-1/K (adjust units to match your E)
T = temperature  # make sure this is set in your code

# Effective dipole operator
D_eff = alpha * np.ones((N_intermediate, N_ground), dtype=complex) / np.sqrt(N_ground)
R_eff = D_eff.conj().T @ D_eff

n_levels = N_ground
all_intensities = np.zeros((len(fields), n_levels, n_levels))
all_energies = np.zeros((len(fields), n_levels, n_levels))

for idx_H, H in enumerate(fields):
    ionObj = diagonalizeC(ionObj, ion, Jz, H, temperature)
    E = np.array(ionObj.eigenvalues)
    V = np.array(ionObj.eigenvectors)

    # Boltzmann weights for initial states
    boltz = np.exp(-E / (kB * T))
    boltz /= boltz.sum()

    for i in range(n_levels):
        psi_i = V[i]
        w_i = boltz[i]
        for f in range(i+1, n_levels):
            psi_f = V[f]
            amp = np.vdot(psi_f, R_eff @ psi_i)
            I = np.abs(amp)**2 * w_i  # apply Boltzmann weight
            all_intensities[idx_H, i, f] = I
            all_energies[idx_H, i, f] = (E[f] - E[i])*meVTocCmInv

# === Plot with opacity ===
cmap = mcolors.ListedColormap(['cyan'])

# Background contour plot
fig, ax = plt.subplots(figsize=(8,6))
plt.contourf(ramanField, ramanWavenums, ramanData, 100, cmap='Oranges')
plt.xlim(0,14)
plt.ylim(0,120)
plt.colorbar()
plt.xlabel('Field (T)')
plt.ylabel('Energy (cm$^{-1}$)')

# plt.title('CsErSe2 H||C with overlayed calc lines\n'
#           f'B20: {B20} B40: {B40} B43: {B43}\n'
#           f'B60: {B60} B63: {B63} B66: {B66}\n'
#           f'temp = {temperature}')

# Build calc_field, arrB, ampB arrays from your Raman data
calc_field = fields

arrB = []
ampB = []

for i in range(n_levels):
    for f in range(n_levels):
        # if i == f:
            # continue
        arrB.append(all_energies[:, i, f])
        ampB.append(all_intensities[:, i, f] / all_intensities.max())  # normalize 0..1

# Overlay calculated lines
for idx in range(len(arrB)):
    # plt.plot(calc_field, arrB[idx]*meVToCm, 'c', alpha=1, linewidth=0.7)
    alphas = ampB[idx]
    points = np.array([calc_field, arrB[idx]]).T.reshape(-1, 1, 2)
    segments = np.concatenate([points[:-1], points[1:]], axis=1)
    lc = LineCollection(segments, cmap=cmap, alpha=alphas, linewidth=0.7)
    ax.add_collection(lc)

ax.autoscale()
plt.show()


## go over all manifolds
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
import matplotlib.colors as mcolors

alpha = 0.1
N_ground = 16
fields = np.linspace(0, 18, 180)

# --- Define intermediate manifold dimensions ---
J_manifolds = {
    13/2: 14,
    11/2: 12,
    9/2: 10,
    3/2: 4
}

# --- Build effective dipole matrices for each manifold ---
def make_D_eff(n_intermediate, n_ground, alpha):
    return alpha * np.ones((n_intermediate, n_ground), dtype=complex) / np.sqrt(n_ground)

# Sum Raman operators from each manifold
R_eff_total = np.zeros((N_ground, N_ground), dtype=complex)
for J, dim in J_manifolds.items():
    D_eff = make_D_eff(dim, N_ground, alpha)
    R_eff_total += D_eff.conj().T @ D_eff

# === Pre-allocate arrays ===
n_levels = N_ground
all_intensities = np.zeros((len(fields), n_levels, n_levels))
all_energies = np.zeros((len(fields), n_levels, n_levels))

for idx_H, H in enumerate(fields):
    ionObj = diagonalizeC(ionObj, ion, Jz, H, temperature)
    E = np.array(ionObj.eigenvalues)
    V = np.array(ionObj.eigenvectors)

    # Boltzmann weighting
    T = temperature
    kB  = 8.617e-2 
    boltz = np.exp(-E / (kB * T))
    boltz /= boltz.sum()

    # Compute all i→f transitions
    for i in range(n_levels):
        psi_i = V[i]
        w_i = boltz[i]
        for f in range(i+1, n_levels):
            psi_f = V[f]
            amp = np.vdot(psi_f, R_eff_total @ psi_i)
            I = np.abs(amp)**2 * w_i
            all_intensities[idx_H, i, f] = I
            all_energies[idx_H, i, f] = (E[f] - E[i])*meVToCm

# === Build arrays for plotting ===
calc_field = fields
arrB = []
ampB = []
for i in range(n_levels):
    for f in range(n_levels):
        if i == f:
            continue
        arrB.append(all_energies[:, i, f])
        ampB.append(all_intensities[:, i, f] / all_intensities.max())

# === Plot ===
cmap = mcolors.ListedColormap(['cyan'])
fig, ax = plt.subplots(figsize=(8,6))

# Background experimental map (replace with your data arrays)
plt.contourf(ramanField, ramanWavenums, ramanData,100, cmap = 'Oranges')
plt.xlim(0,14)
plt.ylim(0,120)
plt.colorbar()
plt.xlabel('Field (T)')
plt.ylabel('Energy (cm$^{-1}$)')
plt.title('CsErSe2 H||C with multiple intermediate manifolds')

# Overlay calculated lines
for idx in range(len(arrB)):
    alphas = ampB[idx]
    points = np.array([calc_field, arrB[idx]]).T.reshape(-1, 1, 2)
    segments = np.concatenate([points[:-1], points[1:]], axis=1)
    lc = LineCollection(segments, cmap=cmap, alpha=alphas, linewidth=0.7)
    ax.add_collection(lc)

ax.autoscale()
plt.show()

##
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors

sigma = .5  # Gaussian width in cm^-1
nF = 200     # field resolution
nE = 400     # energy resolution

# Field and energy grid
F_grid = np.linspace(fields.min(), fields.max(), nF)
E_grid = np.linspace(0, all_energies.max(), nE)
F_mesh, E_mesh = np.meshgrid(F_grid, E_grid, indexing='xy')

# Initialize calculated intensity map
calc_map = np.zeros_like(F_mesh)

# Loop through all transitions and add Gaussian contributions
for i in range(n_levels):
    for f in range(n_levels):
        if i == f:
            continue
        for idx_H, H in enumerate(fields):
            E_trans = all_energies[idx_H, i, f]
            I = all_intensities[idx_H, i, f]

            # Find nearest grid indices for field
            F_idx = np.argmin(np.abs(F_grid - H))

            # Add a Gaussian in the energy dimension
            gauss = np.exp(-0.5 * ((E_grid - E_trans) / sigma)**2)
            gauss /= (sigma * np.sqrt(2*np.pi))  # normalize
            calc_map[:, F_idx] += I * gauss

# === Plot ===
fig, ax = plt.subplots(figsize=(8,6))

# Background experimental data
# plt.contourf(ramanField, ramanWavenums, ramanData, 100, cmap='Oranges')

# Overlay Gaussian-broadened calculated map
extent = [F_grid.min(), F_grid.max(), E_grid.min(), E_grid.max()]
plt.imshow(calc_map, norm=mcolors.LogNorm(vmin=calc_map.min()+.000001, vmax=calc_map.max()), extent=extent, origin='lower', aspect='auto',
           cmap='coolwarm')

plt.xlim(0,14)
plt.ylim(0,120)
plt.colorbar(label='Calculated Raman Intensity')
plt.xlabel('Field (T)')
plt.ylabel('Energy (cm$^{-1}$)')
plt.title('CsErSe2 H||C: Gaussian-Broadened Raman Spectrum')
plt.tight_layout()
plt.show()

##
def gtensorzeeman(ionObj, field=0.1, Temp=0.1):
    '''Returns g tensor computed numerically from zeeman splitting'''
    Jx = cef.Operator.Jx(ionObj.J)
    Jy = cef.Operator.Jy(ionObj.J)
    Jz = cef.Operator.Jz(ionObj.J)

    #print(Jx)
    #print(Jy)
    muB = 5.7883818012e-2  # meV/T
    #mu0 = np.pi*4e-7       # T*m/A

    gg = np.zeros(3)
    #loop through x,y,z
    for i,Field in enumerate([[field,0,0], [0,field,0], [0,0,field]]):
        JdotB = muB*(Field[0]*Jx + Field[1]*Jy + Field[2]*Jz)

        # B) Diagonalize full Hamiltonian
        FieldHam = ionObj.H + JdotB.O
        diagonalH = LA.eigh(FieldHam)

        minE = np.amin(diagonalH[0])
        evals = diagonalH[0] - minE
        evecs = diagonalH[1].T

        DeltaZeeman = evals[1]-evals[0]
        print(DeltaZeeman)

        # Now find the expectation value of J
        JexpVals = np.zeros((len(evals),3))
        for ii, ev in enumerate(evecs):
            kev = cef.Ket(ev)
            JexpVals[ii] =[np.real(kev*kev.Jx()),
                        np.real(kev*kev.Jy()),
                        np.real(kev*kev.Jz())]
        k_B = 8.6173303e-2  # meV/K

        Zz = np.sum(np.exp(-evals/(k_B*Temp)))
        JexpVal = np.dot(np.exp(-evals/(k_B*Temp)),JexpVals)/Zz

        expectationJ = JexpVal[i]

        # calculate g values
        gg[i] = DeltaZeeman/(muB*field*expectationJ)
    
    return gg


## calc just the high energy lines
def zeemanSplitC(field, wavenum, B20, B40, B43, B60, B63, B66, Jz,temperature):    
    dEphonon = 49.3
    phononAmp = 0.499
    phononSig = 0.95
    fun = []
    amp = []
    energies = []
    Bparams =  {'B20': B20, 'B40':B40,'B43': B43, 'B60': B60, 'B63':B63,'B66':B66}
    ionObj = cef.CFLevels.Bdict(ion,Bparams)
    for b in field: 
        ionObj = diagonalizeC(ionObj, ion, Jz, b, temperature)
        dE =[eval for eval in ionObj.eigenvalues] # this is the spitting if everything is in the GS -> not necessarily true for finite temp
        tempAmp = [0 if i == 0 else transition(ionObj, 0,i,temperature)[0]+transition(ionObj, 0,i,temperature)[1] for i in range(len(dE))]
        numLines = len(dE)
        dE_temp = dE
        dE_temp = dE[10:numLines]
        tempAmp = tempAmp[10:numLines]
        # tempAmp = [ionObj.transitionIntensity(0,i, temperature) for i in range(len(dE))]
        # tempAmp[0] = 0
        # numLines = len(dE)
        for i in range(10,numLines): 
            # print(i)
            for j in range(i+1, numLines):
                # print(j)
                temp = dE[j]-dE[i]
                dE_temp.append(temp)
                # tempAmp.append(ionObj.transitionIntensity(i,j,temperature))
                tempAmp.append(transition(ionObj, i,j,temperature)[0]+transition(ionObj, i,j,temperature)[1])
        wid = 1
        centers = dE_temp
        tempfun = lorentzian(wavenum, 0, dEphonon, phononSig)
        for i in range(len(centers)):
            a = tempAmp[i]
            tempfun += lorentzian(wavenum, a, centers[i]*meVToCm, wid)
        fun.append(tempfun)
        amp.append(tempAmp)
        energies.append(centers)
    return fun, amp, energies

# now, generate lines for each, simulated spectrum for each
calc_field = np.arange(0,20, .05)
calc_wavenums = np.arange(0,500, .05)
temperature = 10
kBT = kB*temperature
Jperp = 0

simulated_spec_C, ampC, arrC = zeemanSplitC(calc_field, calc_wavenums, B20, B40, B43, B60, B63, B66, Jz, temperature)
simulated_spec_C  = np.array(simulated_spec_C)
ampC = np.array(ampC).T/np.array(ampC).max(axis=None)
arrC = np.array(arrC).T

# # calc_wavenums = np.linspace(0,400, 400)
# simulated_spec_IR_B, ampB, arrB = zeemanSplitAB(calc_field, calc_wavenums, B20, B40, B43, B60, B63, B66, Jperp, temperature)
# simulated_spec_IR_B  = np.array(simulated_spec_IR_B)
# ampB = np.array(ampB).T
# arrB = np.array(arrB).T

# spectroscopy_calc_fname
spec_calc_fname = 'spectroscopy_fgr_high_energy_2025Aug07.h5'
with h5py.File(spec_calc_fname, 'w') as hdf:
    hdf.create_dataset('linesC', data = arrC)
    # hdf.create_dataset('linesB', data = arrB)
    hdf.create_dataset('ampC', data = ampC)
    # hdf.create_dataset('ampB', data = ampB)
    # hdf.create_dataset('linesC_nomft', data = arrC_nomft)
    # hdf.create_dataset('linesB_nomft', data = arrB_nomft)
    # hdf.create_dataset('ampC_nomft', data = ampC_nomft)
    # hdf.create_dataset('ampB_nomft', data = ampB_nomft)
    hdf.create_dataset('calc_field', data = calc_field)
    hdf.create_dataset('calc_wavenums', data= calc_wavenums)
    hdf.create_dataset('simulated_C', data = simulated_spec_C)
    # hdf.create_dataset('simulated_IR_B', data = simulated_spec_IR_B)
    # hdf.create_dataset('simulated_IR_C', data = simulated_spec_IR_C)
    # hdf.create_dataset('ramanData', data=ramanData.values)
    # hdf.create_dataset('IR_dataB', data=dataB.values)
    # hdf.create_dataset('IR_dataC', data=dataC.values)
    # hdf.create_dataset('IR_B_wavenums', data = IRBwavenums)
    # hdf.create_dataset('IR_C_wavenums', data = IRCwavenums)
    # hdf.create_dataset('raman_wavenums', data = ramanWavenums)
    # hdf.create_dataset('IR_B_field', data = IRBfield)
    # hdf.create_dataset('IR_C_field', data = IRCfield)
    # hdf.create_dataset('raman_field', data = ramanField)
    hdf.attrs['B20'] = B20
    hdf.attrs['B40'] = B40
    hdf.attrs['B43'] = B43
    hdf.attrs['B60'] = B60
    hdf.attrs['B63'] = B63
    hdf.attrs['B66'] = B66
    hdf.attrs['Jz'] = Jz
    hdf.attrs['Jperp'] = Jperp
    hdf.attrs['notes'] = 'recalculated with fgr Jul 30 2025 '
