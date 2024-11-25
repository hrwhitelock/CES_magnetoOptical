# reorganization of magnetic calcs 
import numpy as np
import matplotlib.pyplot as plt
from ipywidgets import interact
import PyCrystalField as cef
import scipy
from scipy.misc import derivative
import lmfit
import pandas as pd
from numba import njit
import lmfit

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


Jperp = -0.9e-3 #meV
Jz = -0.48e-3 #meV

JperpAllen = -0.2 # meV
JzAllen = -2.4e-3 # meV

q= 6

# init objects
B20 = -0.03265325 # init = -0.03559)
B40 = -0.0003849 # fixed)
B43 = -0.01393 # fixed)
B60 =  3.03e-6# fixed)
B63 = -8.4011e-07 # init = -4.695e-06)
B66 =  3.3815e-05 # fixed)

g = cef.LandeGFactor(ion)
myBparams =  {'B20': B20, 'B40':B40,'B43': B43, 'B60': B60, 'B63':B63,'B66':B66}
MyErObj = cef.CFLevels.Bdict(ion,myBparams)


# neutron fit vals
B20 = -3.559e-2
B40 = -3.849e-4
B43 = -1.393e-2
B60 = 3.154e-6
B63 = -4.695e-6
B66 = 3.3815e-5

g = cef.LandeGFactor(ion)
AllenBparams =  {'B20': B20, 'B40':B40,'B43': B43, 'B60': B60, 'B63':B63,'B66':B66}
AllenErObj = cef.CFLevels.Bdict(ion,AllenBparams)

## first, calc c axis mag
magF = np.linspace(0,10,10000)
MFTField = np.linspace(0,8,10000)
# magF = np.concatenate((np.linspace(0,1,100000), np.linspace(1,4.8,1000), np.linspace(4.8,5.8, 10000), np.linspace(5.8,12, 1000)))
field = [[0,0,b] for b in magF]
myCaxisMagnetization = np.array([MyErObj.magnetization(ion, temperature, f) for f in field]).T
myCaxisMagnetization = myCaxisMagnetization[2]
myCaxisMagnetization = [m*-1 for m in myCaxisMagnetization] # make the magnetization the correct sign

allenCaxisMagnetization = np.array([AllenErObj.magnetization(ion, temperature, f) for f in field]).T
allenCaxisMagnetization = allenCaxisMagnetization[2]
allenCaxisMagnetization = [m*-1 for m in allenCaxisMagnetization]

## now to MFT correction
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


myMFTCaxis2K = MolecularFieldTheory(MFTField, magF, myCaxisMagnetization, Jz) # allens def of J is lamb unfortunately
allenMFTCaxis = MolecularFieldTheory(MFTField, magF, allenCaxisMagnetization, JzAllen)

# now we load the data
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

# now plot data with curves
plt.figure()
plt.plot(CESMHdata[6]/1e4,CESMHdata[7],'b.', label='data ($H \\parallel c$)')
plt.plot(MFTField, myMFTCaxis, '-', label = 'Raman fit B params')
plt.plot(MFTField, allenMFTCaxis, '-',label = 'Neutron fit B params')
plt.plot(magF, myCaxisMagnetization2K, '--', label = 'Raman fit Bparams, no MFT')
plt.plot(magF, allenCaxisMagnetization, '--', label = 'Neutron fit Bparams, no MFT')
# plt.xlim(0,6)
# plt.ylim(0,10)
plt.legend()
plt.title('C magnetization')
plt.xlabel('Field (T)')
plt.ylabel('Magnetization (uB/Er)')

###################################################################################
# lets load the MPMS data
fname20K = '/Users/hopeless/Desktop/LeeLab/data/CsErSe2_data/MPMS/CsErSe2/CsErSe2_HParC_MvsH_20K.txt'
fname6K = '/Users/hopeless/Desktop/LeeLab/data/CsErSe2_data/MPMS/CsErSe2/CsErSe2_HParC_MvsH_6K.txt'
fname2K = '/Users/hopeless/Desktop/LeeLab/data/CsErSe2_data/MPMS/CsErSe2/CsErSe2_HParC_MvsH_2K.txt'

Mdata20K = np.genfromtxt(fname20K, delimiter=',',  unpack=True, skip_header=1)
Mdata6K = np.genfromtxt(fname6K, delimiter=',',  unpack=True, skip_header=1)
Mdata2K = np.genfromtxt(fname2K, delimiter=',',  unpack=True, skip_header=1)


# let's make some temperature dependent data
magF_mpms = np.linspace(-8,8,1000)
field = [[0,0,b] for b in magF_mpms]
temperature = 2
magnetization2K = np.array([MyErObj.magnetization(ion, temperature, f) for f in field]).T
magnetization2K = magnetization2K[2]
magnetization2K = [m*-1 for m in magnetization2K] # make the magnetization the correct sign

temperature = 6
magnetization6K = np.array([MyErObj.magnetization(ion, temperature, f) for f in field]).T
magnetization6K = magnetization6K[2]
magnetization6K = [m*-1 for m in magnetization6K] # make the magnetization the correct sign

temperature = 20
magnetization20K = np.array([MyErObj.magnetization(ion, temperature, f) for f in field]).T
magnetization20K = magnetization20K[2]
magnetization20K = [m*-1 for m in magnetization20K] # make the magnetization the correct sign

MFTField_mpms = np.linspace(-8,8,10000)
MFT2K = MolecularFieldTheory(MFTField_mpms, magF_mpms, magnetization2K, Jz)
MFT6K = MolecularFieldTheory(MFTField_mpms, magF_mpms, magnetization6K, Jz)
MFT20K = MolecularFieldTheory(MFTField_mpms, magF_mpms, magnetization20K, Jz)

plt.figure()
plt.plot(Mdata20K[0], Mdata20K[1], 'o', label = '20K MPMS data')
plt.plot(Mdata6K[0], Mdata6K[1], 'o', label = '6K MPMS data')
plt.plot(Mdata2K[0], Mdata2K[1], 'o', label = '2K MPMS data')
plt.plot(CESMHdata[6]/1e4,CESMHdata[7],'b.', label='from Allens paper')
plt.plot(MFTField_mpms, MFT2K, '-', label = 'MFT 2K')
plt.plot(MFTField_mpms, MFT6K, '-', label = 'MFT 6K')
plt.plot(MFTField_mpms, MFT20K, '-', label = 'MFT 20K')
plt.plot(magF_mpms, magnetization2K, '--', label = '2K, no MFT')
plt.plot(magF_mpms, magnetization6K, '--', label = '6K, no MFT')
plt.plot(magF_mpms, magnetization20K, '--', label = '20K, no MFT')
# plt.xlim(0,6)
# plt.ylim(0,10)
plt.legend()
plt.title('C magnetization')
plt.xlabel('Field (T)')
plt.ylabel('Magnetization (uB/Er)')
plt.show()


###################################################################################
# Now let's nail down B60
# init objects
B20 = -0.03265325 # init = -0.03559)
B40 = -0.0003849 # fixed)
B43 = -0.01393 # fixed)
B60 =  3.03e-6# fixed)
B63 = -8.4011e-07 # init = -4.695e-06)
B66 =  3.3815e-05 # fixed)

g = cef.LandeGFactor(ion)
myBparams =  {'B20': B20, 'B40':B40,'B43': B43, 'B60': B60, 'B63':B63,'B66':B66}
MyErObj = cef.CFLevels.Bdict(ion,myBparams)

temperature = 0.1 # I know this is where we'll see a defined peak

magMe = [m*-1 for m in np.array([MyErObj.magnetization(ion, temperature, f) for f in field]).T[2]]
m = MolecularFieldTheory(MFTField, magF, magMe, Jz)
dmdh = np.gradient(m, MFTField)

plt.figure()
plt.plot(MFTField, dmdh)
plt.vlines(x = 5.4, ymin=0, ymax = 300)
plt.xlabel('Field (T)')
plt.ylabel('dM/dH')


###################################################################################
# okay, so now we plot some temp dependance

temps = [0.025, 0.045,0.1, 0.171, .25, .35, .45, .543, .827, 2, 6, 20]
n = len(temps)
colors = plt.cm.jet(np.linspace(0,1,n))

field = [[0,0,b] for b in magF]
tempMag =[]

for temperature in temps: 
    magMe = np.array([MyErObj.magnetization(ion, temperature, f) for f in field]).T
    temp = MolecularFieldTheory(MFTField, magF, -1*magMe[2], Jz)
    tempMag.append(temp) # what the actual fuck is this naming holy shit

plt.figure()

for i in range(0, len(tempMag)): 
    plt.plot(MFTField, tempMag[i], label = str(temps[i])+'K', color = colors[i])

plt.legend()
# plt.ylim(0,3)
plt.xlim(0,9)
plt.title('C axis magnetization MFT \n calculated from Raman fit B params \n test B60 = '+str(B60))
plt.xlabel('Field (T)')
plt.ylabel('Magnetization')

###################################################################################
# now temp dependent dmdh
n = len(temps)
colors = plt.cm.jet(np.linspace(0,1,n))
i = 0

# calculate gradient
dmdH = []
for mag in tempMag: 
    temp = np.gradient(mag)
    dmdH.append(temp)

for i in range(len(temps)): 
    plt.plot(MFTField, dmdH[i], label = str(temps[i])+'K', color = colors[i])

plt.legend()
plt.xlabel(' applied field')
plt.ylabel('dm/dH')
# plt.vlines(x = 5.4, ymin=0, ymax = 12)
# plt.xlim(0,7)
# plt.ylim(-1,12)
plt.title('C axis dM/dH \n calculated from Raman fit B params')

# I really think I should fit B60, J, because they aren't seperate..... I move B60 a different way depending on 
# the sign of B60..................
# this weekend I want to write a custom fit function to fit params simultaneously