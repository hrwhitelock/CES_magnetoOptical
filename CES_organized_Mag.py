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

# define some constants
temperature = 2 # in K

muB = 5.7883818012e-2  # meV/T
mu0 = np.pi*4e-7       # T*m/A
kB  = 8.617e-2         # [meV/K]
meVToCm =meVTocCmInv= 8.066 
ion = 'Er3+'

kBT = kB*temperature
gJ = cef.LandeGFactor('Er3+')


Jperp = -0.2e-3*-2 #meV
Jz = -2.4e-3*-0.2 #meV
q= 6

# init objects
B20 = -0.03265325 # init = -0.03559)
B40 = -0.0003849 # fixed)
B43 = -0.01393 # fixed)
B60 =  3.035e-6# fixed)
B63 = -8.4011e-07 # init = -4.695e-06)
B66 =  3.3815e-05 # fixed)

g = cef.LandeGFactor(ion)
myBparams =  {'B20': B20, 'B40':B40,'B43': B43, 'B60': B60, 'B63':B63,'B66':B66}
MyErObj = cef.CFLevels.Bdict(ion,myBparams)


# neutron fit vals
B20 = -4.73e-2
B40 = -3.7037e-4
B43 = -1.44431e-2
B60 = 3.1605e-6
B63 = 6.5259e-6
B66 = 3.9314e-5

g = cef.LandeGFactor(ion)
AllenBparams =  {'B20': B20, 'B40':B40,'B43': B43, 'B60': B60, 'B63':B63,'B66':B66}
AllenErObj = cef.CFLevels.Bdict(ion,AllenBparams)

## first, calc c axis mag
f = np.linspace(0,10, 1000)
ffine = np.concatenate((np.linspace(0,1,100000), np.linspace(0,12,1000000)))
field = [[0,0,b] for b in ffine]
magMe_C = np.array([MyErObj.magnetization(ion, temperature, f) for f in field]).T


magAllen_C = np.array([AllenErObj.magnetization(ion, temperature, f) for f in field]).T
field = np.array(field).T

plt.plot(field[2], magMe_C[2])
plt.plot(field[2], magAllen_C[2])
plt.title('magnetization')
plt.xlabel('applied field')
plt.ylabel('mag')

## now to MFT correction
def MolecularFieldTheory(H, Hth, Mth, lamb):
    '''Use mean field theory to correct magnetization
    for an exchange strength lamb. H is the number of fields to calculate,
    Hth and Mth are the theoretical single-ion magnetization curves to correct.'''
    n = 10
    # colors = plt.cm.gnuplot(np.linspace(0,1,n))
    # plt.figure()
    # plt.plot(Hth, Mth, 'r', label = 'orginal data')
    newM = np.interp(H, Hth, Mth)
    for i in range(n):
        newH = H - 6*lamb*newM/muB/(gJ)**2
        newM = np.interp(newH,Hth,Mth)
        # for testing
    #     plt.plot(newH, newM, label = str(i), color = colors[i])
    # plt.legend()
    return newM


mme_C = magMe_C[2]
mAllen_C = magAllen_C[2]
f = np.linspace(0,10,1000)
ffine = np.linspace(0,10,len(mme_C))

mft_mZ_me_C = MolecularFieldTheory(f, ffine, mme_C, Jz) # allens def of J is lamb unfortunately
mft_mZ_Allen_C = MolecularFieldTheory(f, ffine, mAllen_C, Jz)

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
plt.plot(f, -1*mft_mZ_me_C, '--m', label = 'Raman fit B params')
plt.plot(f, -1*mft_mZ_Allen_C, '--c',label = 'Neutron fit B params')
plt.plot(ffine, -1*mme_C, '--g', label = 'Raman fit Bparams, no MFT')
plt.plot(ffine, -1*mAllen_C, '--p', label = 'Neutron fit Bparams, no MFT')
plt.xlim(0,6)
plt.ylim(0,10)
plt.legend()
plt.title('C magnetization')
plt.xlabel('Field (T)')
plt.ylabel('Magnetization (uB/Er)')