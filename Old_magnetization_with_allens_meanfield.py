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
Jz = 0.48e-3 #meV

JperpAllen = -0.2 # meV
JzAllen = -2.4e-3 # meV

q= 6

# init objects
B20 = -0.03265325 # init = -0.03559)
B40 = -0.0003849 # fixed)
B43 = -0.01393 # fixed)
B60 =  3.079e-6# fixed)
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
magF = np.linspace(0,10,100)
MFTField = np.linspace(0,8,10000)
temperature = 2
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
    n = 11 # try odd for reasons
    newM = np.interp(H, Hth, Mth)
    for i in range(n):
        newH = H + 6*lamb*newM/muB/(gJ)**2
        newM = np.interp(newH,Hth,Mth)
    return newM


myMFTCaxis = MolecularFieldTheory(MFTField, magF, myCaxisMagnetization, Jz) # allens def of J is lamb unfortunately
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
plt.plot(magF, myCaxisMagnetization, '--', label = 'Raman fit Bparams, no MFT')
plt.plot(magF, allenCaxisMagnetization, '--', label = 'Neutron fit Bparams, no MFT')
# plt.xlim(0,6)
# plt.ylim(0,10)
plt.legend()
plt.title('C magnetization')
plt.xlabel('Field (T)')
plt.ylabel('Magnetization (uB/Er)')


# Save data to HDF5
with h5py.File('magnetization_data_calculation_and_allen_c_axis.h5', 'w') as f:
    # Save the data arrays
    f.create_dataset('magF', data=magF)
    f.create_dataset('MFTField', data=MFTField)
    f.create_dataset('temperature', data=temperature)
    f.create_dataset('field', data=field)
    f.create_dataset('myCaxisMagnetization', data=myCaxisMagnetization)
    f.create_dataset('allenCaxisMagnetization', data=allenCaxisMagnetization)
    f.create_dataset('myMFTCaxis', data=myMFTCaxis)
    f.create_dataset('allenMFTCaxis', data=allenMFTCaxis)
    f.create_dataset('CESMHdata', data=CESMHdata)
    f.create_dataset('CESMTdata', data=CESMTdata)

print("Data saved to 'magnetization_data_caluclation_and_allen_c_axis.h5'.")

###################################################################################
###################################################################################
###################################################################################
###################################################################################

# lets load the MPMS data
fname20K = '/Users/hopeless/Desktop/LeeLab/data/CsErSe2_data/MPMS/CsErSe2/CsErSe2_HParC_MvsH_20K.txt'
fname6K = '/Users/hopeless/Desktop/LeeLab/data/CsErSe2_data/MPMS/CsErSe2/CsErSe2_HParC_MvsH_6K.txt'
fname2K = '/Users/hopeless/Desktop/LeeLab/data/CsErSe2_data/MPMS/CsErSe2/CsErSe2_HParC_MvsH_2K.txt'

Mdata20K = np.genfromtxt(fname20K, delimiter=',',  unpack=True, skip_header=1)
Mdata6K = np.genfromtxt(fname6K, delimiter=',',  unpack=True, skip_header=1)
Mdata2K = np.genfromtxt(fname2K, delimiter=',',  unpack=True, skip_header=1)

testJz = -2.4e-4

# let's make some temperature dependent data
magF_mpms = np.linspace(-8,8,100)
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
plt.grid(True)
plt.plot(Mdata20K[0], Mdata20K[1]/1.35, 'o', label = '20K MPMS data')
plt.plot(Mdata6K[0], Mdata6K[1]/1.35, 'o', label = '6K MPMS data')
plt.plot(Mdata2K[0], Mdata2K[1]/1.35, 'o', label = '2K MPMS data')
plt.plot(CESMHdata[6]/1e4,CESMHdata[7],'b.', label='from Allens paper')
plt.plot(MFTField_mpms, MFT2K, '-', label = 'MFT 2K')
plt.plot(MFTField_mpms, MFT6K, '-', label = 'MFT 6K')
plt.plot(MFTField_mpms, MFT20K, '-', label = 'MFT 20K')
plt.plot(magF_mpms, magnetization2K, '--', label = '2K, no MFT')
plt.plot(magF_mpms, magnetization6K, '--', label = '6K, no MFT')
plt.plot(magF_mpms, magnetization20K, '--', label = '20K, no MFT')
plt.xlim(0,7)
plt.ylim(0,8)
plt.legend()
plt.title('C magnetization \n B20 ='+str(MyErObj.B[0])+ 'B40 = '+str(MyErObj.B[1])+'B43 =' +str(MyErObj.B[2])+'B60 =' +str(MyErObj.B[3]) + 'B63 = '+str(MyErObj.B[4])+ 'B66 = '+str(MyErObj.B[5]) +'\n Jperp = '+ str(Jperp)+'\n JperpAllen = '+ str(JperpAllen)+'\n Jz = '+ str(Jz)+'\n JzAllen = '+ str(JzAllen))
# plt.title('c-axis magnetization with test Jz = ' + str(testJz))

plt.xlabel('Field (T)')
plt.ylabel('Magnetization (uB/Er)')
plt.show()


# Save all the data to an HDF5 file
with h5py.File('magnetization_data_caluclation_mpms_data_c_axis.h5', 'w') as f:
    f.create_dataset('magF_mpms', data=magF_mpms)
    f.create_dataset('MFTField_mpms', data=MFTField_mpms)
    f.create_dataset('magnetization2K', data=magnetization2K)
    f.create_dataset('magnetization6K', data=magnetization6K)
    f.create_dataset('magnetization20K', data=magnetization20K)
    f.create_dataset('MFT2K', data=MFT2K)
    f.create_dataset('MFT6K', data=MFT6K)
    f.create_dataset('MFT20K', data=MFT20K)
    f.create_dataset('Mdata2K', data=Mdata2K)
    f.create_dataset('Mdata6K', data=Mdata6K)
    f.create_dataset('Mdata20K', data=Mdata20K)
    f.create_dataset('CESMHdata', data=CESMHdata)
    hdf.attrs['B20'] = B20
    hdf.attrs['B40'] = B40
    hdf.attrs['B43'] = B43
    hdf.attrs['B60'] = B60
    hdf.attrs['B63'] = B63
    hdf.attrs['B66'] = B66
    hdf.attrs['Jperp'] = Jperp
    hdf.attrs['Jz'] = Jz
    hdf.attrs['JperpAllen'] = JperpAllen
    hdf.attrs['JzAllen'] = JzAllen

print("Data saved to 'magnetization_data.h5'.")

###################################################################################
###################################################################################
###################################################################################
###################################################################################
###################################################################################
# lets show a couple values of Jz

fname2K = '/Users/hopeless/Desktop/LeeLab/data/CsErSe2_data/MPMS/CsErSe2/CsErSe2_HParC_MvsH_2K.txt'

Mdata2K = np.genfromtxt(fname2K, delimiter=',',  unpack=True, skip_header=1)

testJz = [.48e-3, -.48e-3, 1e-3, -1e-3, -2e-3, 2e-3]

# let's make some temperature dependent data
magF_mpms = np.linspace(0,8,100)
field = [[0,0,b] for b in magF_mpms]
temperature = 2
magnetization2K = np.array([MyErObj.magnetization(ion, temperature, f) for f in field]).T
magnetization2K = magnetization2K[2]
magnetization2K = [m*-1 for m in magnetization2K] # make the magnetization the correct sign

MFTField_mpms = np.linspace(0,8,10000)
plt.figure()
mftArr = []
for j in testJz: 
    mft = MolecularFieldTheory(MFTField_mpms, magF_mpms, magnetization2K, j)
    plt.plot(MFTField_mpms, mft, '--', label = 'Jz = ' + str(j))
    mftArr.append(mft)

plt.grid(True)
plt.plot(magF_mpms, magnetization2K, label = 'no mft calculation')
plt.plot(Mdata2K[0], Mdata2K[1]/1.35, 'o', label = '2K MPMS data')
plt.plot(CESMHdata[6]/1e4,CESMHdata[7],'b.', label='from Allens paper')

plt.xlim(0,7)
plt.ylim(0,8)
plt.legend()
plt.title('C magnetization \n B20 ='+str(MyErObj.B[0])+ 'B40 = '+str(MyErObj.B[1])+'B43 =' +str(MyErObj.B[2])+'B60 =' +str(MyErObj.B[3]) + 'B63 = '+str(MyErObj.B[4])+ 'B66 = '+str(MyErObj.B[5]) +'\n Jperp = '+ str(Jperp)+'\n JperpAllen = '+ str(JperpAllen)+'\n Jz = '+ str(Jz)+'\n JzAllen = '+ str(JzAllen))
# plt.title('c-axis magnetization with test Jz = ' + str(testJz))

plt.xlabel('Field (T)')
plt.ylabel('Magnetization (uB/Er)')
plt.show()
with h5py.File('checking_J_B60_3079.h5', 'w') as f:
    f.create_dataset('testJz', data=testJz)
    f.create_dataset('mftArr', data=mftArr)
    f.create_dataset('MFTField', data=MFTField_mpms)
    f.create_dataset('magF_mpms', data=magF_mpms)
    f.create_dataset('magnetization2K', data=magnetization2K)
    f.create_dataset('MData2K', data = Mdata2K)
    f.attrs['B60'] = MyErObj.B[3]




###################################################################################
###################################################################################
###################################################################################
###################################################################################
###################################################################################
# Now let's nail down B60
# init objects
B20 = -0.03265325 # init = -0.03559)
B40 = -0.0003849 # fixed)
B43 = -0.01393 # fixed)
B60 =  3.079e-6# fixed)
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
plt.grid(True)
plt.plot(MFTField, dmdh)
plt.vlines(x = 5.4, ymin=0, ymax = 300)
plt.xlabel('Field (T)')
plt.ylabel('dM/dH')

###################################################################################
###################################################################################
###################################################################################
###################################################################################
###################################################################################

# okay, so now we plot some temp dependance

temps = [2,2.5,3,3.5,4,5, 6]#[0.005, 0.01,0.15, 0.02,  0.025, 0.05,0.075, 0.1,.125, 0.15, 0.171, .25, .35, .45, .543, .827, 1, 2, 6, 20]
n = len(temps)
colors = plt.cm.jet(np.linspace(0,1,n))


field = [[0,0,b] for b in magF]
tempMag =[]

for temperature in temps: 
    magMe = np.array([MyErObj.magnetization(ion, temperature, f) for f in field]).T
    temp = MolecularFieldTheory(MFTField, magF, -1*magMe[2], Jz)
    tempMag.append(temp) # what the actual fuck is this naming holy shit

# Save data to an HDF5 file
with h5py.File('M_vs_H_temperature_dependence_myParams_allensCode.h5', 'w') as hdf:
    hdf.create_dataset('temps', data=temps)
    hdf.create_dataset('MFTField', data=MFTField)
    hdf.create_dataset('tempMag', data=np.array(tempMag))

plt.figure()
plt.grid(True)

for i in range(0, len(tempMag)): 
    plt.plot(MFTField, tempMag[i], label = str(temps[i])+'K', color = colors[i])

plt.legend()
# plt.ylim(0,3)
plt.xlim(0,9)
# plt.title('C axis magnetization MFT \n calculated from Raman fit B params \n test B60 = '+str(B60))
plt.title('C-axi magnetization \n high temp to see when mag curvature turns off')
plt.xlabel('H [T]')
plt.ylabel('M \mu_B/Er')

# do the same for Allen's params

temps = [0.005, 0.01,0.15, 0.02,  0.025, 0.05,0.075, 0.1,.125, 0.15, 0.171, .25, .35, .45, .543, .827, 1, 2, 6, 20]
n = len(temps)
colors = plt.cm.jet(np.linspace(0,1,n))


field = [[0,0,b] for b in magF]
tempMagAllen =[]

for temperature in temps: 
    magAllen = np.array([AllenErObj.magnetization(ion, temperature, f) for f in field]).T
    magAllen = magAllen[2]
    magAllen = [m*-1 for m in magAllen]
    temp = MolecularFieldTheory(MFTField, magF, magAllen, JzAllen)
    tempMagAllen.append(temp) # what the actual fuck is this naming holy shit

# Save data to an HDF5 file
with h5py.File('M_vs_H_temperature_dependence_allens_code_allens_params.h5', 'w') as hdf:
    hdf.create_dataset('temps', data=temps)
    hdf.create_dataset('MFTField', data=MFTField)
    hdf.create_dataset('tempMagAllen', data=np.array(tempMagAllen))

plt.figure()
plt.grid(True)

for i in range(0, len(tempMagAllen)): 
    plt.plot(MFTField, tempMagAllen[i], label = str(temps[i])+'K', color = colors[i])

plt.legend()
# plt.ylim(0,3)
plt.xlim(0,9)
plt.title('C axis magnetization MFT \n calculated from Raman fit B params \n test B60 = '+str(B60))
plt.xlabel('Field (T)')
plt.ylabel('Magnetization')

###################################################################################
###################################################################################
###################################################################################
###################################################################################
###################################################################################
# now temp dependent dmdh
n = len(temps)
colors = plt.cm.jet(np.linspace(0,1,n))
i = 0

# calculate gradient
dmdH = []
for mag in tempMag: 
    temp = np.gradient(mag, MFTField)
    dmdH.append(temp)

plt.figure()
plt.grid(True)
for i in range(len(temps)): 
    plt.plot(MFTField, dmdH[i], label = str(temps[i])+'K', color = colors[i])

plt.legend()
plt.xlabel('Field (T)')
plt.ylabel('dm/dH')
plt.yscale('log')
# plt.vlines(x = 5.4, ymin=0, ymax = 12)
# plt.xlim(0,7)
# plt.ylim(-1,12)
plt.title('C axis dM/dH \n calculated from Raman fit B params')

with h5py.File('dMdH_temperature_dependence.h5', 'w') as hdf:
    hdf.create_dataset('temps', data=temps)
    hdf.create_dataset('MFTField', data=MFTField)
    hdf.create_dataset('dmdH', data=np.array(dmdH))

###################################################################################
###################################################################################
###################################################################################
###################################################################################
###################################################################################
# lets do susceptibility
def susceptibility(ionObj, fieldVal, temps):
    chi = []
    for temp in temps: 
        # f = np.arange(fieldVal-.5*fieldVal, .5*fieldVal+fieldVal, .05*fieldVal) 
        f = np.arange(fieldVal-0.5, fieldVal+0.5, .0012)
        field = [[0,0,b] for b in f]
        mag= np.array([ionObj.magnetization(ion, temp, f) for f in field]).T
        m = MolecularFieldTheory(f, f, -mag[2], Jz)
        m = np.array(m).T
        x = np.gradient(m, f) 
        valIdx = findIdx(f, fieldVal)
        chi.append(x[valIdx])
    return chi
def findIdx(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx

### grab C axis susceptibility
fnames = ['/Users/hopeless/Desktop/LeeLab/data/CsErSe2_data/scm1_chi_T_01T.txt', 
            '/Users/hopeless/Desktop/LeeLab/data/CsErSe2_data/scm1_chi_T_1T.txt', 
            '/Users/hopeless/Desktop/LeeLab/data/CsErSe2_data/scm1_chi_T_3T.txt']
labels = ['0.1T', '1T', '3T']
xArrs = []
yArrs =[]


for fname in fnames: 
    temp = np.genfromtxt(fname, delimiter=',',  unpack=True, skip_header=1)
    xArrs.append(temp[0])
    yArrs.append(temp[1])

fieldVal = 0.0
temps = np.concatenate((np.linspace(.01, 2, 50), np.linspace(2,20,50), np.linspace(20,300, 279)))
mysus0T = susceptibility(MyErObj, fieldVal, temps) 
neutronSus0T = susceptibility(AllenErObj, fieldVal, temps)

fieldVal = 0.1
temps = np.concatenate((np.linspace(.01, 2, 50), np.linspace(2,20,50), np.linspace(20,300, 279)))
mysus01T = susceptibility(MyErObj, fieldVal, temps) 
neutronSus01T = susceptibility(AllenErObj, fieldVal, temps)

fieldVal = 1.0
mysus1T = susceptibility(MyErObj, fieldVal, temps) 
neutronSus1T = susceptibility(AllenErObj, fieldVal, temps)

fieldVal = 3.0
mysus3T = susceptibility(MyErObj, fieldVal, temps) 
neutronSus3T = susceptibility(AllenErObj, fieldVal, temps)


field = [0,0,.1]
df = 0.0012 # same as other spacing
susPCF = MyErObj.susceptibility(ion, temps, field, df)
susPCF = np.array(susPCF).T[2]
len(susPCF)

myinv0T = [1/x for x in mysus0T]
neutroninv0T = [1/x for x in neutronSus0T]
myinv01T = [1/x for x in mysus01T]
neutroninv01T = [1/x for x in neutronSus01T]
myinv1T = [1/x for x in mysus1T]
neutroninv1T = [1/x for x in neutronSus1T]
myinv3T = [1/x for x in mysus3T]
neutroninv3T = [1/x for x in neutronSus3T]



susinvPCF = [-1/x for x in susPCF]
plt.figure()
plt.grid(True)
for i in range(len(labels)): 
    plt.plot(xArrs[i], yArrs[i]*1.35, label = labels[i])

plt.plot(CESMTdata[12], 1/CESMTdata[13]*SCF, label='c-axis data from Allens paper')
plt.plot(temps, myinv0T, '--', label = 'Raman B params MFT 0T')
plt.plot(temps, neutroninv0T, '-.', label = 'neutrons B params MFT 0T' )

plt.plot(temps, myinv01T, '--', label = 'Raman B params MFT 0.1T')
plt.plot(temps, neutroninv01T, '-.', label = 'neutrons B params MFT 0.1T' )

plt.plot(temps, myinv1T, '--', label = 'Raman B params MFT 1T')
plt.plot(temps, neutroninv1T, '-.', label = 'neutrons B params MFT 1T' )

plt.plot(temps, myinv3T, '--', label = 'Raman B params MFT 3T')
plt.plot(temps, neutroninv3T, '-.', label = 'neutrons B params MFT 3T' )

plt.plot(temps, susinvPCF, '--', label = 'Raman B params no MFT, 0.1T')
plt.title('calculated MFT susceptibility at 0.1T ')
plt.xlabel('temperature (K)')
plt.ylabel('1/chi')
plt.legend()
plt.xlim(0, 200)
plt.ylim(0,10)

# Save data to HDF5
with h5py.File('susceptibility_wide_temp_range.h5', 'w') as hdf:
    hdf.create_dataset('xArrs', data=np.array(xArrs, dtype=object), dtype=h5py.vlen_dtype(float))
    hdf.create_dataset('yArrs', data=np.array(yArrs, dtype=object), dtype=h5py.vlen_dtype(float))
    hdf.create_dataset('labels', data=np.array(labels, dtype='S'))
    hdf.create_dataset('temps', data=temps)
    hdf.create_dataset('myinv01T', data=myinv01T)
    hdf.create_dataset('neutroninv01T', data=neutroninv01T)
    hdf.create_dataset('myinv0T', data=myinv0T)
    hdf.create_dataset('neutroninv0T', data=neutroninv0T)
    hdf.create_dataset('myinv1T', data=myinv1T)
    hdf.create_dataset('neutroninv1T', data=neutroninv1T)
    hdf.create_dataset('myinv3T', data=myinv3T)
    hdf.create_dataset('neutroninv3T', data=neutroninv3T)
    hdf.create_dataset('susinvPCF', data=susinvPCF)
    hdf.create_dataset('CESMTdata', data=CESMTdata)

###################################################################################
###################################################################################
###################################################################################
###################################################################################
###################################################################################
###################################################################################
# now i need to import all the dmdh data 
temps = [0.025, 0.1, 0.174, .25, .35, .422, .543, .658, .827, 1.008] 
magF = np.linspace(0,13, 10000)
MFTField = np.linspace(0,12,10000)
field = [[0,0,b] for b in magF]
tempMag =[]

for temperature in temps: 
    magMe = np.array([MyErObj.magnetization(ion, temperature, f) for f in field]).T
    temp = MolecularFieldTheory(MFTField, magF, -1*magMe[2], Jz)
    tempMag.append(temp) # what the actual fuck is this naming holy shi
dmdH = []
for mag in tempMag: 
    temp = np.gradient(mag, MFTField) # MFTField
    dmdH.append(temp)


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
labels = ['25mK', '100mK', '174mK', '249mK', '345mK', '422mK', '543mK', '658mK', '827mK', '1008mK']
xArrs = []
yArrs =[]

for fname in fnames: 
    temp = np.genfromtxt(fname, delimiter=',',  unpack=True, skip_header=1)
    xArrs.append(temp[0])
    yArrs.append(temp[1])

n = len(labels)
colors = plt.cm.inferno(np.linspace(0,0.8,n))

sm = plt.cm.ScalarMappable(cmap='inferno', norm=plt.Normalize(vmin=0.025, vmax=1))


plt.figure()
plt.grid(True)
# cbar = plt.colorbar(plt.cm.ScalarMappable(cmap='inferno', norm=plt.Normalize(vmin=0.025, vmax=1),ax=plt.gca()))
i = 0
for x,y,label, color, dm in zip(xArrs, yArrs, labels, colors,dmdH ): 
    y = np.log(y)
    y = y + min(y)*-1
    y = y/max(y)
    dm = np.log(dm)
    dm = dm+min(dm)*-1
    dm = dm/max(dm)
    plt.plot(x,y+i*.4, label = label, color = color)
    plt.annotate(labels[i], xy = (11, i*.4 +.1), fontsize = 9)
    plt.plot(MFTField, dm +i*.4, '--', label = label, color = color)
    i+=1
plt.title('dM/dH from SCM1 \n calculated dM/dH in dotted line')
plt.ylabel('dM/dH (arb)')
plt.xlabel('Field (T)')

# Save data to HDF5
with h5py.File('dmdh_with_scm1_data.h5', 'w') as hdf:
    hdf.create_dataset('temps', data=temps)
    hdf.create_dataset('MFTField', data=MFTField)
    hdf.create_dataset('magF', data=magF)
    hdf.create_dataset('tempMag', data=np.array(tempMag))
    hdf.create_dataset('dmdH', data=np.array(dmdH))
    hdf.create_dataset('xArrs', data=np.array(xArrs, dtype=object), dtype=h5py.vlen_dtype(float))
    hdf.create_dataset('yArrs', data=np.array(yArrs, dtype=object), dtype=h5py.vlen_dtype(float))
    hdf.create_dataset('labels', data=np.array(labels, dtype='S'))

###################################################################################
###################################################################################
###################################################################################
###################################################################################
###################################################################################
# okay, so now let's integrate
# first we have to sort gross

integratedMag = []
for x,y in zip(xArrs, yArrs): 
    x = np.array(x)
    y = np.array(y)
    y = y/y[len(y)-100] # this is kinda bad, but it's so arbitrary i don't care
    inds = x.argsort()
    x = x[inds[::1]]
    y = y[inds[::1]]
    temp = np.cumsum(y)
    integratedMag.append(temp)

n = len(labels)
colors = plt.cm.inferno(np.linspace(0,0.8,n))

plt.figure()
plt.grid(True)
i = 0
for x,y,label, color in zip(xArrs, integratedMag, labels, colors): 
    x = np.array(x)
    inds = x.argsort()
    x = x[inds[::1]]
    y = y/max(y)
    plt.plot(x, y+i*.4, label = label, color = color)
    plt.annotate(labels[i], xy = (11, i*.4 +1.05), fontsize = 9)
    plt.plot(magF, tempMag[i]/max(tempMag[i])+i*.4, '--', label = label, color = color) 
    i +=1

plt.title('integrated chi(H) \n numerically integrated from SCM1 data \n calculated curve in dotted')
plt.xlabel('Field (T)')
plt.ylabel('Magnetization (arb)')

with h5py.File('integrated_magnetization.h5', 'w') as hdf:
    hdf.create_dataset('temps', data=temps)
    hdf.create_dataset('MFTField', data=MFTField)
    hdf.create_dataset('magF', data=magF)
    hdf.create_dataset('tempMag', data=np.array(tempMag))
    hdf.create_dataset('integratedMag', data=integratedMag)
    hdf.create_dataset('xArrs', data=np.array(xArrs, dtype=object), dtype=h5py.vlen_dtype(float))
    hdf.create_dataset('yArrs', data=np.array(yArrs, dtype=object), dtype=h5py.vlen_dtype(float))
    hdf.create_dataset('labels', data=np.array(labels, dtype='S'))

###################################################################################
###################################################################################
###################################################################################
###################################################################################

# now we load the low temp chi(T)
fnames = ['/Users/hopeless/Desktop/LeeLab/data/CsErSe2_data/00358T0warm24.txt', 
            '/Users/hopeless/Desktop/LeeLab/data/CsErSe2_data/00358T0cool25.txt', 
            '/Users/hopeless/Desktop/LeeLab/data/CsErSe2_data/0116T0warm26.txt',
            '/Users/hopeless/Desktop/LeeLab/data/CsErSe2_data/0116T0cool27.txt', 
            '/Users/hopeless/Desktop/LeeLab/data/CsErSe2_data/037T0warm28.txt', 
            '/Users/hopeless/Desktop/LeeLab/data/CsErSe2_data/037T0cool29.txt',
            '/Users/hopeless/Desktop/LeeLab/data/CsErSe2_data/54T0warm32.txt', 
            '/Users/hopeless/Desktop/LeeLab/data/CsErSe2_data/064T0warm30.txt', 
            '/Users/hopeless/Desktop/LeeLab/data/CsErSe2_data/064T0cool31.txt']
labels = ['.0358T', '.0358T', '.116T', '.116T', ,'5.4T', '0.64T', '0.64T']
xArrs = []
yArrs =[]


for fname in fnames: 
    temp = np.genfromtxt(fname, delimiter=',',  unpack=True, skip_header=1)
    xArrs.append(temp[0])
    yArrs.append(temp[1])


# simulate 
fields = [0.358, 1.16, 3.7, 5.4, 6.4]
# fields = [[0,0,b] for b in fields]
tempArr = np.linspace(0.03,1.5, 100)
susArr = []
for f in fields: 
    sus = susceptibility(MyErObj, f, tempArr)
    susArr.append(sus)

n = len(labels)
colors = plt.cm.viridis(np.linspace(0,.7,n))
fig, ax1 = plt.subplots()
for i in range(len(labels)): 
    # sort first
    x = np.array(xArrs[i])
    inds = x.argsort()
    x = x[inds[::1]]
    y = yArrs[i]
    y = y[inds[::1]]
    color = colors[i]
    if i>0: 
        if labels[i]== labels[i-1]: 
            color = colors[i-1]
    ax1.plot(x, y/max(y), label = labels[i], color = color)

ax1.set_ylabel('data Chi (arbitrarily scaled)')
ax1.set_xlabel('Field (T)')

ax2 = ax1.twinx()

n = len(fields)
colors = plt.cm.viridis(np.linspace(0,0.7,n))

for i  in range(len(susArr)): 
    sus = susArr[i]
    ax2.plot(tempArr, sus/max(sus), '--',  label = str(fields[i]), color = colors[i] )
ax2.set_ylabel('calculated Chi')




# Save data to HDF5
with h5py.File('susceptibility_low_temp.h5', 'w') as hdf:
    hdf.create_dataset('xArrs', data=np.array(xArrs, dtype=object), dtype=h5py.vlen_dtype(float))
    hdf.create_dataset('yArrs', data=np.array(yArrs, dtype=object), dtype=h5py.vlen_dtype(float))
    hdf.create_dataset('labels', data=np.array(labels, dtype='S'))
    hdf.create_dataset('fields', data=fields)
    hdf.create_dataset('tempArr', data=tempArr)
    hdf.create_dataset('susArr', data=np.array(susArr))




###################################################################################
# lets do ab plane chi(T)
def susceptibilityAB(ionObj, fieldVal, temps):
    chi = []
    for temp in temps: 
        f = np.linspace(fieldVal-0.1, 1+fieldVal, 200) 
        field = [[0,b,0] for b in f]
        mag= np.array([ionObj.magnetization(ion, temp, f) for f in field]).T
        m = MolecularFieldTheory(f, f, -mag[1], Jperp)
        m = np.array(m).T
        x = np.gradient(m, f) 
        # now we've gotta access the very low field value
        valIdx = findIdx(f, fieldVal)
        chi.append(x[valIdx])
    return chi



### grab AB plane susceptibility
fnames = ['/Users/hopeless/Desktop/LeeLab/data/CsErSe2_data/scm1_chi_T_AB_-6T.txt', 
            '/Users/hopeless/Desktop/LeeLab/data/CsErSe2_data/scm1_chi_T_AB_01T.txt', 
            '/Users/hopeless/Desktop/LeeLab/data/CsErSe2_data/scm1_chi_T_AB_1T.txt']
labels = ['-6T', '0.1', '1T']
xArrs = []
yArrs =[]


for fname in fnames: 
    temp = np.genfromtxt(fname, delimiter=',',  unpack=True, skip_header=1)
    xArrs.append(temp[0])
    yArrs.append(temp[1])


fieldVal = .1
df = 0.001
temps = np.concatenate((np.linspace(1, 2, 20), np.linspace(2,20,50), np.linspace(20,300, 279)))
mysus01T = susceptibilityAB(MyErObj, fieldVal, temps) #myErObj.susceptibility(ion, temps, field, df)
# sus = ErObj.susceptibility(ion, temps, field, df)
neutronSus01T = susceptibilityAB(AllenErObj, fieldVal, temps)

fieldVal = 1
mysus1T = susceptibilityAB(MyErObj, fieldVal, temps)

fieldVal = 0
mysus0T = susceptibilityAB(MyErObj, fieldVal, temps)

fieldVal = -6
mysus6T = susceptibilityAB(MyErObj, fieldVal, temps)

field = [0,0.1,0]
susPCF = MyErObj.susceptibility(ion, temps, field, df)
susPCF = np.array(susPCF).T[1]
len(susPCF)

myinv01T = [1/x for x in mysus01T]
myinv0T = [1/x for x in mysus0T]
neutroninv01T = [1/x for x in neutronSus01T]
myinv1T = [1/x for x in mysus1T]
myinv6T = [1/x for x in mysus6T]
# sus =  [ErObj.susceptibility(ion, t, field, df) for t in temps]
susinvPCF = [-1/x for x in susPCF]
plt.figure()
plt.grid(True)
for i in range(len(labels)): 
    plt.plot(xArrs[i], yArrs[i]*1.35, label = labels[i])


# plt.plot(CESMTdata[12], 1/CESMTdata[13]*SCF, label='c-axis data from Allens paper')
plt.plot(temps, myinv01T, '--', label = 'Raman B params MFT, 0.1T')
plt.plot(temps, myinv0T, '--', label = 'Raman B params MFT, 0T')
plt.plot(temps, myinv1T, '--', label = 'Raman B params MFT, 1T')
plt.plot(temps, myinv6T, '--', label = 'Raman B params MFT, -6T')
plt.plot(temps, neutroninv01T, '-.', label = 'neutrons B params MFT' )
plt.plot(temps, susinvPCF, '--', label = 'Raman B params no MFT')
plt.title('calculated MFT susceptibility at 0.1T \n AB plane')
plt.xlabel('temperature (K)')
plt.ylabel('1/chi')
plt.legend()
# plt.xlim(0, 200)
# plt.ylim(0,10)

# Save data to HDF5
with h5py.File('susceptibility_AB_plane.h5', 'w') as hdf:
    hdf.create_dataset('xArrs', data=np.array(xArrs, dtype=object), dtype=h5py.vlen_dtype(float))
    hdf.create_dataset('yArrs', data=np.array(yArrs, dtype=object), dtype=h5py.vlen_dtype(float))
    hdf.create_dataset('labels', data=np.array(labels, dtype='S'))
    hdf.create_dataset('temps', data=temps)
    hdf.create_dataset('mysus0T', data=mysus0T)
    hdf.create_dataset('mysus1T', data=mysus1T)
    hdf.create_dataset('mysus01T', data=mysus01T)
    hdf.create_dataset('mysus6T', data=mysus6T)
    hdf.create_dataset('neutronsus01T', data=neutronSus01T)
    hdf.create_dataset('sussusPCF', data=susPCF)


###################################################################################
# make AB plane magnetization
fname20K = '/Users/hopeless/Desktop/LeeLab/data/CsErSe2_data/MPMS/CsErSe2/CsErSe2_HParAB_MvsH_20K.txt'
fname6K = '/Users/hopeless/Desktop/LeeLab/data/CsErSe2_data/MPMS/CsErSe2/CsErSe2_HParAB_MvsH_6K.txt'
fname2K = '/Users/hopeless/Desktop/LeeLab/data/CsErSe2_data/MPMS/CsErSe2/CsErSe2_HParAB_MvsH_2K.txt'

Mdata20K = np.genfromtxt(fname20K, delimiter=',',  unpack=True, skip_header=1)
Mdata6K = np.genfromtxt(fname6K, delimiter=',',  unpack=True, skip_header=1)
Mdata2K = np.genfromtxt(fname2K, delimiter=',',  unpack=True, skip_header=1)



# let's make some temperature dependent data
magF_mpms = np.linspace(-8,8,100)
field = [[0,b,0] for b in magF_mpms]
temperature = 2
magnetization2K = np.array([MyErObj.magnetization(ion, temperature, f) for f in field]).T
magnetization2K = magnetization2K[1]
magnetization2K = [m*-1 for m in magnetization2K] # make the magnetization the correct sign

temperature = 6
magnetization6K = np.array([MyErObj.magnetization(ion, temperature, f) for f in field]).T
magnetization6K = magnetization6K[1]
magnetization6K = [m*-1 for m in magnetization6K] # make the magnetization the correct sign

temperature = 20
magnetization20K = np.array([MyErObj.magnetization(ion, temperature, f) for f in field]).T
magnetization20K = magnetization20K[1]
magnetization20K = [m*-1 for m in magnetization20K] # make the magnetization the correct sign

testJperp = -1e-3
MFTField_mpms = np.linspace(-8,8,10000)
MFT2K = MolecularFieldTheory(MFTField_mpms, magF_mpms, magnetization2K, testJperp)
MFT6K = MolecularFieldTheory(MFTField_mpms, magF_mpms, magnetization6K, testJperp)
MFT20K = MolecularFieldTheory(MFTField_mpms, magF_mpms, magnetization20K, testJperp)

plt.figure()
plt.grid(True)
plt.plot(Mdata20K[0], Mdata20K[1]/1.35, 'o', label = '20K MPMS data')
plt.plot(Mdata6K[0], Mdata6K[1]/1.35, 'o', label = '6K MPMS data')
plt.plot(Mdata2K[0], Mdata2K[1]/1.35, 'o', label = '2K MPMS data')
plt.plot(CESMHdata[0]/1e4,CESMHdata[1],'b.', label='from Allens paper')
plt.plot(MFTField_mpms, MFT2K, '-', label = 'MFT 2K')
plt.plot(MFTField_mpms, MFT6K, '-', label = 'MFT 6K')
plt.plot(MFTField_mpms, MFT20K, '-', label = 'MFT 20K')
plt.plot(magF_mpms, magnetization2K, '--', label = '2K, no MFT')
plt.plot(magF_mpms, magnetization6K, '--', label = '6K, no MFT')
plt.plot(magF_mpms, magnetization20K, '--', label = '20K, no MFT')
plt.xlim(0,8)
plt.ylim(0,7)
plt.legend()
# plt.title('AB magnetization \n B20 ='+str(MyErObj.B[0])+ 'B40 = '+str(MyErObj.B[1])+'B43 =' +str(MyErObj.B[2])+'B60 =' +str(MyErObj.B[3]) + 'B63 = '+str(MyErObj.B[4])+ 'B66 = '+str(MyErObj.B[5]) +'\n Jperp = '+ str(Jperp)+' JperpAllen = '+ str(JperpAllen)+'\n Jz = '+ str(Jz)+' JzAllen = '+ str(JzAllen))
plt.title('AB magnetization, test Jperp = '+str(testJperp))
plt.xlabel('Field (T)')
plt.ylabel('Magnetization ab plane(uB/Er)')
plt.show()


# Save all the data to an HDF5 file
with h5py.File('magnetization_data_caluclation_mpms_data_AB_plane.h5', 'w') as f:
    f.create_dataset('magF_mpms', data=magF_mpms)
    f.create_dataset('MFTField_mpms', data=MFTField_mpms)
    f.create_dataset('magnetization2K', data=magnetization2K)
    f.create_dataset('magnetization6K', data=magnetization6K)
    f.create_dataset('magnetization20K', data=magnetization20K)
    f.create_dataset('MFT2K', data=MFT2K)
    f.create_dataset('MFT6K', data=MFT6K)
    f.create_dataset('MFT20K', data=MFT20K)
    f.create_dataset('Mdata2K', data=Mdata2K)
    f.create_dataset('Mdata6K', data=Mdata6K)
    f.create_dataset('Mdata20K', data=Mdata20K)
    f.create_dataset('CESMHdata', data=CESMHdata)



###################################################################################
###################################################################################
# okay, so now let's do M vs H for the AB plane

temps = [0.025, 0.045,0.1, 0.171, .25, .35, .45, .543, .827, 2, 6, 20]
n = len(temps)
colors = plt.cm.jet(np.linspace(0,1,n))
magF = np.linspace(-1,10,100)
MFTField = np.linspace(-.5,8, 1000)
field = [[0,b,0] for b in magF]
tempMagB =[]
nomftB = []

for temperature in temps: 
    magMe = np.array([MyErObj.magnetization(ion, temperature, f) for f in field]).T
    nomftB.append(-1*magMe[1])
    tempB = MolecularFieldTheory(MFTField, magF, -1*magMe[1], Jperp) 
    tempMagB.append(tempB)


# Save data to an HDF5 file
with h5py.File('M_vs_H_temperature_dependence_AB_plane_myParams_Allens_code.h5', 'w') as hdf:
    hdf.create_dataset('temps', data=temps)
    hdf.create_dataset('MFTField', data=MFTField)
    hdf.create_dataset('tempMag', data=np.array(tempMagB))

plt.figure()

for i in range(0, len(tempMagB)): 
    plt.plot(MFTField, tempMagB[i], label = str(temps[i])+'K', color = colors[i])

plt.legend()
plt.xlim(0,8)
plt.ylim(0,7)
plt.title('B axis magnetization MFT \n calculated from Raman fit B params \n test B60 = '+str(MyErObj.B[3]))
plt.xlabel('Field (T)')
plt.ylabel('Magnetization')
plt.grid(True)

########
# do again w allen's params

temps = [0.025, 0.045,0.1, 0.171, .25, .35, .45, .543, .827, 2, 6, 20]
n = len(temps)
colors = plt.cm.jet(np.linspace(0,1,n))
magF = np.linspace(0,10,100)
MFTField = np.linspace(0,8, 1000)
field = [[0,b,0] for b in magF]
tempMagB =[]
nomftB = []

for temperature in temps: 
    magAllen = np.array([AllenErObj.magnetization(ion, temperature, f) for f in field]).T
    nomftB.append(-1*magAllen[1])
    tempB = MolecularFieldTheory(MFTField, magF, -1*magAllen[1], JperpAllen) 
    tempMagB.append(tempB)


# Save data to an HDF5 file
with h5py.File('M_vs_H_temperature_dependence_AB_plane_allenParams_Allens_code.h5', 'w') as hdf:
    hdf.create_dataset('temps', data=temps)
    hdf.create_dataset('MFTField', data=MFTField)
    hdf.create_dataset('tempMag', data=np.array(tempMagB))

plt.figure()

for i in range(0, len(tempMagB)): 
    plt.plot(MFTField, tempMagB[i], label = str(temps[i])+'K', color = colors[i])

plt.legend()
plt.xlim(0,8)
plt.ylim(0,7)
plt.title('B axis magnetization MFT \n calculated from Raman fit B params \n test B60 = '+str(AllenErObj.B[3]))
plt.xlabel('Field (T)')
plt.ylabel('Magnetization')
plt.grid(True)
###################################################################################
# now let's do dM/dH for the AB plane
dMdH = []
for mag in tempMagB: 
    temporary = np.gradient(mag, MFTField)
    dMdH.append(temporary)

with h5py.File('dMdH_temperature_dependence_AB_plane.h5', 'w') as hdf:
    hdf.create_dataset('temps', data=temps)
    hdf.create_dataset('MFTField', data=MFTField)
    hdf.create_dataset('dMdH', data=np.array(dMdH))

plt.figure()

for dm, temp, color in zip(dMdH, temps, colors): 
    plt.plot(MFTField, dm, label = str(temp)+'K', color = color)

plt.legend()
plt.xlim(0,8)
plt.ylim(0,7)
plt.title('B axis dM/dH MFT \n calculated from Raman fit B params \n test B60 = '+str(MyErObj.B[3]))
plt.xlabel('Field (T)')
plt.ylabel('dM/dH (\muB/T)')
plt.grid(True)
###################################################################################
# okay, so now we plot some temp dependance

temps = np.arange(0.1, 1, .01)
n = len(temps)
colors = plt.cm.jet(np.linspace(0,1,n))

f = np.linspace(0,0.5, 10000)
field = [[0,0,b] for b in f]
tempMag =[]

for temperature in temps: 
    magMe = np.array([MyErObj.magnetization(ion, temperature, f) for f in field]).T
    temp = MolecularFieldTheory(f, f, -1*magMe[2], Jz)
    tempMag.append(temp) # what the actual fuck is this naming holy shit

# Save data to an HDF5 file
with h5py.File('M_vs_H_temperature_dependence_fine_spacing.h5', 'w') as hdf:
    hdf.create_dataset('temps', data=temps)
    hdf.create_dataset('f', data=f)
    hdf.create_dataset('tempMag', data=np.array(tempMag))

plt.figure()

for i in range(0, len(tempMag)): 
    plt.plot(f, tempMag[i], label = str(temps[i])+'K', color = colors[i])

plt.legend()
# plt.ylim(0,3)
plt.xlim(0,9)
plt.title('C axis magnetization MFT \n calculated from Raman fit B params \n test B60 = '+str(B60))
plt.xlabel('Field (T)')
plt.ylabel('Magnetization')


###################################################################################
###################################################################################
###################################################################################
###################################################################################

# okay, so now we're going to scan through a bunch of B60 vals, try a bunch of J vals
# for each B60, and plot the dm/dh to just like see whats going on


# init objects
B20 = -0.03265325 # init = -0.03559)
B40 = -0.0003849 # fixed)
B43 = -0.01393 # fixed)
# B60 =  3.079e-6# fixed)
B63 = -8.4011e-07 # init = -4.695e-06)
B66 =  3.3815e-05 # fixed)

B60Arr = [2.75e-6]#[3.0e-6, 2.9e-6, 2.8e-6]#[3.25e-6, 3.2e-6, 3.154e-6, 3.1e-6, 3.08e-6, 3.06e-6, 3.05e-6, 3.03e-6]
testJz = [.48e-3, -.48e-3, 1e-3, -1e-3, -2e-3, 2e-3, -4e-3, 4e-3, 5e-3, -5e-3]
H = np.linspace(0,8, 50)
for B60 in B60Arr: 
    g = cef.LandeGFactor(ion)
    myBparams =  {'B20': B20, 'B40':B40,'B43': B43, 'B60': B60, 'B63':B63,'B66':B66}
    MyErObj = cef.CFLevels.Bdict(ion,myBparams)
    plt.figure()
    plt.grid(True)
    temperature = 2 # I know this is where we'll see a defined peak
    # temp is low enought that we need to use my mft code
    mArr = []
    dmdhArr = []
    plt.plot(Mdata2K[0], Mdata2K[1]/1.35, 'o', label = '2K MPMS data')
    for Jz in testJz: 
        m = MFTmagC(MyErObj, H, Jz, temperature)
        mArr.append(m)
        plt.plot(H, m, label = 'Jz = '+str(Jz))
        plt.vlines(x = 5.4, ymin=0, ymax = 300)
        plt.xlabel('H [T]')
        plt.ylabel('M [\mu_B/Er')
        plt.legend()
        plt.xlim(3,8)
    plt.title('dmdh for test B60 = '+str(B60))
    fname = 'm_B60_test'+str(B60)+'.h5'
    with h5py.File(fname, 'w') as hdf:
        hdf.create_dataset('testJz', data=testJz)
        hdf.create_dataset('H', data=H)
        hdf.create_dataset('mag', data=np.array(mArr))
        hdf.create_dataset('Mdata2K', data=Mdata2K)
        hdf.attrs['B60'] = B60