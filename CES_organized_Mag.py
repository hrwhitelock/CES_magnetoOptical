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
magF = np.linspace(0,10,1000)
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
    n = 10
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
# plt.xlim(0,6)
# plt.ylim(0,10)
plt.legend()
plt.title('C magnetization')
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

print("Data saved to 'magnetization_data.h5'.")

###################################################################################
# Now let's nail down B60
# init objects
B20 = -0.03265325 # init = -0.03559)
B40 = -0.0003849 # fixed)
B43 = -0.01393 # fixed)
B60 =  3.06e-6# fixed)
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

# this just really isn't fucking good enough
# you can't just iteratively change variables by hand and seriously believe that
# I do not believe this and will NOT put my name on a paper that does that

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

# Save data to an HDF5 file
with h5py.File('M_vs_H_temperature_dependence.h5', 'w') as hdf:
    hdf.create_dataset('temps', data=temps)
    hdf.create_dataset('MFTField', data=MFTField)
    hdf.create_dataset('tempMag', data=np.array(tempMag))

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

plt.figure()
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
# lets do susceptibility
def susceptibility(ionObj, fieldVal, temps):
    chi = []
    for temp in temps: 
        f = np.linspace(fieldVal-0.1, 0.1+fieldVal, 100) 
        field = [[0,0,b] for b in f]
        mag= np.array([ionObj.magnetization(ion, temp, f) for f in field]).T
        m = MolecularFieldTheory(f, f, mag[2], Jz)
        m = np.array(m).T
        x = np.gradient(m, f) 
        # now we've gotta access the very low field value
        valIdx = findIdx(field, [0,0,fieldVal])
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

fieldVal = .1
df = 0.001
temps = np.concatenate((np.linspace(.001, 2, 100), np.linspace(2,20,100), np.linspace(20,300, 279)))
mysus = susceptibility(MyErObj, fieldVal, temps) #myErObj.susceptibility(ion, temps, field, df)
# sus = ErObj.susceptibility(ion, temps, field, df)
neutronSus = susceptibility(AllenErObj, fieldVal, temps)


field = [0,0,.1]
susPCF = MyErObj.susceptibility(ion, temps, field, df)
susPCF = np.array(susPCF).T[2]
len(susPCF)

myinv = [-1/x for x in mysus]
neutroninv = [-1/x for x in neutronSus]
# sus =  [ErObj.susceptibility(ion, t, field, df) for t in temps]
susinvPCF = [-1/x for x in susPCF]
plt.figure()
for i in range(len(labels)): 
    plt.plot(xArrs[i], yArrs[i], label = labels[i])

plt.plot(CESMTdata[12], 1/CESMTdata[13]*SCF, label='c-axis data from Allens paper')
plt.plot(temps, myinv, '--', label = 'Raman B params MFT')
plt.plot(temps, neutroninv, '-.', label = 'neutrons B params MFT' )
plt.plot(temps, susinvPCF, '--', label = 'Raman B params no MFT')
plt.title('calculated MFT susceptibility at 0.1T ')
plt.xlabel('temperature (K)')
plt.ylabel('1/chi')
plt.legend()
plt.xlim(0, 200)
plt.ylim(0,10)


###################################################################################
# now i need to import all the fucking dmdh data because heave for fucking bid that be in an easy format
temps = [0.025, 0.045,0.1, 0.174, .25, .35, .422, .543, .658, .827, 1.008] 
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
    temp = np.gradient(mag)
    dmdH.append(temp)


fnames = ['/Users/hopeless/Desktop/LeeLab/data/CsErSe2_data/25mKdn022.txt', 
            '/Users/hopeless/Desktop/LeeLab/data/CsErSe2_data/45mKdn035.txt', 
            '/Users/hopeless/Desktop/LeeLab/data/CsErSe2_data/100mKUp021.txt',
            '/Users/hopeless/Desktop/LeeLab/data/CsErSe2_data/174mKDn020.txt', 
            '/Users/hopeless/Desktop/LeeLab/data/CsErSe2_data/249mKUp019.txt', 
            '/Users/hopeless/Desktop/LeeLab/data/CsErSe2_data/345mKDn018.txt', 
            '/Users/hopeless/Desktop/LeeLab/data/CsErSe2_data/422mKUp017.txt', 
            '/Users/hopeless/Desktop/LeeLab/data/CsErSe2_data/543mKDn016.txt',
            '/Users/hopeless/Desktop/LeeLab/data/CsErSe2_data/658mKUp015.txt', 
            '/Users/hopeless/Desktop/LeeLab/data/CsErSe2_data/827mKDn013-14.txt', 
            '/Users/hopeless/Desktop/LeeLab/data/CsErSe2_data/1008mKdn033.txt']
labels = ['25mK', '45mK', '100mK', '174mK', '249mK', '345mK', '422mK', '543mK', '658mK', '827mK', '1008mK']
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
# cbar = plt.colorbar(plt.cm.ScalarMappable(cmap='inferno', norm=plt.Normalize(vmin=0.025, vmax=1),ax=plt.gca()))
i = 0
for x,y,label, color, dm in zip(xArrs, yArrs, labels, colors,dmdH ): 
    y = y/max(y)
    dm = np.log(dm)
    dm = dm+min(dm)*-1
    dm = dm/max(dm)
    plt.plot(x,y+i*.5, label = label, color = color)
    plt.plot(MFTField, dm +i*.5, '--', label = label, color = color)
    i+=1
plt.title('dM/dH from SCM1 \n calculated dM/dH in dotted line')
plt.ylabel('dM/dH (arb)')
plt.xlabel('Field (T)')

# Save data to HDF5
with h5py.File('magnetization_analysis.h5', 'w') as hdf:
    hdf.create_dataset('temps', data=temps)
    hdf.create_dataset('MFTField', data=MFTField)
    hdf.create_dataset('magF', data=magF)
    hdf.create_dataset('tempMag', data=np.array(tempMag))
    hdf.create_dataset('dmdH', data=np.array(dmdH))
    hdf.create_dataset('xArrs', data=np.array(xArrs, dtype=object), dtype=h5py.vlen_dtype(float))
    hdf.create_dataset('yArrs', data=np.array(yArrs, dtype=object), dtype=h5py.vlen_dtype(float))
    hdf.create_dataset('labels', data=np.array(labels, dtype='S'))

###################################################################################
# okay, so now let's integrate
# first we have to sort gross

integratedMag = []
for x,y in zip(xArrs, yArrs): 
    x = np.array(x)
    y = np.array(y)
    y = y/max(y)
    inds = x.argsort()
    x = x[inds[::1]]
    y = y[inds[::1]]
    temp = np.cumsum(y)
    integratedMag.append(temp)

n = len(labels)
colors = plt.cm.inferno(np.linspace(0,0.8,n))

plt.figure()
i = 0
for x,y,label, color in zip(xArrs, integratedMag, labels, colors): 
    x = np.array(x)
    inds = x.argsort()
    x = x[inds[::1]]
    y = y/max(y)
    plt.plot(x, y+i*.5, label = label, color = color)
    plt.plot(magF, tempMag[i]/max(tempMag[i])+i*.5, '--', label = label, color = color) 
    i +=1

plt.title('integrated chi(H) \n numerically integrated from SCM1 data \n calculated curve in dotted')
plt.xlabel('Field (T)')
plt.ylabel('Magnetization (arb)')

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
labels = ['.358T', '.358T', '1.16T', '1.16T', '3.7T', '3.7T','5.4T', '6.4T', '6.4T']
xArrs = []
yArrs =[]


for fname in fnames: 
    temp = np.genfromtxt(fname, delimiter=',',  unpack=True, skip_header=1)
    xArrs.append(temp[0])
    yArrs.append(temp[1])


# simulate 
fields = [0.3, 1, 4, 5.5, 6]
# fields = [[0,0,b] for b in fields]
tempArr = np.linspace(0.01,1, 100 )
susArr = []
for f in fields: 
    sus = susceptibility(MyErObj, f, tempArr)
    susArr.append(sus)

n = len(labels)
colors = plt.cm.cool(np.linspace(0,0.8,n))
plt.figure()
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
    plt.plot(x, y/y[-1], label = labels[i], color = color)

n = len(fields)
colors = plt.cm.cool(np.linspace(0,0.8,n))

for i  in range(len(susArr)): 
    plt.plot(tempArr, susArr[i], label = str(fields[i]), color = colors[i] )





###################################################################################
# lets do ab plane chi(T)
def susceptibilityAB(ionObj, fieldVal, temps):
    chi = []
    for temp in temps: 
        f = np.linspace(fieldVal-0.1, 0.1+fieldVal, 100) 
        field = [[0,b,0] for b in f]
        mag= np.array([ionObj.magnetization(ion, temp, f) for f in field]).T
        m = MolecularFieldTheory(f, f, mag[1], Jperp)
        m = np.array(m).T
        x = np.gradient(m, f) 
        # now we've gotta access the very low field value
        valIdx = findIdx(field, [0,fieldVal,0])
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
temps = np.concatenate((np.linspace(.001, 2, 100), np.linspace(2,20,100), np.linspace(20,300, 279)))
mysus = susceptibilityAB(MyErObj, fieldVal, temps) #myErObj.susceptibility(ion, temps, field, df)
# sus = ErObj.susceptibility(ion, temps, field, df)
neutronSus = susceptibilityAB(AllenErObj, fieldVal, temps)


field = [0,0.1,0]
susPCF = MyErObj.susceptibility(ion, temps, field, df)
susPCF = np.array(susPCF).T[1]
len(susPCF)

myinv = [-1/x for x in mysus]
neutroninv = [-1/x for x in neutronSus]
# sus =  [ErObj.susceptibility(ion, t, field, df) for t in temps]
susinvPCF = [-1/x for x in susPCF]
plt.figure()
for i in range(len(labels)): 
    plt.plot(xArrs[i], yArrs[i], label = labels[i])


# plt.plot(CESMTdata[12], 1/CESMTdata[13]*SCF, label='c-axis data from Allens paper')
plt.plot(temps, myinv, '--', label = 'Raman B params MFT')
plt.plot(temps, neutroninv, '-.', label = 'neutrons B params MFT' )
plt.plot(temps, susinvPCF, '--', label = 'Raman B params no MFT')
plt.title('calculated MFT susceptibility at 0.1T \n AB plane')
plt.xlabel('temperature (K)')
plt.ylabel('1/chi')
plt.legend()
plt.xlim(0, 200)
plt.ylim(0,10)

# now save as csv so I can load into matlab 
# this is so fucking stupid 

mysus = pd.DataFrame(data = np.array(myinv).T, index = temps)
mysus.to_csv('/Users/hopeless/Desktop/LeeLab/data/CsErSe2_data/calculated_data/mysus.csv')

neutronsus = pd.DataFrame(data = np.array(neutroninv).T, index = temps)
neutronsus.to_csv('/Users/hopeless/Desktop/LeeLab/data/CsErSe2_data/calculated_data/neutronsus.csv')

