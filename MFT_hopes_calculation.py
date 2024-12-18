## scratch MFT  code

def bmag(mag, J,  h, temperature):
    newh = h +q*J*mag/muB/1.2/1.2
    mag = -1*MyErObj.magnetization(ion, temperature, [0,newh, 0]).T[1]-mag
    return mag

def cmag(mag, J,  h, temperature):
    newh = h +q*J*mag/muB/1.2/1.2
    mag = -1*AllenErObj.magnetization(ion, temperature, [0,0, newh]).T[2]-mag
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
            mag = fsolve(cmag, mag, args = (J, h, temperature))[0] # fsolve spits out an array - this time its one val
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
            mag = fsolve(bmag, mag, args = (J, h, temperature))[0] # fsolve spits out an array - this time its one val
        magArr.append(mag) # fso
    return magArr



temps = [0.025, 0.045,0.1, 0.171, .25, .35, .45, .543, .827, 2, 6, 20]
labels = [str(T)+ 'K' for T in temps]
n = len(temps)
colors = plt.cm.jet(np.linspace(0,1,n))
tempMagB = []
tempMagC = []
H = np.linspace(0,10, 100)
Jz = -2.53e-3#0.48e-3
Jperp = -0.9e-3

tempMagB =[]
for temperature in temps: 
    # temp = MFTmagB(MyErObj, H, Jperp, temperature)
    # tempMagB.append(temp)
    temp = MFTmagC(AllenErObj, H, Jz, temperature)
    tempMagC.append(temp)

B20 = -3.559e-2
B40 = -3.849e-4
B43 = -1.393e-2
B60 =  3.154e-6
B63 = -4.695e-6
B66 =  3.3815e-5
Jz = -2.53e-3


plt.figure()
for mag, color, label in zip(tempMagB, colors, labels): 
    plt.plot(H, mag, '--o', color = color, label = label)
plt.title('B axis MFT magnetization \n calcualted from my mft function \n Jperp = '+str(Jperp))
plt.ylabel('B axis MFT magnetization')
plt.xlabel('Field (T)')
plt.legend()
plt.grid(True)

plt.figure()
for mag, color, label in zip(tempMagC, colors, labels): 
    plt.plot(H, mag, '--o', color = color, label = label)
plt.title('C axis MFT magnetization \n calcualted from my mft function\n Jz = '+str(Jz))
plt.ylabel('C axis MFT magnetization')
plt.xlabel('Field (T)')
plt.legend()
plt.grid(True)

##########################################################################################################
##########################################################################################################
##########################################################################################################
##########################################################################################################
# okay so now let's add the data
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

# calculate Allen's data
magF = np.linspace(0,10,1000)
MFTField = np.linspace(0,8,10000)
temperature = 2
field = [[0,0,h] for h in magF]
allenCaxisMagnetization = np.array([AllenErObj.magnetization(ion, temperature, f) for f in field]).T
allenCaxisMagnetization = allenCaxisMagnetization[2]
allenCaxisMagnetization = [m*-1 for m in allenCaxisMagnetization]

allenMFTCaxis = MolecularFieldTheory(MFTField, magF, allenCaxisMagnetization, JzAllen)

# load MPMS data
fname20K = '/Users/hopeless/Desktop/LeeLab/data/CsErSe2_data/MPMS/CsErSe2/CsErSe2_HParC_MvsH_20K.txt'
fname6K = '/Users/hopeless/Desktop/LeeLab/data/CsErSe2_data/MPMS/CsErSe2/CsErSe2_HParC_MvsH_6K.txt'
fname2K = '/Users/hopeless/Desktop/LeeLab/data/CsErSe2_data/MPMS/CsErSe2/CsErSe2_HParC_MvsH_2K.txt'

Mdata20K = np.genfromtxt(fname20K, delimiter=',',  unpack=True, skip_header=1)
Mdata6K = np.genfromtxt(fname6K, delimiter=',',  unpack=True, skip_header=1)
Mdata2K = np.genfromtxt(fname2K, delimiter=',',  unpack=True, skip_header=1)

# plot everything
plt.figure()
plt.grid(True)
plt.plot(Mdata20K[0], Mdata20K[1]/1.35, 'o', label = '20K MPMS data')
plt.plot(Mdata6K[0], Mdata6K[1]/1.35, 'o', label = '6K MPMS data')
plt.plot(Mdata2K[0], Mdata2K[1]/1.35, 'o', label = '2K MPMS data')
plt.plot(CESMHdata[6]/1e4,CESMHdata[7],'b.', label='from Allens paper')
plt.plot(MFTField, allenMFTCaxis, 'b--', label = 'Allens 2K MFT calculation')
for mag, color, label in zip(tempMagC, colors, labels): 
    plt.plot(H, mag, '-', color = color, label = label)
plt.xlim(0,7)
plt.ylim(0,8)
plt.legend()
plt.title('C magnetization \n B20 ='+str(MyErObj.B[0])+ 'B40 = '+str(MyErObj.B[1])+'B43 =' +str(MyErObj.B[2])+'B60 =' +str(MyErObj.B[3]) + 'B63 = '+str(MyErObj.B[4])+ 'B66 = '+str(MyErObj.B[5]) +'\n Jz = '+ str(Jz)+'\n JzAllen = '+ str(JzAllen))
# plt.title('c-axis magnetization with test Jz = ' + str(testJz))
plt.xlabel('Field (T)')
plt.ylabel('Magnetization (uB/Er)')
plt.show()

with h5py.File('hopes_fitted_Jz_MFT_c_axis_calculation.h5', 'w') as f:
    f.create_dataset('MFTField', data=MFTField)
    f.create_dataset('magnetization2K', data=Mdata2K)
    f.create_dataset('magnetization6K', data=Mdata6K)
    f.create_dataset('magnetization20K', data=Mdata20K)
    f.create_dataset('allenMFTCaxis', data= allenMFTCaxis)
    f.create_dataset('CESMHdata', data=CESMHdata)
    f.create_dataset('tempMagC', data = tempMagC)
    f.create_dataset('calculationH', data = H)
    f.create_dataset('labels', data=np.array(labels, dtype='S'))
    f.attrs['B20'] = B20
    f.attrs['B40'] = B40
    f.attrs['B43'] = B43
    f.attrs['B60'] = B60
    f.attrs['B63'] = B63
    f.attrs['B66'] = B66
    f.attrs['Jperp'] = Jperp
    f.attrs['Jz'] = Jz
    f.attrs['JperpAllen'] = JperpAllen
    f.attrs['JzAllen'] = JzAllen

##########################################################################################################
##########################################################################################################
##########################################################################################################
##########################################################################################################
# now let's do the same thing for the AB plane (b-axis)

# import mpms data
fname20K = '/Users/hopeless/Desktop/LeeLab/data/CsErSe2_data/MPMS/CsErSe2/CsErSe2_HParAB_MvsH_20K.txt'
fname6K = '/Users/hopeless/Desktop/LeeLab/data/CsErSe2_data/MPMS/CsErSe2/CsErSe2_HParAB_MvsH_6K.txt'
fname2K = '/Users/hopeless/Desktop/LeeLab/data/CsErSe2_data/MPMS/CsErSe2/CsErSe2_HParAB_MvsH_2K.txt'

Mdata20K = np.genfromtxt(fname20K, delimiter=',',  unpack=True, skip_header=1)
Mdata6K = np.genfromtxt(fname6K, delimiter=',',  unpack=True, skip_header=1)
Mdata2K = np.genfromtxt(fname2K, delimiter=',',  unpack=True, skip_header=1)

# make Allen's calculation
magF = np.linspace(0,10,1000)
MFTField = np.linspace(0,8,10000)
temperature = 2
field = [[0,h,0] for h in magF]
allenBaxisMagnetization = np.array([AllenErObj.magnetization(ion, temperature, f) for f in field]).T
allenBaxisMagnetization = allenBaxisMagnetization[1]
allenBaxisMagnetization = [m*-1 for m in allenBaxisMagnetization]

allenMFTBaxis = MolecularFieldTheory(MFTField, magF, allenBaxisMagnetization, JperpAllen)

# plot everything
plt.figure()
plt.grid(True)
plt.plot(Mdata20K[0], Mdata20K[1]/1.35, 'o', label = '20K MPMS data')
plt.plot(Mdata6K[0], Mdata6K[1]/1.35, 'o', label = '6K MPMS data')
plt.plot(Mdata2K[0], Mdata2K[1]/1.35, 'o', label = '2K MPMS data')
plt.plot(CESMHdata[0]/1e4,CESMHdata[1],'b.', label='from Allens paper')
plt.plot(MFTField, allenMFTBaxis, 'b--', label = 'Allens 2K MFT calculation')
for mag, color, label in zip(tempMagB, colors, labels): 
    plt.plot(H, mag, '-', color = color, label = label)
plt.xlim(0,7)
plt.ylim(0,8)
plt.legend()
plt.title('B magnetization \n B20 ='+str(MyErObj.B[0])+ 'B40 = '+str(MyErObj.B[1])+'B43 =' +str(MyErObj.B[2])+'B60 =' +str(MyErObj.B[3]) + 'B63 = '+str(MyErObj.B[4])+ 'B66 = '+str(MyErObj.B[5]) +'\n Jperp = '+ str(Jperp)+'\n JperpAllen = '+ str(JperpAllen))
plt.xlabel('Field (T)')
plt.ylabel('Magnetization (uB/Er)')
plt.show()

with h5py.File('hopes_MFT_b_axis_calculation.h5', 'w') as f:
    f.create_dataset('MFTField', data=MFTField)
    f.create_dataset('magnetization2K', data=Mdata2K)
    f.create_dataset('magnetization6K', data=Mdata6K)
    f.create_dataset('magnetization20K', data=Mdata20K)
    f.create_dataset('allenMFTBaxis', data= allenMFTBaxis)
    f.create_dataset('CESMHdata', data=CESMHdata)
    f.create_dataset('tempMagB', data = tempMagB)
    f.create_dataset('calculationH', data = H)
    f.create_dataset('labels', data=np.array(labels, dtype='S'))
    f.attrs['B20'] = B20
    f.attrs['B40'] = B40
    f.attrs['B43'] = B43
    f.attrs['B60'] = B60
    f.attrs['B63'] = B63
    f.attrs['B66'] = B66
    f.attrs['Jperp'] = Jperp
    f.attrs['Jz'] = Jz
    f.attrs['JperpAllen'] = JperpAllen
    f.attrs['JzAllen'] = JzAllen


    ##########################################################################################################
##########################################################################################################
##########################################################################################################
##########################################################################################################
# now, at 0.5K, let's look at the new MFT fn and Allen's MFT fn on AB plane

# first, calculation with Allen's code
temperature = 0.5
magF = np.linspace(-1,10,100)
MFTField = np.linspace(-.5,8, 1000)
field = [[0,b,0] for b in magF]
tempMagB =[]
nomftB = []

mag_noMFT = np.array([MyErObj.magnetization(ion, temperature, f) for f in field]).T
tempB = MolecularFieldTheory(MFTField, magF, -1*mag_noMFT[1], Jperp) 


# now do my mft
H = MFTField
myMFT = MFTmagB(MyErObj, H, Jperp, temperature)


plt.figure()
plt.plot(magF, -mag_noMFT[1], label = 'PCF calculation, no MFT')
plt.plot(MFTField, tempB, label = 'MFT calculation, Allens code')
plt.plot(H, myMFT, label = 'MFT calculation, my code')
plt.grid(True)
plt.title('B axis magnetization at 0.5K,\n B20 ='+str(MyErObj.B[0])+ 'B40 = '+str(MyErObj.B[1])+'B43 =' +str(MyErObj.B[2])+'B60 =' +str(MyErObj.B[3]) + 'B63 = '+str(MyErObj.B[4])+ 'B66 = '+str(MyErObj.B[5]) +'\n Jperp = '+ str(Jperp)+'\n JperpAllen = '+ str(JperpAllen))
plt.xlabel('Field (T)')
plt.ylabel('Magnetization ab-plane')

# regenerate AB plane data
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
with h5py.File('M_vs_H_temperature_dependence_AB_plane_low_temp.h5', 'w') as hdf:
    hdf.create_dataset('temps', data=temps)
    hdf.create_dataset('MFTField', data=MFTField)
    hdf.create_dataset('tempMag', data=np.array(tempMagB))


        ##########################################################################################################
##########################################################################################################
##########################################################################################################
##########################################################################################################

# okay, so now let's make the susceptibility
def susceptibility(fieldVal, temps):
    chi = []
    for temp in temps: 
        # f = np.arange(fieldVal-.5*fieldVal, .5*fieldVal+fieldVal, .05*fieldVal) 
        f = np.arange(fieldVal-0.2, fieldVal+0.2, .012)
        m = MFTmagB(MyErObj, f, Jz, temperature)
        m = np.array(m).T
        x = np.gradient(m, f) 
        valIdx = findIdx(f, fieldVal)
        chi.append(x[valIdx])
    return chi
def findIdx(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx

temps = np.linspace(0.1, 300, 300)
fieldVal = 0
susC = susceptibility(fieldVal, temps)
