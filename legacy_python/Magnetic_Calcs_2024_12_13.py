
# okay, so here I'm doing a bit more organized MFT calculations
# from fitted J params, all J fits done from spectroscopy

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
MyErObj = cef.CFLevels.Bdict(ion,myBparams)

B20 = -3.559e-2
B40 = -3.849e-4
B43 = -1.393e-2
B60 = 3.154e-6
B63 = -4.695e-6
B66 = 3.3815e-5
JperpAllen = -0.2e-3
JzAllen = -2.4e-3
g = cef.LandeGFactor(ion)
AllenBparams =  {'B20': B20, 'B40':B40,'B43': B43, 'B60': B60, 'B63':B63,'B66':B66}
AllenErObj = cef.CFLevels.Bdict(ion,AllenBparams)


temps = [0.025, 0.05,0.1, 0.171, .25, .35, .45, .543, .827,1, 2, 6, 20]
labels = [str(T)+ 'K' for T in temps]
n = len(temps)
colors = plt.cm.jet(np.linspace(0,1,n))
tempMagB_me = []
tempMagC_me = []
tempMagB_Allen = []
tempMagC_Allen = []
H = np.linspace(0,8, 100)


for temperature in temps: 
    temp = MFTmagB(MyErObj, H, Jperp, temperature)
    tempMagB_me.append(temp)
    temp = MFTmagC(MyErObj, H, Jz, temperature)
    tempMagC_me.append(temp)
    temp = MFTmagB(AllenErObj, H, JperpAllen, temperature)
    tempMagB_Allen.append(temp)
    temp = MFTmagC(AllenErObj, H, JzAllen, temperature)
    tempMagC_Allen.append(temp)

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
for mag, color, label in zip(tempMagC_me, colors, labels): 
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
    f.create_dataset('tempMagC_me', data = tempMagC_me)
    f.create_dataset('tempMagC_Allen', data = tempMagC_Allen)
    f.create_dataset('calculationH', data = H)
    f.create_dataset('labels', data=np.array(labels, dtype='S'))
    f.attrs['B20'] = MyErObj.B[0]
    f.attrs['B40'] = MyErObj.B[1]
    f.attrs['B43'] = MyErObj.B[2]
    f.attrs['B60'] = MyErObj.B[3]
    f.attrs['B63'] = MyErObj.B[4]
    f.attrs['B66'] = MyErObj.B[5]
    f.attrs['allen_B20'] = AllenErObj.B[0]
    f.attrs['allen_B40'] = AllenErObj.B[1]
    f.attrs['allen_B43'] = AllenErObj.B[2]
    f.attrs['allen_B60'] = AllenErObj.B[3]
    f.attrs['allen_B63'] = AllenErObj.B[4]
    f.attrs['allen_B66'] = AllenErObj.B[5]
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
for mag, color, label in zip(tempMagB_Allen, colors, labels): 
    plt.plot(H, mag, '-', color = color, label = label)
plt.xlim(0,7)
plt.ylim(0,8)
plt.legend()
plt.title('B magnetization \n B20 ='+str(MyErObj.B[0])+ 'B40 = '+str(MyErObj.B[1])+'B43 =' +str(MyErObj.B[2])+'B60 =' +str(MyErObj.B[3]) + 'B63 = '+str(MyErObj.B[4])+ 'B66 = '+str(MyErObj.B[5]) +'\n Jperp = '+ str(Jperp)+'\n JperpAllen = '+ str(JperpAllen))
plt.xlabel('Field (T)')
plt.ylabel('Magnetization (uB/Er)')
plt.show()

##########################################################################################################
##########################################################################################################
##########################################################################################################
##########################################################################################################
# let's calculate dmdh
temps = [0.025, 0.1, 0.174, .25, .35, .422, .543, .658, .827, 1.008]
# recalculatign to match data perfectly :'(
labels = [str(T)+ 'K' for T in temps]
n = len(temps)
colors = plt.cm.jet(np.linspace(0,1,n))
tempMagC_Allen = []
tempMagC_me = []

# temps = [0.0001]

H = np.concatenate((np.linspace(0,1,30), np.linspace(1.01,10, 100)))
# H = np.linspace(0,1,100)

for temperature in temps: 
    # temp = MFTmagC(MyErObj, H, Jz, temperature)
    # tempMagC_me.append(temp)
    temp = MFTmagC(AllenErObj, H, JzAllen, temperature)
    tempMagC_Allen.append(temp)


dmdhArr_me = []
dmdhArr_Allen = []

for mag_me, mag_Allen in zip(tempMagC_me,tempMagC_Allen): 
    # dmdhArr_me.append(np.gradient(mag_me, H))
    dmdhArr_Allen.append(np.gradient(mag_Allen, H))

# load dmdh data
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
for x,y,label, color, dm in zip(xArrs, yArrs, labels, colors,dmdhArr_me ): 
    y = np.log(y)
    y = y + min(y)*-1
    y = y/max(y)
    dm = np.log(dm)
    dm = dm+min(dm)*-1
    dm = dm/max(dm)
    plt.plot(x,y+i*.4, label = label, color = color)
    plt.annotate(labels[i], xy = (11, i*.4 +.1), fontsize = 9)
    plt.plot(H, dm +i*.4, '--', label = label, color = color)
    i+=1
plt.title('dM/dH from SCM1 \n my calculated dM/dH in dotted line')
plt.ylabel('dM/dH (arb)')
plt.xlabel('Field (T)')

plt.figure()
plt.grid(True)
# cbar = plt.colorbar(plt.cm.ScalarMappable(cmap='inferno', norm=plt.Normalize(vmin=0.025, vmax=1),ax=plt.gca()))
i = 0
for x,y,label, color, dm in zip(xArrs, yArrs, labels, colors,dmdhArr_me ): 
    y = np.log(y)
    y = y + min(y)*-1
    y = y/max(y)
    dm = np.log(dm)
    dm = dm+min(dm)*-1
    dm = dm/max(dm)
    plt.plot(x,y+i*.4, label = label, color = color)
    plt.annotate(labels[i], xy = (11, i*.4 +.1), fontsize = 9)
    plt.plot(H, dm +i*.4, '--', label = label, color = color)
    i+=1
plt.title('dM/dH from SCM1 \n Allens calculated dM/dH in dotted line')
plt.ylabel('dM/dH (arb)')
plt.xlabel('Field (T)')


# save data
with h5py.File('dmdh_with_scm1_data.h5', 'w') as hdf:
    hdf.create_dataset('temps', data=temps)
    hdf.create_dataset('MFTField', data=MFTField)
    hdf.create_dataset('magF', data=magF)
    hdf.create_dataset('tempMagC_me', data=np.array(tempMagC_me))
    hdf.create_dataset('tempMagC_Allen', data=np.array(tempMagC_Allen))
    hdf.create_dataset('dmdH', data=np.array(dmdH))
    hdf.create_dataset('xArrs', data=np.array(xArrs, dtype=object), dtype=h5py.vlen_dtype(float))
    hdf.create_dataset('yArrs', data=np.array(yArrs, dtype=object), dtype=h5py.vlen_dtype(float))
    hdf.create_dataset('labels', data=np.array(labels, dtype='S'))