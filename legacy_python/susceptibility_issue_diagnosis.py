# lets do susceptibility
def susceptibility(ionObj, fieldVal, temps):
    chi = []
    for temp in temps: 
        # f = np.arange(fieldVal-.5*fieldVal, .5*fieldVal+fieldVal, .05*fieldVal) 
        f = np.arange(fieldVal-1, fieldVal+1, .0012)
        field = [[0,0,b] for b in f]
        mag= np.array([ionObj.magnetization(ion, temp, f) for f in field]).T
        m = MolecularFieldTheory(f, f, -mag[2], Jz)
        m = np.array(m).T
        # now let's interpolate for larger spacing
        # interpF = np.linspace()
        x = np.gradient(m, f) 
        # now we've gotta access the very low field value
        valIdx = findIdx(f, fieldVal)
        chi.append(x[valIdx])
    return chi
def findIdx(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx

fieldVal = 1e-1
temps = np.linspace(.01, 2, 10)
mysus01T = susceptibility(MyErObj, fieldVal, temps)
neutronSus01T = susceptibility(AllenErObj, fieldVal, temps)


myinv01T = [1/x for x in mysus01T]
neutroninv01T = [1/x for x in neutronSus01T]
# myinv1T = [-1/x for x in mysus1T]
# neutroninv1T = [-1/x for x in neutronSus1T]
# myinv3T = [-1/x for x in mysus3T]
# neutroninv3T = [-1/x for x in neutronSus3T]



susinvPCF = [-1/x for x in susPCF]
plt.figure()
# for i in range(len(labels)): 
#     plt.plot(xArrs[i], yArrs[i]*1.35, label = labels[i])

# add manual data
xarr = [0.01, 0.025, 0.045, 0.1, .17, 0.25, .35, 0.45, 0.543, 0.827,2]
yarr = [.4,.425,.425,.575,2.35,6,10.25, 11.25,10.325,7.586,3.27]
yinv = [1/y for y in yarr]

plt.plot(xarr, yarr, 'o', label = 'chi extracted manually from Mvs H')
# plt.plot(CESMTdata[12], 1/CESMTdata[13]*SCF, label='c-axis data from Allens paper')
plt.plot(temps, mysus01T, '--', label = 'Raman B params MFT 0.1T')
plt.plot(temps, neutronSus01T, '-.', label = 'neutrons B params MFT 0.1T' )


plt.title('calculated MFT susceptibility at'+str(fieldVal)+'T')
plt.xlabel('temperature (K)')
plt.ylabel('chi')
plt.legend()



with h5py.File('sus_diagnosis.h5', 'w') as hdf:
    hdf.create_dataset('field', data=f)
    hdf.create_dataset('mag', data=-mag[2])
    hdf.create_dataset('m', data=m)


tempArr = [0.025, .1, .175, .249, .345, .422, .54, .65, .8, 1.008]
dmArr = [5e-4, 6e-4, .0028, .0072, .0119, .0135, .0121, .0114, .0094, .0077]
dmFixed = [d/.0012 for d in dmArr]

plt.figure()
# for i in range(len(labels)): 
#     plt.plot(xArrs[i], yArrs[i]*1.35, label = labels[i])

# add manual data
xarr = [0.01, 0.025, 0.045, 0.1, .17, 0.25, .35, 0.45, 0.543, 0.827,2]
yarr = [.4,.425,.425,.575,2.35,6,10.25, 11.25,10.325,7.586,3.27]
yinv = [1/y for y in yarr]

plt.plot(xarr, yarr, 'o', label = 'chi extracted manually from Mvs H')
# plt.plot(CESMTdata[12], 1/CESMTdata[13]*SCF, label='c-axis data from Allens paper')
plt.plot(temps, mysus01T, '--', label = 'Raman B params MFT 0.1T')
plt.plot(temps, neutronSus01T, '-.', label = 'neutrons B params MFT 0.1T' )
plt.plot(tempArr, dmFixed,'x', label = 'from dmdh curve')

plt.title('calculated MFT susceptibility at'+str(fieldVal)+'T')
plt.xlabel('temperature (K)')
plt.ylabel('chi')
plt.legend()

