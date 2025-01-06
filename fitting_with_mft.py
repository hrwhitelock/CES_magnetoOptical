# okay, so I really want to do better fits. I just feel like these fits are like, fine
# I know it's fine, but I also know that I can do better


## let me just add mft to the c axis zeeman and see what that does
def cmag(mag, J,  h, temperature, ionObj):
    newh = h +q*J*mag/muB/1.2/1.2
    mag = -1*ionObj.magnetization(ion, temperature, [0,0, newh]).T[2]-mag
    return mag

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

def zeemanSplitC(field, wavenum, B20, B40, B43, B60, B63, B66, Jz, temperature):    
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
        wid = 1
        centers = dE
        tempfun = lorentzian(wavenum, 0, dEphonon, phononSig)
        for i in range(len(centers)):
            a = tempAmp[i]
            tempfun += lorentzian(wavenum, a, centers[i]*meVToCm, wid)
        fun.append(tempfun)
    return fun
#model 1 from Allen's paper
B20 = -3.559e-2
B40 = -3.849e-4
B43 = -1.393e-2
B60 =  3.154e-6
B63 = -4.695e-6
B66 =  3.3815e-5
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
params['Jz'].set(value = -3e-3)

z = np.array(fitData.to_numpy()) # gotta do it twice with tuples :((((
z = z.T

result = model.fit(z, field=field, wavenum=wavenum, params =params)

print(result.fit_report())

''' result!!! 
[[Model]]
    Model(zeemanSplitC)
[[Fit Statistics]]
    # fitting method   = leastsq
    # function evals   = 616
    # data points      = 95140
    # variables        = 4
    chi-square         = 2584.04379
    reduced chi-square = 0.02716158
    Akaike info crit   = -343066.265
    Bayesian info crit = -343028.413
    R-squared          = -2.98799518
[[Variables]]
    B20: -0.03689049 +/- 1.2695e-04 (0.34%) (init = -0.03559)
    B40: -0.0003849 (fixed)
    B43: -0.01393 (fixed)
    B60:  3.1631e-06 +/- 4.1188e-09 (0.13%) (init = 3.154e-06)
    B63: -3.3295e-06 +/- 1.1926e-07 (3.58%) (init = -4.695e-06)
    B66:  3.3815e-05 (fixed)
    Jz:  -0.00256847 +/- 1.6205e-05 (0.63%) (init = -0.001)
[[Correlations]] (unreported correlations are < 0.100)
    C(B20, B63) = -0.8695
    C(B20, B60) = -0.5207
    C(B20, Jz)  = +0.3671
    C(B60, B63) = +0.3495
    C(B60, Jz)  = +0.1684
'''

def zeemanSplitLinesC(field, B20, B40, B43, B60, B63, B66, Jz, temperature):     
    # assuming only H||B rn
    # assuming that x is an array
    amp = []#[1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1,]#[0, .15, .15, .2, 0.15,0.15,0.15,0.15,0.07,0.07, .1,.1,.1,.1,.1]
    Bparams =  {'B20': B20, 'B40':B40,'B43': B43, 'B60': B60, 'B63':B63,'B66':B66}
    ionObj = cef.CFLevels.Bdict(ion,Bparams)
    # kBT = 10*kB
    # temperature = 10
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


# params from full fit
B20 = -0.03721092 # 3.3632e-04 (0.90%) (init = -0.03559)
B40 = -3.8796e-04 # 8.4839e-07 (0.22%) (init = -0.0003849)
B43 = -0.01406804 # 2.5391e-05 (0.18%) (init = -0.01393)
B60 =  3.1865e-06 # 5.3881e-09 (0.17%) (init = 3.154e-06)
B63 = -3.5930e-06 # 2.4015e-07 (6.68%) (init = -4.695e-06)
B66 =  3.4913e-05 # 3.1308e-07 (0.90%) (init = 3.3815e-05)
Jz =  -0.00263253 # 1.7886e-05 (0.68%) (init = -0.003)

waveArr= np.linspace(0,120, 480)
fieldArr = np.linspace(0,18, 100)

# make simulated data
arrC = zeemanSplitC(fieldArr, waveArr, B20, B40, B43, B60, B63, B66, Jz)

arrC = np.array(arrC)



plt.figure()
plt.contourf(fieldArr,waveArr,arrC.T, 100, cmap = 'Reds')

plt.title('simulated fitted CES H||c data \n B20: '+ str(B20)+' B40: '+str(B40)+' B43: ' +str(B43)+ '\n B60: ' +str(B60) + ' B63: ' + str(B63)+ ' B66: ' + str(B66) + 'Jz = -2.53ueV')
plt.xlabel('Field (T)')
plt.ylabel('Wavenumber (cm$^{-1}$)')
plt.colorbar()



############################################################
# now let's plot jsut the raman data (c-axis)
field = [float(b) for b in ramanData.columns.values]
wavenums = [float(i) for i in ramanData.index.values]

plt.figure()
plt.contourf(field, wavenums, ramanData,100, cmap='Reds')
plt.xlim(0,14)
plt.ylim(0,120)
# plt.clim(0, 1)
plt.colorbar()
plt.title('CsErSe2 H||c raman data')

############################################################
# lets put both on the same plot for easy compare
field = [float(b) for b in ramanData.columns.values]
wavenums = [float(i) for i in ramanData.index.values]
waveArr= np.linspace(0,120, 480)
fieldArr = np.linspace(0,18, 100)
fig = plt.figure()
gs = fig.add_gridspec(1,2, hspace=0, wspace=0)
axs = gs.subplots(sharex=True, sharey=True)
plt1 = axs[0].contourf(field, wavenums, ramanData,100, cmap='viridis')
axs[0].set_title('Raman Data H||c')
plt2 = axs[1].contourf(fieldArr,waveArr,arrC.T, 100, cmap = 'Reds')
axs[1].set_title('Simulated Data H||c')
axs[1].set_xlabel('Field (T)')
axs[0].set_ylabel('Energy (cm$^{-1}$)')
axs[0].set_xlabel('Field (T)')

# Hide x labels and tick labels for all but bottom plot.
for ax in axs:
    ax.label_outer()
plt.xlim(0,14)
plt.ylim(0,120)

fig.subplots_adjust(right=0.8)
cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.7])
fig.colorbar(plt2, cax=cbar_ax)

simulatedRaman = arrC


# now we plot data w/ overlay lines
# params from full fit
B20 = -0.03721092 # +/- 3.3632e-04 (0.90%) (init = -0.03559)
B40 = -3.8796e-04 # +/- 8.4839e-07 (0.22%) (init = -0.0003849)
B43 = -0.01406804 # +/- 2.5391e-05 (0.18%) (init = -0.01393)
B60 =  3.1865e-06 # +/- 5.3881e-09 (0.17%) (init = 3.154e-06)
B63 = -3.5930e-06 # +/- 2.4015e-07 (6.68%) (init = -4.695e-06)
B66 =  3.4913e-05 # +/- 3.1308e-07 (0.90%) (init = 3.3815e-05)
Jz =  -0.00263253 # +/- 1.7886e-05 (0.68%) (init = -0.003)



field = [float(b) for b in ramanData.columns.values]
wavenums = [float(i) for i in ramanData.index.values]
waveArr= wavenums
fieldArr = field 


ampC, arrC = zeemanSplitLinesC(fieldArr, B20, B40, B43, B60, B63, B66, Jz)
arrC = np.array(arrC)
arrC = arrC*meVToCm
arrC = arrC.T

ampC = np.array(ampC)
ampC = ampC.T

field = [float(b) for b in ramanData.columns.values]
wavenums = [float(i) for i in ramanData.index.values]

plt.figure()
plt.contourf(field, wavenums, ramanData,50)
plt.xlim(0,17.5)
plt.ylim(0,120)

plt.title('CsErSe2 H||c with overlayed  calclines\n B20: '+ str(B20)+' B40: '+str(B40)+' B43: ' +str(B43)+ '\n B60: ' +str(B60) + ' B63: ' + str(B63)+ ' B66: ' + str(B66))
for i in range(40):
    if i<16: 
        plt.plot(fieldArr, arrC_neg[i], 'r', alpha=0.7)
    if i>=16:
        plt.plot(fieldArr, arrC_neg[i], 'r--', alpha=0.7)


# IR data
field = [float(b) for b in dataC.columns.values]
wavenums = [float(i) for i in dataC.index.values]

plt.figure()
plt.contourf(field, wavenums, dataC,50)
plt.xlim(0,14)
plt.ylim(0,120)
# plt.clim(0, 1)
# plt.colorbar()
plt.title('CsErSe2 H||c with overlayed  calclines\n B20: '+ str(B20)+' B40: '+str(B40)+' B43: ' +str(B43)+ '\n B60: ' +str(B60) + ' B63: ' + str(B63)+ ' B66: ' + str(B66)+ 'JJz = ' + str(Jz))
for i in range(40):
    if i<16: 
        plt.plot(fieldArr, arrC_neg[i], 'r-', alpha=0.7) 
        plt.plot(fieldArr, arrC_zero[i], 'r--', alpha=0.7) 


# calc and save values no plots: 
# allen's Bparams
# B20 = -3.559e-2
# B40 = -3.849e-4
# B43 = -1.393e-2
# B60 =  3.154e-6
# B63 = -4.695e-6
# B66 =  3.3815e-5
# Jz = -2.4e-3

# my fit params
B20 =   -0.03721092
B40 =   -0.00038796
B43 =   -0.01406804
B60 =    3.1865e-06
B63 =   -3.593e-06
B66 =    3.4913e-05
Jperp = -.53070e-03 # +/- 2.6332e-06 (0.50%) (init = 0)
Jz =    -2.63253e-03
ampC_neg, arrC_neg = zeemanSplitLinesC(fieldArr, B20, B40, B43, B60, B63, B66, Jz)
arrC_neg = np.array(arrC_neg) * meVToCm
arrC_neg = arrC_neg.T
ampC_neg = np.array(ampC_neg).T

ampC_zero, arrC_zero = zeemanSplitLinesC(fieldArr, B20, B40, B43, B60, B63, B66, 0)
arrC_zero = np.array(arrC_zero) * meVToCm
arrC_zero = arrC_zero.T
ampC_zero = np.array(ampC_zero).T

# Save results to HDF5
with h5py.File('zeeman_split_lines_my_parmas_my_fit.h5', 'w') as hdf:
    hdf.create_dataset('fieldArr', data=fieldArr)
    hdf.create_dataset('arrC_neg', data=arrC_neg)
    hdf.create_dataset('ampC_neg', data=ampC_neg)
    hdf.create_dataset('arrC_zero', data=arrC_zero)
    hdf.create_dataset('ampC_zero', data=ampC_zero)
    hdf.create_dataset('dataC', data = dataC.to_numpy())
    hdf.create_dataset('dataField', data = field)
    hdf.create_dataset('dataWavenums', data = wavenums)
    hdf.attrs['B20'] = B20
    hdf.attrs['B40'] = B40
    hdf.attrs['B43'] = B43
    hdf.attrs['B60'] = B60
    hdf.attrs['B63'] = B63
    hdf.attrs['B66'] = B66
    hdf.attrs['Jz'] = Jz
    hdf.attrs['Jperp'] = Jperp

# let's do the same fit for the b-axis
# first import b data
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

# now build ab functions
def zeemanSplitLinesAB(field, B20, B40, B43, B60, B63, B66, Jperp):     
    # assuming only H||B rn
    # assuming that x is an array
    amp = []#[1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1,]#[0, .15, .15, .2, 0.15,0.15,0.15,0.15,0.07,0.07, .1,.1,.1,.1,.1]
    Bparams =  {'B20': B20, 'B40':B40,'B43': B43, 'B60': B60, 'B63':B63,'B66':B66}
    ionObj = cef.CFLevels.Bdict(ion,Bparams)
    # kBT = 10*kB
    temperature = 10
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

def diagonalizeAB(ionObj, ion, J, h, temperature): 
    h = newHAB(ionObj, h, J, temperature)
    JdotB = muB*(h*cef.Operator.Jy(ionObj.J))*cef.LandeGFactor(ion)
    H = np.sum([a*b for a,b in zip(ionObj.O, ionObj.B)], axis=0)
    ionObj.diagonalize(H + JdotB.O) # this is just H = Hcef + Hmag
    evals = ionObj.eigenvalues
    return evals

def bmag(mag, J,  h, temperature, ionObj):
    newh = h +q*J*mag/muB/1.2/1.2
    mag = -1*ionObj.magnetization(ion, temperature, [0,newh, 0]).T[1]-mag
    return mag

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
    temperature = 10
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

'''
## -- End pasted text --
## -- End pasted text --
[[Model]]
    Model(zeemanSplitAB)
[[Fit Statistics]]
    # fitting method   = leastsq
    # function evals   = 35
    # data points      = 116493
    # variables        = 1
    chi-square         = 1296.46434
    reduced chi-square = 0.01112921
    Akaike info crit   = -524005.690
    Bayesian info crit = -523996.025
    R-squared          = 0.24997391
[[Variables]]
    B20:   -0.03721092 (fixed)
    B40:   -0.00038796 (fixed)
    B43:   -0.01406804 (fixed)
    B60:    3.1865e-06 (fixed)
    B63:   -3.593e-06 (fixed)
    B66:    3.4913e-05 (fixed)
    Jperp: -5.3070e-04 +/- 2.6332e-06 (0.50%) (init = 0)

'''

# now let's calculate everything!!
B20 =   -0.03721092
B40 =   -0.00038796
B43 =   -0.01406804
B60 =    3.1865e-06
B63 =   -3.593e-06
B66 =    3.4913e-05
Jperp = -5.3070e-04 # +/- 2.6332e-06 (0.50%) (init = 0)

waveArr= np.linspace(0,120, 480)
fieldArr = np.linspace(0,18, 100)

ampB_neg, arrB_neg = zeemanSplitLinesAB(fieldArr, B20, B40, B43, B60, B63, B66, Jperp)
arrB_neg = np.array(arrB_neg) * meVToCm
arrB_neg = arrB_neg.T
ampB_neg = np.array(ampB_neg).T

ampB_zero, arrB_zero = zeemanSplitLinesAB(fieldArr, B20, B40, B43, B60, B63, B66, 0)
arrB_zero = np.array(arrB_zero) * meVToCm
arrB_zero = arrB_zero.T
ampB_zero = np.array(ampB_zero).T

field = [float(b) for b in dataB.columns.values]
wavenums = [float(i) for i in dataB.index.values]
plt.figure()
plt.contourf(field, wavenums, dataB,50)
plt.xlim(0,14)
plt.ylim(0,120)
# plt.clim(0, 1)
# plt.colorbar()
plt.xlabel('H [T]')
plt.ylabel('Energy [cm^-1]')
plt.title('CsErSe2 H||b with overlayed  calclines\n B20: '+ str(B20)+' B40: '+str(B40)+' B43: ' +str(B43)+ '\n B60: ' +str(B60) + ' B63: ' + str(B63)+ ' B66: ' + str(B66) + 'JJperp :' + str(Jperp) +'\n solid line = mft broken line = no mft')
for i in range(40):
    if i<16: 
        plt.plot(fieldArr, arrB_neg[i], 'r-', alpha=0.7)
        plt.plot(fieldArr, arrB_zero[i], 'r--', alpha=0.7)

# Save results to HDF5
with h5py.File('zeeman_split_lines_AB.h5', 'w') as hdf:
    hdf.create_dataset('fieldArr', data=fieldArr)
    hdf.create_dataset('arrB_neg', data=arrB_neg)
    hdf.create_dataset('ampB_neg', data=ampB_neg)
    hdf.create_dataset('arrB_zero', data=arrB_zero)
    hdf.create_dataset('ampB_zero', data=ampB_zero)
    hdf.create_dataset('dataB', data = dataB.values)
    hdf.create_dataset('dataField', data = field)
    hdf.create_dataset('dataWavenums', data = wavenums)
    hdf.attrs['B20'] = B20
    hdf.attrs['B40'] = B40
    hdf.attrs['B43'] = B43
    hdf.attrs['B60'] = B60
    hdf.attrs['B63'] = B63
    hdf.attrs['B66'] = B66

# show again with Allen's params
B20 = -3.559e-2
B40 = -3.849e-4
B43 = -1.393e-2
B60 = 3.154e-6
B63 = -4.695e-6
B66 = 3.3815e-5
JperpAllen = -0.2e-3
JzAllen = -2.4e-3

waveArr= np.linspace(0,120, 480)
fieldArr = np.linspace(0,18, 100)

ampB_neg, arrB_neg = zeemanSplitLinesAB(fieldArr, B20, B40, B43, B60, B63, B66, JperpAllen)
arrB_neg = np.array(arrB_neg) * meVToCm
arrB_neg = arrB_neg.T
ampB_neg = np.array(ampB_neg).T

ampB_zero, arrB_zero = zeemanSplitLinesAB(fieldArr, B20, B40, B43, B60, B63, B66, 0)
arrB_zero = np.array(arrB_zero) * meVToCm
arrB_zero = arrB_zero.T
ampB_zero = np.array(ampB_zero).T

field = [float(b) for b in dataB.columns.values]
wavenums = [float(i) for i in dataB.index.values]
plt.figure()
plt.contourf(field, wavenums, dataB,50)
plt.xlim(0,14)
plt.ylim(0,120)
# plt.clim(0, 1)
# plt.colorbar()
plt.xlabel('H [T]')
plt.ylabel('Energy [cm^-1]')
plt.title('CsErSe2 H||b with overlayed  calclines\n B20: '+ str(B20)+' B40: '+str(B40)+' B43: ' +str(B43)+ '\n B60: ' +str(B60) + ' B63: ' + str(B63)+ ' B66: ' + str(B66) + 'JJperp :' + str(JperpAllen) +'\n solid line = mft broken line = no mft')
for i in range(40):
    if i<16: 
        plt.plot(fieldArr, arrB_neg[i], 'r-', alpha=0.7)
        plt.plot(fieldArr, arrB_zero[i], 'r--', alpha=0.7)

# Save results to HDF5
with h5py.File('zeeman_split_lines_AB_allens_params.h5', 'w') as hdf:
    hdf.create_dataset('fieldArr', data=fieldArr)
    hdf.create_dataset('arrB_neg', data=arrB_neg)
    hdf.create_dataset('ampB_neg', data=ampB_neg)
    hdf.create_dataset('arrB_zero', data=arrB_zero)
    hdf.create_dataset('ampB_zero', data=ampB_zero)
    hdf.create_dataset('dataB', data = dataB.values)
    hdf.create_dataset('dataField', data = field)
    hdf.create_dataset('dataWavenums', data = wavenums)
    hdf.attrs['B20'] = B20
    hdf.attrs['B40'] = B40
    hdf.attrs['B43'] = B43
    hdf.attrs['B60'] = B60
    hdf.attrs['B63'] = B63
    hdf.attrs['B66'] = B66