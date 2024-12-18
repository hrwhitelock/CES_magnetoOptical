# 2nd iteration of simul fits


# still just using explicit function that fits on b,c axis
# theres a better way to do this, but oh well


# first let's load the data in 2D, AB plane from IR and C axis from raman

bfname = '/Users/hopeless/Desktop/LeeLab/data/CsErSe2_data/Maglab_IRMeasurement_June2022/ProcessingCode/Load1_TrimData/P2_CsEr_100-FIR_RAWAVG.dat'
cfname = '/Users/hopeless/Desktop/LeeLab/data/CsErSe2_data/CsErSe2_raman.csv'

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

# lets drop all the data that's not relevant
dropidx = []
for index in dataB.index:
    if index<20: 
        dropidx.append(index)
    if index>100:
        dropidx.append(index)

dataB = dataB.drop(labels=dropidx , axis = 0)

# now we invert
# dataB = (dataB - dataB.max(axis = None))*-1

# drop everything below 0 
for idx in dataB.index: 
    for col in dataB.columns: 
        if dataB.loc[idx, col]<0: 
            dataB.loc[idx,col] = 0

dataB = dataB/dataB.max(axis=None)

temp = dataB
# dataB = np.log(dataB+1)
# now plot B axis to make sure we did it right
field = [float(b) for b in dataB.columns.values]
wavenums = [float(i) for i in dataB.index.values]

plt.figure()
plt.contourf(field,wavenums, dataB, 100)


## now we clean up raman data
ramanData = pd.read_csv(cfname, index_col=0, skiprows=0, header=1, delimiter=',')

# take the log - make sure to do this FIRST
ramanData = np.log(ramanData)
ramanData = ramanData-ramanData.min(axis=None)
# normalize 
ramanData = ramanData/ramanData.max(axis=None)

dropidx =[]
# here I'm cutting out all the phonons I can - the phonon spec is boring and idgaf
for idx in ramanData.index: 
    if idx>100: 
        dropidx.append(idx)

ramanData = ramanData.drop(labels = dropidx, axis = 0)

field = [float(b) for b in ramanData.columns.values]
wavenums = [float(i) for i in ramanData.index.values]

plt.figure()
plt.contourf(field, wavenums, ramanData, 100)
plt.colorbar()

# put them together -> arrays are jagged so we need to fix that first
# b axis index is much bigger -> easiest way to do is to pad our functions
padidxC = []
for idx in dataB.index: 
    if idx not in ramanData.index: 
        padidxC.append(idx)

padidxB = []
for idx in ramanData.index: 
    if idx not in dataB.index: 
        padidxB.append(idx)

padcolC = []
for col in dataB.columns: 
    if col not in ramanData.columns: 
        padcolC.append(col)

padcolB = []
for col in ramanData.columns: 
    if col not in dataB.columns: 
        padcolB.append(col)


padArrB = np.empty((len(padidxB), len(padcolB)))
padArrB[:] = np.nan
padArrC = np.empty((len(padidxC), len(padcolC)))
padArrC[:] = np.nan

paddfB = pd.DataFrame(data = padArrB, index = padidxB, columns = padcolB)
paddfC = pd.DataFrame(data = padArrC, index = padidxC, columns = padcolC)

cres = pd.concat([ramanData, paddfC], axis=1, join='outer')
bres = pd.concat([dataB, paddfB], axis=1, join='outer')

########
#plot to see if we did it right

field = [float(b) for b in cres.columns.values]
wavenums = [float(i) for i in cres.index.values]

plt.figure()
plt.contourf(field, wavenums, cres, 100)
plt.colorbar()

# put them together

numpyDataC = np.array(cres.to_numpy()).tolist()
numpyDataB = np.array(cres.to_numpy()).tolist()

data = [numpyDataB, numpyDataC]


# now we try fitting and maybe cry, maybe be fine we'll see
field = [float(b) for b in bres.columns.values]
wavenums = [float(i) for i in bres.index.values]


B20 = -3.559e-2
B40 = -3.849e-4
B43 = -1.393e-2
B60 = 3.154e-6
B63 = -4.695e-6
B66 = 3.3815e-5



model = lmfit.Model(zeemanFitModel, independent_vars=['field', 'wavenum'])
params = model.make_params()

params['B20'].set(value= B20, min = -.06, max = 0.06)
params['B40'].set(value= B40, min = -.06, max = 0.06)
params['B43'].set(value= B43, min = -.06, max = 0.06)
params['B60'].set(value= B60, min = -.06, max = 0.06)
params['B63'].set(value= B63, min = -.06, max = 0.06)
params['B66'].set(value= B66, min = -.06, max = 0.06)
params['amp1'].set(value = 0.5, vary = False) # these don't vary because they are outside the bounds of our data
params['amp2'].set(value = 0.3, vary = False)
params['amp3'].set(value = 0.7, min = 0.2, max = 1)
params['amp4'].set(value = 0.3, min = 0.15, max = 1)
params['amp5'].set(value = 0.5, min = 0.2, max = 1)
params['amp6'].set(value = 0.5, min = 0.2, max = 1)
params['amp7'].set(value = 0.5, min = 0.2, max = 1)
params['amp8'].set(value = 0.5, min = 0.2, max = 1)
params['amp9'].set(value = 0.5, min = 0.2, max = 1)
params['amp10'].set(value = 0.5, min = 0.2, max = 1)
params['wid'].set(value = 0.5, min = 0.0, max = 2.0)

result = model.fit(data, field=field, wavenum=wavenums, params =params)
print(result.fit_report())
