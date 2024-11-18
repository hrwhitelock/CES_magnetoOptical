
import numpy as np
import matplotlib.pyplot as plt
import PyCrystalField as cef
import scipy
from scipy.optimize import minimize
from scipy.optimize import leastsq
import pandas as pd
import lmfit

plt.ion()

#model 2 form Allens paper

B20 = -3.559e-2
B40 = -3.849e-4
B43 = -1.393e-2
B60 = 3.154e-6
B63 = 4.695e-6
B66 = 3.3815e-5

# B20 = -2.773e-2
# B40 = -3.987e-4
# B43 = -1.416e-2
# B60 = 3.152e-6
# B63 = -7.616e-6
# B66 = 3.275e-5
''' AB plane fit to start with
[[Model]]
    Model(zeemanSplit)
[[Fit Statistics]]
    # fitting method   = leastsq
    # function evals   = 1499
    # data points      = 116493
    # variables        = 6
    chi-square         = 7136.64640
    reduced chi-square = 0.06126560
    Akaike info crit   = -325304.978
    Bayesian info crit = -325246.985
    R-squared          = -5.46889547
[[Variables]]
    B20: -0.04040883 +/- 3.8025e-05 (0.09%) (init = -0.03559)
    B40: -3.8692e-04 +/- 9.3965e-08 (0.02%) (init = -0.0003849)
    B43: -0.01433959 +/- 3.0376e-06 (0.02%) (init = -0.01393)
    B60:  3.2030e-06 +/- 6.7337e-10 (0.02%) (init = 3.154e-06)
    B63: -1.9732e-06 +/- 3.0535e-08 (1.55%) (init = -4.695e-06)
    B66:  3.7301e-05 +/- 3.4077e-08 (0.09%) (init = 3.3815e-05)
[[Correlations]] (unreported correlations are < 0.100)
    C(B20, B63) = -0.9648
    C(B43, B66) = -0.9576
    C(B20, B43) = +0.9069
    C(B20, B66) = -0.8619
    C(B43, B63) = -0.8107
    C(B63, B66) = +0.7741
    C(B40, B63) = +0.4559
    C(B20, B40) = -0.4012
    C(B20, B60) = -0.3850
    C(B43, B60) = -0.3782
    C(B60, B63) = +0.3693
    C(B40, B66) = +0.2513
    C(B40, B60) = -0.2425
    C(B60, B66) = +0.2009
    C(B40, B43) = -0.1680
'''


temperature = 10 # in K

muB = 5.7883818012e-2  # meV/T
mu0 = np.pi*4e-7       # T*m/A
kB  = 8.617e-2         # [meV/K];
meVToCm =meVTocCmInv= 8.066 
ion = 'Er3+'
kBT =kB*temperature

# let's import the data and do a ZF subtraction
# this ZF subtraction is important because I want to remove all **unchanged** raman line
# there's one field dependent raman line and that's really cool!!!!

fname = '/Users/hopeless/Desktop/LeeLab/data/CsErSe2_raman.csv'
rawData = pd.read_csv(fname, index_col=0, skiprows=0, header=1, delimiter=',')
ramanData = rawData

# take the log - make sure to do this FIRST
ramanData = np.log(ramanData)
ramanData = ramanData-ramanData.min(axis=None)
# normalize 
ramanData = ramanData/ramanData.max(axis=None)


fitData =ramanData


dropidx =[]
# here I'm cutting out all the phonons I can - the phonon spec is boring and idgaf
for idx in fitData.index: 
    if idx>100: 
        dropidx.append(idx)
    # elif idx<190 and idx>100: 
    #     for col in fitData.columns:
    #         fitData.loc[idx,col] = fitData.min(axis=None)
    # elif idx>240 and idx<260: 
    #     dropidx.append(idx)
    # elif idx>285 and idx<330: 
    #     dropidx.append(idx)


# dropcol = []
# for col in fitData.columns: 
#     if float(col) < 3: 
#         dropcol.append(col)

fitData = fitData.drop(labels = dropidx, axis = 0)

# fitData = fitData-fitData.min(axis=None)
# fitData = fitData/fitData.max()

# quickly plot to make sure we're doing this well
field = [float(b) for b in ramanData.columns.values]
wavenums = [float(i) for i in ramanData.index.values]

plt.figure()
plt.contourf(field, wavenums, ramanData, 100)
plt.colorbar()

# plot fitData
field = [float(b) for b in fitData.columns.values]
wavenums = [float(i) for i in fitData.index.values]

plt.figure()
plt.contourf(field, wavenums, fitData, 100)
plt.colorbar()

# let's try adding fake high energy data, say 0 above 285 ev

fakeEnergies = np.linspace(348, 500000, 200000) # should be enough to enforce hard boundary :(
fakeField = ramanData.columns.values

fakeSpec = np.zeros((len(fakeEnergies), len(fakeField)))

fakeData = pd.DataFrame(fakeSpec, index = fakeEnergies, columns = fakeField)

tempFitData = np.array(fitData.to_numpy())


tempDF = pd.concat([fitData, fakeData], axis =0)
fitData = tempDF

# # remove our 0 field column, since we've subtracted it out
# ramanData = ramanData.drop(labels='0.001', axis=1)
# field = [float(b) for b in ramanData.columns.values]
# wavenums = [float(i) for i in ramanData.index.values]

# okay, so yee haw, now we get to actually fit this data, and we have 
# a high energy line to fit to which is super, duper exciting

# I've defined functions in another file for cleanliness

# # from AB plane fit
B20= -0.04040883
B40= -3.8692e-04
B43= -0.01433959
B60=  3.2030e-06
B63= -1.9732e-06
B66=  3.7301e-05

field = [float(b) for b in fitData.columns.values]
wavenums = [float(i) for i in fitData.index.values]

plt.figure()
plt.contourf(field, wavenums, fitData, 100)
plt.colorbar()

# B20= 0
# B40= 0
# B43= 0
# B60= 0
# B63= 0
# B66= 0


model = lmfit.Model(zeemanSplitC, independent_vars=['field', 'wavenum'])
params = model.make_params()

params['B20'].set(value= B20, min = -.06, max = 0.06)
params['B40'].set(value= B40, min = -.06, max = 0.06)
params['B43'].set(value= B43, min = -.06, max = 0.06)
params['B60'].set(value= B60, min = -.06, max = 0.06)
params['B63'].set(value= B63, min = -.06, max = 0.06)
params['B66'].set(value= B66, min = -.06, max = 0.06)
# params['amp1'].set(value = 0.5, vary = False)
# params['amp2'].set(value = 0.5, min = 0.0, max = 1)
# params['amp3'].set(value = 1, min = 0.0, max = 1)
# params['amp4'].set(value = 0.3, min = 0.0, max = 1)
# params['amp5'].set(value = 0.5, min = 0.0, max = 1)
# params['amp6'].set(value = 0.5, min = 0.0, max = 1)
# params['amp7'].set(value = 0.5, min = 0.0, max = 1)
# params['amp8'].set(value = 0.5, min = 0.0, max = 1)
# params['amp9'].set(value = 0.5, min = 0.0, max = 1)
# params['amp10'].set(value = 0.5, min = 0.0, max = 1)
# params['amp11'].set(value = 0.5, min = 0.0, max = 1)
# params['amp12'].set(value = 0.5, min = 0.0, max = 1)
# params['amp13'].set(value = 0.5, min = 0.0, max = 1)
# params['amp14'].set(value = 0.5, min = 0.0, max = 1)
# params['amp15'].set(value = 0.5, min = 0.0, max = 1)
# params['amp16'].set(value = 0.5, min = 0.0, max = 1)

z = np.array(fitData.to_numpy()) # gotta do it twice with tuples :((((
z = z.T

result = model.fit(z, field=field, wavenum=wavenums, params =params)

print(result.fit_report())

'''
[[Model]]
    Model(zeemanSplitC)
[[Fit Statistics]]
    # fitting method   = leastsq
    # function evals   = 29833
    # data points      = 46434
    # variables        = 21
    chi-square         = 689.568990
    reduced chi-square = 0.01485724
    Akaike info crit   = -195432.161
    Bayesian info crit = -195248.499
    R-squared          = -4.99544513
[[Variables]]
    B20:   -0.03553308 +/- 8.5004e-05 (0.24%) (init = -0.04040883)
    B40:   -3.7973e-04 +/- 4.2073e-07 (0.11%) (init = -0.00038692)
    B43:   -0.01424079 +/- 2.2208e-06 (0.02%) (init = -0.01433959)
    B60:    3.0728e-06 +/- 2.9139e-09 (0.09%) (init = 3.203e-06)
    B63:   -2.7710e-07 +/- 7.4939e-08 (27.04%) (init = -1.9732e-06)
    B66:    3.8956e-05 +/- 6.0802e-08 (0.16%) (init = 3.7301e-05)
    amp1:   0.5 (fixed)
    amp2:   0.52534606 +/- 0.07017554 (13.36%) (init = 0.5)
    amp3:   0.36320159 +/- 0.01342228 (3.70%) (init = 1)
    amp4:   0.22092630 +/- 0.00779686 (3.53%) (init = 0.3)
    amp5:   0.20981666 +/- 0.00776056 (3.70%) (init = 0.5)
    amp6:   0.26793413 +/- 0.00769700 (2.87%) (init = 0.5)
    amp7:   0.35832102 +/- 0.00784340 (2.19%) (init = 0.5)
    amp8:   0.22404343 +/- 0.00784676 (3.50%) (init = 0.5)
    amp9:   0.18615851 +/- 0.00778119 (4.18%) (init = 0.5)
    amp10:  0.16854994 +/- 0.00776591 (4.61%) (init = 0.5)
    amp11:  0.19001187 +/- 0.00768504 (4.04%) (init = 0.5)
    amp12:  0.18901404 +/- 0.00768771 (4.07%) (init = 0.5)
    amp13:  0.20604854 +/- 0.00758656 (3.68%) (init = 0.5)
    amp14:  0.17930535 +/- 0.00865938 (4.83%) (init = 0.5)
    amp15:  0.12373527 +/- 0.00997736 (8.06%) (init = 0.5)
    amp16:  0.14891783 +/- 0.00990352 (6.65%) (init = 0.5)
[[Correlations]] (unreported correlations are < 0.100)
    C(B20, B63)     = -0.9217
    C(B20, B40)     = -0.7998
    C(B40, B66)     = +0.7484
    C(B60, B66)     = -0.6852
    C(B40, B60)     = -0.6753
    C(B40, B63)     = +0.6548
    C(B20, B66)     = -0.5231
    C(amp15, amp16) = -0.4936
    C(B63, B66)     = +0.3857
    C(B20, B43)     = -0.3101
    C(B43, B63)     = +0.2889
    C(B20, B60)     = +0.2720
    C(B60, B63)     = -0.2537
    C(B43, B60)     = +0.2498
    C(amp7, amp8)   = -0.1992
    C(amp11, amp12) = -0.1678
    C(amp9, amp10)  = -0.1404
    C(amp4, amp5)   = -0.1190
    C(amp14, amp15) = -0.1086
    C(B20, amp2)    = -0.1069
    C(B40, amp2)    = +0.1058
    C(B43, amp16)   = +0.1023
    C(amp5, amp6)   = -0.1014
'''

# lets plot this bitch
field = [float(b) for b in ramanData.columns.values]
wavenums = [float(i) for i in ramanData.index.values]
waveArr= np.linspace(0,500, 1000)
fieldArr = np.linspace(0,14, 100)

ampC, arrC = zeemanSplitLinesC(fieldArr,result.params['B20'], result.params['B40'], 
                   result.params['B43'], result.params['B60'], 
                   result.params['B63'], result.params['B66'])

# ampC, arrC = zeemanSplitLinesC(fieldArr, .01, .001, .0001, .00001, .001, .0001)

# ampC, arrC = zeemanSplitLinesC(fieldArr, B20, B40, B43, B60, B63, B66)
arrC = np.array(arrC)
arrC = arrC*meVToCm
arrC = arrC.T

ampC = np.array(ampC)
ampC = ampC.T

plt.figure()
plt.contourf(field, wavenums, ramanData,50)
# plt.xlim(0,17.5)
# plt.ylim(0,110)
# plt.clim(0, .7)
plt.colorbar()
plt.title('CsErSe2 H||c with overlayed  calclines')
for i in range(40):
    if i<10: 
        plt.plot(fieldArr, arrC[i], 'r', alpha=0.4)
    # if i>=16: 
    #     plt.plot(fieldArr, arrC[i], 'r--', alpha=0.3)

waveArr= np.linspace(0,500, 1000)
fieldArr = np.linspace(0,14, 100)

arrC = zeemanSplitC(fieldArr, waveArr, result.params['B20'], result.params['B40'],result.params['B43']
                    ,result.params['B60'],result.params['B63'],result.params['B66']) 
#                     # result.params['amp1'],result.params['amp2'],result.params['amp3'],result.params['amp4'],result.params['amp5'],
#                     # result.params['amp6'],result.params['amp7'],result.params['amp8'],result.params['amp9'],result.params['amp10'],
#                     # result.params['amp11'],result.params['amp12'],result.params['amp13'],result.params['amp14'],result.params['amp15'],
#                     # result.params['amp16'] ) 
# arr = zeemanSplit(fieldArr, waveArr,result.params['B20'], result.params['B40'], result.params['B43'], result.params['B60'], result.params['B63'], result.params['B66']) 
# arrC = zeemanSplitC(fieldArr, waveArr, B20, B40, B43, B60, B63, B66)
arrC = np.array(arrC)



plt.figure()
plt.contourf(fieldArr,waveArr,arrC.T, 100)
# plt.xlim(0,17.5)
# plt.ylim(20,100)
# plt.clim(.1,2.5)
plt.title('simulated fitted CES H||c data')
plt.xlabel('appliead field')
plt.ylabel('wavenumber')
plt.colorbar()

''' ampgo notes
first pass after all inputs 0: 
[[Model]]
    Model(zeemanSplitC)
[[Fit Statistics]]
    # fitting method   = ampgo, with L-BFGS-B as local solver
    # function evals   = 12474
    # data points      = 46434
    # variables        = 6
    chi-square         = 1631.86736
    reduced chi-square = 0.03514835
    Akaike info crit   = -155463.286
    Bayesian info crit = -155410.811
    R-squared          = -0.91010329
##  Warning: uncertainties could not be estimated:
[[Variables]]
    B20:  0.05858630 (init = 0)
    B40:  0.00879018 (init = 0)
    B43: -0.09554336 (init = 0)
    B60:  0.00158043 (init = 0)
    B63: -0.00567974 (init = 0)
    B66: -0.01362483 (init = 0)

now i refit with AB plane params and see what happens

'''