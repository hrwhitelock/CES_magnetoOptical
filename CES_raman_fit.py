
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

fname = '/Users/hopeless/Desktop/LeeLab/data/CsErSe2_data/CsErSe2_raman.csv'
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

fakeEnergies = np.linspace(300, 5000000, 100000) # should be enough to enforce hard boundary :(
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

#model 1 from Allen's paper
B20 = -3.559e-2
B40 = -3.849e-4
B43 = -1.393e-2
B60 = 3.154e-6
B63 = -4.695e-6
B66 = 3.3815e-5


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
params['amp1'].set(value = 0.5, vary = False)
params['amp2'].set(value = 0.3, vary = False)
params['amp3'].set(value = 0.7, min = 0.2, max = 1)
params['amp4'].set(value = 0.3, min = 0.15, max = 1)
params['amp5'].set(value = 0.5, min = 0.2, max = 1)
params['amp6'].set(value = 0.5, min = 0.2, max = 1)
params['amp7'].set(value = 0.5, min = 0.2, max = 1)
params['amp8'].set(value = 0.5, min = 0.2, max = 1)
params['amp9'].set(value = 0.5, min = 0.05, max = 1)
params['amp10'].set(value = 0.5, min = 0.05, max = 1)
params['width'].set(value = 0.5, min = 0.0, max = 2.0)
params['phononCen'].set(value = 49.2)
params['phononAmp'].set(value = 0.5)
params['phononWid'].set(value = 1)
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
result 11/11 -> trying maybe one more thing
[[Model]]
    Model(zeemanSplitC)
[[Fit Statistics]]
    # fitting method   = leastsq
    # function evals   = 19613
    # data points      = 15265
    # variables        = 15
    chi-square         = 230.067002
    reduced chi-square = 0.01508636
    Akaike info crit   = -64005.8709
    Bayesian info crit = -63891.3712
    R-squared          = 0.13034021
##  Warning: uncertainties could not be estimated:
    amp2:   at boundary
    amp4:   at boundary
    amp9:   at boundary
    amp10:  at boundary
[[Variables]]
[[Model]]
    Model(zeemanSplitC)
[[Fit Statistics]]
    # fitting method   = leastsq
    # function evals   = 38000
    # data points      = 15265
    # variables        = 18
    chi-square         = 109.449089
    reduced chi-square = 0.00717840
    Akaike info crit   = -75340.4086
    Bayesian info crit = -75203.0089
    R-squared          = 0.58627934
##  Warning: uncertainties could not be estimated:
    amp4:       at boundary
    amp5:       at boundary
    amp6:       at boundary
    amp8:       at boundary
[[Variables]]
    B20:       -0.02926696 (init = -0.03559)
    B40:       -3.9097e-04 (init = -0.0003849)
    B43:       -0.01391860 (init = -0.01393)
    B60:        3.0584e-06 (init = 3.154e-06)
    B63:       -4.2840e-06 (init = -4.695e-06)
    B66:        3.3645e-05 (init = 3.3815e-05)
    amp1:       0.5 (fixed)
    amp2:       0.3 (fixed)
    amp3:       0.30406726 (init = 0.7)
    amp4:       0.15000000 (init = 0.3)
    amp5:       0.20000103 (init = 0.5)
    amp6:       0.20000000 (init = 0.5)
    amp7:       0.28652781 (init = 0.5)
    amp8:       0.20000000 (init = 0.5)
    width:      1.24715342 (init = 0.5)
    phononCen:  49.3020750 (init = 49.2)
    phononAmp:  0.49902960 (init = 0.5)
    phononWid:  0.95028048 (init = 1)
    amp9:       0.13605417 (init = 0.5)
    amp10:      0.09724615 (init = 0.5)
'''

# lets plot this bitch
field = [float(b) for b in ramanData.columns.values]
wavenums = [float(i) for i in ramanData.index.values]
waveArr= wavenums #np.linspace(0,500, 1000)
fieldArr = field #np.linspace(0,14, 100)

# Allen's params
# B20 = -3.559e-2
# B40 = -3.849e-4
# B43 = -1.393e-2
# B60 = 3.154e-6
# B63 = -4.695e-6
# B66 = 3.3815e-5

# my params
# B20 = -0.02926696 # -0.03559)
# B40 = -3.9097e-04 # -0.0003849)
# B43 = -0.01391860 # -0.01393)
# B60 =  3.0584e-06 # 3.154e-06)
# B63 = -4.2840e-06 # -4.695e-06)
# B66 =  3.3645e-05 # 3.3815e-05)

# my params averaged
# B20=-0.02926696
# B40=(-3.9097e-04-3.849e-4)/2
# B43=(-0.01391860-1.393e-2)/2
# B60= 3.0584e-06
# B63=-4.2840e-06
# B66= (3.3645e-05+3.3815e-5)/2

# Allens params averaged
# B20=-3.559e-2
# B40=(-3.9097e-04-3.849e-4)/2
# B43=(-0.01391860-1.393e-2)/2
# B60= (3.0584e-06 +3.154e-6)/2
# B63=-4.695e-6
# B66= (3.3645e-05+3.3815e-5)/2

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
    if i<16: 
        plt.plot(fieldArr, arrC[i], 'r', alpha=0.4)
    # if i>=16: 
    #     plt.plot(fieldArr, arrC[i], 'r--', alpha=0.3)

waveArr= np.linspace(0,500, 1000)
fieldArr = np.linspace(0,14, 100)

arrC = zeemanSplitC(fieldArr, waveArr, result.params['B20'], result.params['B40'],result.params['B43']
                    ,result.params['B60'],result.params['B63'],result.params['B66'],
                    result.params['amp1'],result.params['amp2'],result.params['amp3'],result.params['amp4'],result.params['amp5'],
                    result.params['amp6'],result.params['amp7'],result.params['amp8'], 0.5, 49.2, 1,1,.1,.1)#result.params['amp9'],result.params['amp10'], 


# arr = zeemanSplit(fieldArr, waveArr,result.params['B20'], result.params['B40'], result.params['B43'], result.params['B60'], result.params['B63'], result.params['B66']) 
# arrC = zeemanSplitC(fieldArr, waveArr, B20, B40, B43, B60, B63, B66, 1,1,1,1,1,1,1,1,1,1)
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


####################################################################################################################
####################################################################################################################
####################################################################################################################
####################################################################################################################
# okay, so now I'm redoing the fit, but only letting 3 B params actually do the fit: B20, B63, B60

B20 =  -0.03559
B40 =  -0.0003849 
B43 =  -0.01393
B60 =   3.035e-06
B63 =  -4.695e-06
B66 =   3.3815e-05


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
params['B40'].set(value= B40, vary = False)
params['B43'].set(value= B43, vary = False)
params['B60'].set(value= B60, vary = False)
params['B63'].set(value= B63, min = -1e5, max = 1e5)
params['B66'].set(value= B66, vary = False)
params['amp1'].set(value = 0.5, vary = False)
params['amp2'].set(value = 0.3, vary = False)
params['amp3'].set(value = 0.7, vary = False)
params['amp4'].set(value = 0.4, vary = False)
params['amp5'].set(value = 0.4, vary = False)
params['amp6'].set(value = 0.4, vary = False)
params['amp7'].set(value = 0.4, vary = False)
params['amp8'].set(value = 0.4, vary = False)
params['amp9'].set(value = 0.2, vary = False)
params['amp10'].set(value = 0.05, vary = False)
params['width'].set(value = 1, vary=False)
# params['phononCen'].set(value = 49.2, vary = False)
# params['phononAmp'].set(value = 0.6, vary = False)
# params['phononWid'].set(value = 1, vary = False)
# params['amp11'].set(value = 0.5, min = 0.0, max = 1)
# params['amp12'].set(value = 0.5, min = 0.0, max = 1)
# params['amp13'].set(value = 0.5, min = 0.0, max = 1)
# params['amp14'].set(value = 0.5, min = 0.0, max = 1)
# params['amp15'].set(value = 0.5, min = 0.0, max = 1)
# params['amp16'].set(value = 0.5, min = 0.0, max = 1)

z = np.array(fitData.to_numpy()) # gotta do it twice with tuples :((((
z = z.T

result = model.fit(z, field=field, wavenum=wavenums, params =params, method ='ampgo')

print(result.fit_report())


'''
[[Model]]
    Model(zeemanSplitC)
[[Fit Statistics]]
    # fitting method   = leastsq
    # function evals   = 250
    # data points      = 15265
    # variables        = 2
    chi-square         = 263.383418
    reduced chi-square = 0.01725633
    Akaike info crit   = -61967.4285
    Bayesian info crit = -61952.1618
    R-squared          = 0.00440321
##  Warning: uncertainties could not be estimated:
[[Variables]]
    B20:       -0.03265325 (init = -0.03559)
    B40:       -0.0003849 (fixed)
    B43:       -0.01393 (fixed)
    B60:        3.054e-06 (fixed)
    B63:       -8.4011e-07 (init = -4.695e-06)
    B66:        3.3815e-05 (fixed)
    amp1:       0.5 (fixed)
    amp2:       0.3 (fixed)
    amp3:       0.7 (fixed)
    amp4:       0.4 (fixed)
    amp5:       0.4 (fixed)
    amp6:       0.4 (fixed)
    amp7:       0.4 (fixed)
    amp8:       0.4 (fixed)
    width:      1 (fixed)
    phononCen:  49.2 (fixed)
    phononAmp:  0.6 (fixed)
    phononWid:  1 (fixed)
    amp9:       0.2 (fixed)
    amp10:      0.05 (fixed)
'''