import numpy as np
import matplotlib.pyplot as plt
import PyCrystalField as cef
import scipy
from scipy.optimize import minimize
from scipy.optimize import leastsq
import pandas as pd
import lmfit

plt.ion()

################################################################
# working B params -> from C-axis raman fit 11/13
B20 = -0.02926696 # init = -0.03559
B40 = -3.9097e-04 # init = -0.0003849
B43 = -0.01391860 # init = -0.01393
B60 =  3.0584e-06 # init = 3.154e-06
B63 = -4.2840e-06 # init = -4.695e-06
B66 =  3.3645e-05 # init = 3.3815e-05


################################################################

# first, lets do the waterfall raman plots

# load data

fname = '/Users/hopeless/Desktop/LeeLab/data/CsErSe2_raman.csv'
rawData = pd.read_csv(fname, index_col=0, skiprows=0, header=1, delimiter=',')
ramanData = rawData

# take the log - make sure to do this FIRST
ramanData = np.log(ramanData)
ramanData = ramanData-ramanData.min(axis=None)
# normalize 
ramanData = ramanData/ramanData.max(axis=None)

fields = ['0.001','1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14'] # just want to plot every tesla

n = len(fields)
colors = plt.cm.gnuplot(np.linspace(0,.8,n))

plt.figure()
i = 0
for f in fields: 
    spec = ramanData[f]+i/3
    plt.plot(ramanData.index, spec, label = f, color= colors[i])
    i+=1
plt.legend(draggable = True)
plt.xlim(0,100)
plt.title('Raman Data, 0-100cm')
plt.xlabel('wavenumber (cm^-1)')
plt.ylabel('Intensity (arb)')

# now plot the raman lines
plt.figure()
i = 0
for f in fields: 
    spec = ramanData[f]+i/3
    plt.plot(ramanData.index, spec, label = f, color= colors[i])
    i+=1
# plt.legend()
plt.xlim(100,200)
plt.title('Raman Data, phonon lines, 100-200 cm')
plt.xlabel('wavenumber (cm^-1)')
plt.ylabel('Intensity (arb)')

# now plot high energy lines
plt.figure()
i = 0
for f in fields: 
    spec = ramanData[f]+i/5
    plt.plot(ramanData.index, spec, label = f, color= colors[i])
    i+=1
# plt.legend()
plt.xlim(200,560)
plt.title('Raman Data, high energy, 200-560 cm, possible 13/2 lines')
plt.xlabel('wavenumber (cm^-1)')
plt.ylabel('Intensity (arb)')


################################################################
# now we make some fit plots
# start with the c-axis
# assume we've run the healper functions

# start with lines
field = [float(b) for b in ramanData.columns.values]
wavenums = [float(i) for i in ramanData.index.values]
waveArr= wavenums
fieldArr = field 


ampC, arrC = zeemanSplitLinesC(fieldArr, B20, B40, B43, B60, B63, B66)
arrC = np.array(arrC)
arrC = arrC*meVToCm
arrC = arrC.T

ampC = np.array(ampC)
ampC = ampC.T

plt.figure()
plt.contourf(field, wavenums, ramanData,50)
plt.xlim(0,14)
plt.ylim(0,120)
plt.clim(0, 1)
plt.colorbar()
plt.title('CsErSe2 H||c with overlayed  calclines\n B20: '+ str(B20)+' B40: '+str(B40)+' B43: ' +str(B43)+ '\n B60: ' +str(B60) + ' B63: ' + str(B63)+ ' B66: ' + str(B66))
for i in range(40):
    if i<16: 
        plt.plot(fieldArr, arrC[i], 'r', alpha=0.7)
    if i>=16: 
        plt.plot(fieldArr, arrC[i], 'r--', alpha=0.3)


# plot simulated data

waveArr= np.linspace(0,120, 480)
fieldArr = np.linspace(0,14, 100)
 # field, wavenum, B20, B40, B43, B60, B63, B66, amp1, amp2, amp3, amp4, amp5, amp6, amp7, amp8,width, phononCen, phononAmp, phononWid, amp9, amp10)
amp1 = 0.1
amp2 = 0.3
amp3 = 0.3
amp4 = 0.15
amp5 = 0.2
amp6 = 0.2
amp7 = 0.287
amp8 = 0.2
phononCen = 49.3
phononAmp = 0.499
phononWid = 0.95
amp9 = 0.136
amp10 = 0.097
width = 1.247
arrC = zeemanSplitC(fieldArr, waveArr, B20, B40, B43, B60, B63, B66,
                    amp1, amp2, amp3, amp4, amp5, amp6, amp7, amp8,width, 
                    phononCen, phononAmp, phononWid, amp9, amp10)

arrC = np.array(arrC)



plt.figure()
plt.contourf(fieldArr,waveArr,arrC.T, 100)
# plt.xlim(0,17.5)
# plt.ylim(20,100)
# plt.clim(.1,2.5)
plt.title('simulated fitted CES H||c data \n B20: '+ str(B20)+' B40: '+str(B40)+' B43: ' +str(B43)+ '\n B60: ' +str(B60) + ' B63: ' + str(B63)+ ' B66: ' + str(B66))
plt.xlabel('Field (T)')
plt.ylabel('Wavenumber (cm$^{-1}$)')
plt.colorbar()

############################################################
# now let's plot jsut the raman data (c-axis)
field = [float(b) for b in ramanData.columns.values]
wavenums = [float(i) for i in ramanData.index.values]

plt.figure()
plt.contourf(field, wavenums, ramanData,50)
plt.xlim(0,14)
plt.ylim(0,120)
# plt.clim(0, 1)
plt.colorbar()
plt.title('CsErSe2 H||c raman data')


############################################################
############################################################
############################################################

# okay so now we want to take a look at the H||b data
# first load

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

# lets drop all the data that's not relevant
dropidx = []
for index in dataB.index:
    if index<20: 
        dropidx.append(index)
    if index>100:
        dropidx.append(index)

dataB = dataB.drop(labels=dropidx , axis = 0)

dataB = dataB/dataB.max(axis=None)

####################################################################################
# now plot just the data

field = [float(b) for b in dataB.columns.values]
wavenums = [float(i) for i in dataB.index.values]

plt.figure()
plt.contourf(field,wavenums, dataB, 100)
plt.title('CsErSe$_2$ H||b IR absorption')
plt.xlabel('Field (T)')
plt.ylabel('Wavenumber (cm$^{-1}$)')
plt.colorbar()

####################################################################################

# now we want to plot with overlaid lines
# test B params

# B20 = -0.02926696 # init = -0.03559
# B40 = (-3.9097e-04 + -0.0003849)/2.
# B43 = (-0.01391860  + -0.01393)/2.
# B60 =  3.0584e-06  # 3.154e-06)/2.
# B63 = -4.2840e-06  # init -4.695e-06
# B66 =  (3.3645e-05  + 3.3815e-05)/2.

# B20 = -0.02965287 # init = -0.02926696)
# B40 = -0.000387935 # fixed)
# B43 = -0.0139243 # fixed)
# B60 =  3.0671e-06 # init = 3.0584e-06)
# B63 = -4.3545e-06 # init = -4.284e-06)

# B params from two var fit
B20 = -0.03265325 # init = -0.03559)
B40 = -0.0003849 # fixed)
B43 = -0.01393 # fixed)
B60 =  3.054e-06 # fixed)
B63 = -8.4011e-07 # init = -4.695e-06)
B66 =  3.3815e-05 # fixed)

wavenums= [float(i) for i in dataB.index.values]
field = [float(b) for b in dataB.columns.values]


ampAB, arrAB = zeemanSplitLinesAB(field, B20, B40, B43, B60, B63, B66)
arrAB = np.array(arrAB)
arrAB = arrAB*meVToCm
arrAB = arrAB.T

ampAB = np.array(ampAB)
ampAB = ampAB.T

plt.figure()
plt.contourf(field, wavenums, dataB,100)
# plt.xlim(0,14)
plt.ylim(0,120)
# plt.clim(0, 1)
plt.colorbar()
plt.title('CsErSe$_2$ H||b with overlaid calculated CEF lines \n B20: '+ str(B20)+' B40: '+str(B40)+' B43: ' +str(B43)+ '\n B60: ' +str(B60) + ' B63: ' + str(B63)+ ' B66: ' + str(B66))
for i in range(40):
    if i<16: 
        plt.plot(field, arrAB[i], 'r', alpha=0.7)
    if i>=16: 
        plt.plot(field, arrAB[i], 'r--', alpha=0.3)
