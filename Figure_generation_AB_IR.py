import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import PyCrystalField as cef
import scipy
from scipy.optimize import minimize
from scipy.optimize import leastsq
import pandas as pd
import lmfit
from matplotlib import font_manager
font_manager.fontManager.addfont('/Users/hopeless/Library/Fonts/cmunrm.ttf')
plt.rcParams['mathtext.fontset'] = 'stix'
plt.rcParams['font.family'] = 'CMU Serif'
plt.rcParams['axes.unicode_minus'] = False
plt.rcParams.update({'font.size': 16})

plt.ion()


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


dataB = dataB/dataB.max(axis=None)



################################################################
# now we make some fit plots

# assume we've run the healper functions

# start with lines
field = [float(b) for b in dataB.columns.values]
wavenums = [float(i) for i in dataB.index.values]
waveArr= wavenums
fieldArr = np.linspace(0,20, 100)

B20 = -0.03265325 # init = -0.03559)
B40 = -0.0003849 # fixed)
B43 = -0.01393 # fixed)
B60 =  3.035e-6 #3.054e-06 # fixed) # this B60 param determined from field induced transition
B63 = -8.4011e-07 # init = -4.695e-06)
B66 =  3.3815e-05 # fixed)


ampAB, arrAB = zeemanSplitLinesAB(fieldArr, B20, B40, B43, B60, B63, B66)
arrAB = np.array(arrAB)
arrAB = arrAB*meVToCm
arrAB = arrAB.T

ampAB = np.array(ampAB)
ampAB = ampAB.T

plt.figure()
plt.contourf(field, wavenums, dataB,100, cmap = 'jet')
plt.xlim(0,14)
plt.ylim(0,120)
plt.clim(0, 1)
plt.xlabel('Field (T)')
plt.ylabel('Energy (cm$^{-1}$)')
plt.colorbar()
plt.title('CsErSe2 H||B with overlayed  calclines\n B20: '+ str(B20)+' B40: '+str(B40)+' B43: ' +str(B43)+ '\n B60: ' +str(B60) + ' B63: ' + str(B63)+ ' B66: ' + str(B66))
for i in range(40):
    if i<16: 
        plt.plot(fieldArr, arrAB[i], 'r', alpha=1, linewidth= 1)
    if i>=16: 
        plt.plot(fieldArr, arrAB[i], 'r--', alpha=0.7, linewidth = 0.5)


##################################################################################################
# lets try this plotting lines over a waterfall plot
fields = ['0.25', '0.999', '2', '3', '3.999', '5', '5.999', '7','8','9','10', '11',  '11.999', 
       '13', '14',  '14.999',  '16', '17',]
n = len(fields)
colors = plt.cm.inferno(np.linspace(0,.8,n))

fig,ax1 =  plt.subplots()
i = 0
for f in fields: 
    spec = dataB[f]+i/1.05
    ax1.plot(dataB.index, spec, label = f, color= colors[i])
    if i == 0: 
        ax1.annotate('H = 0T', xy = (90, i +.15), fontsize = 9)
    elif i ==7 or i==14:
        ax1.annotate('H = ' + f + 'T', xy = (90, i +.15), fontsize = 9)
    i+=1
plt.ylim(0,18)
ax2 = ax1.twinx() 
for i in range(40):
    if i<16: 
        ax2.plot(arrAB[i], fieldArr, 'r', alpha=1, linewidth= 1)
    if i>=16: 
        ax2.plot(arrAB[i], fieldArr,'r--', alpha=0.3)
# plt.legend(draggable = True)
ax2.set_ylabel('Field [T]')
plt.xlim(0,100)
plt.title('IR Data, 0-100cm')
ax1.set_xlabel('Energy (cm$^{-1}$)')
ax1.set_ylabel('Intensity (arb)')
plt.ylim(0,18)

################################################
# plot line from Allen's params
B20 =  -0.03559
B40 =  -0.0003849 
B43 =  -0.01393
B60 =   3.054e-06
B63 =  -4.695e-06
B66 =   3.3815e-05

field = [float(b) for b in dataB.columns.values]
wavenums = [float(i) for i in dataB.index.values]
waveArr= wavenums
fieldArr = field 


ampAB, arrAB = zeemanSplitLinesAB(fieldArr, B20, B40, B43, B60, B63, B66)
arrAB = np.array(arrAB)
arrAB = arrAB*meVToCm
arrAB = arrAB.T

ampAB = np.array(ampAB)
ampAB = ampAB.T

plt.figure()
plt.contourf(field, wavenums, dataB,50)
plt.xlim(0,14)
plt.ylim(0,120)
plt.clim(0, 1)
plt.colorbar()
plt.title('CsErSe2 H||AB with overlayed  calclines\n B20: '+ str(B20)+' B40: '+str(B40)+' B43: ' +str(B43)+ '\n B60: ' +str(B60) + ' B63: ' + str(B63)+ ' B66: ' + str(B66))
for i in range(40):
    if i<16: 
        plt.plot(fieldArr, arrAB[i], 'r', alpha=0.7) 

################################################
# plot simulated data
B20 = -0.03265325 # init = -0.03559)
B40 = -0.0003849 # fixed)
B43 = -0.01393 # fixed)
B60 =  3.035e-6 #3.054e-06 # fixed) # this B60 param determined from field induced transition
B63 = -8.4011e-07 # init = -4.695e-06)
B66 =  3.3815e-05 # fixed)

waveArr= np.linspace(0,120, 480)
fieldArr = np.linspace(0,18, 100)
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
width = 0.9
#field, wavenum, B20, B40, B43, B60, B63, B66, amp1, amp2, amp3, amp4, amp5, amp6, amp7, amp8,amp9, amp10,width
arrAB = zeemanSplitAB(fieldArr, waveArr, B20, B40, B43, B60, B63, B66,
                    amp1, amp2, amp3, amp4, amp5, amp6, amp7, amp8,amp9, amp10, 
                    width)

arrAB = np.array(arrAB)



plt.figure()
plt.contourf(fieldArr,waveArr,arrAB.T, 100, cmap = 'jet')
# plt.xlim(0,17.5)
# plt.ylim(20,100)
# plt.clim(.1,2.5)
plt.title('simulated fitted CES H||AB data \n B20: '+ str(B20)+' B40: '+str(B40)+' B43: ' +str(B43)+ '\n B60: ' +str(B60) + ' B63: ' + str(B63)+ ' B66: ' + str(B66))
plt.xlabel('Field (T)')
plt.ylabel('Wavenumber (cm$^{-1}$)')
plt.colorbar()

############################################################
# now let's plot jsut the IR data (H||b)
field = [float(b) for b in dataB.columns.values]
wavenums = [float(i) for i in dataB.index.values]

plt.figure()
plt.contourf(field, wavenums, dataB,100, cmap='jet')
plt.xlim(0,14)
plt.ylim(0,120)
# plt.clim(0, 1)
plt.colorbar()
plt.title('CsErSe2 H||AB IR data')

############################################################
# lets put both on the same plot for easy compare
fig = plt.figure()
gs = fig.add_gridspec(1,2, hspace=0, wspace=0)
axs = gs.subplots(sharex=True, sharey=True)
plt1 = axs[0].contourf(field, wavenums, dataB,100, cmap='jet')
axs[0].set_title('IR Data')
plt2 = axs[1].contourf(fieldArr,waveArr,arrAB.T, 100, cmap = 'jet')

axs[1].set_title('Simulated Data')
axs[1].set_xlabel('Field (T)')
axs[0].set_ylabel('Energy (cm$^{-1}$')
axs[0].set_xlabel('Field (T)')

# Hide x labels and tick labels for all but bottom plot.
for ax in axs:
    ax.label_outer()
plt.xlim(0,17.5)
plt.ylim(0,120)
plt.clim(0,1)

fig.subplots_adjust(right=0.8)
cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.7])
fig.colorbar(plt2, cax=cbar_ax)