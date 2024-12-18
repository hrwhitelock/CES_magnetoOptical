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

################################################################
# working B params -> from C-axis raman fit 11/13
# B20 = -0.02926696 # init = -0.03559
# B40 = -3.9097e-04 # init = -0.0003849
# B43 = -0.01391860 # init = -0.01393
# B60 =  3.0584e-06 # init = 3.154e-06
# B63 = -4.2840e-06 # init = -4.695e-06
# B66 =  3.3645e-05 # init = 3.3815e-05




################################################################

# first, lets do the waterfall raman plots

# load data

fname = '/Users/hopeless/Desktop/LeeLab/data/CsErSe2_data/CsErSe2_raman.csv'
rawData = pd.read_csv(fname, index_col=0, skiprows=0, header=1, delimiter=',')
ramanData = rawData

# take the log - make sure to do this FIRST
ramanData = np.log(ramanData)
ramanData = ramanData-ramanData.min(axis=None)
# normalize 
ramanData = ramanData/ramanData.max(axis=None)

# zfsub = ramanData
# for col in ramanData.columns: 
#     zfsub[col] = ramanData[col]-ramanData['0.001']

fields = ['0.001','1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14'] # just want to plot every tesla

n = len(fields)
colors = plt.cm.inferno(np.linspace(0,.8,n))

plt.figure()
i = 0
for f in fields: 
    spec = ramanData[f]+i/3
    plt.plot(ramanData.index, spec, label = f, color= colors[i])
    if i == 0: 
        plt.annotate('H = 0T', xy = (90, i/3 +.15), fontsize = 9)
    elif i ==7 or i==14:
        plt.annotate('H = ' + f + 'T', xy = (90, i/3 +.15), fontsize = 9)
    i+=1
# plt.legend(draggable = True)
plt.xlim(0,100)
plt.title('Raman Data, 0-100cm')
plt.xlabel('Energy (cm$^{-1}$)')
plt.ylabel('Intensity (arb)')

# now plot the raman lines
plt.figure()
i = 0
for f in fields: 
    spec = ramanData[f]+i/3
    plt.plot(ramanData.index, spec, label = f, color= colors[i])
    if i == 0: 
        plt.annotate('H = 0T', xy = (176, i/3 +.05), fontsize = 9)
    elif i ==7 or i==14:
        plt.annotate('H = ' + f + 'T', xy = (176, i/3 +.05), fontsize = 9)
    i+=1
# plt.legend()
plt.xlim(100,500)
plt.title('Raman Data, phonon lines, 100-200 cm')
plt.xlabel('Energy (cm$^{-1}$)')
plt.ylabel('Intensity (arb)')

# now plot high energy lines
plt.figure()
i = 0
for f in fields: 
    spec = ramanData[f]+i/5
    plt.plot(ramanData.index, spec, label = f, color= colors[i])
    if i == 0: 
        plt.annotate('H = 0T', xy = (176, i/3 +.05), fontsize = 9)
    elif i ==7 or i==14:
        plt.annotate('H = ' + f + 'T', xy = (176, i/3 +.05), fontsize = 9)
    i+=1
# plt.legend()
plt.xlim(200,560)
plt.title('Raman Data, high energy, 200-560 cm \n possible 13/2 lines')
plt.xlabel('Energy (cm$^{-1}$)')
plt.ylabel('Intensity (arb)')


################################################################
# now we make some fit plots
# start with the c-axis
# assume we've run the healper functions

# start with lines
field = [float(b) for b in ramanData.columns.values]
wavenums = [float(i) for i in ramanData.index.values]
waveArr= wavenums
fieldArr = np.linspace(0,20, 100)

B20 = -0.03265325 # init = -0.03559)
B40 = -0.0003849 # fixed)
B43 = -0.01393 # fixed)
B60 =  2.75e-6#3.079e-6 #3.054e-06 # fixed) # this B60 param determined from field induced transition
B63 = -8.4011e-07 # init = -4.695e-06)
B66 =  3.3815e-05 # fixed)


ampC, arrC = zeemanSplitLinesC(fieldArr, B20, B40, B43, B60, B63, B66)
arrC = np.array(arrC)
arrC = arrC*meVToCm
arrC = arrC.T

ampC = np.array(ampC)
ampC = ampC.T

cmap = mcolors.ListedColormap(['magenta'])

fig, ax2 = plt.subplots()
plt.contourf(field, wavenums, ramanData,100, cmap = 'jet')
plt.xlim(0,14)
plt.ylim(0,120)
plt.clim(0, 1)
plt.colorbar()
plt.title('CsErSe2 H||c with overlayed  calclines\n B20: '+ str(B20)+' B40: '+str(B40)+' B43: ' +str(B43)+ '\n B60: ' +str(B60) + ' B63: ' + str(B63)+ ' B66: ' + str(B66))
for i in range(len(arrC)):
    if i<16: 
        plt.plot(fieldArr, arrC[i], 'r', alpha=1, linewidth= .7)
    if i>=16:  
        alphas = ampC[i]
        # Create a LineCollection
        points = np.array([fieldArr, arrC[i]]).T.reshape(-1, 1, 2)
        segments = np.concatenate([points[:-1], points[1:]], axis=1)
        lc = LineCollection(segments, cmap=cmap, alpha=alphas)
        ax2.add_collection(lc)
        ax2.autoscale()


##################################################################################################
# lets try this plotting lines over a waterfall plot

n = len(fields)
colors = plt.cm.inferno(np.linspace(0,.8,n))
cmap = mcolors.ListedColormap(['red'])
fig,ax1 =  plt.subplots()
i = 0
for f in fields: 
    spec = ramanData[f]*1.2+i/1.05
    ax1.plot(ramanData.index, spec, label = f, color= colors[i])
    if i == 0: 
        ax1.annotate('H = 0T', xy = (90, i +.15), fontsize = 9)
    elif i ==7 or i==14:
        ax1.annotate('H = ' + f + 'T', xy = (90, i +.15), fontsize = 9)
    i+=1
plt.ylim(0,15)
ax2 = ax1.twinx() 
for i in range(len(arrC)):
    if i<16: 
        plt.plot(arrC[i], fieldArr, 'r', alpha=1, linewidth= .7)
    if i>=16:  
        alphas = ampC[i]
        # Create a LineCollection
        points = np.array([arrC[i], fieldArr]).T.reshape(-1, 1, 2)
        segments = np.concatenate([points[:-1], points[1:]], axis=1)
        lc = LineCollection(segments, cmap=cmap, alpha=alphas)
        ax2.add_collection(lc)
        ax2.autoscale()
# plt.legend(draggable = True)
ax2.set_ylabel('Field [T]')
plt.xlim(0,100)
plt.title('Raman Data, 0-100cm')
ax1.set_xlabel('Energy (cm$^{-1}$)')
ax1.set_ylabel('Intensity (arb)')
plt.ylim(0,15)

################################################
# plot line from Allen's params
B20 =  -0.03559
B40 =  -0.0003849 
B43 =  -0.01393
B60 =   3.154e-06
B63 =  -4.695e-06
B66 =   3.3815e-05

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

################################################
# plot simulated data
B20 = -0.03265325 # init = -0.03559)
B40 = -0.0003849 # fixed)
B43 = -0.01393 # fixed)
B60 =  3.079e-6 #3.054e-06 # fixed) # this B60 param determined from field induced transition
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
width = 0.5
#field, wavenum, B20, B40, B43, B60, B63, B66, amp1, amp2, amp3, amp4, amp5, amp6, amp7, amp8,amp9, amp10,width
arrC = zeemanSplitC(fieldArr, waveArr, B20, B40, B43, B60, B63, B66,
                    amp1, amp2, amp3, amp4, amp5, amp6, amp7, amp8,amp9, amp10, 
                    width)

arrC = np.array(arrC)



plt.figure()
plt.contourf(fieldArr,waveArr,arrC.T, 100, cmap = 'Reds')
# plt.xlim(0,17.5)
# plt.ylim(20,100)
# plt.clim(.1,2.5)
plt.title('simulated fitted CES H||c data \n B20: '+ str(B20)+' B40: '+str(B40)+' B43: ' +str(B43)+ '\n B60: ' +str(B60) + ' B63: ' + str(B63)+ ' B66: ' + str(B66))
plt.xlabel('Field (T)')
plt.ylabel('Wavenumber (cm$^{-1}$)')
plt.colorbar()
with h5py.File('CsErSe2_raman_simulation.h5', 'w') as hdf:
    hdf.create_dataset('ramanData', data=ramanData.values)
    hdf.create_dataset('fields', data=np.array(field, dtype='S'))
    hdf.create_dataset('wavenums', data=wavenums)
    hdf.create_dataset('waveArr', data=waveArr)
    hdf.create_dataset('fieldArr', data=fieldArr)
    hdf.create_dataset('simulatedData', data=arrC)
    hdf.attrs['B20'] = B20
    hdf.attrs['B40'] = B40
    hdf.attrs['B43'] = B43
    hdf.attrs['B60'] = B60
    hdf.attrs['B63'] = B63
    hdf.attrs['B66'] = B66


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

############################################################
# without the data, just a plot of my lines in red, Allen's in blue
# me first!
B20 = -0.03265325 # init = -0.03559)
B40 = -0.0003849 # fixed)
B43 = -0.01393 # fixed)
B60 =  2.75e-6#3.079e-6 # fixed)
B63 = -8.4011e-07 # init = -4.695e-06)
B66 =  3.3815e-05 # fixed)

# field = [float(b) for b in ramanData.columns.values]
# wavenums = [float(i) for i in ramanData.index.values]
# waveArr= wavenums
# fieldArr = field 
fieldArr = np.linspace(0,18,100)
wavenums = np.linspace(0,120,200)

ampC, arrC = zeemanSplitLinesC(fieldArr, B20, B40, B43, B60, B63, B66)
arrC = np.array(arrC)
arrC = arrC*meVToCm
arrC = arrC.T

ampC = np.array(ampC)
ampC = ampC.T

plt.figure()
plt.xlabel('Field (T)')
plt.ylabel('Energy (cm$^{-1}$)')
# plt.title('CsErSe2 H||c with overlayed  calclines\n B20: '+ str(B20)+' B40: '+str(B40)+' B43: ' +str(B43)+ '\n B60: ' +str(B60) + ' B63: ' + str(B63)+ ' B66: ' + str(B66))
for i in range(40):
    if i ==1: 
        plt.plot(fieldArr, arrC[i], 'r', alpha=0.7, label = 'calculated from Raman fit')
    elif i<16 and i !=1: 
        plt.plot(fieldArr, arrC[i], 'r', alpha=0.7)

# now Allen
B20 =  -0.03559
B40 =  -0.0003849 
B43 =  -0.01393
B60 =   3.154e-06
B63 =  -4.695e-06
B66 =   3.3815e-05

ampC, arrC = zeemanSplitLinesC(fieldArr, B20, B40, B43, B60, B63, B66)
arrC = np.array(arrC)
arrC = arrC*meVToCm
arrC = arrC.T

ampC = np.array(ampC)
ampC = ampC.T


# plt.figure()
# plt.contourf(field, wavenums, ramanData,50)
# plt.xlim(0,14)
# plt.ylim(0,120)
# plt.clim(0, 1)
# plt.colorbar()

# plt.title('CsErSe2 H||c with overlayed  calclines\n B20: '+ str(B20)+' B40: '+str(B40)+' B43: ' +str(B43)+ '\n B60: ' +str(B60) + ' B63: ' + str(B63)+ ' B66: ' + str(B66))
for i in range(40):
    if i == 1: 
        plt.plot(fieldArr, arrC[i], 'b', alpha=0.7, label = 'calculated from Neutron fit')
    elif i<16 and i !=1: 
        plt.plot(fieldArr, arrC[i], 'b', alpha=0.7)
plt.legend()
# plt.title('Calculated CES lines')


with h5py.File('CsErSe2_test_B60_lines.h5', 'w') as hdf:
    hdf.create_dataset('ramanData', data=ramanData.values)
    hdf.create_dataset('fields', data=np.array(field, dtype='S'))
    hdf.create_dataset('wavenums', data=wavenums)
    hdf.create_dataset('waveArr', data=waveArr)
    hdf.create_dataset('fieldArr', data=fieldArr)
    hdf.create_dataset('ampC', data=ampC)
    hdf.create_dataset('arrC', data=arrC)
    hdf.attrs['B20'] = B20
    hdf.attrs['B40'] = B40
    hdf.attrs['B43'] = B43
    hdf.attrs['B60'] = B60
    hdf.attrs['B63'] = B63
    hdf.attrs['B66'] = B66



############################################################
# without the data, just a plot of my lines in red, Allen's in blue
# No Jz first
B20 =  -0.03559
B40 =  -0.0003849 
B43 =  -0.01393
B60 =   3.154e-06
B63 =  -4.695e-06
B66 =   3.3815e-05
Jz = 0

# field = [float(b) for b in ramanData.columns.values]
# wavenums = [float(i) for i in ramanData.index.values]
# waveArr= wavenums
# fieldArr = field 
fieldArr = np.linspace(0,18,100)
wavenums = np.linspace(0,120,200)

ampC, arrC = zeemanSplitLinesC(fieldArr, B20, B40, B43, B60, B63, B66, Jz)
arrC = np.array(arrC)
arrC = arrC*meVToCm
arrC = arrC.T

ampC = np.array(ampC)
ampC = ampC.T

plt.figure()
plt.xlabel('Field (T)')
plt.ylabel('Energy (cm$^{-1}$)')
# plt.title('CsErSe2 H||c with overlayed  calclines\n B20: '+ str(B20)+' B40: '+str(B40)+' B43: ' +str(B43)+ '\n B60: ' +str(B60) + ' B63: ' + str(B63)+ ' B66: ' + str(B66))
for i in range(40):
    if i ==1: 
        plt.plot(fieldArr, arrC[i], 'r', alpha=0.7, label = 'Jz = 0')
    elif i<16 and i !=1: 
        plt.plot(fieldArr, arrC[i], 'r', alpha=0.7)

# now Allen
B20 =  -0.03559
B40 =  -0.0003849 
B43 =  -0.01393
B60 =   3.154e-06
B63 =  -4.695e-06
B66 =   3.3815e-05
J = -2.53e-3
ampC, arrC = zeemanSplitLinesC(fieldArr, B20, B40, B43, B60, B63, B66, Jz)
arrC = np.array(arrC)
arrC = arrC*meVToCm
arrC = arrC.T

ampC = np.array(ampC)
ampC = ampC.T


# plt.figure()
# plt.contourf(field, wavenums, ramanData,50)
# plt.xlim(0,14)
# plt.ylim(0,120)
# plt.clim(0, 1)
# plt.colorbar()

# plt.title('CsErSe2 H||c with overlayed  calclines\n B20: '+ str(B20)+' B40: '+str(B40)+' B43: ' +str(B43)+ '\n B60: ' +str(B60) + ' B63: ' + str(B63)+ ' B66: ' + str(B66))
for i in range(40):
    if i == 1: 
        plt.plot(fieldArr, arrC[i], 'b', alpha=0.7, label = 'Jz = 2.53\mu eV')
    elif i<16 and i !=1: 
        plt.plot(fieldArr, arrC[i], 'b', alpha=0.7)
plt.legend()
# plt.title('Calculated CES lines')