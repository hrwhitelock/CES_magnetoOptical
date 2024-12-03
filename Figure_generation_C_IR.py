import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import PyCrystalField as cef
import scipy
from scipy.optimize import minimize
from scipy.optimize import leastsq
import pandas as pd
import lmfit
import matplotlib.colors as mcolors
from matplotlib import font_manager
font_manager.fontManager.addfont('/Users/hopeless/Library/Fonts/cmunrm.ttf')
plt.rcParams['mathtext.fontset'] = 'stix'
plt.rcParams['font.family'] = 'CMU Serif'
plt.rcParams['axes.unicode_minus'] = False
plt.rcParams.update({'font.size': 16})

plt.ion()


# okay so now we want to take a look at the H||b data
# first load

cfname = '/Users/hopeless/Desktop/LeeLab/data/CsErSe2_data/Maglab_IRMeasurement_June2022/ProcessingCode/Load2_TrimData/P3_CsEr_100_RAWAVG.dat'

# clean up AB first
dataC = pd.read_csv(cfname, index_col=0, skiprows=0, header=1)
dataC = dataC.dropna(axis = 0)
dataC = dataC.dropna(axis=1)
dataC = dataC.drop(labels = '-1.1', axis=1)
rawData = dataC
rawData = rawData.drop(labels = '-1.1', axis=1)

normSpec = dataC['0.001']/max(dataC['0.001'])*-1
avgSpec = normSpec
for column in dataC.columns: 
    dataC[column] = max(dataC[column]) -dataC[column]
    dataC[column] = dataC[column]/(max(dataC[column])) -normSpec
    avgSpec = avgSpec + dataC[column]

for column in dataC.columns: 
    dataC[column] = dataC[column]-avgSpec/len(dataC.columns)
    dataC[column] = dataC[column]-(sum(dataC[column])/len(dataC[column]))

dataC = dataC.drop(labels='0.001', axis=1) # drop this column because we used it as bg


dataC = dataC/dataC.max(axis=None)



################################################################
# now we make some fit plots

# assume we've run the healper functions

# start with lines
field = [float(b) for b in dataC.columns.values]
wavenums = [float(i) for i in dataC.index.values]
waveArr= wavenums
fieldArr = np.linspace(0,20, 100)

B20 = -0.03265325 # init = -0.03559)
B40 = -0.0003849 # fixed)
B43 = -0.01393 # fixed)
B60 =  3.079e-6 #3.054e-06 # fixed) # this B60 param determined from field induced transition
B63 = -8.4011e-07 # init = -4.695e-06)
B66 =  3.3815e-05 # fixed)


ampC, arrC = zeemanSplitLinesC(fieldArr, B20, B40, B43, B60, B63, B66)
arrC = np.array(arrC)
arrC = arrC*meVToCm
arrC = arrC.T

ampC = np.array(ampC)
ampC = ampC.T

cmap = mcolors.ListedColormap(['cyan'])

fig, ax = plt.subplots()
plt.contourf(field, wavenums, dataC,100, cmap = 'gnuplot')
plt.xlim(0,14)
plt.ylim(0,120)
plt.clim(-0.3, 1)
plt.xlabel('Field (T)')
plt.ylabel('Energy (cm$^{-1}$)')
plt.colorbar()
plt.title('CsErSe2 H||C with overlayed  calclines\n B20: '+ str(B20)+' B40: '+str(B40)+' B43: ' +str(B43)+ '\n B60: ' +str(B60) + ' B63: ' + str(B63)+ ' B66: ' + str(B66))
for i in range(len(arrC)):
    if i<16: 
        plt.plot(fieldArr, arrC[i], 'c', alpha=1, linewidth= .7)
    if i>=16:  
        alphas = ampC[i]
        # Create a LineCollection
        points = np.array([fieldArr, arrC[i]]).T.reshape(-1, 1, 2)
        segments = np.concatenate([points[:-1], points[1:]], axis=1)
        lc = LineCollection(segments, cmap=cmap, alpha=alphas)
        ax.add_collection(lc)
        ax.autoscale()
plt.show()


##################################################################################################
# lets try this plotting lines over a waterfall plot
fields = ['0.249', '1', '1.999', '2.999', '3.999', '5', '6', '7','7.999','8.999','9.999', '10.999',  '12', 
       '12.999', '14',  '15',  '15.999', '17']
n = len(fields)
colors = plt.cm.cool(np.linspace(0,1,n))

fig,ax1 =  plt.subplots()
i = 0
for f in fields: 
    spec = dataC[f]+i
    ax1.plot(dataC.index, spec, label = f, color= colors[i])
    if i == 0: 
        ax1.annotate('H = 0T', xy = (90, i +.15), fontsize = 9)
    elif i ==7 or i==14:
        ax1.annotate('H = ' + f + 'T', xy = (90, i +.15), fontsize = 9)
    i+=1
plt.ylim(0,18)
ax2 = ax1.twinx() 
for i in range(len(arrC)):
    if i<16: 
        plt.plot(arrC[i], fieldArr, 'c', alpha=1, linewidth= .7)
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
plt.title('IR Data, 0-100cm')
ax1.set_xlabel('Energy (cm$^{-1}$)')
ax1.set_ylabel('Intensity (arb)')
plt.ylim(0,17.5)

################################################
# plot line from Allen's params
B20 =  -0.03559
B40 =  -0.0003849 
B43 =  -0.01393
B60 =   3.154e-06
B63 =  -4.695e-06
B66 =   3.3815e-05

field = [float(b) for b in dataC.columns.values]
wavenums = [float(i) for i in dataC.index.values]
waveArr= wavenums
fieldArr = field 


ampC, arrC = zeemanSplitLinesC(fieldArr, B20, B40, B43, B60, B63, B66)
arrC = np.array(arrC)
arrC = arrC*meVToCm
arrC = arrC.T

ampC = np.array(ampC)
ampC = ampC.T

plt.figure()
plt.contourf(field, wavenums, dataC,50)
plt.xlim(0,14)
plt.ylim(0,120)
plt.clim(0, 1)
plt.colorbar()
plt.title('CsErSe2 H||C with overlayed  calclines\n B20: '+ str(B20)+' B40: '+str(B40)+' B43: ' +str(B43)+ '\n B60: ' +str(B60) + ' B63: ' + str(B63)+ ' B66: ' + str(B66))
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
width = 0.9
#field, wavenum, B20, B40, B43, B60, B63, B66, amp1, amp2, amp3, amp4, amp5, amp6, amp7, amp8,amp9, amp10,width
arrC = zeemanSplitC(fieldArr, waveArr, B20, B40, B43, B60, B63, B66,
                    amp1, amp2, amp3, amp4, amp5, amp6, amp7, amp8,amp9, amp10, 
                    width)

arrC = np.array(arrC)



plt.figure()
plt.contourf(fieldArr,waveArr,arrC.T, 100, cmap = 'gnuplot')
# plt.xlim(0,17.5)
# plt.ylim(20,100)
# plt.clim(.1,2.5)
plt.title('simulated fitted CES H||C data \n B20: '+ str(B20)+' B40: '+str(B40)+' B43: ' +str(B43)+ '\n B60: ' +str(B60) + ' B63: ' + str(B63)+ ' B66: ' + str(B66))
plt.xlabel('Field (T)')
plt.ylabel('Wavenumber (cm$^{-1}$)')
plt.colorbar()

############################################################
# now let's plot jsut the IR data (H||b)
field = [float(b) for b in dataC.columns.values]
wavenums = [float(i) for i in dataC.index.values]

plt.figure()
plt.contourf(field, wavenums, dataC,100, cmap='gnuplot')
plt.xlim(0,14)
plt.ylim(0,120)
# plt.clim(0, 1)
plt.colorbar()
plt.title('CsErSe2 H||C IR data')

############################################################
# lets put both on the same plot for easy compare
fig = plt.figure()
gs = fig.add_gridspec(1,2, hspace=0, wspace=0)
axs = gs.subplots(sharex=True, sharey=True)
plt1 = axs[0].contourf(field, wavenums, dataC,100, cmap='gnuplot')
axs[0].set_title('IR Data')
plt2 = axs[1].contourf(fieldArr,waveArr,arrC.T, 100, cmap = 'gnuplot')

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

########################################################################################################################
# let's make some plots like the IR data

fields = ['0.249', '1', '1.999', '2.999', '3.999', '5', '6', '7','7.999','8.999','9.999', '10.999',  '12', 
       '12.999', '14',  '15',  '15.999', '17']
n = len(fields)
colors = plt.cm.cool(np.linspace(0,.8,n))

plt.figure()
i = 0
for f in fields: 
    spec = dataC[f]+i/3
    plt.plot(dataC.index, spec, label = f, color= colors[i])
    if i == 0: 
        plt.annotate('H = 0T', xy = (90, i/3 +.15), fontsize = 9)
    elif i ==7 or i==14:
        plt.annotate('H = ' + f + 'T', xy = (90, i/3 +.15), fontsize = 9)
    i+=1
# plt.legend(draggable = True)
plt.xlim(0,100)
plt.title('IR data, 0-100cm')
plt.xlabel('Energy (cm$^{-1}$)')
plt.ylabel('Intensity (arb)')

# now plot the raman lines
plt.figure()
i = 0
for f in fields: 
    spec = dataC[f]+i/3
    plt.plot(dataC.index, spec, label = f, color= colors[i])
    if i == 0: 
        plt.annotate('H = 0T', xy = (176, i/3 +.05), fontsize = 9)
    elif i ==7 or i==14:
        plt.annotate('H = ' + f + 'T', xy = (176, i/3 +.05), fontsize = 9)
    i+=1
# plt.legend()
plt.xlim(100,200)
plt.title('IR data, IR deadzone, 100-200 cm')
plt.xlabel('Energy (cm$^{-1}$)')
plt.ylabel('Intensity (arb)')

# now plot high energy lines
plt.figure()
i = 0
for f in fields: 
    spec = dataC[f]+i/5
    plt.plot(dataC.index, spec, label = f, color= colors[i])
    if i == 0: 
        plt.annotate('H = 0T', xy = (176, i/3 +.05), fontsize = 9)
    elif i ==7 or i==14:
        plt.annotate('H = ' + f + 'T', xy = (176, i/3 +.05), fontsize = 9)
    i+=1
# plt.legend()
# plt.xlim(200,560)
plt.title('IR data, high energy')
plt.xlabel('Energy (cm$^{-1}$)')
plt.ylabel('Intensity (arb)')

########################################################################################################################
# let's make some raw data plots to see 
field = [float(b) for b in rawData.columns.values]
wavenums = [float(i) for i in rawData.index.values]

from mpl_toolkits.mplot3d import axes3d
from matplotlib.ticker import LinearLocator

fig, ax = plt.subplots(subplot_kw={"projection": "3d"})

# Make data.
[X,Y] = np.meshgrid(field, wavenums)
Z =rawData

# Plot the surface.
surf = ax.plot_surface(X, Y, Z, cmap=cm.cool,
                       linewidth=0, antialiased=False)

# Customize the z axis.
ax.set_zlim(-1.01, 1.01)
ax.zaxis.set_major_locator(LinearLocator(10))
# A StrMethodFormatter is used automatically
ax.zaxis.set_major_formatter('{x:.02f}')

# Add a color bar which maps values to colors.
fig.colorbar(surf, shrink=0.5, aspect=5)

plt.show()


########################################################################################################################
# make the same data vs calculation plot
field = [float(b) for b in dataC.columns.values]
wavenums = [float(i) for i in dataC.index.values]
waveArr= np.linspace(0,120, 480)
fieldArr = np.linspace(0,18, 100)
fig = plt.figure()
gs = fig.add_gridspec(1,2, hspace=0, wspace=0)
axs = gs.subplots(sharex=True, sharey=True)
plt1 = axs[0].contourf(field, wavenums, dataC,100, cmap='jet')
plt1.set_clim(vmin = 0,vmax = 1)
axs[0].set_title('IR Data H||c')
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


########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
 # ookay, let's make some line plots for fig 1
 # panel 1, hline at evals
 # panel 2, ab plane
 # panel 3, c axis

 def 