import h5py





file_path = '/Users/hopeless/Desktop/LeeLab/data/MagnetoOpticalProject_prelimData/KErSe2_FIR_2.h5'

with h5py.File(file_path, 'r') as f:
    # Access a dataset
    wavenums = f['wavenumbers'] 
    wavenums = wavenums[:] 
    wavenums = wavenums[0]

    fields = f['avgField']
    fields = fields[:]
    fields = fields.T[0]

    data = f['avgData0field']
    data = data[:]


# now we invert and normalize
data = (data-np.max(data))*-1
data = data/np.max(data)
data = np.delete(data, [0], axis = 1) # drop zero field (this is easier in a dataframe)
fields = np.delete(fields, [0], axis = 0)

# lets drop anything th
# now we have our data, do the fit as usual
# assume we've loaded the functions (god i need to build a class)
# call the class Magneto Optical Object (MOO)

#model 1 from Allen's paper
B20 = -2.773e-2
B40 = -3.987e-4
B43 = -1.416e-2
B60 = 3.152e-6
B63 = -7.616e-6
B66 =  3.275e-5
Jz =  -1.8e-3*2

# now do fit
model = lmfit.Model(zeemanSplitC_IR, independent_vars=['field', 'wavenum'])
params = model.make_params()

params['B20'].set(value= B20, min = -.06, max = 0.06)
params['B40'].set(value= B40, min = -.06, max = 0.06)
params['B43'].set(value= B43, min = -.06, max = 0.06)
params['B60'].set(value= B60, min = -.06, max = 0.06)
params['B63'].set(value= B63, min = -.06, max = 0.06)
params['B66'].set(value= B66, min = -.06, max = 0.06)
params['Jz'].set(value = Jz)
params['temperature'].set(value= 5, vary = False)


result = model.fit(data.T, field=fields, wavenum=wavenums, params =params)

print(result.fit_report())

# let's plot 
B20 = -0.02341326 # +/- 3.8435e-04 (1.64%) (init = -0.02773)
B40 = -4.1396e-04 # +/- 9.3228e-07 (0.23%) (init = -0.0003987)
B43 = -0.01422081 # +/- 2.8180e-05 (0.20%) (init = -0.01416)
B60 =  3.1573e-06 # +/- 6.3770e-09 (0.20%) (init = 3.152e-06)
B63 = -1.1273e-05 # +/- 2.9517e-07 (2.62%) (init = -7.616e-06)
B66 =  3.1026e-05 # +/- 3.3434e-07 (1.08%) (init = 3.275e-05)
Jz =  -0.00390730 # +/- 2.1525e-05 (0.55%) (init = -0.0036)
temperature =  5 #(fixed)

B20 =-0.02773027 # init = -0.02773)
B40 =-4.6467e-04 # init = -0.0003987)
B43 =-0.01415711 # init = -0.01416)
B60 =-1.1925e-06 # init = 3.152e-06)
B63 =-1.8048e-05 # init = -7.616e-06)
B66 = 3.8798e-05 # init = 3.275e-05)
Jz = -0.00707328 # init = -0.0018)
temperature =  5 # fixed)


waveArr= wavenums
fieldArr = fields 
temperature = 1

ampC, arrC = zeemanSplitLinesC(fieldArr, B20, B40, B43, B60, B63, B66, Jz, temperature)
arrC = np.array(arrC)
arrC = arrC*meVToCm
arrC = arrC.T

ampC = np.array(ampC)
ampC = ampC.T


plt.figure()
plt.contourf(fields, wavenums, data,50)
plt.xlim(0,17.5)
plt.ylim(0,120)

plt.title('KErSe2 H||c with overlayed  calclines\n B20: '+ str(B20)+' B40: '+str(B40)+' B43: ' +str(B43)+ '\n B60: ' +str(B60) + ' B63: ' + str(B63)+ ' B66: ' + str(B66))
for i in range(40):
    if i<16: 
        plt.plot(fieldArr, arrC[i], 'r', alpha=0.7)
    if i>=16:
        plt.plot(fieldArr, arrC[i], 'r--', alpha=0.7)