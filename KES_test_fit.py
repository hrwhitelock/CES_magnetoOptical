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
# now we have our data, do the fit as usual
# assume we've loaded the functions (god i need to build a class)
# call the class Magneto Optical Object (MOO)

#model 1 from Allen's paper
B20 = 0#-3.559e-2
B40 = 0#-3.849e-4
B43 = 0#-1.393e-2
B60 = 0# 3.154e-6
B63 = 0#-4.695e-6
B66 = 0# 3.3815e-5

# now do fit
model = lmfit.Model(zeemanSplitC, independent_vars=['field', 'wavenum'])
params = model.make_params()

params['B20'].set(value= B20, min = -.06, max = 0.06)
params['B40'].set(value= B40, min = -.06, max = 0.06)
params['B43'].set(value= B43, min = -.06, max = 0.06)
params['B60'].set(value= B60, min = -.06, max = 0.06)
params['B63'].set(value= B63, min = -.06, max = 0.06)
params['B66'].set(value= B66, min = -.06, max = 0.06)
params['Jz'].set(value = 0)
params['temperature'].set(value= 10)


result = model.fit(data.T, field=fields, wavenum=wavenums, params =params, method = 'ampgo')

print(result.fit_report())

# let's plot

calc_field = np.arange(0,18, 0.02)
calc_wavenums = np.arange(0,120, 0.1)
temperature = 5
kBT = kB*temperature

ampC,arrC = zeemanSplitLinesC(calc_field, B20, B40, B43, B60, B63, B66, Jz, temperature)
arrC = np.array(arrC)
arrC = arrC*meVToCm
arrC = arrC.T