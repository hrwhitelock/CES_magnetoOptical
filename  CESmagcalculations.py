# CES mag calculations

import numpy as np
import matplotlib.pyplot as plt
import PyCrystalField as cef
import scipy
from scipy.misc import derivative
import lmfit

plt.ion()

temperature = 5 # in K

muB = 5.7883818012e-2  # meV/T
mu0 = np.pi*4e-7       # T*m/A
kB  = 8.617e-2         # [meV/K]
meVToCm =meVTocCmInv= 8.066 
ion = 'Er3+'

# # my values from fit
B20 = -0.04910422
B40 =  -3.6483e-04
B43 = -0.01474077
B60 =  3.1547e-06
B63 =  3.2378e-06
B66 =  4.2797e-05

# alans neutron fit vals
# B20 = -3.559e-2
# B40 = -3.849e-4
# B43 = -1.393e-2
# B60 = 3.154e-6
# B63 = -4.695e-6
# B66 = 3.3815e-5

g = cef.LandeGFactor(ion)
Bparams =  {'B20': B20, 'B40':B40,'B43': B43, 'B60': B60, 'B63':B63,'B66':B66}
ErObj = cef.CFLevels.Bdict(ion,Bparams)
ionObj = ErObj

kBT = kB*temperature

Jperp = -0.2e-3
Jz = -2.5e-3

def zeemanCEF(ionObj, ion, Field): 
    JdotB = muB*(Field[0]*cef.Operator.Jx(ionObj.J) + Field[1]*cef.Operator.Jy(ionObj.J) + Field[2]*cef.Operator.Jz(ionObj.J))*cef.LandeGFactor(ion)
    H = np.sum([a*b for a,b in zip(ionObj.O, ionObj.B)], axis=0)
    ionObj.diagonalize(H + JdotB.O) # this is just H = Hcef + Hmag
    return ionObj

def xxzShift(ionObj): #, Jperp, Jz): 
    # let's assume we already did the zeeman part
    # starting w/ zeeman int so we can get some exp vals
    # then we'll compute an exp val for J, and apply that as a mf
    q=6
    evals = ionObj.eigenvaluesNoNorm
    evecs = ionObj.eigenvectors
    JexpVals = np.zeros((len(evals),3))
    for i, ev in enumerate(evecs):
        kev = cef.Ket(ev)
        # print np.real(np.dot(ev,kev.Jy().ket)), np.real(np.dot(ev,np.dot(Jy.O,ev)))
        # print np.real(kev*kev.Jy()) - np.real(np.dot(ev,np.dot(Jy.O,ev)))
        JexpVals[i] =[np.real(kev*kev.Jx()),
                      np.real(kev*kev.Jy()),
                      np.real(kev*kev.Jz())]
    Jhat = sum(JexpVals)/sum(sum(JexpVals)) #don't think this is correct, but it's fine for now i guess
    JexpVals = JexpVals.T
    perpPart = [[Jperp*j,0,0]  for j in JexpVals[0]]
    zPart = [[0,0,Jz*j]  for j in JexpVals[2]]
    Exxz = q/2*(Jperp*JexpVals[0]**2+Jz*JexpVals[2]) - q*np.dot(np.add(perpPart, zPart), Jhat)
    # temp = q*np.dot(([Jperp*sum(JexpVals[0]),0.,0.]+[0,0,Jz*sum(JexpVals[2])]), Jhat)
    newEvals = np.add(evals, Exxz)
    minE = min(newEvals)
    # newEvals = [e-minE for e in newEvals]
    return newEvals

xxz = []
field = [[0,0,i] for i in np.linspace(0,10, 100)]
for b in field: 
    ionObj = zeemanCEF(ErObj, ion, b)
    evals = ionObj.eigenvalues#NoNorm
    evals = [e*meVToCm for e in evals]
    shifted = xxzShift(ionObj)
    xxz.append(shifted)

xxz = np.array(xxz)

for x in xxz.T: 
    plt.plot(np.array(field).T[2], x)

def newfreeEnergy(field=[0,0,0], temperature=1, ionObj = ErObj): 
    try: 
        kBT= kB*temperature
        ionObj = zeemanCEF(ionObj, ion, field)
        evals = xxzShift(ionObj)
        Emin - min(evals)
        # E =[eval-Emin for eval in evals]
        Z = [np.exp(-Ei/kBT) for Ei in E]
        Z = sum(Z)
        f = -kBT*np.log(Z) 
    except AttributeError: 
        f =[]
        for b in field: 
            kBT= kB*temperature
            if temperature ==0: 
                temperature = 1e-10
                kBT= kB*temperature
            ionObj = zeemanCEF(ionObj, ion, b)
            evals = xxzShift(ionObj)
            E =[eval for eval in evals]
            minE = min(E)
            # E = [e-minE for e in E]
            Z = [np.exp(-Ei/kBT) for Ei in E]
            Z = sum(Z)
            f.append(-kBT*np.log(Z))
            # f.append(Z)
    return f

f = newfreeEnergy(field, temperature=.5)
plt.plot(np.array(field).T[2], f)
plt.title('free Energy')
plt.xlabel('applied field (T)')
# plt.ylim(0,10e180)
plt.grid(True)
# plt.yscale('log')

def newnumericalmagnetization(lowfield, highfield, temperature = 0.5, numpoints = 300): 
    b, dx = np.linspace(lowfield,highfield,numpoints, retstep=True)
    field = [[0,0,i] for i in b]
    f = newfreeEnergy(field, temperature)
    # fuck this, writing my own gradient?
    m = -np.gradient(f, dx)
    return m,dx

def newnumericalsusceptibility(fieldVal, temps):
    chi = []
    for temp in temps: 
        m, dx = newnumericalmagnetization(1-fieldVal,1+fieldVal, temperature = temp, numpoints =100)
        x = np.gradient(m, .1) 
        # now we've gotta access the very low field value
        valIdx = findIdx(field, [0,0,fieldVal])
        #chi.append(x[valIdx])
        chi.append(x[52])
    return chi

m,dx = newnumericalmagnetization(-7,7, temperature = 2)


plt.plot(np.linspace(-7,7,300), m)
plt.grid(True)
# plt.ylim(0,0.6)
# plt.xlim(-1,1)
plt.title('CsErSe2 calculated magnetization @0.5K')
plt.xlabel('applied field (T)')
plt.ylabel('magnetization (uB/Er)')

temps = np.linspace(.5,300,350)
chi = newnumericalsusceptibility(1, temps)

chiinv = [1/x for x in chi]
plt.plot(temps, chiinv)
plt.grid(True)
plt.title('CsErSe2 calculated susceptibility')
plt.ylabel('1/chi (I think the units are wrong')
plt.xlabel('Temperature (K)')
# plt.xlim(0,5)