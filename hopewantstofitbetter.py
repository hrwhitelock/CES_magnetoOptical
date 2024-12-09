# okay, so I really want to do better fits. I just feel like these fits are like, fine
# I know it's fine, but I also know that I can do better


# first step is to make a custom cost function to minimize
# this should take AB plane IR, c axis raman, and 100mK dm/dh
# I really think fitting the dm/dh will work idk i believe
# the ir and raman need a huge penatly 

def costfn(params): 
    B20 = params[0] 
    B40 = params[1]
    B43 = params[2]
    B60 = params[3]
    B63 = params[4]
    B66 = params[5]
    Jz = params[6]
    
    ion = 'Er3+'

    Bparams =  {'B20': B20, 'B40':B40,'B43': B43, 'B60': B60, 'B63':B63,'B66':B66}
    ionObj = cef.CFLevels.Bdict(ion,Bparams)
    evals = diagonalizeC(ionObj, ion, 0)
    if evals[8]*meVToCm >200: 
        cost = 1e12
    else:
        calcC = zeemanSplitC(fieldC, wavenumsC, ionObj)
        calcB = zeemanSplitAB(fieldB, wavenumsB,ionObj)
        calcdmdh = dmdhfn(arr, ionObj, Jz)
        costC = np.sum((dataC - calcC)**2/calcC) # make amp set val, amp1, amp2, amp3, amp4, amp5, amp6, amp7, amp8,amp9, amp10,width)
        costB = np.sum((dataB - calcB)**2/calcB)
        costdmdh = np.sum((dmdhData - calcdmdh)**2/calcdmdh)
    
    
    cost = costC+costB+costdmdh
    # check if line 8 > 200cm^-1 
    # picking this line, because both c and b data show this much lower
    
    return cost

def lorentzian( wave, amp, cen, wid ):
    return [amp * wid**2 / ( wid**2 + ( x - cen )**2) for x in wave]


def diagonalizeC(ionObj, ion, Jz, H, temperature): 
    # first calc effective h
    H = MFTmagC(ionObj, H, Jz, temperature)
    JdotB = muB*(H*cef.Operator.Jz(ionObj.J))*cef.LandeGFactor(ion)
    H = np.sum([a*b for a,b in zip(ionObj.O, ionObj.B)], axis=0)
    ionObj.diagonalize(H + JdotB.O) # this is just H = Hcef + Hmag
    return ionObj.eigenvalues 

def zeemanSplitC(field, wavenum, B20, B40, B43, B60, B63, B66, Jz):    
    # assuming that x is an array
    # amp = [amp1, amp2, amp3, amp4, amp5, amp6, amp7, amp8, amp9, amp10]#, amp11, amp12, amp13, amp14, amp15, amp16]
    amp = [0.1,0.3,0.3,0.15,0.2,0.2,0.287, 0.2, 0.135, 0.097]
    dEphonon = 49.3
    phononAmp = 0.499
    phononSig = 0.95
    fun = []
    Bparams =  {'B20': B20, 'B40':B40,'B43': B43, 'B60': B60, 'B63':B63,'B66':B66}
    ionObj = cef.CFLevels.Bdict(ion,Bparams)
    temperature = 10
    for b in field: 
        diagonalizeC(ionObj, ion, Jz, b, temperature)
        dE =[eval for eval in evals] # this is the spitting if everything is in the GS -> not necessarily true for finite temp
        dE = dE[0:10] # only want to look at the bottome few lines - lets see if this works??
        tempAmp = amp
        Z = [np.exp(-Ei/kBT) for Ei in dE]
        Z = sum(Z)
        p = [1/Z*np.exp(-Ei/kBT) for Ei in dE]
        numLines = len(dE)
        for i in range(1,numLines): 
            for j in range(i+1, numLines):
                temp = dE[j]-dE[i]
                dE.append(temp)
                tempAmp.append(p[i]*amp[j])
        wid = 0.6
        centers = dE
        tempfun = lorentzian(wavenum, phononAmp, dEphonon, phononSig)
        for i in range(len(centers)):
            a = tempAmp[i]
            tempfun += lorentzian(wavenum, a, centers[i]*meVToCm, wid)
        fun.append(tempfun)
    return fun
def diagonalizeAB(ionObj, ion, field): 
    JdotB = muB*(field*cef.Operator.Jy(ionObj.J))*cef.LandeGFactor(ion)
    H = np.sum([a*b for a,b in zip(ionObj.O, ionObj.B)], axis=0)
    ionObj.diagonalize(H + JdotB.O) # this is just H = Hcef + Hmag
    evals = ionObj.eigenvalues
    return evals
def zeemanSplitAB(field, wavenum, ionObj):     
    amp = [0.1,0.3,0.3,0.15,0.2,0.2,0.287, 0.2, 0.135, 0.097]
    fun = []
    wid = 0.7
    for b in field: 
        evals = diagonalizeAB(ionObj, ion, b)
        dE =[eval for eval in evals] # this is the spitting if everything is in the GS -> not necessarily true for finite temp
        dE = dE[0:10] # only want to look at the bottome few lines - lets see if this works??
        tempAmp = amp
        Z = [np.exp(-Ei/kBT) for Ei in dE]
        Z = sum(Z)
        p = [1/Z*np.exp(-Ei/kBT) for Ei in dE]
        numLines = len(dE)
        for i in range(1,numLines): 
            for j in range(i+1, numLines):
                temp = dE[j]-dE[i]
                dE.append(temp)
                tempAmp.append(p[i]*amp[j])
        centers = dE
        tempfun = np.zeros(len(wavenum)).tolist()
        for i in range(len(centers)):
            a = tempAmp[i]
            tempfun += lorentzian(wavenum, a, centers[i]*meVToCm, wid)
        fun.append(tempfun)
    return np.array(fun)

def dmdhfn(field, ionObj, J): 
    #first calc mag
    mag = MFTmagC(ionObj, field, J)
    dmdh = np.gradient(mag, field)*1e-6
    return dmdh
def cmag(mag, J,  h):
    newh = h +q*J*mag/muB/1.2/1.2
    mag = -1*MyErObj.magnetization(ion, 0.1, [0,0, newh]).T[2]-mag
    return mag
def MFTmagC(ionObj, H, J): 
    # okay so let's start by defining the MF ham
    # only doing c direction rn because that is faster
    q = 6
    ion = 'Er3+' # hard coded for now dont' @ me
    n = 3 # iterations
    # first step, call the pcf mag because that's a good place to start
    magArr = []
    for h in H: 
        mag = 0
        for i in range(n): 
            mag = fsolve(cmag, mag, args = (J, h))[0] # fsolve spits out an array - this time its one val
        magArr.append(mag) # fso
    return magArr

fieldC = [float(b) for b in dataC.columns.values]
wavenumsC = [float(i) for i in dataC.index.values]

fieldB = [float(b) for b in dataB.columns.values]
wavenumsB = [float(i) for i in dataB.index.values]


# load the 100mK data
dmdhData = np.genfromtxt('/Users/hopeless/Desktop/LeeLab/data/CsErSe2_data/CES_dmdh_100mK.txt', 
                       delimiter='\t', unpack=True, skip_header=1)
# let's interpolate this to save on computation time
arr = np.linspace(0,7,100)
dmdhData = np.interp(arr.flatten(), dmdhData[0].flatten(),dmdhData[1].flatten()) 
dmdhData[0] = arr

import scipy

B20 = -0.03265325 # init = -0.03559)
B40 = -0.0003849 # fixed)
B43 = -0.01393 # fixed)
B60 =  3.079e-6 #3.054e-06 # fixed) # this B60 param determined from field induced transition
B63 = -8.4011e-07 # init = -4.695e-06)
B66 =  3.3815e-05 # fixed)
Jz=0.48e-3

result = scipy.optimize(costfn, [B20,B40,B43,B60,B63,B66,Jz])
\


## let me just add mft to the c axis zeeman and see what that does
def cmag(mag, J,  h, temperature):
    newh = h +q*J*mag/muB/1.2/1.2
    mag = -1*MyErObj.magnetization(ion, temperature, [0,0, newh]).T[2]-mag
    return mag

def newH(ionObj, H, J, temperature): 
    # okay so let's start by defining the MF ham
    # only doing c direction rn because that is faster
    q = 6
    ion = 'Er3+' # hard coded for now dont' @ me
    n = 4 # iterations
    # first step, call the pcf mag because that's a good place to start
    mag = 0
    for i in range(n): 
        mag = fsolve(cmag, mag, args = (J, H, temperature))[0] # fsolve spits out an array - this time its one val
    newh = H+q*J*mag/muB/1.2/1.2
    return newh

def newHAB(ionObj, H, J, temperature): 
    # okay so let's start by defining the MF ham
    # only doing c direction rn because that is faster
    q = 6
    ion = 'Er3+' # hard coded for now dont' @ me
    n = 4 # iterations
    # first step, call the pcf mag because that's a good place to start
    mag = 0
    for i in range(n): 
        mag = fsolve(cmag, mag, args = (J, H, temperature))[0] # fsolve spits out an array - this time its one val
    newh = H+q*J*mag/muB/1.2/1.2
    return newh

def diagonalizeC(ionObj, ion, Jz, H, temperature): 
    # first calc effective h
    H = MFTmagC(ionObj, H, Jz, temperature)
    JdotB = muB*(H*cef.Operator.Jz(ionObj.J))*cef.LandeGFactor(ion)
    H = np.sum([a*b for a,b in zip(ionObj.O, ionObj.B)], axis=0)
    ionObj.diagonalize(H + JdotB.O) # this is just H = Hcef + Hmag
    return ionObj.eigenvalues 

def zeemanSplitC(field, wavenum, B20, B40, B43, B60, B63, B66, Jz):    
    # assuming that x is an array
    # amp = [amp1, amp2, amp3, amp4, amp5, amp6, amp7, amp8, amp9, amp10]#, amp11, amp12, amp13, amp14, amp15, amp16]
    amp = [0.1,0.3,0.3,0.15,0.2,0.2,0.287, 0.2, 0.135, 0.097]
    dEphonon = 49.3
    phononAmp = 0.499
    phononSig = 0.95
    fun = []
    Bparams =  {'B20': B20, 'B40':B40,'B43': B43, 'B60': B60, 'B63':B63,'B66':B66}
    ionObj = cef.CFLevels.Bdict(ion,Bparams)
    temperature = 10
    for b in field: 
        diagonalizeC(ionObj, ion, Jz, b, temperature)
        dE =[eval for eval in evals] # this is the spitting if everything is in the GS -> not necessarily true for finite temp
        dE = dE[0:10] # only want to look at the bottome few lines - lets see if this works??
        tempAmp = amp
        Z = [np.exp(-Ei/kBT) for Ei in dE]
        Z = sum(Z)
        p = [1/Z*np.exp(-Ei/kBT) for Ei in dE]
        numLines = len(dE)
        for i in range(1,numLines): 
            for j in range(i+1, numLines):
                temp = dE[j]-dE[i]
                dE.append(temp)
                tempAmp.append(p[i]*amp[j])
        wid = 0.6
        centers = dE
        tempfun = lorentzian(wavenum, phononAmp, dEphonon, phononSig)
        for i in range(len(centers)):
            a = tempAmp[i]
            tempfun += lorentzian(wavenum, a, centers[i]*meVToCm, wid)
        fun.append(tempfun)
    return fun
#model 1 from Allen's paper
B20 = -3.559e-2
B40 = -3.849e-4
B43 = -1.393e-2
B60 = 3.154e-6
B63 = -4.695e-6
B66 = 3.3815e-5
field = [float(b) for b in fitData.columns.values]
wavenum = [float(i) for i in fitData.index.values]
# now do fit
model = lmfit.Model(zeemanSplitC, independent_vars=['field', 'wavenum'])
params = model.make_params()

params['B20'].set(value= B20, min = -.06, max = 0.06)
params['B40'].set(value= B40, min = -.06, max = 0.06)
params['B43'].set(value= B43, min = -.06, max = 0.06)
params['B60'].set(value= B60, min = -.06, max = 0.06)
params['B63'].set(value= B63, min = -.06, max = 0.06)
params['B66'].set(value= B66, min = -.06, max = 0.06)
params['Jz'].set(value = 0.5e-3, vary = False)

z = np.array(fitData.to_numpy()) # gotta do it twice with tuples :((((
z = z.T

result = model.fit(z, field=field, wavenum=wavenum, params =params)

print(result.fit_report())

def zeemanSplitLinesC(field, B20, B40, B43, B60, B63, B66, Jz):     
    # assuming only H||B rn
    # assuming that x is an array
    amp = []#[1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1,]#[0, .15, .15, .2, 0.15,0.15,0.15,0.15,0.07,0.07, .1,.1,.1,.1,.1]
    Bparams =  {'B20': B20, 'B40':B40,'B43': B43, 'B60': B60, 'B63':B63,'B66':B66}
    ionObj = cef.CFLevels.Bdict(ion,Bparams)
    # kBT = 10*kB
    temperature = 10
    amp = []
    dE = []
    for b in field: 
        evals = diagonalizeC(ionObj, ion, Jz, b, temperature)
        dE_temp =[eval for eval in evals] # this is the spitting if everything is in the GS -> not necessarily true for finite temp
        # first, lets calculate the partition function - the microstates of the system don't change, jsut the number of abs line
        # we can make between them, so we can make this now
        # Z = sum_i (e^(-beta*E_i))
        Z = [np.exp(-Ei/kBT) for Ei in dE_temp]
        Z = sum(Z)
        # now that we have the partition fn, we can calculate probabilities
        p = [1/Z*np.exp(-Ei/kBT) for Ei in dE_temp]
        # temp_amp = [1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1,]
        # okay, so now we want to add lines between non ground states 
        # we want the amplitude to be the probability -> main lines are already determined
        numLines = len(dE_temp)
        for i in range(1,numLines): 
            # skip GS - already have those dE
            for j in range(i+1, numLines):
                temp = dE_temp[j]-dE_temp[i]
                dE_temp.append(temp)
                p.append(p[i])
        dE.append(dE_temp)
        amp.append(p)
    return amp, dE