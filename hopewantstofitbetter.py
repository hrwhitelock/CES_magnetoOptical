# okay, so I really want to do better fits. I just feel like these fits are like, fine
# I know it's fine, but I also know that I can do better


# first step is to make a custom cost function to minimize
# this should take AB plane IR, c axis raman, and 100mK dm/dh
# I really think fitting the dm/dh will work idk i believe
# the ir and raman need a huge penatly 

def costfn(field, params): 
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

    costC = np.sum((dataC - zeemanSplitC(fieldC, wavenumC, ionObj))**2/zeemanSplitC(fieldC, wavenumC, ionObj)) # make amp set val, amp1, amp2, amp3, amp4, amp5, amp6, amp7, amp8,amp9, amp10,width)
    costB = np.sum((dataB - zeemanSplitAB(fieldB, wavenumB,ionObj))/**2/zeemanSplitAB(fieldB, wavenumB, ionObj))
    costdmdh = np.sum((dmdhData - dmdhFn(field, ionObj, J))**2/dmdhFn(field, ionObj, J))
    
    
    cost = costC+costB+costdmdh
    # check if line 8 > 200cm^-1 
    # picking this line, because both c and b data show this much lower
    evals = diagonalizeC(ionObj, ion, 0)
    if evals[8]*meVToCm >200: 
        cost = cost*1e8
    return cost



def zeemanSplitC(field, wavenum, ionObj):    
    # assuming that x is an array
    # amp = [amp1, amp2, amp3, amp4, amp5, amp6, amp7, amp8, amp9, amp10]#, amp11, amp12, amp13, amp14, amp15, amp16]
    amp = [0.1,0.3,0.3,0.15,0.2,0.2,0.287, 0.2, 0.135, 0.097]
    dEphonon = 49.3
    phononAmp = 0.499
    phononSig = 0.95
    fun = []
    for b in field: 
        evals = diagonalizeC(ionObj, ion, b)
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
    return np.array(fun)

def diagonalizeC(ionObj, ion, Field): 
    JdotB = muB*(Field*cef.Operator.Jz(ionObj.J))*cef.LandeGFactor(ion)
    H = np.sum([a*b for a,b in zip(ionObj.O, ionObj.B)], axis=0)
    ionObj.diagonalize(H + JdotB.O) # this is just H = Hcef + Hmag
    return ionObj.eigenvalues 

def diagonalizeAB(ionObj, ion, field): 
    JdotB = muB*(field*cef.Operator.Jy(ionObj.J))*cef.LandeGFactor(ion)
    H = np.sum([a*b for a,b in zip(ionObj.O, ionObj.B)], axis=0)
    ionObj.diagonalize(H + JdotB.O) # this is just H = Hcef + Hmag
    evals = ionObj.eigenvalues
    return evals
def zeemanSplitAB(field, wavenum, ionObj):     
    amp = [0.1,0.3,0.3,0.15,0.2,0.2,0.287, 0.2, 0.135, 0.097]
    fun = []
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
        tempfun = 0
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
            mag = fsolve(cmag, mag, args = (J, h, 0.1))[0] # fsolve spits out an array - this time its one val
        magArr.append(mag) # fso
    return magArr

fieldC = [float(b) for b in dataC.columns.values]
wavenumsC = [float(i) for i in dataC.index.values]

fieldB = [float(b) for b in dataB.columns.values]
wavenumsB = [float(i) for i in dataB.index.values]


# load the 100mK data
RawMTdata = np.genfromtxt('/Users/hopeless/Desktop/LeeLab/data/CsErSe2_data/CsErSe2_MTall.dat', 
                       delimiter='\t', unpack=True, skip_header=1)