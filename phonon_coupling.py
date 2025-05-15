
# doing some shenanigans w/ phonon coupling

# start w/ calculating gap 

# ass. we've already run all the stuff we gotta

temperature = .1# in K

muB = 5.7883818012e-2  # meV/T
mu0 = np.pi*4e-7       # T*m/A
kB  = 8.617e-2         # [meV/K]
meVToCm =meVTocCmInv= 8.066 
ion = 'Er3+'

kBT = kB*temperature
gJ = cef.LandeGFactor('Er3+')
q = 6

calc_field = np.linspace(0,5,100)
ampC, arrC = zeemanSplitLinesC(calc_field, B20, B40, B43, B60, B63, B66, Jz, temperature)
arrC = np.array(arrC)
def zeemanSplitLinesC(field, B20, B40, B43, B60, B63, B66, Jz, temperature):     
    # assuming only H||B rn
    # assuming that x is an array
    amp = []#[1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1,]#[0, .15, .15, .2, 0.15,0.15,0.15,0.15,0.07,0.07, .1,.1,.1,.1,.1]
    Bparams =  {'B20': B20, 'B40':B40,'B43': B43, 'B60': B60, 'B63':B63,'B66':B66}
    ionObj = cef.CFLevels.Bdict(ion,Bparams)
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
#arrC[1] is our gamma

# constants
hbarvs  = 12.93 # meV A # taken from CYS paper shamelessly
# def e(H, k, eta): 
#     e0 = hbarvs*k # figure this out??
#     na, gaps = zeemanSplitLinesC(H, B20, B40, B43, B60, B63, B66, Jz, temperature)
#     gaps = np.array(gaps)
#     gaps = gaps.T
#     gamma = gaps[1]
#     ep = (e0 + gamma*H)/2 + np.sqrt(((e0-gamma*H)/2)*((e0-gamma*H)/2) + eta*e0*gamma*gamma*H*H)
#     em = (e0 + gamma*H)/2 - np.sqrt(((e0-gamma*H)/2)*((e0-gamma*H)/2) + eta*e0*gamma*gamma*H*H)
#     return ep, em, e0

# okay, so now we wanna plot for k = 0 to pi/a
# say a = 16.64 A 

# kArr = np.linspace(0, np.pi/16.64, 10000)
# epArr = []
# emArr = []
# e0Arr = []
# H = [.5]
# eta = 1
# for k in kArr: 
#     ep, em, e0 = e(H, k, eta)
#     epArr.append(ep)
#     emArr.append(em)
#     e0Arr.append(e0)

# 
# i = 9
# T = 0.08617328149741*temperature
# plt.figure()
# plt.hlines(y = T, xmin = kArr[0], xmax = kArr[-1], colors = 'black')
# plt.plot(kArr, epArr[i], label = '+')
# plt.plot(kArr, emArr[i], label = '-')
# # plt.plot(kArr, e0Arr[0], 'g--')
# plt.title('H = '+str(H[i]))
# plt.legend()
# # plt.xlabel('k')
# plt.ylabel('E')
## okay, so now as a function of field we want to plot available DOS

hbarvs  = 10
eta = 1
def evk(H):
    kArr = np.linspace(0, np.pi/10, 1000)
    emArr = []
    epArr = []
    for h in H: 
        na, gaps = zeemanSplitLinesC([h], B20, B40, B43, B60, B63, B66, Jz, temperature)
        gaps = np.array(gaps)
        gaps = gaps.T
        gamma = gaps[1]
        sump = 0
        summ = 0
        ep = []
        em = []
        for k in kArr: 
            e0 = hbarvs*k
            ep.append((e0 + gamma*h)/2 + np.sqrt(((e0-gamma*h)/2)*((e0-gamma*h)/2) + eta*e0*gamma*gamma*h*h))
            em.append((e0 + gamma*h)/2 - np.sqrt(((e0-gamma*h)/2)*((e0-gamma*h)/2) + eta*e0*gamma*gamma*h*h))
        epArr.append(ep)
        emArr.append(em)
        print(h)
    return emArr, epArr 

H = np.linspace(0,15, 5000)
kArr = np.linspace(0, np.pi/10, 1000)
emArr, epArr = evk(H)


plt.figure()
T = 0.08617328149741*temperature

for h, i in zip(H, np.arange(len(H))): 
    plt.plot(kArr, emArr[45], label = str(h))
    plt.plot(kArr, epArr[45], label = str(h))

plt.hlines(y = T, xmin = kArr[0], xmax = kArr[-1], colors = 'black')

# so now we plot DOS vs H
# have E vs k, want to integrate over k
dosm = []
dosp = []
dos = []
E = []
kmArr = []
kpArr = []
Ep = np.linspace(0, 2*T, 1000)
for h, i in zip(H, np.arange(len(H))):
    m = [l[0] for l in emArr[i]]
    p = [l[0] for l in epArr[i]]

    # we want to integrate over energy, which is currently m, p
    # so let's linearly interp
    newKp = np.interp(Ep, p, kArr)
    newKm = np.interp(Ep, m, kArr)
    kmArr.append(newKm)
    kpArr.append(newKp)
    print(h)

dosm = []
dosp = []
dos = []


for h, i in zip(H, np.arange(len(H))):
    temporaryp =[]
    temporarym = []
    kp = kpArr[i]
    km = kmArr[i]
    for i in range(len(kp)-1): 
        temporaryp.append((Ep[1]-Ep[0])*(kp[i+1]-kp[i]))
        temporarym.append((Ep[1]-Ep[0])*(km[i+1]-km[i]))
    dosm.append(temporarym)
    dosp.append(temporaryp)
    print(h)

plt.figure()
plt.plot(Ep[0:len(Ep)-1], dosm[-1])
plt.plot(Ep[0:len(Ep)-1], dosp[-1])

# make full DOS
dos = []
for i in range(len(H)): 
    mm = dosm[i]
    pp = dosp[i]
    tempdos = [m+p for m,p in zip(mm, pp)]
    dos.append(tempdos)

plt.figure()
plt.hlines(y = T, xmin = kArr[0], xmax = 1e-6, colors = 'black')
n = len(H)
colors = plt.cm.jet(np.linspace(0,1,n))
for i in np.arange(0, len(H), 2): 
    plt.plot( dos[i], Ep[0:len(Ep)-1], label = str(np.round(H[i], decimals=2)), color = colors[i])
plt.legend()

# now we intgrate DOS up to the temperature
available_states = []
for i in range(len(H)): 
    # loop through energy, add if below temperature
    temp = 0
    for e, j  in zip(Ep, range(len(Ep))): 
        if e<= T: 
            temp = temp+dos[i][j]
    available_states.append(temp)

plt.figure()
plt.plot(H, available_states)



plt.figure()
idx = [40, 45, 46, 47] 
for i in idx: 
    plt.plot( dos[i], Ep[0:len(Ep)-1], label = str(np.round(H[i], decimals=2)))
plt.legend()