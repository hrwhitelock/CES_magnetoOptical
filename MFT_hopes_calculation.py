## scratch MFT  code
# rlly annoyed that all of my concerns were like real lol

def bmag(mag, J,  h, temperature):
    newh = h +q*J*mag/muB/1.2/1.2
    mag = -1*MyErObj.magnetization(ion, temperature, [0,newh, 0]).T[1]-mag
    return mag

def cmag(mag, J,  h, temperature):
    newh = h +q*J*mag/muB/1.2/1.2
    mag = -1*MyErObj.magnetization(ion, temperature, [0,0, newh]).T[2]-mag
    return mag

def MFTmagC(ionObj, H, J, temperature): 
    # okay so let's start by defining the MF ham
    # only doing c direction rn because that is faster
    q = 6
    ion = 'Er3+' # hard coded for now dont' @ me
    n = 10 # iterations
    # first step, call the pcf mag because that's a good place to start
    magArr = []
    for h in H: 
        mag = 0
        for i in range(n): 
            mag = fsolve(cmag, mag, args = (J, h, temperature))[0] # fsolve spits out an array - this time its one val
        magArr.append(mag) # fso
    return magArr

def MFTmagB(ionObj, H, J, temperature): 
    # okay so let's start by defining the MF ham
    # only doing c direction rn because that is faster
    q = 6
    ion = 'Er3+' # hard coded for now dont' @ me
    n = 10 # iterations
    # first step, call the pcf mag because that's a good place to start
    magArr = []
    for h in H: 
        mag = 0
        for i in range(n): 
            mag = fsolve(bmag, mag, args = (J, h, temperature))[0] # fsolve spits out an array - this time its one val
        magArr.append(mag) # fso
    return magArr



temps = [0.025, 0.045,0.1, 0.171, .25, .35, .45, .543, .827, 2, 6, 20]
labels = [str(T)+ 'K' for T in temps]
n = len(temps)
colors = plt.cm.jet(np.linspace(0,1,n))
tempMagB = []
tempMagC = []
H = np.linspace(0,10, 100)
Jz = .48e-3
Jperp = -.9e-3

# H = 1
# temperature = 1
# temp = MFTmag(MyErObj, H, Jz, temperature)

temps = [1]
tempMagB =[]
for temperature in temps: 
    temp = MFTmagB(MyErObj, H, -Jperp, temperature)
    tempMagB.append(temp)
    temp = MFTmagC(MyErObj, H, Jz, temperature)
    tempMagC.append(temp)

temp = MFTmagB(MyErObj, H, Jperp, temperature)

plt.figure()
for mag, color, label in zip(tempMagB, colors, labels): 
    plt.plot(H, mag, '--o', color = color, label = label)
plt.plot(H, temp, label = 'J = -.9ueV')
plt.title('B axis MFT magnetization \n calcualted from my mft function \n Jperp = '+str(Jperp))
plt.ylabel('B axis MFT magnetization')
plt.xlabel('Field (T)')
plt.legend()
plt.grid(True)

plt.figure()
for mag, color, label in zip(tempMagC, colors, labels): 
    plt.plot(H, mag, '--o', color = color, label = label)
plt.title('C axis MFT magnetization \n calcualted from my mft function')
plt.ylabel('C axis MFT magnetization')
plt.xlabel('Field (T)')
plt.legend()


##########################################################################################################
##########################################################################################################
##########################################################################################################
##########################################################################################################
# okay so now let's calculate for real