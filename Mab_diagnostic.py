# testing convergence

def MolecularFieldTheory_diagnostic(H, Hth, Mth, lamb):
    '''Use mean field theory to correct magnetization
    for an exchange strength lamb. H is the number of fields to calculate,
    Hth and Mth are the theoretical single-ion magnetization curves to correct.'''
    n = 50
    Marr = []
    newM = np.interp(H, Hth, Mth)
    colors = plt.cm.cool(np.linspace(0,1,n))
    plt.figure()
    plt.plot(Hth, Mth, 'r--', label = 'input curve')
    MvalArr =[]
    HvalArr =[]
    for i in range(n):
        idx = findIdx(H, 0.1)
        temp = newM[idx]
        tempH = H[idx] +6*lamb*newM[idx]/muB/(gJ)**2
        MvalArr.append(temp)
        HvalArr.append(tempH)
        newH = H + 6*lamb*newM/muB/(gJ)**2
        newM = np.interp(newH,Hth,Mth)
        Marr.append(newM)
        plt.plot(H, newM, 'o', label = str(n), color = colors[i])
    plt.grid(True)
    plt.title(' M vs H for n='+str(n))
    plt.xlabel('H (T)')
    plt.ylabel('M (\mu_B /Er')
    plt.figure()
    plt.plot(range(n), MvalArr)
    plt.axhline(y = max(MvalArr), xmin = 0, xmax = n)
    plt.axhline(y = min(MvalArr), xmin = 0, xmax = n)
    plt.title('New M convergence issue at 0.1T, 0.1K')
    plt.xlabel('iteration number')
    plt.ylabel('New M')
    plt.grid(True)
    plt.figure()
    plt.plot(range(n), HvalArr)
    plt.axhline(y = max(HvalArr), xmin = 0, xmax = n)
    plt.axhline(y = min(HvalArr), xmin = 0, xmax = n)
    plt.title('New H convergence issue at 0.1T, 0.1K')
    plt.xlabel('iteration number')
    plt.ylabel('New H')
    plt.grid(True)
    return Marr

temperature =0.1
H = np.concatenate((np.linspace(0,0.2,100), np.linspace(0.2, .4, 100)))
Hth= np.linspace(-0.4,0.4, 1000)#[0,0.1,0.2,0.3,0.4]
field = [[0,b,0] for b in Hth]
magMe = np.array([MyErObj.magnetization(ion, temperature, f) for f in field]).T
Mth = -magMe[1]
tempA = MolecularFieldTheory_diagnostic(H, Hth, -1*magMe[1], Jperp)
# tempB = MolecularFieldTheory(MFTField, magF, -1*magMe[1], Jperp/2)
# tempMagA.append(tempA) 
# tempMagB.append(tempB)