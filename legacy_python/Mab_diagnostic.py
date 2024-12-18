# testing convergence

def MolecularFieldTheory_diagnostic(H, Hth, Mth, lamb):
    '''Use mean field theory to correct magnetization
    for an exchange strength lamb. H is the number of fields to calculate,
    Hth and Mth are the theoretical single-ion magnetization curves to correct.'''
    n = 50
    Marr = []
    newM = np.interp(H, Hth, Mth)
    colors = plt.cm.cool(np.linspace(0,1,n))
    # plt.figure()
    plt.plot(Hth, Mth, 'r--', label = 'input curve')
    MvalArr =[]
    HvalArr =[]
    for i in range(n):
        idx = findIdx(H, 0.1)
        newH = H + 6*lamb*newM[idx]/muB/(gJ)**2
        newM = np.interp(newH,Hth,Mth)
        Marr.append(newM)
        plt.plot(H, newM, 'o', label = str(n), color = colors[i])
    plt.grid(True)
    plt.title(' M vs H for n='+str(n))
    plt.xlabel('H (T)')
    plt.ylabel('M (\mu_B /Er')
    return newM

temperature =0.1
Jperp = -0.9e-3
H = np.linspace(-1, 1)
Hth= np.linspace(-1,1, 1000)#[0,0.1,0.2,0.3,0.4]
field = [[0,b,0] for b in Hth]
magMe = np.array([MyErObj.magnetization(ion, temperature, f) for f in field]).T
Mth = -magMe[1]
tempA = MolecularFieldTheory_diagnostic(H, Hth, -1*magMe[1], Jperp)
# tempB = MolecularFieldTheory(MFTField, magF, -1*magMe[1], Jperp/2)
# tempMagA.append(tempA) 
# tempMagB.append(tempB)