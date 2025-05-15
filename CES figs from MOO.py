## CES figs from MOO

# # to init my terminal
B20 =   -0.03721092
B40 =   -0.00038796
B43 =   -0.01406804
B60 =    3.1865e-06
B63 =   -3.593e-06
B66 =    3.4913e-05
Jperp = -.53070e-03 # +/- 2.6332e-06 (0.50%) (init = 0)
Jz =    -2.63253e-03



# make my er obj
Bparams =  {'B20': B20, 'B40':B40,'B43': B43, 'B60': B60, 'B63':B63,'B66':B66}
moo = MOO('Er3+', Bparams, Jperp, Jz, q =6)


temps = [0.025, 0.045,0.1, 0.171, .25, .35, .45, .543, .827, 1.008, 2, 6, 20]
labels = [str(T)+ 'K' for T in temps]
n = len(temps)
colors = plt.cm.jet(np.linspace(0,1,n))
tempMagA = []
tempMagB = []
tempMagC = []
# H = np.concatenate((np.linspace(0,1,50), np.linspace(1.01,15, 150)))
H = np.linspace(0,8,8*5)

res = []
for T in temps:
    temp_res = moo.calc_over_field(T, H)
    res.append(temp_res)

plt.figure()
n = len(temps)
colors = plt.cm.jet(np.linspace(0,1,n))
dmdh = []
for i, T in enumerate(temps): 
    n = res[i]
    plt.plot(n['H'], n['m'].T[2], label = str(T)+'K', color = colors[i]) 
    dmdh.append(np.gradient(n['m'].T[2], n['H']))

plt.figure()
for i, T in enumerate(temps): 
    n = res[i]
    plt.plot(n['H'], dmdh[i], label = str(T)+'K', color = colors[i]) 

## redo splitting 
H_2 = np.linspace(0,20,40)
# res2 = []
temps_2 = [0.025, 0.5, 1, 4,10]
for T in temps:
    temp_res = moo.calc_over_field(T, H)
    res2.append(temp_res)


plt.figure()
n = len(temps)
colors = plt.cm.jet(np.linspace(0,1,n))
for i, T in enumerate(temps): 
    n = res2[i]
    evals = np.array(n['evals']).T
    for j in range(len(evals)):
        plt.plot(n['H'], evals[j], label = str(T)+'K', color = colors[i]) 

## save as hdf5

with h5py.File('moo_results.h5', 'w') as f:
    # First dataset: res (per-temperature magnetization)
    grp1 = f.create_group('res')
    grp1.create_dataset('temps', data=np.array(temps))
    # For each temperature, save H and m_z (c-axis magnetization) and store gradient
    for i, T in enumerate(temps):
        gi = grp1.create_group(f'T_{i}')
        H = res[i]['H']
        mz = res[i]['m'][:,2]
        dmdh = np.gradient(mz, H)
        gi.create_dataset('H', data=H)
        gi.create_dataset('mz', data=mz)
        gi.create_dataset('dmdh', data=dmdh)

    # Second dataset: res2 (per-temperature CEF splitting)
    grp2 = f.create_group('res2')
    grp2.create_dataset('temps2', data=np.array(temps_2))
    for i, T in enumerate(temps_2):
        gi = grp2.create_group(f'T_{i}')
        H2 = res2[i]['H']
        evals = np.array(res2[i]['evals'])  # shape (Nfields x Nlevels)
        gi.create_dataset('H', data=H2)
        gi.create_dataset('evals', data=evals)



## hmmm, so what now? let's calc at a couple different temps, and then w/ no J
temps = [0.001, 0.025, 0.045,0.1, 0.171, .25, .35, .45, .543, .827, 1.008, 2, 6, 20]
labels = [str(T)+ 'K' for T in temps]
n = len(temps)
colors = plt.cm.jet(np.linspace(0,1,n))
H = np.concatenate((np.linspace(0,1,50), np.linspace(1.01,10, 100)))
# H = np.linspace(0,10,200)

res = []
for T in temps:
    temp_res = moo.calc_over_field(T, H)
    res.append(temp_res)
# # now calc no mft
# moo.JJz = 0
# nomft = moo.calc_over_field(0.1,H)

plt.figure()
n = len(temps)
colors = plt.cm.jet(np.linspace(0,1,n))
for i, T in enumerate(temps): 
    result = res[i]['x']
    # first plot magnetization 
    mx = result['m'].T[0]
    plt.plot(result['H'], mx, label = str(T)+'K', color= colors[i])
plt.title('H||ab')
plt.grid(True)
plt.xlabel('H [T]')
plt.ylabel('M [uB/Er]')

# now plot evals
plt.figure()
n = len(temps)
colors = plt.cm.jet(np.linspace(0,1,n))
for i, T in enumerate(temps): 
    result = res[i]['x']
    evals = np.array(result['evals']).T
    for j in range(len(evals)):
        plt.plot(result['H'], evals[j], label = str(T)+'K', color = colors[i]) 
plt.title('H||ab') 
plt.xlabel('H [T]')
plt.ylabel('E [meV]')

# now plot H vs Heff
plt.figure()
n = len(temps)
colors = plt.cm.jet(np.linspace(0,1,n))
for i, T in enumerate(temps): 
    result = res[i]['x']
    # first plot magnetization 
    mx = result['m'].T[0]
    H_eff = np.array(result['H']) + moo.q*moo.JJperp*mx/(muB*moo.gL**2)
    plt.plot(result['H'], H_eff, label = str(T)+'K', color= colors[i])
plt.title('H||ab')
plt.grid(True)
plt.xlabel('H [T]')
plt.ylabel('H_eff [T]')

## make all figs for c axis
plt.figure()
n = len(temps)
colors = plt.cm.jet(np.linspace(0,1,n))
for i, T in enumerate(temps): 
    result = res[i]['z']
    # first plot magnetization 
    mz = result['m'].T[2]
    plt.plot(result['H'], mz, label = str(T)+'K', color= colors[i])
plt.title('H||c')
plt.grid(True)
plt.xlabel('H [T]')
plt.ylabel('M [uB/Er]')

# now plot evals
plt.figure()
n = len(temps)
colors = plt.cm.jet(np.linspace(0,1,n))
for i, T in enumerate(temps): 
    result = res[i]['z']
    evals = np.array(result['evals']).T
    for j in range(len(evals)):
        plt.plot(result['H'], evals[j], label = str(T)+'K', color = colors[i]) 
plt.title('H||c') 
plt.xlabel('H [T]')
plt.ylabel('E [meV]')

# now plot H vs Heff
plt.figure()
n = len(temps)
colors = plt.cm.jet(np.linspace(0,1,n))
for i, T in enumerate(temps): 
    result = res[i]['z']
    # first plot magnetization 
    mz = result['m'].T[2]
    H_eff = [h + moo.q*moo.JJz*m/(muB*moo.gL**2) for h, m in zip(result['H'], mz)]
    plt.plot(result['H'], H_eff, label = str(T)+'K', color= colors[i])
plt.title('H||c')
plt.grid(True)
plt.xlabel('H [T]')
plt.ylabel('H_eff [T]')


# let's plot the eigenvectors

names = ['15/2', '13/2', '11/2', '9/2', '7/2', '5/2', '3/2', '1/2', 
         '-1/2', '-3/2', '-5/2', '-7/2', '-9/2', '-11/2', '-13/2', '-15/2']
n_temps = len(temps)
colors = plt.cm.jet(np.linspace(0, 1, n_temps))

for state_index in range(16):  # Loop over eigenstates
    fig, axs = plt.subplots(4, 4, figsize=(12, 10), sharex=True, sharey=True)
    axs = axs.flatten()

    for basis_index in range(16):  # Loop over spin basis states
        ax = axs[basis_index]

        for i, T in enumerate(temps):  # Loop over temperatures
            result = res[i]['z']
            Hvals = result['H']
            evecs = result['evecs']  # List of eigenvector matrices (shape: [N_H, 16, 16])
            
            # Pull out the `state_index`-th eigenvector at each field, and get weight on `basis_index`
            weights = [np.abs(evec[:, state_index][basis_index])**2 for evec in evecs]
            ax.plot(Hvals, weights, color=colors[i], label=f'{T}K' if basis_index == 0 else "")

        ax.set_title(names[basis_index], fontsize=8)
        ax.grid(True)
        if basis_index % 4 == 0:
            ax.set_ylabel('Weight')
        if basis_index // 4 == 3:
            ax.set_xlabel('H [T]')

    fig.suptitle(f'H || c: Eigenstate {state_index} composition', fontsize=14)
    fig.tight_layout()
    fig.subplots_adjust(top=0.92, bottom=0.12)
    fig.legend(loc='lower center', ncol=6)
    plt.show()

with h5py.File('moo_3d_result.h5', 'w') as f:
    # Save temperature and field arrays
    f.create_dataset('temps', data=np.array(temps))
    f.create_dataset('H', data=np.array(H))
    
    # For each direction ('x','z'), save subgroups
    for dir_key in ('x','y','z'):
        gdir = f.create_group(dir_key)
        for i, T in enumerate(temps):
            grp = gdir.create_group(f'T_{i}')
            result = res[i][dir_key]
            grp.create_dataset('m', data=result['m'])         # Nx3 magnetizations
            grp.create_dataset('evals', data=result['evals']) # NxK eigenvalues
            grp.create_dataset('evecs', data=result['evecs'])
            # Compute H_eff on the fly if desired, or save true H (same)
            # H_eff depends on direction:
            m = result['m']
            if dir_key=='x':
                mx = m[:,0]
                Heff = H + moo.q*moo.JJperp*mx/(muB*moo.gL**2)
            elif dir_key=='y':
                mx = m[:,1]
                Heff = H + moo.q*moo.JJperp*mx/(muB*moo.gL**2)
            else:
                mz = m[:,2]
                Heff = H + moo.q*moo.JJz  *mz/(muB*moo.gL**2)
            grp.create_dataset('Heff', data=Heff)


# rq let's look at the toy model

def halfJdotB(h): 
    zb = muB * 2 # g factor is like 2 for an electron
    op = cef.Operator(J=1/2)
    JdotB = zb*h*op.Jz(1/2).O
    return JdotB

def _expectation_m(evals, evecs, T):
    beta = 1/(kB*T)
    w = np.exp(-beta*evals)
    Z = w.sum()
    Jexp = np.zeros((len(evals),3))
    for i, vec in enumerate(evecs):
        ket = cef.Ket(vec)
        Jexp[i] = [np.real(ket*ket.Jx()),
                    np.real(ket*ket.Jy()),
                    np.real(ket*ket.Jz())]
    Javg = (w[:,None]*Jexp).sum(axis=0)/Z
    return -2 * Javg


def self_const_diag(T, H):
    hx = hy = 0
    hz = H
    muBg2 = muB*2**2
    # locals for speed
    diag = ionObj.diagonalize
    expect = _expectation_m
    # JJz = -.5

    # root functions
    def fix_m(i, m, h0, JJc):
        def f(m_val):
            h_eff = h0 + q*JJc*m_val/muBg2

            evals, evecs = np.linalg.eig(halfJdotB(h_eff))
            return expect(evals, evecs, T)[i] - m_val
        return f

    # initial guesses
    mx = my = mz = 0.0
    # iterate
    test_mz = []
    for _ in range(10):
        mz = fsolve(fix_m(2, mz, hz, JJz), mz)[0]
        test_mz.append(mz)
    # plt.figure()
    # plt.plot(range(10), test_mz, 'o')
    # final diag
    newH = hz + q*JJz*mz/muBg2
    evals, evecs = np.linalg.eig(halfJdotB(newH))
    return (mx,my,mz), evals, evecs


# temps = [0.001, 0.025, 0.045,0.1, 0.171, .25, .35, .45, .543, .827, 1.008, 2, 6, 20]
half_temps = [.1, .5, 1,2,5,10,100, 500 ]
labels = [str(T)+ 'K' for T in temps]
n = len(half_temps)
colors = plt.cm.jet(np.linspace(0,1,n))
field = np.linspace(0,1000, 1000)
# H = np.linspace(0,10,200)

JJz = -1
m_temp_arr = []
evals_temp_arr = []
plt.figure()
for i, T in enumerate(half_temps):
    m_arr = []
    evals_arr = []
    for h in field: 
        m, evals, evecs = self_const_diag(T, h)
        m_arr.append(m)
        evals_arr.append(evals)
    plt.figure()
    plt.plot(field, m_arr, color = colors[i])
    m_temp_arr.append(m_arr)
    evals_temp_arr.append(evals_arr)



# now let's plot
plt.figure()
for i, T in enumerate(half_temps): 
    m = np.array(m_temp_arr[i]).T
    plt.plot(field, m[2], color = colors[i])


plt.figure()
for i, T in enumerate(half_temps): 
    evals = np.array(evals_temp_arr[i]).T
    plt.plot(field, evals[0], color = colors[i])
    plt.plot(field, evals[1], color = colors[i])







