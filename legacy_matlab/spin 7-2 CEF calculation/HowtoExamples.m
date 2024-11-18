

%% Paramaters
clear 
Bmn = [-0.4233    0.0117    0.5494    0.00035    0.0052   -0.00045];
Jxx = 0.5388; 
Jzz = 0.6145;
X4 = [9.95 Bmn Jxx Jzz];
fac = [1,10,1,1e3,1e2,1e2];

%% Plot the single ion spectrum
CEF_Spectrum_D3d('full spectrum',Bmn)

%% 
CEF_Spectrum_D3d('Finite Field Mag',Bmn)


%% Plot the MF susceptiblity
CEF_Magnetotropic_Modeling('Chi v T', X4);

%% Plot the MF susceptiblity faster (confusingly the parameters are scaled by fac in this version)
CEF_Magnetotropic_Modeling_Fast('DC Susc', [Bmn.*fac Jxx Jzz]);

%% Plot kab, kc
CEF_Magnetotropic_Modeling_Fast('principle k', [Bmn.*fac Jxx Jzz]);