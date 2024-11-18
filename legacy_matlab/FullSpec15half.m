enableFullSpecNorm = true; %Plots a, b, and c field orientations, normalized to E0
enabledoubleplot = false; %Plots two sets of B params in the same plot for comparison

%% Spin 15/2 matrices, {|15/2,15/2>.,...,|15/2,-15/2>| basis

ID = eye(16);

% Create a vector with values from 15/2 down to -15/2 in steps of -1
Jz_val = 15/2:-1:-15/2;

% Create the diagonal matrix
Jz = diag(Jz_val);
    
% Create the vector for the superdiagonal
Jp_val = [sqrt(15), sqrt(28), sqrt(39), sqrt(48), sqrt(55), sqrt(60), sqrt(63), 8, ...
               sqrt(63), sqrt(60), sqrt(55), sqrt(48), sqrt(39), sqrt(28), sqrt(15)];

% Create the matrix Jp
Jp = diag(Jp_val, 1);  % The second argument '1' places the values on the first superdiagonal
    
    
% Define the vector for the subdiagonal
Jm_val = sqrt([1*15, 2*14, 3*13, 4*12, 5*11, 6*10, 7*9, 8*8, 7*9, 6*10, 5*11, 4*12, 3*13, 2*14, 1*15]);

% Create the matrix Jm
Jm = diag(Jm_val, -1);  % The '-1' places the values on the first subdiagonal

Jx = (Jp+Jm)/2;

Jy = (Jp-Jm)/2i;


%% Steven's Operators
J=15/2; X = J*(J+1); A = Jp*Jp*Jp + Jm*Jm*Jm;

O20 = 3*Jz*Jz - X*ID;
O22 = (1/2)*(Jp*Jp + Jm*Jm);
O40 = 35*power(Jz,4) - (30*X - 25)*Jz*Jz + (3*X*X - 6*X)*ID;
O42 = (1/4)*( (Jp*Jp + Jm*Jm)*(7*power(Jz,2) - (X + 5)*ID) + (7*power(Jz,2) - (X + 5)*ID)*(Jp*Jp + Jm*Jm));  
O44 = (1/2)*(Jp*Jp*Jp*Jp + Jm*Jm*Jm*Jm);
O43 = (1/4)*( (A)*Jz + Jz*(A) );
O60 = 231*power(Jz,6) - (315*X-735)*power(Jz,4) + (105*X*X - 525*X +294)*power(Jz,2) - (5*X*X*X - 40*X*X + 60*X)*ID;
O62 = (1/4)*((Jp*Jp + Jm*Jm)*(33*power(Jz,4) - (18*X + 123)*power(Jz,2)+(X*X+10*X+102)*ID) + (33*power(Jz,4) - (18*X + 123)*power(Jz,2)+(X*X+10*X+102)*ID)*(Jp*Jp + Jm*Jm));
O63 = (1/4)*( A*(11*power(Jz,3) - (3*X + 59)*Jz ) + (11*power(Jz,3) -(3*X + 59)*Jz)*A );
O64 = (1/4)*((Jp*Jp*Jp*Jp + Jm*Jm*Jm*Jm)*(11*power(Jz,2) - (X + 38)*ID) + (11*power(Jz,2) - (X + 38)*ID)*(Jp*Jp*Jp*Jp + Jm*Jm*Jm*Jm));
O66 = (1/2)*(Jp*Jp*Jp*Jp*Jp*Jp + Jm*Jm*Jm*Jm*Jm*Jm);

% Change B param values as needed 
%params = [-2.773e-2    -3.987e-4    -1.416e-2    3.152e-6    -7.616e-6   3.275e-5];
params = [6.741e-2  1.363e-3  -8.998e-3   9.565e-6  1.113e-4  1.661e-4]; %ETO Petit
params_1 = [7.50e-2 1.41e-3  1.25-2  1.09e-5 -1.8e-4 1.5e-4]; %ETO Bertin

B20_1 = params_1(1);
B40_1 = params_1(2);
B43_1 = params_1(3);
B60_1 = params_1(4);
B63_1 = params_1(5);
B66_1 = params_1(6);

B20 = params(1);
B40 = params(2);
B43 = params(3);
B60 = params(4);
B63 = params(5);
B66 = params(6);

HCEF_1 = B20_1*O20 + B40_1*O40 + B43_1*O43 + B60_1*O60 + B63_1*O63 + B66_1*O66;
HCEF = B20*O20 + B40*O40 + B43*O43 + B60*O60 + B63*O63 + B66*O66;

%% This is specifically for EGO which has 9 B params
params_EGO = [0.1271, -0.4371, 6.6574e-4, -0.0017, 0.0033, 1.0300e-5, 9.0100e-5, 5.0300e-5, -8.5100e-6];

B20_EGO = params_EGO(1);
B22_EGO = params_EGO(2);
B40_EGO = params_EGO(3);
B42_EGO = params_EGO(4);
B44_EGO = params_EGO(5);
B60_EGO = params_EGO(6);
B62_EGO = params_EGO(7);
B64_EGO = params_EGO(8);
B66_EGO = params_EGO(9);

HCEF_EGO = B20_EGO*O20 + B22_EGO*O22 + B40_EGO*O40 + B42_EGO*O42 + B44_EGO*O44 + B60_EGO*O60 + B62_EGO*O62 + B64_EGO*O64 + B66_EGO*O66;
%% Constants
muB = 5.78838e-2; % [meV/T];
kB  = 8.617e-2  ; % [meV/K];
gJ = 1.2; % L=7, S=1/2, J=15/2 g-lande factor;
%%
if enableFullSpecNorm
    % This first part is used to calculate gxx0 and gzz0 values, which are
    % only used in the plot title
    [P,D] = eig(HCEF_EGO + Jz*1e-10);

    En = real([D(1,1) D(2,2) D(3,3) D(4,4) D(5,5) D(6,6) D(7,7) D(8,8) D(9,9) D(10,10) D(11,11) D(12,12) D(13,13) D(14,14) D(15,15) D(16,16)]);
    Ev = sort(En);

    E0 = Ev(1);

    ind1 = find(En==Ev(1));
    ind2 = find(En==Ev(2));

    ev1 = P(:,ind1);
    ev2 = P(:,ind2);

    gzz0 = 2*gJ*abs(ev1'*Jz*ev1); 
    gxx0 = 2*gJ*abs(ev1'*Jx*ev2);

    H = 0:1:20; %Define the magnetic field range

    n = 1;
    for h = H
        HZeeman = -gJ*muB*h*Jx; %For ||a orientation
        [P,D] = eig(HCEF_EGO + HZeeman);

        for m = 1:16 %16 eigenvalues in spin 15/2
            Eigenval(m) = real(D(m,m));
        end
        Eigenval = sort(Eigenval);
        for m = 1:16
            E1{m}(n) = Eigenval(m);
        end
        
        HZeeman = -gJ*muB*h*Jy; %For ||b orientation
        [P,D] = eig(HCEF_EGO + HZeeman);

        for m = 1:16 %16 eigenvalues in spin 15/2
            Eigenval(m) = real(D(m,m));
        end
        Eigenval = sort(Eigenval);
        for m = 1:16
            E2{m}(n) = Eigenval(m);
        end
        
        HZeeman = -gJ*muB*h*Jz; %for ||c orientation
        [P,D] = eig(HCEF_EGO + HZeeman);

        for m = 1:16
            Eigenval(m) = real(D(m,m));
        end
        Eigenval = sort(Eigenval);
        for m = 1:16
            E3{m}(n) = Eigenval(m);
        end

        n = n + 1;
    end
    
    % Plot figure of ||a field orientation
    figure; hold on;
    colz = {'k','k','r','r','b','b','c','c', 'k','k','r','r','b','b','c','c'};
    for m = 1:16
        E1cm = (E1{m}) * 8.06554;
        E1zcm = (E1{1}) * 8.06554;
      
        plot(H,E1cm-E1zcm,'color',colz{m})
      

    end

    ylabel('Wavenumber (cm^{-1})'); xlabel('\mu_0H [T]'); grid on; box on;
    title(['H||a,   g_{||} = ' num2str(gxx0)])
    
    % Plot figure of ||b orientation
    figure; hold on;
    colz = {'k','k','r','r','b','b','c','c', 'k','k','r','r','b','b','c','c'};
    for m = 1:16
        E2cm = (E2{m}) * 8.06554;
        E2zcm = (E2{1}) * 8.06554;
      
        plot(H,E2cm-E2zcm,'color',colz{m})
      
        % Extract and display the values at H = 0
        E2cm_at_H0 = (E2{m}(1)) * 8.06554;  % Convert to wavenumbers
        E2zcm_at_H0 = (E2{1}(1)) * 8.06554;  % Convert to wavenumbers
        difference_at_H0 = E2cm_at_H0 - E2zcm_at_H0;

       % disp(['E2cm - E2zcm at H = 0: ', num2str(difference_at_H0)]);
    end

    ylabel('Wavenumber (cm^{-1})'); xlabel('\mu_0H [T]'); grid on; box on;
    title('H||b,   g_{||} = ')
    
    % Plot figure of ||c field orientation
    figure; hold on;
    colz = {'k','k','r','r','b','b','c','c', 'k','k','r','r','b','b','c','c'};
    for m = 1:16
        E3cm = (E3{m}) * 8.06554;
        E3zcm = (E3{1}) * 8.06554;
      
        plot(H,E3cm-E3zcm,'color',colz{m})
    end

    ylabel('Wavenumber (cm^{-1})'); xlabel('\mu_0H [T]'); grid on; box on;
    title(['H||c,   g_{||} = ' num2str(gxx0)])
end

% Check if normalization and double plotting are enabled
if enableFullSpecNorm && enabledoubleplot
    % Initialize plot
    figure;
    hold on;  % Hold on to plot multiple datasets on the same plot

    % Define colors for clarity in plotting
    colors = {'r', 'b'};  % 'r' for the first dataset, 'b' for the second

    % Calculate and plot for the first set of parameters
    HCEF = B20*O20 + B40*O40 + B43*O43 + B60*O60 + B63*O63 + B66*O66;
    plotEigenSpectrum(HCEF, Jz, Jx, colors{1});

    % Calculate and plot for the second set of parameters
    HCEF_1 = B20_1*O20 + B40_1*O40 + B43_1*O43 + B60_1*O60 + B63_1*O63 + B66_1*O66;
    plotEigenSpectrum(HCEF_1, Jz, Jx, colors{2});

    % Finalize plot
    ylabel('Wavenumber (cm^{-1})');
    xlabel('\mu_0H [T]');
    grid on;
    box on;
    title(['H||ab, g_{||} comparison']);
    legend('Params set 1', 'Params set 2');
    hold off;
end

function plotEigenSpectrum(HCEF, Jz, Jx, plotColor)
    muB = 5.78838e-2; % [meV/T]
    gJ = 1.2; % Lande g-factor
    H = 0:1:20; % Magnetic field range
    
    E2 = zeros(16, length(H));
   
    for n = 1:length(H)
    HZeeman = -gJ * muB * H(n) * Jx;
    [P, D] = eig(HCEF + HZeeman);
    E2(:, n) = sort(real(diag(D)));
    end

    % Convert energies to wavenumbers and plot
    E2cm = E2 * 8.06554; % Conversion factor
    E2z = E2(1, :) * 8.06554;
    plot(H, E2cm-E2z, plotColor); 
end
