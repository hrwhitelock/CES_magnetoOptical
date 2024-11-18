function output = CEF_Magnetotropic_Modeling_Fast(option, params)


%% Fixed angle/Temp ZF resonance frequencies
amp = 9.95;
T_RTM = [1.5,4,6,12,20, 30,50,70];
w0 = [49.5241, 49.5276, 49.5238, 49.5215, 49.5149, 49.4941, 49.4440, 49.3838];

%% 2nd Deriv coefficients
C1 = [1/12, 	-2/3, 	0,      2/3, 	-1/12];
C2 = [-1/12 	4/3 	-5/2 	4/3 	-1/12];

%% Spin 7/2 matrices, {|7/2,7/2>.,...,|7/2,-7/2>| basis

ID = [1,0,0,0,0,0,0,0;
    0,1,0,0,0,0,0,0;
    0,0,1,0,0,0,0,0;
    0,0,0,1,0,0,0,0;
    0,0,0,0,1,0,0,0;
    0,0,0,0,0,1,0,0;
    0,0,0,0,0,0,1,0;
    0,0,0,0,0,0,0,1];

Jz = [ 7/2,   0,   0,   0,   0,   0,   0,   0;
    0, 5/2,   0,   0,   0,   0,   0,   0;
    0,   0, 3/2,   0,   0,   0,   0,   0;
    0,   0,   0, 1/2,   0,   0,   0,   0;
    0,   0,   0,   0,-1/2,   0,   0,   0;
    0,   0,   0,   0,   0,-3/2,   0,   0;
    0,   0,   0,   0,   0,   0,-5/2,   0;
    0,   0,   0,   0,   0,   0,   0,-7/2];

Jp =  [       0,sqrt( 7),       0,       0,       0,       0,       0,       0;
    0,       0,sqrt(12),       0,       0,       0,       0,       0;
    0,       0,       0,sqrt(15),       0,       0,       0,       0;
    0,       0,       0,       0,       4,       0,       0,       0;
    0,       0,       0,       0,       0,sqrt(15),       0,       0;
    0,       0,       0,       0,       0,       0,sqrt(12),       0;
    0,       0,       0,       0,       0,       0,       0,sqrt( 7);
    0,       0,       0,       0,       0,       0,       0,       0];

Jm =  [       0,       0,       0,       0,       0,       0,       0,       0;
    sqrt( 7),       0,       0,       0,       0,       0,       0,       0;
    0,sqrt(12),       0,       0,       0,       0,       0,       0;
    0,       0,sqrt(15),       0,       0,       0,       0,       0;
    0,       0,       0,       4,       0,       0,       0,       0;
    0,       0,       0,       0,sqrt(15),       0,       0,       0;
    0,       0,       0,       0,       0,sqrt(12),       0,       0;
    0,       0,       0,       0,       0,       0,sqrt( 7),       0];

Jx = (Jp+Jm)/2;

Jy = (Jp-Jm)/2i;

%% Steven's Operators
J=7/2; X = J*(J+1); A = Jp*Jp*Jp + Jm*Jm*Jm;

O20 = 3*Jz*Jz - X*ID;
O40 = 35*power(Jz,4) - (30*X - 25)*Jz*Jz + (3*X*X - 6*X)*ID;
O60 = 231*power(Jz,6) - (315*X-735)*power(Jz,4) + (105*X*X - 525*X +294)*power(Jz,2) - (5*X*X*X + 40*X*X -60*X)*ID;
O43 = (1/4)*( (A)*Jz + Jz*(A) );
O63 = (1/4)*( A*(11*power(Jz,3) - (3*X + 59)*Jz ) + (11*power(Jz,3) -(3*X + 59)*Jz)*A );
O66 = (1/2)*(Jp*Jp*Jp*Jp*Jp*Jp + Jm*Jm*Jm*Jm*Jm*Jm);

fac = [1,10,1,1e3,1e2,1e2];

HCEFf = @(Bmn) Bmn(1)*O20/fac(1) + Bmn(2)*O40/fac(2) + Bmn(3)*O43/fac(3)...
    + Bmn(4)*O60/fac(4) + Bmn(5)*O63/fac(5) + Bmn(6)*O66/fac(6);

muB = 5.78828e-2; % [meV/T];
kB  = 8.617e-2  ; % [meV/K];
gJ = 8/7; % L=3, S=1/2, J=7/2 g-lande factor;
A1 = 6*kB/(muB*gJ);
C0 = 2.0416;

%% Magnetization/Suscceptiblity data at fixed temperature (empirical form);
% Order is [Ma(4K), ... Ma(70K),  Mc(4K), ... Mc(70K),]

M_params = [0.118245775053758   0.118390928075276   3.145609938797031
   0.104837077149736   0.146091228622670   4.623002478131845
   0.054780448161023   0.568467821697383  12.717071808417087
   0.069143916974718   0.028051501174311   4.730178324588485
   0.056265911714515   0.003620474888152   4.140704470550365
   0.040781498724497                   0   1.000000000000000
   0.032655967521624                   0   1.000000000000000
   0.049439466905029   0.132365001450139   3.117775521471259
   0.048877317082171   0.125084661033167   4.172347457632269
   0.052891925660000   0.047475052836257   4.745946312072956
   0.053376191995845   0.003407661479687   1.650687863769174
   [0.049792977847172   0.002366684862856   0.798577722620777]/1.0374
   0.045978434856910                   0   1.000000000000000
   0.041776613518058                   0   1.000000000000000];

f_mag = @(p,x) p(1).*x + p(2).*tanh(x/p(3));

%% Fixed fields/tempratures for the full set of MvH interpolation

dH = .05; H_single = [0:dH:7+2*dH]; N_H = length(H_single);
Hset = kron(ones(14,1),H_single');
Tset = reshape( ([T_RTM(2:8) T_RTM(2:8)]'*ones(1,length(H_single)))', size(Hset));

Mset0 = [];
for i = 1:14
    Mset0 = [Mset0 ; f_mag(M_params(i,:),H_single')];
end

L1 = length(H_single);
N  = length(Hset); 


%% Curve Functions

    function out = Mab_vH_alt(HCEF,E0,Jxx,Jzz,t,H,dH)
        
        B = 1./(kB*t); hs = .01;
        n = 1; s = size(H);
        mx = zeros(s);
        mzT0 = zeros(s);
        mzT1 = zeros(s);
        
        for h = [H(1)-2*dH, H(1)-dH, H, H(end)+dH, H(end)+2*dH]
            
            HpZ = HCEF - gJ*muB*h*Jx;
            F1(n) = -(1./B)*log(trace(expm(-B*(HpZ - E0*ID)) ));
            n = n+1;
        end
        
        M1 = -(C1(1)*F1(1:end-4) + C1(2)*F1(2:end-3) + C1(3)*F1(3:end-2) + C1(4)*F1(4:end-1) + C1(5)*F1(5:end))/dH;
        
        %X1 = H + Jxx*A1*mx;
        out = M1;%Interp1NonUnique([0 X1], [0 gJ*muB.*mx], H);
    end

    function out = Mc_vH_alt(HCEF,E0,Jxx,Jzz,t,H,dH)
        
        B = 1./(kB*t); hs = .01;
        n = 1; s = size(H);
        mx = zeros(s);
        mzT0 = zeros(s);
        mzT1 = zeros(s);
        
        for h = [H(1)-2*dH, H(1)-dH, H, H(end)+dH, H(end)+2*dH]
            
            HpZ = HCEF - gJ*muB*h*Jz;
            F1(n) = -(1./B)*log(trace(expm(-B*(HpZ - E0*ID)) ));
            n = n+1;
        end
        
        M1 = -(C1(1)*F1(1:end-4) + C1(2)*F1(2:end-3) + C1(3)*F1(3:end-2) + C1(4)*F1(4:end-1) + C1(5)*F1(5:end))/dH;
        
        %X1 = H + Jxx*A1*mx;
        out = M1;%Interp1NonUnique([0 X1], [0 gJ*muB.*mx], H);
    end

    function out = Mab_vH(HCEF,E0,Jxx,Jzz,t,H)
        
        B = 1./(kB*t); hs = .01;
        n = 1; s = size(H);
        mx = zeros(s);
        mzT0 = zeros(s);
        mzT1 = zeros(s);
        
        for h = H
            
            HpZ = HCEF - gJ*muB*h*Jx;
            mx(n) = trace( Jx*expm(-B*(HpZ - E0*ID)) )./trace( expm(-B*(HpZ - E0*ID)));
            n = n+1;
        end
        
        X1 = H + Jxx*A1*mx;
        out = Interp1NonUnique([0 X1], [0 gJ*muB.*mx], H);
    end

    function out = Mc_vH(HCEF,E0,Jxx,Jzz,t,H)
        
        B = 1./(kB*t); hs = .01;
        n = 1; s = size(H);
        mz = zeros(s);
        mxT0 = zeros(s);
        mxT1 = zeros(s);
        
        for h = H
            
            HpZ = HCEF - gJ*muB*h*Jz;
            mz(n) = trace( Jz*expm(-B*(HpZ - E0*ID)) )./trace( expm(-B*(HpZ - E0*ID)));
            n = n+1;
        end
        
        X1 = H + Jzz*A1*mz;
        out = Interp1NonUnique([0 X1], [0 gJ*muB.*mz], H);
    end

    function out = kab_vH(HCEF,E0,Jxx,Jzz,t,H)
        
        B = 1./(kB*t); hs = .01;
        n = 1; s = size(H);
        mx = zeros(s);
        mzT0 = zeros(s);
        mzT1 = zeros(s);
        
        for h = H
            
            HpZ = HCEF - gJ*muB*h*Jx;
            mx(n) = trace( Jx*expm(-B*(HpZ - E0*ID)) )./trace( expm(-B*(HpZ - E0*ID)));
            mzT0(n) = trace( Jz*expm(-B*(HpZ - E0*ID)) )./trace( expm(-B*(HpZ - E0*ID)));
            
            HpZ = HpZ - gJ*muB*hs*Jz;
            mzT1(n) = trace( Jz*expm(-B*(HpZ - E0*ID)) )./trace( expm(-B*(HpZ - E0*ID)));
            
            n = n+1;
        end
        
        X1 = H + Jxx*A1*mx; X2 = power(X1,2);
        chiT = (mzT1-mzT0)./(hs + Jzz*A1*(mzT1 - mzT0));
        
        out = Interp1NonUnique([0 X1], [0 gJ*muB*(X1.*mx - X2.*chiT)], H);
    end

    function out = kc_vH(HCEF,E0,Jxx,Jzz,t,H)
        
        B = 1./(kB*t); hs = .01;
        n = 1; s = size(H);
        mz = zeros(s);
        mxT0 = zeros(s);
        mxT1 = zeros(s);
        
        for h = H
            
            HpZ = HCEF - gJ*muB*h*Jz;
            mz(n) = trace( Jz*expm(-B*(HpZ - E0*ID)) )./trace( expm(-B*(HpZ - E0*ID)));
            mxT0(n) = trace( Jx*expm(-B*(HpZ - E0*ID)) )./trace( expm(-B*(HpZ - E0*ID)));
            
            HpZ = HpZ - gJ*muB*hs*Jx;
            mxT1(n) = trace( Jx*expm(-B*(HpZ - E0*ID)) )./trace( expm(-B*(HpZ - E0*ID)));
            
            n = n+1;
        end
        
        X1 = H + Jzz*A1*mz; X2 = power(X1,2);
        chiT = (mxT1-mxT0)./(hs + Jxx*A1*(mxT1 - mxT0));
        
        out = Interp1NonUnique([0 X1], [0 gJ*muB*(X1.*mz - X2.*chiT)], H);
    end

    function out = XTc_vH(HCEF,E0,Jxx,Jzz,t,H)
        
        B = 1./(kB*t); hs = .01;
        n = 1; s = size(H);
        mx = zeros(s);
        mzT0 = zeros(s);
        mzT1 = zeros(s);
        
        for h = H
            
            HpZ = HCEF - gJ*muB*h*Jx;
            mx(n) = trace( Jx*expm(-B*(HpZ - E0*ID)) )./trace( expm(-B*(HpZ - E0*ID)));
            mzT0(n) = trace( Jz*expm(-B*(HpZ - E0*ID)) )./trace( expm(-B*(HpZ - E0*ID)));
            
            HpZ = HpZ - gJ*muB*hs*Jz;
            mzT1(n) = trace( Jz*expm(-B*(HpZ - E0*ID)) )./trace( expm(-B*(HpZ - E0*ID)));
            
            n = n+1;
        end
        
        X1 = H + Jxx*A1*mx; 
        chiT = (mzT1-mzT0)./(hs + Jzz*A1*(mzT1 - mzT0));
        
        out = Interp1NonUnique([0 X1], gJ*muB*[chiT(1) chiT], H);
    end

    function out = XTab_vH(HCEF,E0,Jxx,Jzz,t,H)
        
        B = 1./(kB*t); hs = .01;
        n = 1; s = size(H);
        mz = zeros(s);
        mxT0 = zeros(s);
        mxT1 = zeros(s);
        
        for h = H
            
            HpZ = HCEF - gJ*muB*h*Jz;
            mz(n) = trace( Jz*expm(-B*(HpZ - E0*ID)) )./trace( expm(-B*(HpZ - E0*ID)));
            mxT0(n) = trace( Jx*expm(-B*(HpZ - E0*ID)) )./trace( expm(-B*(HpZ - E0*ID)));
            
            HpZ = HpZ - gJ*muB*hs*Jx;
            mxT1(n) = trace( Jx*expm(-B*(HpZ - E0*ID)) )./trace( expm(-B*(HpZ - E0*ID)));
            
            n = n+1;
        end
        
        X1 = H + Jzz*A1*mz;
        chiT = (mxT1-mxT0)./(hs + Jxx*A1*(mxT1 - mxT0));
        
        out = Interp1NonUnique([0 X1], gJ*muB*[chiT(1) chiT], H);
    end

    function out = Xc_inv_vT(HCEF,E0,Jaa,T)
        
        hs = .01;
        HpZ = HCEF - gJ*muB*hs*Jz;
        
        n = 1; m = zeros(size(T));
        for t1 = T
            B = 1./(kB*t1);
            m(n) = trace( Jz*expm(-B*(HpZ - E0*ID)) )./trace( expm(-B*(HpZ - E0*ID)));
            n = n+1;
        end
        
        out = C0*((gJ*muB/kB)*hs./m + 6*Jaa);
        
    end

    function out = Xab_inv_vT(HCEF,E0,Jaa,T)
        
        hs = .01;
        HpZ = HCEF - gJ*muB*hs*Jx;
        
        n = 1; m = zeros(size(T));
        for t1 = T
            B = 1./(kB*t1);
            m(n) = trace( Jx*expm(-B*(HpZ - E0*ID)) )./trace( expm(-B*(HpZ - E0*ID)));
            n = n+1;
        end
        
        out = C0*((gJ*muB/kB)*hs./m + 6*Jaa);
        
    end

%% Minimization Function
XaDatInt = [10.1772   11.5966   15.7263   21.0341   27.3460   38.4183   47.9779   60.7196   80.8536  100.1618  119.1967];
%[60.7196   80.8536  100.1618  119.1967];
XcDatInt = [17.4148   20.1512   25.1159   28.4515   30.8158   34.0760   37.5031   43.6232   57.1081   71.9755   87.9933];
%[43.6232   57.1081   71.9755   87.9933];

%% best fit parameters
% LTparams = [3.533, 2.724, -0.01153,  3.987, 0.5996, -0.07036];
% LTparams = [3.533, 2.724, -0.01153,  3.958, 0.5931, -0.07042]; % for exchange consts (0.5388, 0.6145)
%LTparams = [3.569, 2.787, -0.009363, 3.991, 0.5981, -0.07052]; % for exchange consts (0.5948, 0.6095)
LTparams = [6*params(7), 2.774, -0.009435, 6*params(8), 0.6173,  -0.07029];

Ax0 = LTparams(2);
Bx0 = LTparams(3);

Az0 = LTparams(5);
Bz0 = LTparams(6);

gxx = 2*gJ*power(Ax0,1/2);
gzz = 2*gJ*power(Az0,1/2);
ax0 = power(muB*gJ,2)*Bx0/(2*kB);
az0 = power(muB*gJ,2)*-0.07052/(2*kB);
%
% gxx = 3.83;
% gzz = 1.748;
% az0 = -0.001796;
% ax0 = -0.00021;

    function out = CummErr(p0)
        
        
        HCEF1 = HCEFf(p0(1:6));
        Jxx1 = p0(7); Jzz1 = p0(8); chi0 = p0(9);

        [P1,D1] = eig(HCEF1 + Jz*1e-10);

        En = real([D1(1,1) D1(2,2) D1(3,3) D1(4,4) D1(5,5) D1(6,6) D1(7,7) D1(8,8)]);
        Ev1 = sort(En);
        E1 = Ev1(1);

        for n = 1:4
            ind1 = find(En==Ev1(2*(n-1) + 1));
            ind2 = find(En==Ev1(2*n));

            ev1{n} = P1(:,ind1(1));
            ev2{n} = P1(:,ind2(end));
        end

        gzz_m = 2*gJ*abs(ev1{1}'*Jz*ev1{1});
        gxx_m = 2*gJ*abs(ev1{1}'*Jx*ev2{1});

        % Ax1_m = power(ev1{2}'*Jx*ev2{2},2);
        % Az1_m = power(ev1{2}'*Jz*ev1{2},2);

        az0_m = -power(muB*gJ,2).*((power(ev1{1}'*Jz*ev1{2},2) + power(ev1{1}'*Jz*ev2{2},2) )./(Ev1(3) - E1)...
                                 +(power(ev1{1}'*Jz*ev1{3},2) + power(ev1{1}'*Jz*ev2{3},2) )./(Ev1(5) - E1)...
                                 +(power(ev1{1}'*Jz*ev1{4},2) + power(ev1{1}'*Jz*ev2{4},2) )./(Ev1(7) - E1));
%
        ax0_m = -power(muB*gJ,2).*((power(ev1{1}'*Jx*ev1{2},2) + power(ev1{1}'*Jx*ev2{2},2) )./(Ev1(3) - E1)...
                                  +(power(ev1{1}'*Jx*ev1{3},2) + power(ev1{1}'*Jx*ev2{3},2) )./(Ev1(5) - E1)...
                                  +(power(ev1{1}'*Jx*ev1{4},2) + power(ev1{1}'*Jx*ev2{4},2) )./(Ev1(7) - E1));

        out = 0;
                              
%         out = 500*(power(gxx_m - gxx,2) + power(gzz_m - gzz,2))...
%                 + 1e7*power(az0_m - az0,2) + 1e7*power(ax0_m - ax0,2);...
                %+ power(Az1_m - Az1,2) + power(Ax1_m - Ax1,2);
% 
         for n3 = 2:8
             out = out + sum(power(kd_int{1}{n3} - kab_vH(HCEF1,E1,Jxx1,Jzz1,T_RTM(n3),Hd_int{1}{n3}), 2));
             out = out + sum(power(kd_int{2}{n3} -  kc_vH(HCEF1,E1,Jxx1,Jzz1,T_RTM(n3),Hd_int{2}{n3}), 2));
         end 

        T1 = [4,6,12,20,30,50,70,100:50:250];
         Xa1 = 1./(1./Xab_inv_vT(HCEF1,E1,Jxx1,T1) + chi0);
         Xc1 = 1./(1./Xc_inv_vT(HCEF1,E1,Jzz1,T1)  + chi0);

        out = out + 10*(sum(power(Xa1 - XaDatInt,2)) + sum(power(Xc1 - XcDatInt,2))) ;
%

    end

    function out = CummErrSuscOnly(p0)
        
        
        HCEF1 = HCEFf(p0(1:6));
        Jxx1 = p0(7); Jzz1 = p0(8); chi0 = p0(9);

        [P1,D1] = eig(HCEF1 + Jz*1e-10);

        En = real([D1(1,1) D1(2,2) D1(3,3) D1(4,4) D1(5,5) D1(6,6) D1(7,7) D1(8,8)]);
        Ev1 = sort(En);
        E1 = Ev1(1);

        T1 = [4,6,12,20,30,50,70,100:50:250];
         Xa1 = 1./(1./Xab_inv_vT(HCEF1,E1,Jxx1,T1) + chi0);
         Xc1 = 1./(1./Xc_inv_vT(HCEF1,E1,Jzz1,T1)  + chi0);

        out = 100*(sum(power(Xa1 - XaDatInt,2)) + sum(power(Xc1 - XcDatInt,2))) ;

    end

    function out = CummErr_MvH(p0)
        
        HCEF1 = HCEFf(p0(1:6));
        Jxx1 = p0(7); Jzz1 = p0(8); chi0 = p0(9);

        [P1,D1] = eig(HCEF1 + Jz*1e-10);

        En = real([D1(1,1) D1(2,2) D1(3,3) D1(4,4) D1(5,5) D1(6,6) D1(7,7) D1(8,8)]);
        Ev1 = sort(En);
        E1 = Ev1(1);

        out = 0;
        
         for n3 = 1:7
             ind1 = (n3-1)*L1+1:n3*L1;
             out = out + sum(power(gJ*Mset0(ind1)' + H_single*chi0 - Mab_vH(HCEF1,E1,Jxx1,Jzz1,T_RTM(n3+1),H_single)/muB, 2));
             ind1 = (n3-1+7)*L1+1:(n3+7)*L1;
             out = out + sum(power(gJ*Mset0(ind1)' + H_single*chi0 - Mc_vH(HCEF1,E1,Jxx1,Jzz1,T_RTM(n3+1),H_single)/muB, 2));
         end 

    end

    
%% Analysis

switch option
    
    case 'DC Susc'
        
        Jxx = params(7);
        Jzz = params(8);
        
        HCEF = HCEFf(params(1:6));
        [P, D] = eig(HCEF);
        
        Ev = sort(real([D(1,1) D(2,2) D(3,3) D(4,4) D(5,5) D(6,6) D(7,7) D(8,8)]));
        E0 = Ev(1);
        
        figure; hold on;
        T = 1:300;
        %plot([4,6,12,20,30,50,70,100:50:250], XaDatInt, 'csq')
        %plot([4,6,12,20,30,50,70,100:50:250], XcDatInt, 'msq')
        plot(T, Xab_inv_vT(HCEF,E0,Jxx,T),'b')
        plot(T,  Xc_inv_vT(HCEF,E0,Jzz,T),'r')
        
    case 'DC Susc model diff'
        
        Jxx = params(7);
        Jzz = params(8);
        
        HCEF = HCEFf(params(1:6));
        [P, D] = eig(HCEF);
        
        Ev = sort(real([D(1,1) D(2,2) D(3,3) D(4,4) D(5,5) D(6,6) D(7,7) D(8,8)]));
        E0 = Ev(1);
        
        figure; hold on;
        T = [4,6,12,20,30,50,70,100:50:250];
        plot(T, XaDatInt - Xab_inv_vT(HCEF,E0,Jxx,T), 'bsq-')
        plot(T, XcDatInt -  Xc_inv_vT(HCEF,E0,Jzz,T), 'rsq-')
    
    case 'principle k'
        
        Jxx = params(7);
        Jzz = params(8);
        
        HCEF = HCEFf(params(1:6));
        [P, D] = eig(HCEF);
        Ev = sort(real([D(1,1) D(2,2) D(3,3) D(4,4) D(5,5) D(6,6) D(7,7) D(8,8)]));
        E0 = Ev(1);
        
        figure; hold on;
        
        H = 0:60;
        for t = [4,6,12,20,30,50,70]
            
            subplot(1,2,2); hold on;
            plot(H,kc_vH(HCEF,E0,Jxx,Jzz,t,H))
            
            subplot(1,2,1); hold on;
            plot(H,kab_vH(HCEF,E0,Jxx,Jzz,t,H))
            
        end
        
        subplot(1,2,1); grid on; box on; xlabel('H [T]'); ylabel('k_a [meV]'); xlim([0,60]); %ylim([0,2.5]);
        subplot(1,2,2); grid on; box on; xlabel('H [T]'); ylabel('k_c [meV]'); xlim([0,60]); %([0,2.5]);
        
    case 'H-Dept. k'
        
        % pi/2
        dat{1}{1} = load('p048_031020_TDO004.dat');
        dat{1}{2} = load('p043_031020_TDO003.dat');
        dat{1}{3} = load('p050_031020_TDO002.dat');
        dat{1}{4} = load('p052_031020_TDO001.dat');
        dat{1}{5} = load('p058_031020_TDO001.dat');
        dat{1}{6} = load('p067_031020_TDO001.dat');
        dat{1}{7} = load('p079_031020_TDO001.dat');
        dat{1}{8} = load('p092_031020_TDO001.dat');
        
        % 0
        dat{2}{1} = load('p059_031120_TDO001.dat');
        dat{2}{2} = load('p042_030920_TDO002.dat');
        dat{2}{3} = load('p057_031120_TDO001.dat');
        dat{2}{4} = load('p055_031120_TDO001.dat');
        dat{2}{5} = load('p050_031120_TDO001.dat');
        dat{2}{6} = load('p044_031120_TDO001.dat');
        dat{2}{7} = load('p037_031120_TDO001.dat');
        dat{2}{8} = load('p097_031020_TDO001.dat');
        
        startind = 251;
        endind = 2749-450;
        
        figure; hold on;
        
        T_RTM = [1.5,4,6,12,20, 30,50,70]; o = .75;
        
        colz = [0,0,1;...
            0.3,0.75,0.93;...
            0.47,0.67,0.19;...
            0,1,0;...
            0.87,0.49,0;...
            0.75,0,0.75;...
            1,0,0;...
            ];
        
        for n1 = [1,2]
            subplot(1,2,n1); hold on;
            for n2 = [2:8]
                x = dat{n1}{n2}(startind:end-endind,1);
                y = dat{n1}{n2}(startind:end-endind,2);
                
                ind = find(x == min(x));
                freq0{n1}(n2) = y(ind);
                
                plot( x, (y-freq0{n1}(n2))/(1e3*amp/w0(n2)) + o*(n2-1),'o','color',colz(n2-1,:),...
                    'displayname', [num2str(T_RTM(n2)) ' K Data'],...
                    'markerfacecolor',[1,1,1],'linewidth',1.5,...
                    'markersize',8);
            end
        end
        
        H = 1:60;
        
        Jxx = params(7);
        Jzz = params(8);
        
        HCEF = HCEFf(params(1:6));
        [P, D] = eig(HCEF);
        Ev = sort(real([D(1,1) D(2,2) D(3,3) D(4,4) D(5,5) D(6,6) D(7,7) D(8,8)]));
        E0 = Ev(1);
        
        n2 = 1;
        for t = [4,6,12,20, 30,50,70]
            subplot(1,2,2); hold on;
            Y = [0 kc_vH(HCEF,E0,Jxx,Jzz,t,H) + o*n2];
            plot([0 H],Y,'k--','displayname' ,[num2str(t) ' K model'],...
                'linewidth',2)
            subplot(1,2,1); hold on;
            Y = [0 kab_vH(HCEF,E0,Jxx,Jzz,t,H) + o*n2];
            plot([0 H],Y,'k--','displayname' ,[num2str(t) ' K  model'],...
                'linewidth',2)
            n2 = n2+1;
        end
        
        subplot(1,2,2); title('\theta = 0, H || c')
        %legend('location','southwest');
        grid on; box on; set(gca,'fontsize',16);
        xlabel('\mu_0H [T]')
        ylabel('\Deltaf [Hz]')
        
        subplot(1,2,1); title('\theta = \pi/2, H || ab')
        grid on; box on; set(gca,'fontsize',16);
        xlabel('\mu_0H [T]')
        ylabel('\Deltaf [Hz]')
        
    case 'H-Dept. k LF limit'
        
        % pi/2
        dat{1}{1} = load('p048_031020_TDO004.dat');
        dat{1}{2} = load('p043_031020_TDO003.dat');
        dat{1}{3} = load('p050_031020_TDO002.dat');
        dat{1}{4} = load('p052_031020_TDO001.dat');
        dat{1}{5} = load('p058_031020_TDO001.dat');
        dat{1}{6} = load('p067_031020_TDO001.dat');
        dat{1}{7} = load('p079_031020_TDO001.dat');
        dat{1}{8} = load('p092_031020_TDO001.dat');
        
        % 0
        dat{2}{1} = load('p059_031120_TDO001.dat');
        dat{2}{2} = load('p042_030920_TDO002.dat');
        dat{2}{3} = load('p057_031120_TDO001.dat');
        dat{2}{4} = load('p055_031120_TDO001.dat');
        dat{2}{5} = load('p050_031120_TDO001.dat');
        dat{2}{6} = load('p044_031120_TDO001.dat');
        dat{2}{7} = load('p037_031120_TDO001.dat');
        dat{2}{8} = load('p097_031020_TDO001.dat');
        
        startind = 251;
        endind = 2749-450;
        
        figure; hold on;
        
        T_RTM = [1.5,4,6,12,20, 30,50,70]; o = .75;
        
        colz = [0,0,1;...
            0.3,0.75,0.93;...
            0.47,0.67,0.19;...
            0,1,0;...
            0.87,0.49,0;...
            0.75,0,0.75;...
            1,0,0;...
            ];
        
        for n1 = [1,2]
            subplot(1,2,n1); hold on;
            for n2 = [2:8]
                x = dat{n1}{n2}(startind:end-endind,1);
                y = dat{n1}{n2}(startind:end-endind,2);
                
                ind = find(x == min(x));
                freq0{n1}(n2) = y(ind);
                
                plot( x, (y-freq0{n1}(n2))/(1e3*amp/w0(n2)) + o*(n2-1),'o','color',colz(n2-1,:),...
                    'displayname', [num2str(T_RTM(n2)) ' K Data'],...
                    'markerfacecolor',[1,1,1],'linewidth',1.5,...
                    'markersize',8);
            end
        end
        
        H = 1:30;
        
        Jxx = params(7);
        Jzz = params(8);
        
        HCEF = HCEFf(params(1:6));
        [P, D] = eig(HCEF);
        Ev = sort(real([D(1,1) D(2,2) D(3,3) D(4,4) D(5,5) D(6,6) D(7,7) D(8,8)]));
        E0 = Ev(1);
        
        n2 = 1;
        for t = [4,6,12,20,30,50,70]
            subplot(1,2,2); hold on;
            Xa = power(gJ*muB,2)*(C0/kB)./XaDatInt(n2);
            Xc = power(gJ*muB,2)*(C0/kB)./XcDatInt(n2);
            Y = [0 -power(H,2)*(Xa-Xc) ] + o*n2;
            plot([0 H],Y,'k--','displayname' ,[num2str(t) ' K model'],...
                'linewidth',2)
            subplot(1,2,1); hold on;
            Y = [0 power(H,2)*(Xa-Xc) ] + o*n2;
            plot([0 H],Y,'k--','displayname' ,[num2str(t) ' K  model'],...
                'linewidth',2)
            n2 = n2+1;
        end
        
        subplot(1,2,2); title('\theta = 0, H || c')
        %legend('location','southwest');
        grid on; box on; set(gca,'fontsize',16);
        xlabel('\mu_0H [T]')
        ylabel('\Deltaf [Hz]')
        
        subplot(1,2,1); title('\theta = \pi/2, H || ab')
        grid on; box on; set(gca,'fontsize',16);
        xlabel('\mu_0H [T]')
        ylabel('\Deltaf [Hz]')
    
    case 'MvH data'
        
        Jxx = params(7); Jzz = params(8);
        %Jxx = LTparams(1)/6; Jzz = LTparams(4)/6;
        
        p0 = params(1:6);
        pf = p0;
        
        HCEF = HCEFf(pf(1:6));
        [P, D] = eig(HCEF);
        Ev = sort(real([D(1,1) D(2,2) D(3,3) D(4,4) D(5,5) D(6,6) D(7,7) D(8,8)]));
        E0 = Ev(1);
        
        figure; hold on; colz = {[.7,.7,.7],'r','b','g','c','m','y','k'};
        
        subplot(1,2,1); hold on;
        for i = 1:7
            ind = (i-1)*L1+1:i*L1;
            Y = gJ*Mset0(ind);
            %plot(H_single,Y, '-', 'color', colz{i+1})
            plot( H_single(2:end-1), diff(diff(Y)) , '-', 'color', colz{i+1})
        end
        
         subplot(1,2,2); hold on;
        for i = 8:14
            ind = (i-1)*L1+1:i*L1;
            Y = gJ*Mset0(ind);
            %plot(H_single,Y, '-', 'color', colz{i+1-7})
            plot( H_single(2:end-1), diff(diff(Y)) , '-', 'color', colz{i-7+1})
        end  
        
        gxx0 = 3.5237; gzz0 = 1.3336;
        
        for t = [4,6,12,20,30,50,70]
            subplot(1,2,1);
            plot([1,1]*.66/(gxx0*muB/(2*kB*t)), [-3,1]*1e-5 , 'k--') 
            
            subplot(1,2,2);
            plot([1,1]*.66/(gzz0*muB/(2*kB*t)), [-3,1]*1e-5 , 'k--') 
        end
         
    case 'MvT model'
        
        Jxx = params(7); Jzz = params(8);
        %Jxx = LTparams(1)/6; Jzz = LTparams(4)/6;
        
        p0 = params(1:6);
        pf = p0;
        
        HCEF = HCEFf(pf(1:6));
        [P, D] = eig(HCEF);
        Ev = sort(real([D(1,1) D(2,2) D(3,3) D(4,4) D(5,5) D(6,6) D(7,7) D(8,8)]));
        E0 = Ev(1);
        
        figure; hold on; colz = {'r','b','g','c','m','y','k'};
%         
        
        T = [.1:.1:5 6:150]; i = 1;
        for H = [1,5]
            
            
            n1 = 1;
            for t = T
                Y1(n1) = Mab_vH(HCEF,E0,Jxx,Jzz,t,H)/muB;
                Y2(n1) =  Mc_vH(HCEF,E0,Jxx,Jzz,t,H)/muB;
                n1 = n1+1;
            end
            subplot(1,2,1); hold on;  set(gca,'fontsize',20);
            plot(T,1./(Y1/H),'color',colz{i},'displayname', ['H = ' num2str(H) ' T']);
            
            subplot(1,2,2); hold on;  set(gca,'fontsize',20);
            plot(T,1./(Y2/H),'color',colz{i},'displayname', ['H = ' num2str(H) ' T']);
            
            i = i + 1;
        end
        subplot(1,2,1); grid on; box on;
        ylabel('m_{ab} [\mu_B]'); xlabel('T [K]'); legend('location','northeast')
        
        subplot(1,2,2); grid on; box on;
        ylabel('m_c [\mu_B]'); xlabel('T [K]'); legend('location','northeast')    
           
    case 'nmag'
        
        Jxx = params(7); Jzz = params(8);
        %Jxx = LTparams(1)/6; Jzz = LTparams(4)/6;
        
        p0 = params(1:6);
        pf = p0;
        
        HCEF = HCEFf(pf(1:6));
        [P, D] = eig(HCEF);
        Ev = sort(real([D(1,1) D(2,2) D(3,3) D(4,4) D(5,5) D(6,6) D(7,7) D(8,8)]));
        E0 = Ev(1);
        
        figure; hold on; colz = {'r','b','g','c','m','y','k'};
         m0 = 4;
        
        T = 1:150; i = 1;
        for H = [1,3,6,9,12,15,18]
            
            
            n1 = 1;
            for t = T
                Y1(n1) = Mab_vH(HCEF,E0,Jxx,Jzz,t,H)/muB;
                Y2(n1) =  Mc_vH(HCEF,E0,Jxx,Jzz,t,H)/muB;
                n1 = n1+1;
            end
            subplot(1,2,1); hold on;  set(gca,'fontsize',20);
            plot(T,1-Y1/m0,'color',colz{i},'displayname', ['H = ' num2str(H) ' T']);
            
            subplot(1,2,2); hold on;  set(gca,'fontsize',20);
            plot(T,1-Y2/m0,'color',colz{i},'displayname', ['H = ' num2str(H) ' T']);
            
            i = i + 1;
        end
        subplot(1,2,1); grid on; box on;
        ylabel('n_{mag} = 1- m_{ab}/m_s'); xlabel('T [K]'); legend('location','northeast')
        
        subplot(1,2,2); grid on; box on;
        ylabel('n_{mag} = 1 - m_{c}/m_s'); xlabel('T [K]'); legend('location','northeast')    
        
    case 'nmag vs H'
        Jxx = params(7); Jzz = params(8);
        %Jxx = LTparams(1)/6; Jzz = LTparams(4)/6;
        
        p0 = params(1:6);
        pf = p0;
        
        HCEF = HCEFf(pf(1:6));
        [P, D] = eig(HCEF);
        Ev = sort(real([D(1,1) D(2,2) D(3,3) D(4,4) D(5,5) D(6,6) D(7,7) D(8,8)]));
        E0 = Ev(1);
        
        figure; hold on; colz = {[.7,.7,.7],'r','b','g','c','m','y','k'};
        
      ms = 4;
%         
        H = 0:.1:60; i = 1;
        for t = [1.5,4,6,12,20,30,50,70]
            
            subplot(1,2,2); hold on;
            
            Y = Mc_vH(HCEF,E0,Jxx,Jzz,t,H)/muB;
            plot(H,1-Y/ms,'-','linewidth',2, 'color', colz{i},'displayname', [num2str(t) ' K'])
            %d2Y = (C2(1)*Y(1:end-4) + C2(2)*Y(2:end-3) + C2(3)*Y(3:end-2) + C2(4)*Y(4:end-1) + C2(5)*Y(5:end));
            %plot(H(2:end-1),diff(diff(smooth(Y,3))),'-','linewidth',2, 'color', colz{i},'displayname', [num2str(t) ' K']);
            
            subplot(1,2,1); hold on;
            
            Y = Mab_vH(HCEF,E0,Jxx,Jzz,t,H)/muB;
            plot(H,1-Y/ms,'-','linewidth',2, 'color', colz{i},'displayname', [num2str(t) ' K'])
            %d2Y = (C2(1)*Y(1:end-4) + C2(2)*Y(2:end-3) + C2(3)*Y(3:end-2) + C2(4)*Y(4:end-1) + C2(5)*Y(5:end));
            %plot(H(2:end-1),diff(diff(smooth(Y,3))),'-','linewidth',2, 'color', colz{i},'displayname', [num2str(t) ' K']);
            
            i = i +1;
        end
        subplot(1,2,1); grid on; box on;
        ylabel('n_{mag} = 1-m_{ab}/m_s'); xlabel('\mu_0H [T]'); legend('location','southeast')
        
        subplot(1,2,2); grid on; box on;
        ylabel('n_{mag} = 1-m_{c}/m_s'); xlabel('\mu_0H [T]'); legend('location','southeast')  
               
    case 'MvH model alt'
        
        Jxx = params(7); Jzz = params(8);
        %Jxx = LTparams(1)/6; Jzz = LTparams(4)/6;
        
        p0 = params(1:6);
        pf = p0;
        
        HCEF = HCEFf(pf(1:6));
        [P, D] = eig(HCEF);
        Ev = sort(real([D(1,1) D(2,2) D(3,3) D(4,4) D(5,5) D(6,6) D(7,7) D(8,8)]));
        E0 = Ev(1);
        
        figure; hold on; colz = {[.7,.7,.7],'r','b','g','c','m','y','k'};
        
        dH = 0.01;
        H = 0:dH:60; i = 1;
        for t = [1.5,4,6,12,20,30,50,70]
            
            
            subplot(1,2,1); hold on;
            
            Y = Mab_vH_alt(HCEF,E0,Jxx,Jzz,t,H,dH)/muB;
            plot(H,Y,'-','linewidth',2, 'color', colz{i},'displayname', [num2str(t) ' K'])
            %d2Y = (C2(1)*Y(1:end-4) + C2(2)*Y(2:end-3) + C2(3)*Y(3:end-2) + C2(4)*Y(4:end-1) + C2(5)*Y(5:end));
            %plot(H(2:end-1),diff(diff(smooth(Y,3))),'-','linewidth',2, 'color', colz{i},'displayname', [num2str(t) ' K']);
            
            
            subplot(1,2,2); hold on;
            
            Y = Mc_vH_alt(HCEF,E0,Jxx,Jzz,t,H,dH)/muB;
            plot(H,Y,'-','linewidth',2, 'color', colz{i},'displayname', [num2str(t) ' K'])
            
            i = i +1;
        end
        subplot(1,2,1); grid on; box on;
        ylabel('m_{ab} [\mu_B]'); xlabel('\mu_0H [T]'); legend('location','southeast')
        
        subplot(1,2,2); grid on; box on;
        ylabel('m_c [\mu_B]'); xlabel('\mu_0H [T]'); legend('location','southeast')
        
    case 'MvH model'
        Jxx = params(7); Jzz = params(8);
        %Jxx = LTparams(1)/6; Jzz = LTparams(4)/6;
        
        p0 = params(1:6);
        pf = p0;
        
        HCEF = HCEFf(pf(1:6));
        [P, D] = eig(HCEF);
        Ev = sort(real([D(1,1) D(2,2) D(3,3) D(4,4) D(5,5) D(6,6) D(7,7) D(8,8)]));
        E0 = Ev(1);
        
        figure; hold on; colz = {[.7,.7,.7],'r','b','g','c','m','y','k'};
        
        subplot(1,2,1); hold on;
        for i = 1:7
            ind = (i-1)*L1+1:i*L1;
            plot( [0 H_single], gJ*[0; Mset0(ind)], 'sq', 'color', colz{i+1})
            
            %plot([0 H_single], power(gJ,2)*power(power(C0*(muB/kB),-1)*XaDatInt(i),-1)*[0 H_single],'k--','displayname',[num2str(T_RTM(i+1))] )
        end
        
         subplot(1,2,2); hold on;
        for i = 8:14
            ind = (i-1)*L1+1:i*L1;
            plot( [0 H_single], gJ*[0; Mset0(ind)], 'sq', 'color', colz{i-7+1})
            %plot([0 H_single], power(gJ,2)*power(power(C0*(muB/kB),-1)*XcDatInt(i-7),-1)*[0 H_single],'k--','displayname',[num2str(T_RTM(i+1-7))] )
        end
%         
        
        H = 0:.1:60; i = 1;
        for t = [2,4,6,8,14,20,30,50];
            
            subplot(1,2,2); hold on;
            
            Y = Mc_vH(HCEF,E0,Jxx,Jzz,t,H)/muB;
            
            output.Mc{i} = Y;
            plot(H,Y,'-','linewidth',2, 'color', colz{i},'displayname', [num2str(t) ' K'])
            %d2Y = (C2(1)*Y(1:end-4) + C2(2)*Y(2:end-3) + C2(3)*Y(3:end-2) + C2(4)*Y(4:end-1) + C2(5)*Y(5:end));
            %plot(H(2:end-1),diff(diff(smooth(Y,3))),'-','linewidth',2, 'color', colz{i},'displayname', [num2str(t) ' K']);
            
            subplot(1,2,1); hold on;
            
            Y = Mab_vH(HCEF,E0,Jxx,Jzz,t,H)/muB;
            
            output.Ma{i} = Y;
            plot(H,Y,'-','linewidth',2, 'color', colz{i},'displayname', [num2str(t) ' K'])
            %d2Y = (C2(1)*Y(1:end-4) + C2(2)*Y(2:end-3) + C2(3)*Y(3:end-2) + C2(4)*Y(4:end-1) + C2(5)*Y(5:end));
            %plot(H(2:end-1),diff(diff(smooth(Y,3))),'-','linewidth',2, 'color', colz{i},'displayname', [num2str(t) ' K']);
            
            i = i +1;
        end
        subplot(1,2,1); grid on; box on;
        ylabel('m_{ab} [\mu_B]'); xlabel('\mu_0H [T]'); legend('location','southeast')
        
        subplot(1,2,2); grid on; box on;
        ylabel('m_c [\mu_B]'); xlabel('\mu_0H [T]'); legend('location','southeast')
        
        
        
    case 'fmin MvH'
        
        Jxx = params(7); Jzz = params(8);
        %Jxx = LTparams(1)/6; Jzz = LTparams(4)/6;
        
        
        f = @(x) CummErr_MvH([x, Jxx, Jzz, 0]);
        
        p0 = params(1:6);
        pf = fminsearchbnd(f,[p0],[-1,-1,-1,-1,-1,-1], [1,1,1,1,1,1]);
        %pf = fminsearchbnd(f,[p0, Jxx, Jzz, 0],[-1,-1,-1,-1,-1,-1,Jxx,0,-10], [1,1,1,1,1,1,Jxx,10,10]);
        Jxx = pf(7);  Jzz = pf(8);
        display(pf); X0 = pf(9);
        
        HCEF = HCEFf(pf(1:6));
        [P, D] = eig(HCEF);
        Ev = sort(real([D(1,1) D(2,2) D(3,3) D(4,4) D(5,5) D(6,6) D(7,7) D(8,8)]));
        E0 = Ev(1);
        
        figure; hold on; colz = {[.7,.7,.7],'r','b','g','c','m','y','k'};
        
        subplot(1,2,1); hold on;
        for i = 1:7
            ind = (i-1)*L1+1:i*L1;
            plot( [0 H_single], gJ*[0; Mset0(ind)], 'sq-', 'color', colz{i+1})
        end
        
         subplot(1,2,2); hold on;
        for i = 8:14
            ind = (i-1)*L1+1:i*L1;
            plot( [0 H_single], gJ*[0; Mset0(ind)], 'sq-', 'color', colz{i-7+1})
        end
%         
        
        H = 0:.1:7; i = 1;
        for t = [1.5,4,6,12,20,30,50,70]
            
            subplot(1,2,2); hold on;
            
            Y = Mc_vH(HCEF,E0,Jxx,Jzz,t,H)/muB - X0*H;
            plot(H,Y,'-','linewidth',2, 'color', colz{i},'displayname', [num2str(t) ' K'])
            %d2Y = (C2(1)*Y(1:end-4) + C2(2)*Y(2:end-3) + C2(3)*Y(3:end-2) + C2(4)*Y(4:end-1) + C2(5)*Y(5:end));
            %plot(H(2:end-1),diff(diff(smooth(Y,3))),'-','linewidth',2, 'color', colz{i},'displayname', [num2str(t) ' K']);
            
            subplot(1,2,1); hold on;
            
            Y = Mab_vH(HCEF,E0,Jxx,Jzz,t,H)/muB - X0*H ;
            plot(H,Y,'-','linewidth',2, 'color', colz{i},'displayname', [num2str(t) ' K'])
            %d2Y = (C2(1)*Y(1:end-4) + C2(2)*Y(2:end-3) + C2(3)*Y(3:end-2) + C2(4)*Y(4:end-1) + C2(5)*Y(5:end));
            %plot(H(2:end-1),diff(diff(smooth(Y,3))),'-','linewidth',2, 'color', colz{i},'displayname', [num2str(t) ' K']);
            
            i = i +1;
        end
        subplot(1,2,1); grid on; box on;
        ylabel('m_{ab} [\mu_B]'); xlabel('\mu_0H [T]'); legend('location','southeast')
        
        subplot(1,2,2); grid on; box on;
        ylabel('m_c [\mu_B]'); xlabel('\mu_0H [T]'); legend('location','southeast')    
        
    case 'fmin MvH Diamag'
        
        %Jxx = params(7); Jzz = params(8);
        Jxx = LTparams(1)/6; Jzz = LTparams(4)/6;
        
        f = @(x) CummErr_MvH_wDiamag([x, 0]);
        
        p0 = params(1:6);
        %pf = fminsearchbnd(f,[p0],[-1,-1,-1,-1,-1,-1], [1,1,1,1,1,1]);
        pf = fminsearchbnd(f,[p0, Jxx, Jzz,0],[-1,-1,-1,-1,-1,-1,Jxx,0,-10], [1,1,1,1,1,1,Jxx,10,10]);
        Jxx = pf(7);  Jzz = pf(8);
        display(pf);
        
        HCEF = HCEFf(pf(1:6));
        [P, D] = eig(HCEF);
        Ev = sort(real([D(1,1) D(2,2) D(3,3) D(4,4) D(5,5) D(6,6) D(7,7) D(8,8)]));
        E0 = Ev(1);
        
        figure; hold on; colz = {[.7,.7,.7],'r','b','g','c','m','y','k'};
        
        subplot(1,2,1); hold on;
        for i = 1:7
            ind = (i-1)*L1+1:i*L1;
            plot( [0 H_single], gJ*[0; Mset0(ind)], 'sq-', 'color', colz{i+1})
        end
        
         subplot(1,2,2); hold on;
        for i = 8:14
            ind = (i-1)*L1+1:i*L1;
            plot( [0 H_single], gJ*[0; Mset0(ind)], 'sq-', 'color', colz{i-7+1})
        end
%         
        
        H = 0:.1:7; i = 1;
        for t = [1.5,4,6,12,20,30,50,70]
            
            subplot(1,2,2); hold on;
            
            Y = (Mc_vH(HCEF,E0,Jxx,Jzz,t,H) + X0*H)/muB;
            plot(H,Y,'-','linewidth',2, 'color', colz{i},'displayname', [num2str(t) ' K'])
            %d2Y = (C2(1)*Y(1:end-4) + C2(2)*Y(2:end-3) + C2(3)*Y(3:end-2) + C2(4)*Y(4:end-1) + C2(5)*Y(5:end));
            %plot(H(2:end-1),diff(diff(smooth(Y,3))),'-','linewidth',2, 'color', colz{i},'displayname', [num2str(t) ' K']);
            
            subplot(1,2,1); hold on;
            
            Y = (Mab_vH(HCEF,E0,Jxx,Jzz,t,H) + X0*H)/muB;
            plot(H,Y,'-','linewidth',2, 'color', colz{i},'displayname', [num2str(t) ' K'])
            %d2Y = (C2(1)*Y(1:end-4) + C2(2)*Y(2:end-3) + C2(3)*Y(3:end-2) + C2(4)*Y(4:end-1) + C2(5)*Y(5:end));
            %plot(H(2:end-1),diff(diff(smooth(Y,3))),'-','linewidth',2, 'color', colz{i},'displayname', [num2str(t) ' K']);
            
            i = i +1;
        end
        subplot(1,2,1); grid on; box on;
        ylabel('m_{ab} [\mu_B]'); xlabel('\mu_0H [T]'); legend('location','southeast')
        
        subplot(1,2,2); grid on; box on;
        ylabel('m_c [\mu_B]'); xlabel('\mu_0H [T]'); legend('location','southeast')        
        
    case 'MvH data tropic response'
        
        figure; hold on; colz = {'r','b','g','c','m','y','k'};
        
        %subplot(1,2,1); hold on;
        for i = 1:7
            ind = (i-1)*L1+1:i*L1;
            Ma{i} = muB*gJ*[Mset0(ind)];
            %plot( [H_single], gJ*[Mset0(ind)], '-', 'color', colz{i})
            
            Xa(i) = muB*gJ*f_mag(M_params(i,:),.1)./.1;
            %plot(H_single, H_single*Xa(i),'k--')
        end
        
        %subplot(1,2,2); hold on;
        for i = 8:14
            ind = (i-1)*L1+1:i*L1;
            Mc{i-7} = muB*gJ*[Mset0(ind)];
            %plot( [H_single], gJ*[Mset0(ind)], '-', 'color', colz{i-7})
            
            Xc(i-7) = muB*gJ*f_mag(M_params(i,:),.1)./.1;
            %plot(H_single, H_single*Xc(i-7),'k--')
        end
        
%         figure; hold on; C = 1.7909;
%         plot(T_RTM(2:8), XaDatInt(1:7),'-')
%         plot(T_RTM(2:8), XcDatInt(1:7),'-')
%         
%         plot(T_RTM(2:8), C./Xa,'bsq')
%         plot(T_RTM(2:8), C./Xc,'rsq')
%         
        
            X_1 = H_single;
            X_2 = power(X_1,2);

        for i = 1:7

            subplot(1,2,1); hold on;
            plot(X_1, X_2.*(Xa(i) - Xc(i)) ,'k--')
            plot(X_1, X_1.*Ma{i}' - X_2.*Xc(i) , colz{i} )
            
            subplot(1,2,2); hold on;
            plot(X_1, -X_2.*(Xa(i) - Xc(i)) ,'k--')
            plot(X_1, X_1.*Mc{i}' - X_2.*Xa(i) , colz{i})
            
        end

    case 'Trans Succ Analytical'
        
        Jxx = LTparams(1)/6; Jzz = LTparams(4)/6;
        
        p0 = params(1:6);
        
        Dat2 = CEF_Spectrum_D3d('Perturbative expression',p0);

        ax0 = Dat2.ax0{1}; ax1 = Dat2.ax0{2}; ax2 = Dat2.ax0{3}; ax3 = Dat2.ax0{4};
        az0 = Dat2.az0{1}; az1 = Dat2.az0{2}; az2 = Dat2.az0{3}; az3 = Dat2.az0{4};
        
        Bx0 = 2*kB*ax0/power(muB*gJ,2);
        Bx1 = 2*kB*ax1/power(muB*gJ,2);
        Bx2 = 2*kB*ax2/power(muB*gJ,2);
        Bx3 = 2*kB*ax3/power(muB*gJ,2);
        
        Bz0 = 2*kB*az0/power(muB*gJ,2);
        Bz1 = 2*kB*az1/power(muB*gJ,2);
        Bz2 = 2*kB*az2/power(muB*gJ,2);
        Bz3 = 2*kB*az3/power(muB*gJ,2);
        
        gxx0 = Dat2.gxx{1}; gxx1 = Dat2.gxx{2}; gxx2 = Dat2.gxx{3}; gxx3 = Dat2.gxx{4};
        gzz0 = Dat2.gzz{1}; gzz1 = Dat2.gzz{2}; gzz2 = Dat2.gzz{3}; gzz3 = Dat2.gzz{4};
        
        Ax0 = power( gxx0/(2*gJ),2);
        Ax1 = power( gxx1/(2*gJ),2);
        Ax2 = power( gxx2/(2*gJ),2);
        Ax3 = power( gxx3/(2*gJ),2);
        
        Az0 = power( gzz0/(2*gJ),2);
        Az1 = power( gzz1/(2*gJ),2);
        Az2 = power( gzz2/(2*gJ),2);
        Az3 = power( gzz3/(2*gJ),2);
        
        D10 = Dat2.D{1}; D20 = Dat2.D{2}; D30 = Dat2.D{3};
        
        T = 1:300;
        
        E10 = exp(-Dat2.D{1}./(kB.*T));
        E20 = exp(-Dat2.D{2}./(kB.*T));
        E30 = exp(-Dat2.D{3}./(kB.*T));
        
        Xinv1 = 2.0416*(T.*(1 + E10 + E20 + E30)./(...
            Ax0 - Bx0*T + (Ax1 - Bx1*T).*E10 + (Ax2 - Bx2*T).*E20 + (Ax3 - Bx3*T).*E30) + 6*Jxx);
        
        Xinv2 = 2.0416*(T.*(1 + E10 + E20 + E30)./(...
            Az0 - Bz0*T + (Az1 - Bz1*T).*E10 + (Az2 - Bz2*T).*E20 + (Az3 - Bz3*T).*E30) + 6*Jzz);
        
        %figure; hold on;
        %plot(T, Xinv1, 'r')
        %plot(T, Xinv2, 'b')
        
        
        D10 = Dat2.D{1}./kB;
        D20 = Dat2.D{2}./kB;
        D30 = Dat2.D{3}./kB;
        
        XT_inv1 = @(t,h) 2.0416*(t.*(  exp(-az0*power(h,2)./(kB*t)).*cosh(muB*h*gzz0./(2*kB*t))...
            + exp(-az1*power(h,2)./(kB*t)).*cosh(muB*h*gzz1./(2*kB*t)).*exp(-D10./t)...
            + exp(-az2*power(h,2)./(kB*t)).*cosh(muB*h*gzz2./(2*kB*t)).*exp(-D20./t)...
            + exp(-az3*power(h,2)./(kB*t)).*cosh(muB*h*gzz3./(2*kB*t)).*exp(-D30./t))./(...
            (Ax0.*sinh(muB*h*gzz0./(2*kB*t))./(muB*h*gzz0./(2*kB*t)) - Bx0*t.*cosh(muB*h*gzz0./(2*kB*t)) ).*exp(-az0*power(h,2)./(kB*t))...
            + (Ax1.*sinh(muB*h*gzz1./(2*kB*t))./(muB*h*gzz1./(2*kB*t)) - Bx1*t.*cosh(muB*h*gzz1./(2*kB*t)) ).*exp(-az1*power(h,2)./(kB*t)).*exp(-D10./t)...
            + (Ax2.*sinh(muB*h*gzz2./(2*kB*t))./(muB*h*gzz2./(2*kB*t)) - Bx2*t.*cosh(muB*h*gzz2./(2*kB*t)) ).*exp(-az2*power(h,2)./(kB*t)).*exp(-D20./t)...
            + (Ax3.*sinh(muB*h*gzz3./(2*kB*t))./(muB*h*gzz3./(2*kB*t)) - Bx3*t.*cosh(muB*h*gzz3./(2*kB*t)) ).*exp(-az3*power(h,2)./(kB*t)).*exp(-D30./t) ) + TCWx);
        
        XT_inv2 = @(t,h) 2.0416*(t.*(  exp(-ax0*power(h,2)./(kB*t)).*cosh(muB*h*gxx0./(2*kB*t))...
            + exp(-ax1*power(h,2)./(kB*t)).*cosh(muB*h*gxx1./(2*kB*t)).*exp(-D10./t)...
            + exp(-ax2*power(h,2)./(kB*t)).*cosh(muB*h*gxx2./(2*kB*t)).*exp(-D20./t)...
            + exp(-ax3*power(h,2)./(kB*t)).*cosh(muB*h*gxx3./(2*kB*t)).*exp(-D30./t))./(...
            (Az0.*sinh(muB*h*gxx0./(2*kB*t))./(muB*h*gxx0./(2*kB*t)) - Bz0*t.*cosh(muB*h*gxx0./(2*kB*t)) ).*exp(-ax0*power(h,2)./(kB*t))...
            + (Az1.*sinh(muB*h*gxx1./(2*kB*t))./(muB*h*gxx1./(2*kB*t)) - Bz1*t.*cosh(muB*h*gxx1./(2*kB*t)) ).*exp(-ax1*power(h,2)./(kB*t)).*exp(-D10./t)...
            + (Az2.*sinh(muB*h*gxx2./(2*kB*t))./(muB*h*gxx2./(2*kB*t)) - Bz2*t.*cosh(muB*h*gxx2./(2*kB*t)) ).*exp(-ax2*power(h,2)./(kB*t)).*exp(-D20./t)...
            + (Az3.*sinh(muB*h*gxx3./(2*kB*t))./(muB*h*gxx3./(2*kB*t)) - Bz3*t.*cosh(muB*h*gxx3./(2*kB*t)) ).*exp(-ax3*power(h,2)./(kB*t)).*exp(-D30./t) ) + TCWz);
        
        
        
        
        h = 0:.5:60;
        
        figure; hold on;
        for t = [4,6,12,20,30,50,70]
            subplot(1,2,1); hold on;
            plot(h, (.0104/.1002)./XT_inv1(t,h))
            
            subplot(1,2,2); hold on;
            plot(h, (.0104/.1002)./XT_inv2(t,h))
        end
        
    case 'Transervse Susc'
        
        Jxx = params(7); Jzz = params(8);
        %Jxx = LTparams(1)/6; Jzz = LTparams(4)/6;
        
        p0 = params(1:6);
        pf = p0;
        
        HCEF = HCEFf(pf(1:6));
        [P, D] = eig(HCEF);
        Ev = sort(real([D(1,1) D(2,2) D(3,3) D(4,4) D(5,5) D(6,6) D(7,7) D(8,8)]));
        E0 = Ev(1);
        
        
        for i = 1:7
            ind = (i-1)*L1+1:i*L1;
            %Ma{i} = muB*gJ*[Mset0(ind)];
            Xa(i) = muB*gJ*f_mag(M_params(i,:),.1)./.1;
        end
        
        %subplot(1,2,2); hold on;
        for i = 8:14
            ind = (i-1)*L1+1:i*L1;
            %Mc{i-7} = muB*gJ*[Mset0(ind)];
            Xc(i-7) = muB*gJ*f_mag(M_params(i,:),.1)./.1;

        end
        
        
        figure; hold on
        H = 0:60; i = 1;
        for t = [4]
            
            subplot(1,2,1); hold on;
            plot([H],[power(H,2).*XTc_vH(HCEF,E0,Jxx,Jzz,t,H)],'linewidth',2)
            
            subplot(1,2,2); hold on;
            plot([H],[power(H,2).*XTab_vH(HCEF,E0,Jxx,Jzz,t,H)],'linewidth',2)
            
            i = i+1;
        end
    
    case 'k Data'
        
        Hset1 = [2:3:29];
        Hset2 = [4:3:58];
        
        Hd_int{1} = {Hset2, Hset2, Hset1, Hset1, Hset1, Hset2, Hset1, Hset2};
        Hd_int{2} = {Hset2, Hset2, Hset1, Hset2, Hset1, Hset2, Hset1, Hset2};
        
        dat{1}{1} = load('p048_031020_TDO004.dat');
        dat{1}{2} = load('p043_031020_TDO003.dat');
        dat{1}{3} = load('p050_031020_TDO002.dat');
        dat{1}{4} = load('p052_031020_TDO001.dat');
        dat{1}{5} = load('p058_031020_TDO001.dat');
        dat{1}{6} = load('p067_031020_TDO001.dat');
        dat{1}{7} = load('p079_031020_TDO001.dat');
        dat{1}{8} = load('p092_031020_TDO001.dat');
        
        % 0
        dat{2}{1} = load('p059_031120_TDO001.dat');
        dat{2}{2} = load('p042_030920_TDO002.dat');
        dat{2}{3} = load('p057_031120_TDO001.dat');
        dat{2}{4} = load('p055_031120_TDO001.dat');
        dat{2}{5} = load('p050_031120_TDO001.dat');
        dat{2}{6} = load('p044_031120_TDO001.dat');
        dat{2}{7} = load('p037_031120_TDO001.dat');
        dat{2}{8} = load('p097_031020_TDO001.dat');
        
        for i = 1:7
            ind = (i-1)*L1+1:i*L1;
            Xa(i) = muB*gJ*f_mag(M_params(i,:),.1)./.1;
        end
        for i = 8:14
            ind = (i-1)*L1+1:i*L1;
            Xc(i-7) = muB*gJ*f_mag(M_params(i,:),.1)./.1;
        end
        X_1 = 3:10;
        X_2 = power(X_1,2);
        
        
        startind = 251;
        endind = 2749-450;
        
        figure; hold on;
        
        o = 0;
        
        colz = [0,0,1;...
            0.3,0.75,0.93;...
            0.47,0.67,0.19;...
            0,1,0;...
            0.87,0.49,0;...
            0.75,0,0.75;...
            1,0,0;...
            ];
        
        shape = {'^','sq'};
        for n1 = [1,2]
            hold on;
            for n2 = [2:8]
                x = dat{n1}{n2}(startind:end-endind,1);
                y = dat{n1}{n2}(startind:end-endind,2);
                
                %ind = find(x == min(x));
                %freq0{n1}(n2) = y(ind);
                
                ind = find(x <= 12);
                k_susc = power(-1,n1-1)*X_2.*(Xa(n2-1) - Xc(n2-1));    
                k_LFint = Interp1NonUnique(x(ind),y(ind),X_1); 
                f = @(x) sum(power( (k_LFint - x*1e3)./(1e3*amp/x) - k_susc,2));
                
                freq0{n1}(n2) = fminsearch(f,w0(n2));
                display([w0(n2), freq0{n1}(n2)])
                k_dat = (y-freq0{n1}(n2)*1e3)/(1e3*amp/freq0{n1}(n2));
                plot( x, k_dat.*power(-1,n1-1),shape{n1},'color',[.5,.5,.5],...
                    'displayname', [num2str(T_RTM(n2)) ' K Data'],...
                    'markerfacecolor',[1,1,1],'linewidth',1.5,...
                    'markersize',8);
                
                kd_int{n1}{n2} = Interp1NonUnique(x, k_dat, Hd_int{n1}{n2} );
                %display(sum(isnan(kd_int{n1}{n2} )))
                
                %plot(Hd_int{n1}{n2}, kd_int{n1}{n2}, 'sq', 'color', colz(n2-1,:))
            end
        end    
        
    case 'model'
        
        Hset1 = [2:3:29];
        Hset2 = [4:3:58];
        
        Hd_int{1} = {Hset2, Hset2, Hset1, Hset1, Hset1, Hset2, Hset1, Hset2};
        Hd_int{2} = {Hset2, Hset2, Hset1, Hset2, Hset1, Hset2, Hset1, Hset2};
        
        dat{1}{1} = load('p048_031020_TDO004.dat');
        dat{1}{2} = load('p043_031020_TDO003.dat');
        dat{1}{3} = load('p050_031020_TDO002.dat');
        dat{1}{4} = load('p052_031020_TDO001.dat');
        dat{1}{5} = load('p058_031020_TDO001.dat');
        dat{1}{6} = load('p067_031020_TDO001.dat');
        dat{1}{7} = load('p079_031020_TDO001.dat');
        dat{1}{8} = load('p092_031020_TDO001.dat');
        
        % 0
        dat{2}{1} = load('p059_031120_TDO001.dat');
        dat{2}{2} = load('p042_030920_TDO002.dat');
        dat{2}{3} = load('p057_031120_TDO001.dat');
        dat{2}{4} = load('p055_031120_TDO001.dat');
        dat{2}{5} = load('p050_031120_TDO001.dat');
        dat{2}{6} = load('p044_031120_TDO001.dat');
        dat{2}{7} = load('p037_031120_TDO001.dat');
        dat{2}{8} = load('p097_031020_TDO001.dat');
        
        startind = 251;
        endind = 2749-450;
        
        figure; hold on;
        
        o = 0;
        
        colz = [0,0,1;...
            0.3,0.75,0.93;...
            0.47,0.67,0.19;...
            0,1,0;...
            0.87,0.49,0;...
            0.75,0,0.75;...
            1,0,0;...
            ];
        
        
        for i = 1:7
            ind = (i-1)*L1+1:i*L1;
            Xa(i) = muB*gJ*f_mag(M_params(i,:),.1)./.1;
        end
        for i = 8:14
            ind = (i-1)*L1+1:i*L1;
            Xc(i-7) = muB*gJ*f_mag(M_params(i,:),.1)./.1;
        end
        X_1 = 3:10;
        X_2 = power(X_1,2);
        
        for n1 = [1,2]
            subplot(1,3,n1+1); hold on;
            for n2 = [2:8]
                x = dat{n1}{n2}(startind:end-endind,1);
                y = dat{n1}{n2}(startind:end-endind,2);
                
                %ind = find(x == min(x));
                %freq0{n1}(n2) = y(ind);
                
                ind = find(x <= 12);
                k_susc = power(-1,n1-1)*X_2.*(Xa(n2-1) - Xc(n2-1));    
                k_LFint = Interp1NonUnique(x(ind),y(ind),X_1); 
                g = @(x) sum(power( (k_LFint - x*1e3)./(1e3*amp/x) - k_susc,2));
                
                freq0{n1}(n2) = fminsearch(g,w0(n2));
                
                k_dat = (y-freq0{n1}(n2)*1e3)/(1e3*amp/w0(n2));
                plot( x, k_dat + o*(n2-1),'o','color',[.5,.5,.5],...
                    'displayname', [num2str(T_RTM(n2)) ' K Data'],...
                    'markerfacecolor',[1,1,1],'linewidth',1.5,...
                    'markersize',8);
                
                kd_int{n1}{n2} = Interp1NonUnique(x, k_dat, Hd_int{n1}{n2} );
                %display(sum(isnan(kd_int{n1}{n2} )))
                
                %plot(Hd_int{n1}{n2}, kd_int{n1}{n2}, 'sq', 'color', colz(n2-1,:))
            end
        end
        
        subplot(1,3,3); title('\theta = 0, H || c')
        %legend('location','southwest');
        grid on; box on; set(gca,'fontsize',16);
        xlabel('\mu_0H [T]')
        ylabel('\Deltaf [Hz]')
        
        subplot(1,3,2); title('\theta = \pi/2, H || ab')
        grid on; box on; set(gca,'fontsize',16);
        xlabel('\mu_0H [T]')
        ylabel('\Deltaf [Hz]')
        
        subplot(1,3,1); hold on
        grid on; box on; set(gca,'fontsize',16);
        xlabel('T [K]')
        ylabel('\chi^{-1}')
        plot([4,6,12,20,30,50,70,100:50:250], XaDatInt, 'bsq')
        plot([4,6,12,20,30,50,70,100:50:250], XcDatInt, 'rsq')
        
        
        Jxx = params(7); Jzz = params(8);
        %Jxx = LTparams(1)/6; Jzz = LTparams(4)/6;
        display([Jxx,Jzz])
        f = @(x) CummErr([x, Jxx, Jzz, 0]);
        
        p0 = params(1:6);
        pf = p0;%fminsearchbnd(f,p0,[-1,-1,-1,-1,-1,-1], [1,1,1,1,1,1]);
        
        %display(pf)
        
        HCEF = HCEFf(pf(1:6));
        [P, D] = eig(HCEF);
        Ev = sort(real([D(1,1) D(2,2) D(3,3) D(4,4) D(5,5) D(6,6) D(7,7) D(8,8)]));
        E0 = Ev(1);
        
        
        H = 0:60;
        for t = [4,6,12,20,30,50,70]
            
            subplot(1,3,3); hold on;
            plot(H,kc_vH(HCEF,E0,Jxx,Jzz,t,H),'linewidth',2)
            
            subplot(1,3,2); hold on;
            plot(H,kab_vH(HCEF,E0,Jxx,Jzz,t,H),'linewidth',2)
            
        end
        
        T = 1:300;
        subplot(1,3,1)
        plot(T, Xab_inv_vT(HCEF,E0,Jxx,T),'b')
        plot(T,  Xc_inv_vT(HCEF,E0,Jzz,T),'r')
    
    case 'fmin susc only'
        
        Hset1 = [2:3:29];
        Hset2 = [4:3:58];
        
        Hd_int{1} = {Hset2, Hset2, Hset1, Hset1, Hset1, Hset2, Hset1, Hset2};
        Hd_int{2} = {Hset2, Hset2, Hset1, Hset2, Hset1, Hset2, Hset1, Hset2};
        
        dat{1}{1} = load('p048_031020_TDO004.dat');
        dat{1}{2} = load('p043_031020_TDO003.dat');
        dat{1}{3} = load('p050_031020_TDO002.dat');
        dat{1}{4} = load('p052_031020_TDO001.dat');
        dat{1}{5} = load('p058_031020_TDO001.dat');
        dat{1}{6} = load('p067_031020_TDO001.dat');
        dat{1}{7} = load('p079_031020_TDO001.dat');
        dat{1}{8} = load('p092_031020_TDO001.dat');
        
        % 0
        dat{2}{1} = load('p059_031120_TDO001.dat');
        dat{2}{2} = load('p042_030920_TDO002.dat');
        dat{2}{3} = load('p057_031120_TDO001.dat');
        dat{2}{4} = load('p055_031120_TDO001.dat');
        dat{2}{5} = load('p050_031120_TDO001.dat');
        dat{2}{6} = load('p044_031120_TDO001.dat');
        dat{2}{7} = load('p037_031120_TDO001.dat');
        dat{2}{8} = load('p097_031020_TDO001.dat');
        
        startind = 251;
        endind = 2749-450;
        
        f= figure; hold on;
        f.Position = [100 100 2000 800];
        
        o = 0;
        
        colz = [0,0,1;...
            0.3,0.75,0.93;...
            0.47,0.67,0.19;...
            0,1,0;...
            0.87,0.49,0;...
            0.75,0,0.75;...
            1,0,0;...
            ];
        
        for n1 = [1,2]
            subplot(1,3,n1+1); hold on;
            for n2 = [2:8]
                x = dat{n1}{n2}(startind:end-endind,1);
                y = dat{n1}{n2}(startind:end-endind,2);
                
                ind = find(x == min(x));
                freq0{n1}(n2) = y(ind);
                
                k_dat = (y-freq0{n1}(n2))/(1e3*amp/w0(n2));
                plot( x, k_dat + o*(n2-1),'o','color',[.5,.5,.5],...
                    'displayname', [num2str(T_RTM(n2)) ' K Data'],...
                    'markerfacecolor',[1,1,1],'linewidth',1.5,...
                    'markersize',8);
                
                kd_int{n1}{n2} = Interp1NonUnique(x, k_dat, Hd_int{n1}{n2} );
                %display(sum(isnan(kd_int{n1}{n2} )))
                
                %plot(Hd_int{n1}{n2}, kd_int{n1}{n2}, 'sq', 'color', colz(n2-1,:))
            end
        end
        
        subplot(1,3,3); title('\theta = 0, H || c')
        %legend('location','southwest');
        grid on; box on; set(gca,'fontsize',16);
        xlabel('\mu_0H [T]')
        ylabel('\Deltaf [Hz]')
        
        subplot(1,3,2); title('\theta = \pi/2, H || ab')
        grid on; box on; set(gca,'fontsize',16);
        xlabel('\mu_0H [T]')
        ylabel('\Deltaf [Hz]')
        
        subplot(1,3,1); hold on
        grid on; box on; set(gca,'fontsize',16);
        xlabel('T [K]')
        ylabel('\chi^{-1}')
        plot([4,6,12,20,30,50,70,100:50:250], XaDatInt, 'bsq')
        plot([4,6,12,20,30,50,70,100:50:250], XcDatInt, 'rsq')
        
        
        
        Jxx = LTparams(1)/6; Jzz = LTparams(4)/6;
        
        f = @(x) CummErrSuscOnly([x, Jxx, Jzz, 0]);
        
        p0 = params(1:6);
        pf = fminsearchbnd(f,p0, 2*[-1,-1,-1,-1,-1,-1], 2*[1,1,1,1,1,1]);
        
        display(pf)
        
        HCEF = HCEFf(pf(1:6));
        [P, D] = eig(HCEF);
        Ev = sort(real([D(1,1) D(2,2) D(3,3) D(4,4) D(5,5) D(6,6) D(7,7) D(8,8)]));
        E0 = Ev(1);
        
        H = 0:60;
        for t = [4,6,12,20,30,50,70]
            
            subplot(1,3,3); hold on;
            plot(H,kc_vH(HCEF,E0,Jxx,Jzz,t,H),'linewidth',2)
            
            subplot(1,3,2); hold on;
            plot(H,kab_vH(HCEF,E0,Jxx,Jzz,t,H),'linewidth',2)
            
        end
        
        T = 1:300;
        subplot(1,3,1)
        plot(T, Xab_inv_vT(HCEF,E0,Jxx,T),'b')
        plot(T,  Xc_inv_vT(HCEF,E0,Jzz,T),'r')    
        
    case 'fmin k'
        
        Hset1 = [2:3:29];
        Hset2 = [4:3:58];
        
        Hd_int{1} = {Hset2, Hset2, Hset1, Hset1, Hset1, Hset2, Hset1, Hset2};
        Hd_int{2} = {Hset2, Hset2, Hset1, Hset2, Hset1, Hset2, Hset1, Hset2};
        
        dat{1}{1} = load('p048_031020_TDO004.dat');
        dat{1}{2} = load('p043_031020_TDO003.dat');
        dat{1}{3} = load('p050_031020_TDO002.dat');
        dat{1}{4} = load('p052_031020_TDO001.dat');
        dat{1}{5} = load('p058_031020_TDO001.dat');
        dat{1}{6} = load('p067_031020_TDO001.dat');
        dat{1}{7} = load('p079_031020_TDO001.dat');
        dat{1}{8} = load('p092_031020_TDO001.dat');
        
        % 0
        dat{2}{1} = load('p059_031120_TDO001.dat');
        dat{2}{2} = load('p042_030920_TDO002.dat');
        dat{2}{3} = load('p057_031120_TDO001.dat');
        dat{2}{4} = load('p055_031120_TDO001.dat');
        dat{2}{5} = load('p050_031120_TDO001.dat');
        dat{2}{6} = load('p044_031120_TDO001.dat');
        dat{2}{7} = load('p037_031120_TDO001.dat');
        dat{2}{8} = load('p097_031020_TDO001.dat');
        
        startind = 251;
        endind = 2749-450;
        
        figure; hold on;
        
        o = 0;
        
        colz = [0,0,1;...
            0.3,0.75,0.93;...
            0.47,0.67,0.19;...
            0,1,0;...
            0.87,0.49,0;...
            0.75,0,0.75;...
            1,0,0;...
            ];
        
        for n1 = [1,2]
            subplot(1,3,n1+1); hold on;
            for n2 = [2:8]
                x = dat{n1}{n2}(startind:end-endind,1);
                y = dat{n1}{n2}(startind:end-endind,2);
                
                ind = find(x == min(x));
                freq0{n1}(n2) = y(ind);
                
                k_dat = (y-freq0{n1}(n2))/(1e3*amp/w0(n2));
                plot( x, k_dat + o*(n2-1),'o','color',[.5,.5,.5],...
                    'displayname', [num2str(T_RTM(n2)) ' K Data'],...
                    'markerfacecolor',[1,1,1],'linewidth',1.5,...
                    'markersize',8);
                
                kd_int{n1}{n2} = Interp1NonUnique(x, k_dat, Hd_int{n1}{n2} );
                %display(sum(isnan(kd_int{n1}{n2} )))
                
                %plot(Hd_int{n1}{n2}, kd_int{n1}{n2}, 'sq', 'color', colz(n2-1,:))
            end
        end
        
        subplot(1,3,3); title('\theta = 0, H || c')
        %legend('location','southwest');
        grid on; box on; set(gca,'fontsize',16);
        xlabel('\mu_0H [T]')
        ylabel('\Deltaf [Hz]')
        
        subplot(1,3,2); title('\theta = \pi/2, H || ab')
        grid on; box on; set(gca,'fontsize',16);
        xlabel('\mu_0H [T]')
        ylabel('\Deltaf [Hz]')
        
        subplot(1,3,1); hold on
        grid on; box on; set(gca,'fontsize',16);
        xlabel('T [K]')
        ylabel('\chi^{-1}')
        plot([4,6,12,20,30,50,70,100:50:250], XaDatInt, 'bsq')
        plot([4,6,12,20,30,50,70,100:50:250], XcDatInt, 'rsq')
        
        
        
        Jxx = params(7);
        Jzz = params(8);
        %LTparams(1)/6; Jzz = LTparams(4)/6;
        
        f = @(x) CummErr([x, Jxx, Jzz, 0]);
        
        p0 = params(1:6);
        pf = fminsearchbnd(f,p0,[-1,-1,-1,-1,-1,-1], [1,1,1,1,1,1]);
        
        display(pf)
        
        HCEF = HCEFf(pf(1:6));
        [P, D] = eig(HCEF);
        Ev = sort(real([D(1,1) D(2,2) D(3,3) D(4,4) D(5,5) D(6,6) D(7,7) D(8,8)]));
        E0 = Ev(1);
        
        
        H = 0:60;
        for t = [4,6,12,20,30,50,70]
            
            subplot(1,3,3); hold on;
            plot(H,kc_vH(HCEF,E0,Jxx,Jzz,t,H),'linewidth',2)
            
            subplot(1,3,2); hold on;
            plot(H,kab_vH(HCEF,E0,Jxx,Jzz,t,H),'linewidth',2)
            
        end
        
        T = 1:300;
        subplot(1,3,1)
        plot(T, Xab_inv_vT(HCEF,E0,Jxx,T),'b')
        plot(T,  Xc_inv_vT(HCEF,E0,Jzz,T),'r')
        
    case 'fmin k no plot'
        
        Hset1 = [2:3:29];
        Hset2 = [4:3:58];
        
        Hd_int{1} = {Hset2, Hset2, Hset1, Hset1, Hset1, Hset2, Hset1, Hset2};
        Hd_int{2} = {Hset2, Hset2, Hset1, Hset2, Hset1, Hset2, Hset1, Hset2};
        
        dat{1}{1} = load('p048_031020_TDO004.dat');
        dat{1}{2} = load('p043_031020_TDO003.dat');
        dat{1}{3} = load('p050_031020_TDO002.dat');
        dat{1}{4} = load('p052_031020_TDO001.dat');
        dat{1}{5} = load('p058_031020_TDO001.dat');
        dat{1}{6} = load('p067_031020_TDO001.dat');
        dat{1}{7} = load('p079_031020_TDO001.dat');
        dat{1}{8} = load('p092_031020_TDO001.dat');
        
        % 0
        dat{2}{1} = load('p059_031120_TDO001.dat');
        dat{2}{2} = load('p042_030920_TDO002.dat');
        dat{2}{3} = load('p057_031120_TDO001.dat');
        dat{2}{4} = load('p055_031120_TDO001.dat');
        dat{2}{5} = load('p050_031120_TDO001.dat');
        dat{2}{6} = load('p044_031120_TDO001.dat');
        dat{2}{7} = load('p037_031120_TDO001.dat');
        dat{2}{8} = load('p097_031020_TDO001.dat');
        
        startind = 251;
        endind = 2749-450;
        
        for n1 = [1,2]
            for n2 = [2:8]
                x = dat{n1}{n2}(startind:end-endind,1);
                y = dat{n1}{n2}(startind:end-endind,2);
                
                ind = find(x == min(x));
                freq0{n1}(n2) = y(ind);
                
                k_dat = (y-freq0{n1}(n2))/(1e3*amp/w0(n2));
                kd_int{n1}{n2} = Interp1NonUnique(x, k_dat, Hd_int{n1}{n2} );

            end
        end

        Jxx = LTparams(1)/6; Jzz = LTparams(4)/6;
        
        f = @(x) CummErr([x, Jxx, Jzz, 0]);
        
        p0 = params(1:6);
        pf = fminsearchbnd(f,p0,2*[-1,-1,-1,-1,-1,-1], 2*[1,1,1,1,1,1]);
        
       % display(pf)
        
        output.p = pf;
        output.Res = f(pf);
        
    case 'get Res'
        
        Hset1 = [2:3:29];
        Hset2 = [4:3:58];
        
        Hd_int{1} = {Hset2, Hset2, Hset1, Hset1, Hset1, Hset2, Hset1, Hset2};
        Hd_int{2} = {Hset2, Hset2, Hset1, Hset2, Hset1, Hset2, Hset1, Hset2};
        
        dat{1}{1} = load('p048_031020_TDO004.dat');
        dat{1}{2} = load('p043_031020_TDO003.dat');
        dat{1}{3} = load('p050_031020_TDO002.dat');
        dat{1}{4} = load('p052_031020_TDO001.dat');
        dat{1}{5} = load('p058_031020_TDO001.dat');
        dat{1}{6} = load('p067_031020_TDO001.dat');
        dat{1}{7} = load('p079_031020_TDO001.dat');
        dat{1}{8} = load('p092_031020_TDO001.dat');
        
        % 0
        dat{2}{1} = load('p059_031120_TDO001.dat');
        dat{2}{2} = load('p042_030920_TDO002.dat');
        dat{2}{3} = load('p057_031120_TDO001.dat');
        dat{2}{4} = load('p055_031120_TDO001.dat');
        dat{2}{5} = load('p050_031120_TDO001.dat');
        dat{2}{6} = load('p044_031120_TDO001.dat');
        dat{2}{7} = load('p037_031120_TDO001.dat');
        dat{2}{8} = load('p097_031020_TDO001.dat');
        
        startind = 251;
        endind = 2749-450;
        
        for n1 = [1,2]
            for n2 = [2:8]
                x = dat{n1}{n2}(startind:end-endind,1);
                y = dat{n1}{n2}(startind:end-endind,2);
                
                ind = find(x == min(x));
                freq0{n1}(n2) = y(ind);
                
                k_dat = (y-freq0{n1}(n2))/(1e3*amp/w0(n2));
                kd_int{n1}{n2} = Interp1NonUnique(x, k_dat, Hd_int{n1}{n2} );
      
            end
        end
        
        f = @(x) CummErr([x 0]);
        
        p0 = params(1:8);
        output = f(p0);
          
    case 'fmin k float Jaa'
        
        Hset1 = [2:3:29];
        Hset2 = [4:3:58];
        
        Hd_int{1} = {Hset2, Hset2, Hset1, Hset1, Hset1, Hset2, Hset1, Hset2};
        Hd_int{2} = {Hset2, Hset2, Hset1, Hset2, Hset1, Hset2, Hset1, Hset2};
        
        dat{1}{1} = load('p048_031020_TDO004.dat');
        dat{1}{2} = load('p043_031020_TDO003.dat');
        dat{1}{3} = load('p050_031020_TDO002.dat');
        dat{1}{4} = load('p052_031020_TDO001.dat');
        dat{1}{5} = load('p058_031020_TDO001.dat');
        dat{1}{6} = load('p067_031020_TDO001.dat');
        dat{1}{7} = load('p079_031020_TDO001.dat');
        dat{1}{8} = load('p092_031020_TDO001.dat');
        
        % 0
        dat{2}{1} = load('p059_031120_TDO001.dat');
        dat{2}{2} = load('p042_030920_TDO002.dat');
        dat{2}{3} = load('p057_031120_TDO001.dat');
        dat{2}{4} = load('p055_031120_TDO001.dat');
        dat{2}{5} = load('p050_031120_TDO001.dat');
        dat{2}{6} = load('p044_031120_TDO001.dat');
        dat{2}{7} = load('p037_031120_TDO001.dat');
        dat{2}{8} = load('p097_031020_TDO001.dat');
        
        startind = 251;
        endind = 2749-450;
        
        figure; hold on;
        
        o = 0;
        
        colz = [0,0,1;...
            0.3,0.75,0.93;...
            0.47,0.67,0.19;...
            0,1,0;...
            0.87,0.49,0;...
            0.75,0,0.75;...
            1,0,0;...
            ];
        
        for n1 = [1,2]
            subplot(1,3,n1+1); hold on;
            for n2 = [2:8]
                x = dat{n1}{n2}(startind:end-endind,1);
                y = dat{n1}{n2}(startind:end-endind,2);
                
                ind = find(x == min(x));
                freq0{n1}(n2) = y(ind);
                
                k_dat = (y-freq0{n1}(n2))/(1e3*amp/w0(n2));
                plot( x, k_dat + o*(n2-1),'o','color',[.5,.5,.5],...
                    'displayname', [num2str(T_RTM(n2)) ' K Data'],...
                    'markerfacecolor',[1,1,1],'linewidth',1.5,...
                    'markersize',8);
                
                kd_int{n1}{n2} = Interp1NonUnique(x, k_dat, Hd_int{n1}{n2} );
                %display(sum(isnan(kd_int{n1}{n2} )))
                
                %plot(Hd_int{n1}{n2}, kd_int{n1}{n2}, 'sq', 'color', colz(n2-1,:))
            end
        end
        
        subplot(1,3,3); title('\theta = 0, H || c')
        %legend('location','southwest');
        grid on; box on; set(gca,'fontsize',16);
        xlabel('\mu_0H [T]')
        ylabel('\Deltaf [Hz]')
        
        subplot(1,3,2); title('\theta = \pi/2, H || ab')
        grid on; box on; set(gca,'fontsize',16);
        xlabel('\mu_0H [T]')
        ylabel('\Deltaf [Hz]')
        
        subplot(1,3,1); hold on
        grid on; box on; set(gca,'fontsize',16);
        xlabel('T [K]')
        ylabel('\chi^{-1}')
        plot([4,6,12,20,30,50,70,100:50:250], XaDatInt, 'bsq')
        plot([4,6,12,20,30,50,70,100:50:250], XcDatInt, 'rsq')
        
        
        f = @(x) CummErr([x 0]);
        
        p0 = params(1:8);
        pf = fminsearchbnd(f,p0,[-1,-1,-1,-1,-1,-1,0,0], [1,1,1,1,1,1,10,10]);
        
        display(pf)
        
        HCEF = HCEFf(pf(1:6));
        [P, D] = eig(HCEF);
        Ev = sort(real([D(1,1) D(2,2) D(3,3) D(4,4) D(5,5) D(6,6) D(7,7) D(8,8)]));
        E0 = Ev(1);
        
        Jxx = pf(7); Jzz = pf(8);
        
        H = 0:60;
        for t = [4,6,12,20,30,50,70]
            
            subplot(1,3,3); hold on;
            plot(H,kc_vH(HCEF,E0,Jxx,Jzz,t,H),'linewidth',2)
            
            subplot(1,3,2); hold on;
            plot(H,kab_vH(HCEF,E0,Jxx,Jzz,t,H),'linewidth',2)
            
        end
        
        T = 1:300;
        subplot(1,3,1)
        plot(T, Xab_inv_vT(HCEF,E0,Jxx,T),'b')
        plot(T,  Xc_inv_vT(HCEF,E0,Jzz,T),'r')
        
        output = pf;
        
    case 'fmin k float Jaa, include diamagnetic'
        
        Hset1 = [2:3:29];
        Hset2 = [4:3:58];
        
        Hd_int{1} = {Hset2, Hset2, Hset1, Hset1, Hset1, Hset2, Hset1, Hset2};
        Hd_int{2} = {Hset2, Hset2, Hset1, Hset2, Hset1, Hset2, Hset1, Hset2};
        
        dat{1}{1} = load('p048_031020_TDO004.dat');
        dat{1}{2} = load('p043_031020_TDO003.dat');
        dat{1}{3} = load('p050_031020_TDO002.dat');
        dat{1}{4} = load('p052_031020_TDO001.dat');
        dat{1}{5} = load('p058_031020_TDO001.dat');
        dat{1}{6} = load('p067_031020_TDO001.dat');
        dat{1}{7} = load('p079_031020_TDO001.dat');
        dat{1}{8} = load('p092_031020_TDO001.dat');
        
        % 0
        dat{2}{1} = load('p059_031120_TDO001.dat');
        dat{2}{2} = load('p042_030920_TDO002.dat');
        dat{2}{3} = load('p057_031120_TDO001.dat');
        dat{2}{4} = load('p055_031120_TDO001.dat');
        dat{2}{5} = load('p050_031120_TDO001.dat');
        dat{2}{6} = load('p044_031120_TDO001.dat');
        dat{2}{7} = load('p037_031120_TDO001.dat');
        dat{2}{8} = load('p097_031020_TDO001.dat');
        
        startind = 251;
        endind = 2749-450;
        
        figure; hold on;
        
        o = 0;
        
        colz = [0,0,1;...
            0.3,0.75,0.93;...
            0.47,0.67,0.19;...
            0,1,0;...
            0.87,0.49,0;...
            0.75,0,0.75;...
            1,0,0;...
            ];
        
        for n1 = [1,2]
            subplot(1,3,n1+1); hold on;
            for n2 = [2:8]
                x = dat{n1}{n2}(startind:end-endind,1);
                y = dat{n1}{n2}(startind:end-endind,2);
                
                ind = find(x == min(x));
                freq0{n1}(n2) = y(ind);
                
                k_dat = (y-freq0{n1}(n2))/(1e3*amp/w0(n2));
                plot( x, k_dat + o*(n2-1),'o','color',[.5,.5,.5],...
                    'displayname', [num2str(T_RTM(n2)) ' K Data'],...
                    'markerfacecolor',[1,1,1],'linewidth',1.5,...
                    'markersize',8);
                
                kd_int{n1}{n2} = Interp1NonUnique(x, k_dat, Hd_int{n1}{n2} );
                %display(sum(isnan(kd_int{n1}{n2} )))
                
                %plot(Hd_int{n1}{n2}, kd_int{n1}{n2}, 'sq', 'color', colz(n2-1,:))
            end
        end
        
        subplot(1,3,3); title('\theta = 0, H || c')
        %legend('location','southwest');
        grid on; box on; set(gca,'fontsize',16);
        xlabel('\mu_0H [T]')
        ylabel('\Deltaf [Hz]')
        
        subplot(1,3,2); title('\theta = \pi/2, H || ab')
        grid on; box on; set(gca,'fontsize',16);
        xlabel('\mu_0H [T]')
        ylabel('\Deltaf [Hz]')
        
        subplot(1,3,1); hold on
        grid on; box on; set(gca,'fontsize',16);
        xlabel('T [K]')
        ylabel('\chi^{-1}')
        plot([4,6,12,20,30,50,70,100:50:250], XaDatInt, 'bsq')
        plot([4,6,12,20,30,50,70,100:50:250], XcDatInt, 'rsq')
        
        
        f = @(x) CummErr(x);
        
        p0 = params(1:9);
        pf = fminsearchbnd(f,p0,[-1,-1,-1,-1,-1,-1,0,0,10], [1,1,1,1,1,1,10,10,10]);
        
        display(pf)
        
        HCEF = HCEFf(pf(1:6));
        [P, D] = eig(HCEF);
        Ev = sort(real([D(1,1) D(2,2) D(3,3) D(4,4) D(5,5) D(6,6) D(7,7) D(8,8)]));
        E0 = Ev(1);
        
        Jxx = pf(7); Jzz = pf(8); X0 = pf(9);
        
        H = 0:60;
        for t = [4,6,12,20,30,50,70]
            
            subplot(1,3,3); hold on;
            plot(H,kc_vH(HCEF,E0,Jxx,Jzz,t,H),'linewidth',2)
            
            subplot(1,3,2); hold on;
            plot(H,kab_vH(HCEF,E0,Jxx,Jzz,t,H),'linewidth',2)
            
        end
        
        T = 1:300;
        subplot(1,3,1)
        plot(T, 1./(1./Xab_inv_vT(HCEF,E0,Jxx,T) + X0),'b')
        plot(T, 1./(1./Xc_inv_vT(HCEF,E0,Jzz,T)  + X0),'r')
          
    case 'fmin k fix Jaa, include diamagnetic'
        
        Hset1 = [2:3:29];
        Hset2 = [4:3:58];
        
        Hd_int{1} = {Hset2, Hset2, Hset1, Hset1, Hset1, Hset2, Hset1, Hset2};
        Hd_int{2} = {Hset2, Hset2, Hset1, Hset2, Hset1, Hset2, Hset1, Hset2};
        
        dat{1}{1} = load('p048_031020_TDO004.dat');
        dat{1}{2} = load('p043_031020_TDO003.dat');
        dat{1}{3} = load('p050_031020_TDO002.dat');
        dat{1}{4} = load('p052_031020_TDO001.dat');
        dat{1}{5} = load('p058_031020_TDO001.dat');
        dat{1}{6} = load('p067_031020_TDO001.dat');
        dat{1}{7} = load('p079_031020_TDO001.dat');
        dat{1}{8} = load('p092_031020_TDO001.dat');
        
        % 0
        dat{2}{1} = load('p059_031120_TDO001.dat');
        dat{2}{2} = load('p042_030920_TDO002.dat');
        dat{2}{3} = load('p057_031120_TDO001.dat');
        dat{2}{4} = load('p055_031120_TDO001.dat');
        dat{2}{5} = load('p050_031120_TDO001.dat');
        dat{2}{6} = load('p044_031120_TDO001.dat');
        dat{2}{7} = load('p037_031120_TDO001.dat');
        dat{2}{8} = load('p097_031020_TDO001.dat');
        
        startind = 251;
        endind = 2749-450;
        
        figure; hold on;
        
        o = 0;
        
        colz = [0,0,1;...
            0.3,0.75,0.93;...
            0.47,0.67,0.19;...
            0,1,0;...
            0.87,0.49,0;...
            0.75,0,0.75;...
            1,0,0;...
            ];
        
        for n1 = [1,2]
            subplot(1,3,n1+1); hold on;
            for n2 = [2:8]
                x = dat{n1}{n2}(startind:end-endind,1);
                y = dat{n1}{n2}(startind:end-endind,2);
                
                ind = find(x == min(x));
                freq0{n1}(n2) = y(ind);
                
                k_dat = (y-freq0{n1}(n2))/(1e3*amp/w0(n2));
                plot( x, k_dat + o*(n2-1),'o','color',[.5,.5,.5],...
                    'displayname', [num2str(T_RTM(n2)) ' K Data'],...
                    'markerfacecolor',[1,1,1],'linewidth',1.5,...
                    'markersize',8);
                
                kd_int{n1}{n2} = Interp1NonUnique(x, k_dat, Hd_int{n1}{n2} );
                %display(sum(isnan(kd_int{n1}{n2} )))
                
                %plot(Hd_int{n1}{n2}, kd_int{n1}{n2}, 'sq', 'color', colz(n2-1,:))
            end
        end
        
        subplot(1,3,3); title('\theta = 0, H || c')
        %legend('location','southwest');
        grid on; box on; set(gca,'fontsize',16);
        xlabel('\mu_0H [T]')
        ylabel('\Deltaf [Hz]')
        
        subplot(1,3,2); title('\theta = \pi/2, H || ab')
        grid on; box on; set(gca,'fontsize',16);
        xlabel('\mu_0H [T]')
        ylabel('\Deltaf [Hz]')
        
        subplot(1,3,1); hold on
        grid on; box on; set(gca,'fontsize',16);
        xlabel('T [K]')
        ylabel('\chi^{-1}')
        plot([4,6,12,20,30,50,70,100:50:250], XaDatInt, 'bsq')
        plot([4,6,12,20,30,50,70,100:50:250], XcDatInt, 'rsq')
        
        Jxx = LTparams(1)/6; Jzz = LTparams(4)/6;
        
        f = @(x) CummErr([x(1:6), Jxx, Jzz, x(7)]);
        
        p0 = [params(1:6) params(9)];
        pf = fminsearchbnd(f,p0,[-1,-1,-1,-1,-1,-1,-100], [1,1,1,1,1,1,0]);
        
        display(pf)
        
        HCEF = HCEFf(pf(1:6));
        [P, D] = eig(HCEF);
        Ev = sort(real([D(1,1) D(2,2) D(3,3) D(4,4) D(5,5) D(6,6) D(7,7) D(8,8)]));
        E0 = Ev(1);
        
        %Jxx = pf(7); Jzz = pf(8); 
        X0 = pf(7);
        
        H = 0:60;
        for t = [4,6,12,20,30,50,70]
            
            subplot(1,3,3); hold on;
            plot(H,kc_vH(HCEF,E0,Jxx,Jzz,t,H),'linewidth',2)
            
            subplot(1,3,2); hold on;
            plot(H,kab_vH(HCEF,E0,Jxx,Jzz,t,H),'linewidth',2)
            
        end
        
        T = 1:300;
        subplot(1,3,1)
        plot(T, 1./(1./Xab_inv_vT(HCEF,E0,Jxx,T) + X0),'b')
        plot(T, 1./(1./Xc_inv_vT(HCEF,E0,Jzz,T)  + X0),'r')
           
end


end