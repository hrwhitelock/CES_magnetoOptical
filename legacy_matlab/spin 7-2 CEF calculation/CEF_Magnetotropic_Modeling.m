function output = CEF_Magnetotropic_Modeling(option,params,opt2)

amp = params(1);
paramsB = params(2:7);
paramsJJ = params(8:9);

%% Fixed angle/Temp ZF resonance frequencies
w0 = [49.5241, 49.5276, 49.5238, 49.5215, 49.5149, 49.4941, 49.4440, 49.3838];

%% Spin 7/2 matrices, {|7/2,7/2>.,...,|7/2,-7/2>} basis

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

% B20 = params(1);
% B40 = params(2);
% B43 = params(3);
% B60 = params(4);
% B63 = params(5);
% B66 = params(6);

HCEFf = @(Bnm) Bnm(1)*O20 + Bnm(2)*O40 + Bnm(3)*O43 + Bnm(4)*O60 + Bnm(5)*O63 + Bnm(6)*O66;

%% Numerical Physical Consants
muB = 5.78828e-2; % [meV/T];
kB  = 8.617e-2  ; % [meV/K];
gJ = 8/7; % L=3, S=1/2, J=7/2 g-lande factor;
z = 6; % #NN
A1 = 6*kB/(muB*gJ);

%% 2nd derivative finite difference coeficients
C1 = [1/12, 	-2/3,      0,   2/3, 	-1/12];
C2 = [-1/12, 	 4/3, 	-5/2, 	4/3,	-1/12];
C3 = [-1/2, 	   1, 	   0, 	 -1, 	  1/2];

%% Define the J-J coupling matrix
% JJxx = paramsJJ; % [K]
% JJzz = paramsJJ; % [K]
% JJ = [JJxx,0,0; 0,JJxx,0; 0,0,JJzz]*kB;
JJf = @(JJab) [JJab(1),0,0; 0,JJab(1),0; 0,0,JJab(2)]*kB;

%% Define the MF hamiltonian for a single particle,
% given a 8x8 HCEF matrix (HCEFm)
% given a 3x3 J-J coupling constant matrix (JJab)
% at finite field hv = (hx,hy,hz)
% in terms of MF magnetizations mv = (mx,my,mz)

HMF = @(HCEFm, JJab, hv, mv) HCEFm - (z/2)*(mv*JJab*mv')*ID...
    - (muB*gJ*hv(1) - (z/2)*(JJab(:,1)'+JJab(1,:))*mv')*Jx...
    - (muB*gJ*hv(2) - (z/2)*(JJab(:,2)'+JJab(2,:))*mv')*Jy...
    - (muB*gJ*hv(3) - (z/2)*(JJab(:,3)'+JJab(3,:))*mv')*Jz;

% Edipdip = -4.0707e-4; %[meV]
%
% HMF = @(HCEFm, JJab, hv, mv) HCEFm - (z/2)*(mv*JJab*mv')*ID...
%     - (muB*gJ*hv(1) - (z/2)*(JJab(:,1)'+JJab(1,:))*mv' + 3*Edipdip*mv(1))*Jx...  - (muB*gJ*hv(2) - (z/2)*(JJab(:,2)'+JJab(2,:))*mv' + 3*Edipdip*mv(2))*Jy...
%     - (muB*gJ*hv(3) - (z/2)*(JJab(:,3)'+JJab(3,:))*mv' - 6*Edipdip*mv(3))*Jz...
%     + (Edipdip*mv*diag([3,3,-6])*mv'/2)*ID;

% The only diagonal version
% HMF = @(HCEFm, JJab, hv, mv) HCEFm - (z/2)*(JJxx*mv(1)*mv(1) + JJxx*mv(2)*mv(2) + JJzz*mv(3)*mv(3))*ID...
%                               - (muB*gJ*hv(1) - (z)*(JJab(1,1))*mv(1))*Jx...
%                               - (muB*gJ*hv(2) - (z)*(JJab(2,2))*mv(2))*Jy...
%                               - (muB*gJ*hv(3) - (z)*(JJab(3,3))*mv(3))*Jz;

%% function for determining a self-consitent magnetization
% given a 8x8 HCEF matrix (HCEFm)
% given a 3x3 J-J coupling constant matrix (JJab)
% at finite field hv = (hx,hy,hz)
% at a temperature t
    function out = M(HCEFm, JJab, hv, t)
        
        B =  1./(kB*t);
        
        [P,D] = eig(HCEFm);
        Ev = sort(real([D(1,1) D(2,2) D(3,3) D(4,4) D(5,5) D(6,6) D(7,7) D(8,8)]));
        E0 = Ev(1);
        %         if E0 == 0
        %             E0 = min(D(D>0));
        %         end
        
        f1 = @(mx,my,mz) trace( Jx*expm(-B*(HMF(HCEFm, JJab, hv, [mx,my,mz]) - E0*ID)) )/trace( expm(-B*(HMF(HCEFm, JJab, hv, [mx,my,mz]) - E0*ID)) ) - mx;
        %f2 = @(mx,my,mz) trace( Jy*expm(-B*(HMF(HCEFm, JJab, hv, [mx,my,mz]) - E0*ID)) )/trace( expm(-B*(HMF(HCEFm, JJab, hv, [mx,my,mz]) - E0*ID)) ) - my;
        f3 = @(mx,my,mz) trace( Jz*expm(-B*(HMF(HCEFm, JJab, hv, [mx,my,mz]) - E0*ID)) )/trace( expm(-B*(HMF(HCEFm, JJab, hv, [mx,my,mz]) - E0*ID)) ) - mz;
        
        mx_g = 0; my_g = 0; mz_g = 0;
        for m = 1:7
            f4 = @(x) real(f1(  x ,my_g,mz_g)); mx_g = fzero(f4,mx_g);
            %f5 = @(x) real(f2(mx_g,  x ,mz_g)); my_g = fzero(f5,my_g);
            f6 = @(x) real(f3(mx_g,my_g,  x )); mz_g = fzero(f6,mz_g);
        end
        out = [mx_g, my_g, mz_g];
    end

%     function out = M_aph(HCEFm, JJab, hv, t)
%
%         B =  1./(kB*t);
%
%         [P,D] = eig(HCEFm);
%         Ev = sort(real([D(1,1) D(2,2) D(3,3) D(4,4) D(5,5) D(6,6) D(7,7) D(8,8)]));
%         E0 = Ev(1);
%         %         if E0 == 0
%         %             E0 = min(D(D>0));
%         %         end
%
%         f1 = @(mx,my,mz) trace( Jx*expm(-B*(HMF(HCEFm, JJab, hv, [mx,my,mz]) - E0*ID)) )/trace( expm(-B*(HMF(HCEFm, JJab, hv, [mx,my,mz]) - E0*ID)) ) - mx;
%         f2 = @(mx,my,mz) trace( Jy*expm(-B*(HMF(HCEFm, JJab, hv, [mx,my,mz]) - E0*ID)) )/trace( expm(-B*(HMF(HCEFm, JJab, hv, [mx,my,mz]) - E0*ID)) ) - my;
%         f3 = @(mx,my,mz) trace( Jz*expm(-B*(HMF(HCEFm, JJab, hv, [mx,my,mz]) - E0*ID)) )/trace( expm(-B*(HMF(HCEFm, JJab, hv, [mx,my,mz]) - E0*ID)) ) - mz;
%
%         mx_g = 1e-8; my_g = 0; mz_g = 1e-8;
%         for m = 1:3
%             f4 = @(x) real(f1(  x ,my_g,mz_g)); mx_g = fzero(f4,mx_g);
%             f5 = @(x) real(f2(mx_g,  x ,mz_g)); my_g = fzero(f5,my_g);
%             f6 = @(x) real(f3(mx_g,my_g,  x )); mz_g = fzero(f6,mz_g);
%         end
%         out = [mx_g, my_g, mz_g];
%     end

%% function for k(Th) given equispaced Th
    function out = k_vTh_equispc(HCEF,JJ,t,h,Th,dth)
        
        m1 = 1; b1 = 1./(kB*t);
        for th1 = Th
            
            hx1 = h*sin(th1);
            hz1 = h*cos(th1);
            hv = [hx1,0,hz1];
            
            mv = M(HCEF, JJ, hv, t);
            HMF1 = HMF(HCEF, JJ, hv, mv);
            [P,D] = eig(HMF1);
            Ev = sort(real([D(1,1) D(2,2) D(3,3) D(4,4) D(5,5) D(6,6) D(7,7) D(8,8)]));
            E0(m1) = Ev(1);
            
            ZZ(m1) = trace(expm(-b1*(HMF1 - E0(m1)*ID)));
            
            m1=m1+1;
        end
        FF = -kB*t*log(ZZ) + E0 ;
        
        out = (C2(1)*FF(1:end-4) + C2(2)*FF(2:end-3) + C2(3)*FF(3:end-2)...
            + C2(4)*FF(4:end-1) + C2(5)*FF(5:end))/dth/dth;
        
    end

    function out = F_vTh_equispc(HCEF,JJ,t,h,Th)
        
        m1 = 1; b1 = 1./(kB*t);
        for th1 = Th
            
            hx1 = h*sin(th1);
            hz1 = h*cos(th1);
            hv = [hx1,0,hz1];
            
            mv = M(HCEF, JJ, hv, t);
            HMF1 = HMF(HCEF, JJ, hv, mv);
            [P,D] = eig(HMF1);
            Ev = sort(real([D(1,1) D(2,2) D(3,3) D(4,4) D(5,5) D(6,6) D(7,7) D(8,8)]));
            E0(m1) = Ev(1);
            
            ZZ(m1) = trace(expm(-b1*(HMF1 - E0(m1)*ID)));
            
            m1=m1+1;
        end
        FF = -kB*t*log(ZZ) + E0 ;
        out = FF;
        
    end

    function out = FiniteDiff5Point(C,Y)
        
        out = (C(1)*Y(1:end-4) + C(2)*Y(2:end-3) + C(3)*Y(3:end-2)...
            + C(4)*Y(4:end-1) + C(5)*Y(5:end));
    end

%% function for k(Th) at a single point (h,Th)
    function out = k_1pt(HCEF,JJ,t,h,th,dth)
        
        m1 = 1; b1 = 1./(kB*t);
        
        for th1 = (th-2*dth):dth:(th+2*dth)
            hx1 = h*sin(th1);
            hz1 = h*cos(th1);
            hv = [hx1,0,hz1];
            
            mv = M(HCEF, JJ, hv, t);
            HMF1 = HMF(HCEF, JJ, hv, mv);
            [P,D] = eig(HMF1);
            Ev = sort(real([D(1,1) D(2,2) D(3,3) D(4,4) D(5,5) D(6,6) D(7,7) D(8,8)]));
            E0(m1) = Ev(1);
            
            ZZ(m1) = trace(expm(-b1*(HMF1 - E0(m1)*ID)));
            
            m1=m1+1;
        end
        FF = -kB*t*log(ZZ) + E0;
        
        out = (C2(1)*FF(1) + C2(2)*FF(2) + C2(3)*FF(3)...
            + C2(4)*FF(4) + C2(5)*FF(5))/dth/dth;
    end

%% function for single point at arb field (not defined by angle)
    function out = F_1pt_hv(HCEF,JJ,t,hv)
        
        m1 = 1; b1 = 1./(kB*t);
        
        mv = M(HCEF, JJ, hv, t);
        HMF1 = HMF(HCEF, JJ, hv, mv);
        [P,D] = eig(HMF1);
        Ev = sort(real([D(1,1) D(2,2) D(3,3) D(4,4) D(5,5) D(6,6) D(7,7) D(8,8)]));
        E0 = Ev(1);
        
        ZZ = trace(expm(-b1*(HMF1 - E0*ID)));
        FF = -kB*t*log(ZZ) + E0;
        
        out = FF;
    end


    function out = k_1pt_aph(HCEF,JJ,t,h,th,ph,dth)
        
        m1 = 1; b1 = 1./(kB*t);
        
        for th1 = (th-2*dth):dth:(th+2*dth)
            hx1 = h*sin(th1)*cos(ph);
            hy1 = h*sin(th1)*sin(ph);
            hz1 = h*cos(th1);
            hv = [hx1,hy1,hz1];
            
            mv = M(HCEF, JJ, hv, t);
            HMF1 = HMF(HCEF, JJ, hv, mv);
            [P,D] = eig(HMF1);
            Ev = sort(real([D(1,1) D(2,2) D(3,3) D(4,4) D(5,5) D(6,6) D(7,7) D(8,8)]));
            E0(m1) = Ev(1);
            
            ZZ(m1) = trace(expm(-b1*(HMF1 - E0(m1)*ID)));
            
            m1=m1+1;
        end
        FF = -kB*t*log(ZZ) + E0;
        
        out = (C2(1)*FF(1) + C2(2)*FF(2) + C2(3)*FF(3)...
            + C2(4)*FF(4) + C2(5)*FF(5))/dth/dth;
    end

%% function for k(Th) at an aribrary set of Th
    function out = k_vTh(HCEF,JJ,t,h,thset,dth)
        ni = 1;
        for th1 = thset
            out(ni) = k_1pt(HCEF,JJ,t,h,th1,dth);
            ni = ni + 1;
        end
    end

%% function for k(H) at an aribrary set of H
    function out = k_vH(HCEF,JJ,t,hset,th,dth)
        ni = 1;
        for h1 = hset
            out(ni) = k_1pt(HCEF,JJ,t,h1,th,dth);
            ni = ni + 1;
        end
    end

%% function for k(H) at an aribrary set of H, arb phi;
    function out = k_vH_aph(HCEF,JJ,t,hset,th,ph,dth)
        ni = 1;
        for h1 = hset
            out(ni) = k_1pt_aph(HCEF,JJ,t,h1,th,ph,dth);
            ni = ni + 1;
        end
    end

%% function for Th-fminsearch at 4 K
    function out = CummErr_ThDept0(n0,HCEF,JJ)
        
        Err = zeros([1,length(Hset)]);
        for ni = 1:length(Hset)
            h1 = Hset(ni);
            nj = 1;
            K = n0*k_vTh(HCEF,JJ,4,h1,Thset,.001);
            %             for th1 = Thset
            %                 K(nj) = n0*k_1pt(HCEF,JJ,4,h1,th1,.001);
            %                 nj = nj+1;
            %             end
            Err(ni) = sum(power( (K - Fv1{ni}), 2));
        end
        out = sum(Err);
    end

%% function for nonlinlsq of 4K Th-dept data
    function out = Func1(x0,xdata)
        
        t1 = 4;
        n0 = x0(1);
        HCEF1 = HCEFf(x0(2:7));
        JJ1   = JJf(x0(8:9));
        
        
        out = n0*[k_vTh(HCEF1,JJ1,t1,Hset(1),Thset,.001)...
            k_vTh(HCEF1,JJ1,t1,Hset(2),Thset,.001)...
            k_vTh(HCEF1,JJ1,t1,Hset(3),Thset,.001)...
            k_vTh(HCEF1,JJ1,t1,Hset(4),Thset,.001)...
            k_vTh(HCEF1,JJ1,t1,Hset(5),Thset,.001)...
            k_vTh(HCEF1,JJ1,t1,Hset(6),Thset,.001)...
            k_vTh(HCEF1,JJ1,t1,Hset(7),Thset,.001)...
            ];
    end

%% function for nonlinlsq of(4<=T<=70,H||ab,c)-dept data

    function out = Func2(x0,xdata)
        
        t1 = 70;
        n0 = 1e3*amp;
        HCEF1 = HCEFf(x0(2:7));
        JJ1   = JJf(x0(8:9));
        
        out = n0*[k_vH(HCEF1,JJ1,T(2),Hd_int{1}{2},pi/2,.001)/w0(2)...
            k_vH(HCEF1,JJ1,T(3),Hd_int{1}{3},pi/2,.001)/w0(3)...
            k_vH(HCEF1,JJ1,T(4),Hd_int{1}{4},pi/2,.001)/w0(4)...
            k_vH(HCEF1,JJ1,T(5),Hd_int{1}{5},pi/2,.001)/w0(5)...
            k_vH(HCEF1,JJ1,T(6),Hd_int{1}{6},pi/2,.001)/w0(6)...
            k_vH(HCEF1,JJ1,T(7),Hd_int{1}{7},pi/2,.001)/w0(7)...
            k_vH(HCEF1,JJ1,T(8),Hd_int{1}{8},pi/2,.001)/w0(8)...
            k_vH(HCEF1,JJ1,T(2),Hd_int{2}{2},0,.001)/w0(2)...
            k_vH(HCEF1,JJ1,T(3),Hd_int{2}{3},0,.001)/w0(3)...
            k_vH(HCEF1,JJ1,T(4),Hd_int{2}{4},0,.001)/w0(4)...
            k_vH(HCEF1,JJ1,T(5),Hd_int{2}{5},0,.001)/w0(5)...
            k_vH(HCEF1,JJ1,T(6),Hd_int{2}{6},0,.001)/w0(6)...
            k_vH(HCEF1,JJ1,T(7),Hd_int{2}{7},0,.001)/w0(7)...
            k_vH(HCEF1,JJ1,T(8),Hd_int{2}{8},0,.001)/w0(8)...
            ];
    end

    function out = Func3(x0,xdata)
        
        n0 = 1e3*amp;
        HCEF1 = HCEFf(x0(2:7));
        JJ1   = JJf(x0(8:9));
        
        out = n0*[k_vH(HCEF1,JJ1,T(4),Hd_int{1}{4},pi/2,.001)/w0(4)...
            k_vH(HCEF1,JJ1,T(5),Hd_int{1}{5},pi/2,.001)/w0(5)...
            k_vH(HCEF1,JJ1,T(6),Hd_int{1}{6},pi/2,.001)/w0(6)...
            k_vH(HCEF1,JJ1,T(7),Hd_int{1}{7},pi/2,.001)/w0(7)...
            k_vH(HCEF1,JJ1,T(8),Hd_int{1}{8},pi/2,.001)/w0(8)...
            k_vH(HCEF1,JJ1,T(4),Hd_int{2}{4},0,.001)/w0(4)...
            k_vH(HCEF1,JJ1,T(5),Hd_int{2}{5},0,.001)/w0(5)...
            k_vH(HCEF1,JJ1,T(6),Hd_int{2}{6},0,.001)/w0(6)...
            k_vH(HCEF1,JJ1,T(7),Hd_int{2}{7},0,.001)/w0(7)...
            k_vH(HCEF1,JJ1,T(8),Hd_int{2}{8},0,.001)/w0(8)...
            ];
        
    end

    function out = Func4(x0,xdata)
        
        t1 = 4;
        n0 = x0(1);
        HCEF1 = HCEFf(x0(2:7));
        JJ1   = JJf(x0(8:9));
        
        out = n0*[k_vH(HCEF1,JJ1,T(2),Hd_int{2}{2},0,.001)...
            k_vH(HCEF1,JJ1,T(3),Hd_int{2}{3},0,.001)...
            k_vH(HCEF1,JJ1,T(4),Hd_int{2}{4},0,.001)...
            k_vH(HCEF1,JJ1,T(5),Hd_int{2}{5},0,.001)...
            k_vH(HCEF1,JJ1,T(6),Hd_int{2}{6},0,.001)...
            k_vH(HCEF1,JJ1,T(7),Hd_int{2}{7},0,.001)...
            k_vH(HCEF1,JJ1,T(8),Hd_int{2}{8},0,.001)...
            ];
    end

    function out = Func8(x0,xdata)
        
        t1 = 4;
        n0 = x0(1)*1e3;
        HCEF1 = HCEFf(x0(2:7));
        JJ1   = JJf(x0(8:9));
        
        out = n0*[k_vH(HCEF1,JJ1,T(2),Hd_int{1}{2},pi/2,.001)/w0(2)...
            k_vH(HCEF1,JJ1,T(8),Hd_int{1}{8},pi/2,.001)/w0(8)...
            k_vH(HCEF1,JJ1,T(2),Hd_int{2}{2},0,.001)/w0(2)...
            k_vH(HCEF1,JJ1,T(8),Hd_int{2}{8},0,.001)/w0(8)...
            k_vTh(HCEF1,JJ1,4,28,Thset,.001)/w0(2)...
            ];
    end

    function out = Func9(x0,xdata)
        
        t1 = 4;
        n0 = x0(1)*1e3;
        HCEF1 = HCEFf(x0(2:7));
        JJ1   = JJf(x0(8:9));
        
        out = n0*[k_vH(HCEF1,JJ1,T(8),Hd_int{1}{8},pi/2,.001)/w0(8)...
            k_vH(HCEF1,JJ1,T(8),Hd_int{2}{8},0,.001)/w0(8)...
            ];
    end

    function out = Func_test(x0,xdata)
        
        n0 = x0(1);
        HCEF1 = HCEFf(x0(2:7));
        JJ1   = JJf(x0(8:9));
        
        out = n0*[k_vTh(HCEF1,JJ1,70,60,xdata,.001)];
    end

    function out = Func_test2(x0,xdata)
        
        n0 = x0(1)*1e3;
        HCEF1 = HCEFf(x0(2:7));
        JJ1   = JJf(x0(8:9));
        
        out = n0*[k_vTh(HCEF1,JJ1,4,28,X,.001)/w0(2)...
            k_vTh(HCEF1,JJ1,4,56,X,.001)/w0(2)...
            k_vTh(HCEF1,JJ1,70,28,X,.001)/w0(8)...
            k_vTh(HCEF1,JJ1,70,56,X,.001)/w0(8)];
        
    end

    function out = Func10(x0,xdata)
        
        t1 = 70;
        n0 = x0(1);
        HCEF1 = HCEFf(x0(2:7));
        JJ1   = JJf(x0(8:9));
        
        out = n0*[k_vH(HCEF1,JJ1,T(5),Hd_int{1}{5},pi/2,.001)...
            k_vH(HCEF1,JJ1,T(6),Hd_int{1}{6},pi/2,.001)...
            k_vH(HCEF1,JJ1,T(7),Hd_int{1}{7},pi/2,.001)...
            k_vH(HCEF1,JJ1,T(8),Hd_int{1}{8},pi/2,.001)...
            k_vH(HCEF1,JJ1,T(5),Hd_int{2}{5},0,.001)...
            k_vH(HCEF1,JJ1,T(6),Hd_int{2}{6},0,.001)...
            k_vH(HCEF1,JJ1,T(7),Hd_int{2}{7},0,.001)...
            k_vH(HCEF1,JJ1,T(8),Hd_int{2}{8},0,.001)...
            ];
    end

%% Analysis
switch option
    
    
    case 'CEF Spec alt'
        
        
        
        
        
        H = 0:.2:18;
        
        HCEF = HCEFf(paramsB);
        JJ = JJf(paramsJJ);
        Jxx = paramsJJ(1);
        Jzz = paramsJJ(2);
        
        [P1,D1] = eig(HCEF + Jz*1e-10);
        
        En = real([D1(1,1) D1(2,2) D1(3,3) D1(4,4) D1(5,5) D1(6,6) D1(7,7) D1(8,8)]);
        Ev1 = sort(En);
        
        E_0 = Ev1(1);
        
        ind1 = find(En==Ev1(1));
        ind2 = find(En==Ev1(2));
        
        ev1 = P1(:,ind1);
        ev2 = P1(:,ind2);
        
        
        gzz0 = 2*gJ*abs(ev1'*Jz*ev1);
        gxx0 = 2*gJ*abs(ev1'*Jx*ev2);
        
        
        Tvals = [2,4,6,8,10,14,20,30,50,70];%
        colz = varycolor(length(Tvals));
        
        figure; hold on
        
        
        j = 1;
        for t = Tvals%
            n = 1;
            for h = H
                
                hv1 = [h,0,0];
                mv0 = M(HCEF, JJ, hv1, t);
                %HMFa = HMF(HCEF, JJ, hv1, mv0);
                
                Mx(n) = mv0(1);
                HZeeman = -gJ*muB*h*Jx;
                [P1,D1] = eig(HCEF + HZeeman);
                
                for q = 1:8
                    Eigenval(q) = real(D1(q,q));
                end
                Eigenval = sort(Eigenval);
                for q = 1:8
                    E1{q}(n) = Eigenval(q);
                end
                
                hv1 = [0,0,h];
                mv0 = M(HCEF, JJ, hv1, t);
                %HMFc = HMF(HCEF, JJ, hv1, mv0);
                
                Mz(n) = mv0(3);
                HZeeman = -gJ*muB*h*Jz;
                [P1,D1] = eig(HCEF + HZeeman);
                
                for q = 1:8
                    Eigenval(q) = real(D1(q,q));
                end
                Eigenval = sort(Eigenval);
                for q = 1:8
                    E2{q}(n) = Eigenval(q);
                end
                
                n = n + 1;
            end
            
            for q = 1:2
                subplot(1,2,1); hold on;
                plot(H + Jxx*A1*Mx,E1{q}-E_0,'color',colz(j,:))
                
                subplot(1,2,2); hold on;
                plot(H + Jzz*A1*Mz,E2{q}-E_0,'color',colz(j,:))
            end
            
            Dpm_a{j} = E1{2}-E1{1};
            Dpm_c{j} = E2{2}-E2{1};
            
            output.Dpm_a{j} = Interp1NonUnique(H + Jxx*A1*Mx, Dpm_a{j}, H);
            output.Dpm_c{j} = Interp1NonUnique(H + Jzz*A1*Mz, Dpm_c{j}, H);
            
            j = j+1;
        end
        
        
        
        
        
        
        subplot(1,2,1);
        xlabel('\mu_0H [T]'); ylabel('\omega [meV]'); grid on; box on;
        title(['H||ab,   g_{\perp} = ' num2str(gxx0)])
        
        subplot(1,2,2);
        xlabel('\mu_0H [T]'); ylabel('\omega [meV]'); grid on; box on;
        title(['H||c,   g_{||} = ' num2str(gzz0)])
        
    case 'CEF Spectrum'
        
        t = 4;%
        
        H = 0:.2:18;
        
        HCEF = HCEFf(paramsB);
        JJ = JJf(paramsJJ);
        
        [P1,D1] = eig(HCEF + Jz*1e-10);
        
        En = real([D1(1,1) D1(2,2) D1(3,3) D1(4,4) D1(5,5) D1(6,6) D1(7,7) D1(8,8)]);
        Ev1 = sort(En);
        
        E_0 = Ev1(1);
        
        ind1 = find(En==Ev1(1));
        ind2 = find(En==Ev1(2));
        
        ev1 = P1(:,ind1);
        ev2 = P1(:,ind2);
        
        
        gzz0 = 2*gJ*abs(ev1'*Jz*ev1);
        gxx0 = 2*gJ*abs(ev1'*Jx*ev2);
        
        
        n = 1;
        for h = H
            
            hv1 = [h,0,0];
            mv0 = M(HCEF, JJ, hv1, t);
            HMFa = HMF(HCEF, JJ, hv1, mv0);
            
            [P1,D1] = eig(HMFa);
            
            for q = 1:8
                Eigenval(q) = real(D1(q,q));
            end
            Eigenval = sort(Eigenval);
            for q = 1:8
                E1{q}(n) = Eigenval(q);
            end
            
            hv1 = [0,0,h];
            mv0 = M(HCEF, JJ, hv1, t);
            HMFc = HMF(HCEF, JJ, hv1, mv0);
            
            [P1,D1] = eig(HMFc);
            
            for q = 1:8
                Eigenval(q) = real(D1(q,q));
            end
            Eigenval = sort(Eigenval);
            for q = 1:8
                E2{q}(n) = Eigenval(q);
            end
            
            n = n + 1;
        end
        
        figure; hold on
        colz = {'k','k','r','r','b','b','c','c'};
        for q = 1:8
            subplot(1,2,1); hold on;
            plot(H,E1{q}-E_0,'color',colz{q})
            
            subplot(1,2,2); hold on;
            plot(H,E2{q}-E_0,'color',colz{q})
        end
        
        
        
        subplot(1,2,1);
        xlabel('\mu_0H [T]'); ylabel('\omega [meV]'); grid on; box on;
        title(['H||ab,   g_{\perp} = ' num2str(gxx0)])
        
        subplot(1,2,2);
        xlabel('\mu_0H [T]'); ylabel('\omega [meV]'); grid on; box on;
        title(['H||c,   g_{||} = ' num2str(gzz0)])
        
    case 'GS CEF Spectrum finite T'
        
        figure; hold on
        
        Tvals = [1.5,3,4,5,6,8,10,12,14,17,20,30,50,70];
        %Tvals = [2,4,6,8,10,14,20,30,50,70];%
        colz = varycolor(length(Tvals));
        
        H = 0:.2:18;
        
        HCEF = HCEFf(paramsB);
        JJ = JJf(paramsJJ);
        
        [P1,D1] = eig(HCEF + Jz*1e-10);
        
        En = real([D1(1,1) D1(2,2) D1(3,3) D1(4,4) D1(5,5) D1(6,6) D1(7,7) D1(8,8)]);
        Ev1 = sort(En);
        
        E_0 = Ev1(1);
        
        ind1 = find(En==Ev1(1));
        ind2 = find(En==Ev1(2));
        
        ev1 = P1(:,ind1);
        ev2 = P1(:,ind2);
        
        
        gzz0 = 2*gJ*abs(ev1'*Jz*ev1);
        gxx0 = 2*gJ*abs(ev1'*Jx*ev2);
        
        j = 1;
        for t = Tvals
            n = 1;
            for h = H
                
                hv1 = [h,0,0];
                mv0 = M(HCEF, JJ, hv1, t);
                HMFa = HMF(HCEF, JJ, hv1, mv0);
                
                [P1,D1] = eig(HMFa);
                
                for q = 1:8
                    Eigenval(q) = real(D1(q,q));
                end
                Eigenval = sort(Eigenval);
                for q = 1:8
                    E1{q}(n) = Eigenval(q);
                end
                
                hv1 = [0,0,h];
                mv0 = M(HCEF, JJ, hv1, t);
                HMFc = HMF(HCEF, JJ, hv1, mv0);
                
                [P1,D1] = eig(HMFc);
                
                for q = 1:8
                    Eigenval(q) = real(D1(q,q));
                end
                Eigenval = sort(Eigenval);
                for q = 1:8
                    E2{q}(n) = Eigenval(q);
                end
                
                n = n + 1;
            end
            
%             for q = 1:2
%                 subplot(1,2,1); hold on;
%                 plot(H,E1{q}-E_0,'color',colz(j,:))
%                 
%                 subplot(1,2,2); hold on;
%                 plot(H,E2{q}-E_0,'color',colz(j,:))
%             end
            subplot(1,2,1); hold on;
            D_1 = E1{2}-E1{1};
            plot(H,D_1)
            output.D1{j} = D_1;
            
            subplot(1,2,2); hold on;
            D_2 = E2{2}-E2{1};
            plot(H,D_2)
            output.D2{j} = D_2;
            

            
            j = j + 1;
        end
        
        output.H = H;
        Z3 = zeros(3,3);
        
        
        n = 1;
        for h = H
            hv1 = [h,0,0];
            mv0 = M(HCEF, Z3, hv1, t);
            HMFa = HMF(HCEF, Z3, hv1, mv0);
            
            [P1,D1] = eig(HMFa);
            
            for q = 1:8
                Eigenval(q) = real(D1(q,q));
            end
            Eigenval = sort(Eigenval);
            for q = 1:8
                E1{q}(n) = Eigenval(q);
            end
            
            hv1 = [0,0,h];
            mv0 = M(HCEF, Z3, hv1, t);
            HMFc = HMF(HCEF, Z3, hv1, mv0);
            
            [P1,D1] = eig(HMFc);
            
            for q = 1:8
                Eigenval(q) = real(D1(q,q));
            end
            Eigenval = sort(Eigenval);
            for q = 1:8
                E2{q}(n) = Eigenval(q);
            end
        end
        
        
%         
%         for q = 1
%             subplot(1,2,1); hold on;
%             plot(H,E1{2}-E1{1})
%             
%             subplot(1,2,2); hold on;
%             plot(H,E2{2}-E2{1})
%         end
        
        
        subplot(1,2,1);
        xlabel('\mu_0H [T]'); ylabel('\omega [meV]'); grid on; box on; xlim([0,18]);
        title(['H||ab,   g_{\perp} = ' num2str(gxx0)])
        set(gca,'fontsize',20);
        
        subplot(1,2,2);
        xlabel('\mu_0H [T]'); ylabel('\omega [meV]'); grid on; box on; xlim([0,18]);
        title(['H||c,   g_{||} = ' num2str(gzz0)])
        set(gca,'fontsize',20);
        
    case 'GS CEF gap finite T'
        
        figure; hold on
        
        Tvals = [2,4,6,8,10,14,20,30,50,70];%
        colz = varycolor(length(Tvals));
        
        H = 0:.2:18;
        
        HCEF = HCEFf(paramsB);
        JJ = JJf(paramsJJ);
        
        [P1,D1] = eig(HCEF + Jz*1e-10);
        
        En = real([D1(1,1) D1(2,2) D1(3,3) D1(4,4) D1(5,5) D1(6,6) D1(7,7) D1(8,8)]);
        Ev1 = sort(En);
        
        E_0 = Ev1(1);
        
        ind1 = find(En==Ev1(1));
        ind2 = find(En==Ev1(2));
        
        ev1 = P1(:,ind1);
        ev2 = P1(:,ind2);
        
        
        gzz0 = 2*gJ*abs(ev1'*Jz*ev1);
        gxx0 = 2*gJ*abs(ev1'*Jx*ev2);
        
        j = 1;
        for t = Tvals
            n = 1;
            for h = H
                
                hv1 = [h,0,0];
                mv0 = M(HCEF, JJ, hv1, t);
                HMFa = HMF(HCEF, JJ, hv1, mv0);
                
                [P1,D1] = eig(HMFa);
                
                for q = 1:8
                    Eigenval(q) = real(D1(q,q));
                end
                Eigenval = sort(Eigenval);
                for q = 1:8
                    E1{q}(n) = Eigenval(q);
                end
                
                hv1 = [0,0,h];
                mv0 = M(HCEF, JJ, hv1, t);
                HMFc = HMF(HCEF, JJ, hv1, mv0);
                
                [P1,D1] = eig(HMFc);
                
                for q = 1:8
                    Eigenval(q) = real(D1(q,q));
                end
                Eigenval = sort(Eigenval);
                for q = 1:8
                    E2{q}(n) = Eigenval(q);
                end
                
                n = n + 1;
            end
            
            Dpm_a{j} = E1{2}-E1{1};
            subplot(1,2,1); hold on;
            plot(H, Dpm_a{j}, 'color',colz(j,:))
            
            Dpm_c{j} = E2{2}-E2{1};
            subplot(1,2,2); hold on;
            plot(H, Dpm_c{j}, 'color',colz(j,:))
            
            j = j + 1;
        end
        
        output.Dpm_a = Dpm_a;
        output.Dpm_c = Dpm_c;
        
        Z3 = zeros(3,3);
        
        
        n = 1;
        for h = H
            hv1 = [h,0,0];
            mv0 = M(HCEF, Z3, hv1, t);
            HMFa = HMF(HCEF, Z3, hv1, mv0);
            
            [P1,D1] = eig(HMFa);
            
            for q = 1:8
                Eigenval(q) = real(D1(q,q));
            end
            Eigenval = sort(Eigenval);
            for q = 1:8
                E1{q}(n) = Eigenval(q);
            end
            
            hv1 = [0,0,h];
            mv0 = M(HCEF, Z3, hv1, t);
            HMFc = HMF(HCEF, Z3, hv1, mv0);
            
            [P1,D1] = eig(HMFc);
            
            for q = 1:8
                Eigenval(q) = real(D1(q,q));
            end
            Eigenval = sort(Eigenval);
            for q = 1:8
                E2{q}(n) = Eigenval(q);
            end
        end
        
        
        
        
        subplot(1,2,1); hold on;
        plot(H,E1{2}-E1{1},'y--')
        
        subplot(1,2,2); hold on;
        plot(H,E2{2}-E2{1},'y--')
        
        
        subplot(1,2,1);
        xlabel('\mu_0H [T]'); ylabel('\omega [meV]'); grid on; box on; xlim([0,18]);
        title(['H||ab,   g_{\perp} = ' num2str(gxx0)])
        set(gca,'fontsize',20);
        
        subplot(1,2,2);
        xlabel('\mu_0H [T]'); ylabel('\omega [meV]'); grid on; box on; xlim([0,18]);
        title(['H||c,   g_{||} = ' num2str(gzz0)])
        set(gca,'fontsize',20);
        
    case 'CEF Spectrum ac plane'
        
        t = 4;%
        
        H = 0:1:60;
        
        HCEF = HCEFf(paramsB);
        JJ = JJf(paramsJJ);
        
        
        
        
        
        [P1,D1] = eig(HCEF + Jz*1e-10);
        
        En = real([D1(1,1) D1(2,2) D1(3,3) D1(4,4) D1(5,5) D1(6,6) D1(7,7) D1(8,8)]);
        Ev1 = sort(En);
        
        E_0 = Ev1(1);
        
        ind1 = find(En==Ev1(1));
        ind2 = find(En==Ev1(2));
        
        ev1 = P1(:,ind1);
        ev2 = P1(:,ind2);
        
        
        gzz0 = 2*gJ*abs(ev1'*Jz*ev1);
        gxx0 = 2*gJ*abs(ev1'*Jx*ev2);
        
        
        n = 1; Th = linspace(0,pi,1e2); h = 28;
        for th = Th
            
            hv1 = [h*sin(th),0,h*cos(th)];
            mv0 = M(HCEF, JJ, hv1, t);
            HMFa = HMF(HCEF, JJ, hv1, mv0);
            
            [P1,D1] = eig(HMFa);
            
            for q = 1:8
                Eigenval(q) = real(D1(q,q));
            end
            Eigenval = sort(Eigenval);
            for q = 1:8
                E1{q}(n) = Eigenval(q);
            end
            
            
            n = n + 1;
        end
        
        figure; hold on
        colz = {'k','k','r','r','b','b','c','c'};
        for q = 1:8
            plot(Th,E1{q}-E_0,'color',colz{q})
        end
        
    case 'Chi inv v T'
        
        HCEF = HCEFf(paramsB);
        JJ = JJf(paramsJJ);
        
        [P1,D1] = eig(HCEF);
        
        En = real([D1(1,1) D1(2,2) D1(3,3) D1(4,4) D1(5,5) D1(6,6) D1(7,7) D1(8,8)]);
        
        Ev_1 = sort(En);
        
        for k = 0:3
            ind1 = find(En==Ev_1(2*k+1));
            ind2 = find(En==Ev_1(2*k+2));
            
            ev1 = P1(:,ind1(1));
            ev2 = P1(:,ind2(end));
            
            paramsJME(2*k+1) = abs(ev1'*Jx*ev2);
            paramsJME(2*k+2) = abs(ev1'*Jz*ev1);
            
        end
        
        Dm(1) = Ev_1(3)-Ev_1(1);
        Dm(2) = Ev_1(5)-Ev_1(3);
        Dm(3) = Ev_1(7)-Ev_1(5);
        
        Gx10 = power(paramsJME(3)/paramsJME(1),2);
        Gx20 = power(paramsJME(5)/paramsJME(1),2);
        Gx30 = power(paramsJME(7)/paramsJME(1),2);
        
        Gz10 = power(paramsJME(4)/paramsJME(2),2);
        Gz20 = power(paramsJME(6)/paramsJME(2),2);
        Gz30 = power(paramsJME(8)/paramsJME(2),2);
        
        n1 = 1; Tv = 4:300;
        for t = Tv
            b = 1./(kB*t);
            P1 = 1 + exp(-b*Dm(1)) + exp(-b*(Dm(1) + Dm(2))) + exp(-b*(Dm(1) + Dm(2) + Dm(3)));
            Q1 = power(paramsJME(1),2)*(1 + Gx10*exp(-b*Dm(1)) + Gx20*exp(-b*(Dm(1) + Dm(2))) + Gx30*exp(-b*(Dm(1) + Dm(2) + Dm(3))));
            
            Tp1(n1) = t*P1./Q1;
            
            P2 = 1 + exp(-b*Dm(1)) + exp(-b*(Dm(1) + Dm(2))) + exp(-b*(Dm(1) + Dm(2) + Dm(3)));
            Q2 = power(paramsJME(2),2)*(1 + Gz10*exp(-b*Dm(1)) + Gz20*exp(-b*(Dm(1) + Dm(2))) + Gz30*exp(-b*(Dm(1) + Dm(2) + Dm(3))));
            
            Tp2(n1) = t*P2./Q2;
            
            n1 = n1+1;
        end
        
        C = 3*.125*gJ*gJ;
        
        figure; hold on;
        plot(Tv, C./(Tp1 + 6*paramsJJ(1)) , 'r')
        plot(Tv, C./(Tp2 + 6*paramsJJ(2)) , 'b')
        
        
        plot(Tv, 1./(C./(Tp1 + 6*paramsJJ(1))) , 'r')
        plot(Tv,1./( C./(Tp2 + 6*paramsJJ(2))) , 'b')
        
    case 'Chi v T alt'
        
        HCEF = HCEFf(paramsB);
        JJ = JJf(paramsJJ);
        
        C1 = 3*.125*gJ/(muB./kB);
        
        h = 1;
        %Tv = 2:300;
        Tv = 100:100:1e4;
        
        % Calculate Partion Functions : expm (not exp) is crucial for proper
        % exponentiation of matrices.
        n = 1;
        for t = Tv
            %M(HCEFm, JJab, hv, t)
            mV1 = M(HCEF, JJ, [h,0,0], t);
            Xa(n) = C1*mV1(1)./h;
            %Z1(n) = ZMF(t, [h,0,0], mV1, JJ, paramsD, paramsJME);
            
            mV2 = M(HCEF, JJ, [0,0,h], t);
            Xc(n) = C1*mV2(3)./h;
            %Z2(n) = ZMF(t, [0,0,h], mV2, JJ, paramsD, paramsJME);
            n = n + 1;
        end
        
        
        figure; hold on;
        plot(Tv,1./Xa,'r')
        plot(Tv,1./Xc,'b')
        
    case 'Chi v T'
        
        HCEF = HCEFf(paramsB);
        JJ = JJf(paramsJJ);
        
        h = .1;
        T = 2:300;
        
        HZeeman_x = -muB*h*Jx;
        HZeeman_z = -muB*h*Jz;
        
        % Calculate Partion Functions : expm (not exp) is crucial for proper
        % exponentiation of matrices.
        n = 1;
        for t = T
            b = 1./(kB*t);
            
            mv2 = M(HCEF, JJ, [0,0,0], t);
            HMF2 = HMF(HCEF, JJ, [0,0,0], mv2);
            [P2,D2] = eig(HMF2);
            Ev2 = sort(real([D2(1,1) D2(2,2) D2(3,3) D2(4,4) D2(5,5) D2(6,6) D2(7,7) D2(8,8)]));
            E02(n) = Ev2(1);
            
            Z0(n) = trace(expm(-b*(HMF2 - E02(n)*ID)));
            
            mv2 = M(HCEF, JJ, [h,0,0], t);
            HMF2 = HMF(HCEF, JJ, [h,0,0], mv2);
            
            Z1(n) = trace(expm(-b*(HMF2 - E02(n)*ID)));
            
            mv2 = M(HCEF, JJ, [0,0,h], t);
            HMF2 = HMF(HCEF, JJ, [0,0,h], mv2);
            
            Z2(n) = trace(expm(-b*(HMF2 - E02(n)*ID)));
            n = n + 1;
        end
        
        
        F0 = (kB*T).*log(Z0);
        F1 = (kB*T).*log(Z1);
        F2 = (kB*T).*log(Z2);
        
        Xa = (F1-F0)/h/h;
        Xc = (F2-F0)/h/h;
        
        C = 0.0422*108.3/88.14;
        
        figure; hold on;
        plot(T,C./Xa,'r')
        plot(T,C./Xc,'b')
        
    case 'Transverse Susceptibility'
        
        HCEF = HCEFf(paramsB);
        JJ = JJf(paramsJJ);
        
        h = .01;
        T = 2:300;
        
        
        % Calculate Partion Functions : expm (not exp) is crucial for proper
        % exponentiation of matrices.
        HT= [0:60]; Temps = [4,6,12,20,30,50,70];
        
        q = 1;
        for t = Temps
            n = 1;
            b = 1./(kB*t);
            for hT = HT
                
                mv2 = M(HCEF, JJ, [0,0,hT], t);
                HMF2 = HMF(HCEF, JJ, [0,0,hT], mv2);
                [P2,D2] = eig(HMF2);
                Ev2 = sort(real([D2(1,1) D2(2,2) D2(3,3) D2(4,4) D2(5,5) D2(6,6) D2(7,7) D2(8,8)]));
                E02(n) = Ev2(1);
                
                Z0a(n) = trace(expm(-b*(HMF2 - E02(n)*ID)));
                
                mv2 = M(HCEF, JJ, [h,0,hT], t);
                HMF2 = HMF(HCEF, JJ, [h,0,hT], mv2);
                
                Z1a(n) = trace(expm(-b*(HMF2 - E02(n)*ID)));
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                mv2 = M(HCEF, JJ, [hT,0,0], t);
                HMF2 = HMF(HCEF, JJ, [hT,0,0], mv2);
                [P2,D2] = eig(HMF2);
                Ev2 = sort(real([D2(1,1) D2(2,2) D2(3,3) D2(4,4) D2(5,5) D2(6,6) D2(7,7) D2(8,8)]));
                E02(n) = Ev2(1);
                
                Z0c(n) = trace(expm(-b*(HMF2 - E02(n)*ID)));
                
                mv2 = M(HCEF, JJ, [hT,0,h], t);
                HMF2 = HMF(HCEF, JJ, [hT,0,h], mv2);
                
                Z1c(n) = trace(expm(-b*(HMF2 - E02(n)*ID)));
                
                n = n + 1;
            end
            F0a = (kB*t).*log(Z0a);
            F1a = (kB*t).*log(Z1a);
            
            F0c = (kB*t).*log(Z0c);
            F1c = (kB*t).*log(Z1c);
            
            Xc_Tr(q,:) = (F1c-F0c)/h/h;
            Xa_Tr(q,:) = (F1a-F0a)/h/h;
            
            q = q + 1;
        end
        
        C = 0.0422*108.3/88.14;
        
        figure; hold on;
        subplot(1,2,1); grid on; box on; hold on; set(gca,'fontsize',20)
        xlabel('H_c [T]'); ylabel('\chi_{ab}')
        subplot(1,2,2); grid on; box on; hold on; set(gca,'fontsize',20)
        xlabel('H_a [T]'); ylabel('\chi_{c}')
        
        for q = 1:length(Temps)
            
            subplot(1,2,1);
            plot(HT,Xa_Tr(q,:)/C,'displayname',['T = ' num2str(Temps(q)) ' K'])
            
            subplot(1,2,2);
            plot(HT,Xc_Tr(q,:)/C,'displayname',['T = ' num2str(Temps(q)) ' K'])
        end
        
        output.Xa = Xa_Tr;
        output.Xc = Xa_Tr;
        
    case 'Finite Field Susceptibility'
        
        HCEF = HCEFf(paramsB);
        JJ = JJf(paramsJJ);
        
        h = .01;
        T = 2:300;
        
        
        % Calculate Partion Functions : expm (not exp) is crucial for proper
        % exponentiation of matrices.
        HT= [0:1:60]; Temps = [4,6,12,20,30,50,70];
        
        q = 1;
        for t = Temps
            n = 1;
            b = 1./(kB*t);
            for hT = HT
                
                mv2_00(n,:) = M(HCEF, JJ, [hT,0,0], t);
                mv_a0(n) = mv2_00(n,1);
                
                mv2_01(n,:) = M(HCEF, JJ, [hT+h,0,0], t);
                mv_a1(n) = mv2_01(n,1);
                
                mv2_10(n,:) = M(HCEF, JJ, [0,0,hT], t);
                mv_c0(n) = mv2_10(n,3);
                
                mv2_11(n,:) = M(HCEF, JJ, [0,0,hT+h], t);
                mv_c1(n) = mv2_11(n,3);
                
                n = n + 1;
            end
            %display(sum(mv_a0==0));
            %display(sum(mv_c0==0));
            Xa_Tr(q,:) = (mv_a1 - mv_a0)/h;
            Xc_Tr(q,:) = (mv_c1 - mv_c0)/h;
            
            q = q + 1;
        end
        
        %C = 0.0422*108.3/88.14;
        
        figure; hold on;
        subplot(1,2,1); grid on; box on; hold on; set(gca,'fontsize',20)
        xlabel('H_c [T]'); ylabel('\chi_{ab}')
        subplot(1,2,2); grid on; box on; hold on; set(gca,'fontsize',20)
        xlabel('H_a [T]'); ylabel('\chi_{c}')
        
        for q = 1:length(Temps)
            
            subplot(1,2,1);
            plot(HT,Xa_Tr(q,:),'displayname',['T = ' num2str(Temps(q)) ' K'])
            
            subplot(1,2,2);
            plot(HT,Xc_Tr(q,:),'displayname',['T = ' num2str(Temps(q)) ' K'])
        end
        
        output.Xa = Xa_Tr;
        output.Xc = Xc_Tr;
        
    case 'Finite Field Mag'
        
        HCEF = HCEFf(paramsB);
        JJ = JJf(paramsJJ);
        
        h = .01;
        T = 2:300;
        
        
        % Calculate Partion Functions : expm (not exp) is crucial for proper
        % exponentiation of matrices.
        HT= [0:1:60]; Temps = [4,6,12,20,30,50,70];
        
        q = 1;
        for t = Temps
            n = 1;
            b = 1./(kB*t);
            for hT = HT
                
                mv2_00(n,:) = M(HCEF, JJ, [hT,0,0], t);
                mv_a0(n) = mv2_00(n,1);
                
                
                mv2_10(n,:) = M(HCEF, JJ, [0,0,hT], t);
                mv_c0(n) = mv2_10(n,3);
                
                
                n = n + 1;
            end
            %display(sum(mv_a0==0));
            %display(sum(mv_c0==0));
            Xa_Tr(q,:) = mv_a0;
            Xc_Tr(q,:) = mv_c0;
            
            q = q + 1;
        end
        
        %C = 0.0422*108.3/88.14;
        
        figure; hold on;
        subplot(1,2,1); grid on; box on; hold on; set(gca,'fontsize',20)
        xlabel('H_c [T]'); ylabel('\chi_{ab}')
        subplot(1,2,2); grid on; box on; hold on; set(gca,'fontsize',20)
        xlabel('H_a [T]'); ylabel('\chi_{c}')
        
        for q = 1:length(Temps)
            
            subplot(1,2,1);
            plot(HT,Xa_Tr(q,:),'displayname',['T = ' num2str(Temps(q)) ' K'])
            
            subplot(1,2,2);
            plot(HT,Xc_Tr(q,:),'displayname',['T = ' num2str(Temps(q)) ' K'])
        end
        
        output.Ma = Xa_Tr;
        output.Mc = Xc_Tr;
        
    case 'Transverse Susceptibility alt'
        
        HCEF = HCEFf(paramsB);
        JJ = JJf(paramsJJ);
        
        h = .01;
        T = 2:300;
        
        
        % Calculate Partion Functions : expm (not exp) is crucial for proper
        % exponentiation of matrices.
        HT= [0:1:60]; Temps = [4,6,12,20,30,50,70];
        
        q = 1;
        for t = Temps
            n = 1;
            b = 1./(kB*t);
            for hT = HT
                
                mv2_00(n,:) = M(HCEF, JJ, [0,0,hT], t);
                mv_a0(n) = mv2_00(n,1);
                
                mv2_01(n,:) = M(HCEF, JJ, [h,0,hT], t);
                mv_a1(n) = mv2_01(n,1);
                
                mv2_10(n,:) = M(HCEF, JJ, [hT,0,0], t);
                mv_c0(n) = mv2_10(n,3);
                
                mv2_11(n,:) = M(HCEF, JJ, [hT,0,h], t);
                mv_c1(n) = mv2_11(n,3);
                
                n = n + 1;
            end
            %display(sum(mv_a0==0));
            %display(sum(mv_c0==0));
            Xa_Tr(q,:) = (mv_a1 - mv_a0)/h;
            Xc_Tr(q,:) = (mv_c1 - mv_c0)/h;
            
            q = q + 1;
        end
        
        %C = 0.0422*108.3/88.14;
        
        figure; hold on;
        subplot(1,2,1); grid on; box on; hold on; set(gca,'fontsize',20)
        xlabel('H_c [T]'); ylabel('\chi_{ab}')
        subplot(1,2,2); grid on; box on; hold on; set(gca,'fontsize',20)
        xlabel('H_a [T]'); ylabel('\chi_{c}')
        
        for q = 1:length(Temps)
            
            subplot(1,2,1);
            plot(HT,Xa_Tr(q,:),'displayname',['T = ' num2str(Temps(q)) ' K'])
            
            subplot(1,2,2);
            plot(HT,Xc_Tr(q,:),'displayname',['T = ' num2str(Temps(q)) ' K'])
        end
        
        output.Xa = Xa_Tr;
        output.Xc = Xc_Tr;
        
    case '(M dot H) calc'
        
        HCEF = HCEFf(paramsB);
        JJ = JJf(paramsJJ);
        
        HT= [4:4:60]; Temps = [4,6,12,20,30,50,70];
        
        dth = .001;
        figure; hold on;
        q = 1; k = zeros(size(HT)); FF1 = zeros(1,5);
        for Th = [pi/2,0]
            subplot(1,2,q); hold on;
            for t = Temps
                n1 = 1;
                for h = HT
                    n2 = 1;
                    for th = [(Th-2*dth):dth:(Th+2*dth)]
                        hv1 = [h*sin(th),0,h*cos(th)];
                        mv1 = M(HCEF, JJ, hv1, t);
                        tau(n2) = -muB*gJ*(mv1(1)*hv1(3) - mv1(3)*hv1(1));
                        n2 = n2+1;
                    end
                    k(n1) = (C1(1)*tau(1) + C1(2)*tau(2) + C1(3)*tau(3) ...
                        + C1(4)*tau(4) + C1(5)*tau(5) )/dth;
                    
                    n1 = n1 + 1;
                end
                plot([0,HT], [0,k]);
            end
            q = q+1;
        end
        
    case 'k alt'
        
        HCEF = HCEFf(paramsB);
        JJ = JJf(paramsJJ);
        
        h = .01;
        
        HT= [4:4:60]; Temps = [4,6,12,20,30,50,70];
        
        init1N = zeros(size(HT));
        initN3 = zeros(length(HT),3);
        initMN = zeros(length(Temps),length(HT));
        
        mv2_00 = initN3;
        mv2_01 = initN3;
        mv2_10 = initN3;
        mv2_11 = initN3;
        mv_00  = initN3;
        mv_01  = initN3;
        
        mv_a0 = init1N;
        mv_a1 = init1N;
        mv_c0 = init1N;
        mv_c1 = init1N;
        mv_a2 = init1N;
        mv_c2 = init1N;
        
        Xa_Tr = initMN;
        Xc_Tr = initMN;
        Ma    = initMN;
        Mc    = initMN;
        
        q = 1;
        for t = Temps
            n = 1;
            b = 1./(kB*t);
            for hT = HT
                
                mv2_00(n,:) = M(HCEF, JJ, [0,0,hT], t);
                mv_a0(n) = mv2_00(n,1);
                
                mv2_01(n,:) = M(HCEF, JJ, [h,0,hT], t);
                mv_a1(n) = mv2_01(n,1);
                
                mv2_10(n,:) = M(HCEF, JJ, [hT,0,0], t);
                mv_c0(n) = mv2_10(n,3);
                
                mv2_11(n,:) = M(HCEF, JJ, [hT,0,h], t);
                mv_c1(n) = mv2_11(n,3);
                
                
                mv_00(n,:) = M(HCEF, JJ, [hT,0,0], t);
                mv_a2(n) = mv_00(n,1);
                
                mv_01(n,:) = M(HCEF, JJ, [0,0,hT], t);
                mv_c2(n) = mv_01(n,3);
                
                n = n + 1;
            end
            
            Xa_Tr(q,:) = (mv_a1 - mv_a0)/h;
            Xc_Tr(q,:) = (mv_c1 - mv_c0)/h;
            
            Ma(q,:) = mv_a2;
            Mc(q,:) = mv_c2;
            
            q = q + 1;
        end
        
        %C = 0.0422*108.3/88.14;
        
        figure; hold on;
        subplot(1,2,1); grid on; box on; hold on; set(gca,'fontsize',20)
        xlabel('H [T]'); %ylabel('\chi_{ab}')
        subplot(1,2,2); grid on; box on; hold on; set(gca,'fontsize',20)
        xlabel('H [T]'); %ylabel('\chi_{c}')
        
        H2 = power(HT,2);
        
        for q = 1:length(Temps)
            
            subplot(1,2,2);
            Y = HT.*Mc(q,:) - H2.*Xa_Tr(q,:);
            plot([0, HT], [0, Y]*muB*gJ ,'displayname',['T = ' num2str(Temps(q)) ' K'])
            
            subplot(1,2,1);
            Y = HT.*Ma(q,:) - H2.*Xc_Tr(q,:);
            plot([0, HT], [0, Y]*muB*gJ ,'displayname',['T = ' num2str(Temps(q)) ' K'])
        end
        
    case 'M_a(H), M_c(H)'
        
        HCEF = HCEFf(paramsB);
        JJ = JJf(paramsJJ);
        
        
        T = 4;
        N_T = 1;
        
        N_Th = 2;
        Th = [0,pi/2]; ind = [3,1];
        
        H = 2:2:60;
        N_H = length(H);
        
        n1 = 1; Z = ones(N_H,N_Th);
        for t = T
            
            n2 = 1;
            b = 1./(kB*t);
            for h = H
                
                n3 = 1;
                for th = Th
                    
                    hx = h*sin(th);
                    hz = h*cos(th);
                    hV = [hx,0,hz];
                    
                    mV = M(HCEF, JJ, hV, t);
                    %HMF0 = HMF(HCEF, JJ, hV, mV);
                    
                    mp(n2,n3) = mV(ind(n3));
                    %mb(n2,n3) = mV(2);
                    
                    n3=n3+1;
                end
                
                n2=n2+1;
            end
            n1=n1+1;
        end
        
        figure; hold on;
        plot([0 H],[0 mp(:,2)'],'r');
        plot([0 H],[0 mp(:,1)'],'b');
        
        %plot([0 H],[0 mb(:,2)'],'c--');
        %plot([0 H],[0 mb(:,1)'],'m--');
        legend({'M_{ab}','M_c'})
        
    case 'H-Dept. k 1 temp'
        
        
        
        
        
        startind = 251;
        endind = 2749-450;
        
        figure; hold on;
        
        T_RTM = [1.5,4,6,12,20, 30,50,70];
        
        
        
        H = 2:2:60; t = 4;
        
        HCEF = HCEFf(paramsB);
        JJ = JJf(paramsJJ);
        
        n2 = 1;
        for t = [4, 10, 30, 70]
            subplot(1,2,2); hold on;
            Y = [0 1e3*amp*k_vH(HCEF,JJ,t,H,0   ,.001)]/w0(n2);
            plot([0 H],Y,'displayname' ,[num2str(t) ' K model'])
            %display(sum(isnan(Y)));
            subplot(1,2,1); hold on;
            Y = [0 1e3*amp*k_vH(HCEF,JJ,t,H,pi/2,.001)]/w0(n2);
            plot([0 H],Y,'displayname' ,[num2str(t) ' K  model'])
            %display(sum(isnan(Y)));
            n2 = n2+1;
        end
        
        subplot(1,2,2); title('\theta = 0, H || c')
        legend('location','northwest'); grid on; box on;
        xlabel('\mu_0H [T]')
        ylabel('\Deltaf [Hz]')
        
        subplot(1,2,1); title('\theta = \pi/2, H || ab')
        grid on; box on;
        xlabel('\mu_0H [T]')
        ylabel('\Deltaf [Hz]')
        
    case 'H-Dept. k '
        
        
        
        dat{1}{1} = load('p048_031020_TDO004.dat');
        dat{1}{2} = load('p043_031020_TDO003.dat');
        dat{1}{3} = load('p050_031020_TDO002.dat');
        dat{1}{4} = load('p052_031020_TDO001.dat');
        dat{1}{5} = load('p058_031020_TDO001.dat');
        dat{1}{6} = load('p067_031020_TDO001.dat');
        dat{1}{7} = load('p079_031020_TDO001.dat');
        dat{1}{8} = load('p092_031020_TDO001.dat');
        
        % 0 deg
        %cd /Users/iCAPOS/Google Drive/CsYbSe2_RTM
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
        
        T_RTM = [1.5,4,6,12,20, 30,50,70];
        
        for n1 = [1,2]
            subplot(1,2,n1); hold on;
            for n2 = [2:8]
                x = dat{n1}{n2}(startind:end-endind,1);
                y = dat{n1}{n2}(startind:end-endind,2);
                
                ind = find(x == min(x));
                freq0{n1}(n2) = y(ind);
                
                plot( x, y-freq0{n1}(n2) ,'sq','color',[.5,.5,.5], 'displayname', [num2str(T_RTM(n2)) ' K Data']);
            end
        end
        
        H = 4:8:60; t = 4;
        
        HCEF = HCEFf(paramsB);
        JJ = JJf(paramsJJ);
        
        n2 = 1;
        for t = [1.5,4,6,12,20, 30,50,70]
            subplot(1,2,2); hold on;
            Y = [0 1e3*amp*k_vH(HCEF,JJ,t,H,0   ,.001)]/w0(n2);
            plot([0 H],Y,'displayname' ,[num2str(t) ' K model'])
            %display(sum(isnan(Y)));
            subplot(1,2,1); hold on;
            Y = [0 1e3*amp*k_vH(HCEF,JJ,t,H,pi/2,.001)]/w0(n2);
            plot([0 H],Y,'displayname' ,[num2str(t) ' K  model'])
            %display(sum(isnan(Y)));
            n2 = n2+1;
        end
        
        subplot(1,2,2); title('\theta = 0, H || c')
        legend('location','northwest'); grid on; box on;
        xlabel('\mu_0H [T]')
        ylabel('\Deltaf [Hz]')
        
        subplot(1,2,1); title('\theta = \pi/2, H || ab')
        grid on; box on;
        xlabel('\mu_0H [T]')
        ylabel('\Deltaf [Hz]')
        
    case 'H-Dept. k real units'
        
        
        
        dat{1}{1} = load('p048_031020_TDO004.dat');
        dat{1}{2} = load('p043_031020_TDO003.dat');
        dat{1}{3} = load('p050_031020_TDO002.dat');
        dat{1}{4} = load('p052_031020_TDO001.dat');
        dat{1}{5} = load('p058_031020_TDO001.dat');
        dat{1}{6} = load('p067_031020_TDO001.dat');
        dat{1}{7} = load('p079_031020_TDO001.dat');
        dat{1}{8} = load('p092_031020_TDO001.dat');
        
        % 0 deg
        %cd /Users/iCAPOS/Google Drive/CsYbSe2_RTM
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
        
        T_RTM = [1.5,4,6,12,20, 30,50,70];
        
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
                
                plot( x, (y-freq0{n1}(n2))/(1e3*amp/w0(n2)),'o','color',colz(n2-1,:),...
                    'displayname', [num2str(T_RTM(n2)) ' K Data'],...
                    'markerfacecolor',[1,1,1],'linewidth',1.5,...
                    'markersize',8);
            end
        end
        
        H = 4:8:60; t = 4;
        
        HCEF = HCEFf(paramsB);
        JJ = JJf(paramsJJ);
        
        n2 = 1;
        for t = [1.5,4,6,12,20, 30,50,70]
            subplot(1,2,2); hold on;
            Y = [0 k_vH(HCEF,JJ,t,H,0   ,.001)] ;
            plot([0 H],Y,'k','displayname' ,[num2str(t) ' K model'])
            %display(sum(isnan(Y)));
            subplot(1,2,1); hold on;
            Y = [0 k_vH(HCEF,JJ,t,H,pi/2,.001)];
            plot([0 H],Y,'k','displayname' ,[num2str(t) ' K  model'])
            %display(sum(isnan(Y)));
            n2 = n2+1;
        end
        
        subplot(1,2,2); set(gca,'fontsize',24); title('H || c')
        %        legend('location','northwest'); grid on;
        box on;
        xlabel('\mu_0H [T]','fontsize',26.4)
        ylabel('k [meV per Yb^{3+}]','fontsize',26.4)
        
        subplot(1,2,1); set(gca,'fontsize',24); title('H || ab','fontsize',26.4)
        box on;
        xlabel('\mu_0H [T]')
        ylabel('k [meV per Yb^{3+}]','fontsize',26.4)
        
    case 'H-Dept. k real units just calc'
        
        figure; hold on;
        
        H = 4:4:60;
        
        HCEF = HCEFf(paramsB);
        JJ = JJf(paramsJJ);
        
        n2 = 1;
        for t = [4,6,12,20, 30,50,70]
            subplot(1,2,2); hold on;
            Y = [0 k_vH(HCEF,JJ,t,H,0   ,.001)] ;
            plot([0 H],Y,'k','displayname' ,[num2str(t) ' K model'])
            %display(sum(isnan(Y)));
            subplot(1,2,1); hold on;
            Y = [0 k_vH(HCEF,JJ,t,H,pi/2,.001)];
            plot([0 H],Y,'k','displayname' ,[num2str(t) ' K  model'])
            %display(sum(isnan(Y)));
            n2 = n2+1;
        end
        
    case 'H-Dept. k real units w offset'
        
        
        
        dat{1}{1} = load('p048_031020_TDO004.dat');
        dat{1}{2} = load('p043_031020_TDO003.dat');
        dat{1}{3} = load('p050_031020_TDO002.dat');
        dat{1}{4} = load('p052_031020_TDO001.dat');
        dat{1}{5} = load('p058_031020_TDO001.dat');
        dat{1}{6} = load('p067_031020_TDO001.dat');
        dat{1}{7} = load('p079_031020_TDO001.dat');
        dat{1}{8} = load('p092_031020_TDO001.dat');
        
        % 0 deg
        %cd /Users/iCAPOS/Google Drive/CsYbSe2_RTM
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
        
        T_RTM = [1.5,4,6,12,20, 30,50,70];
        
        colz = [0,0,1;...
            0.3,0.75,0.93;...
            0.47,0.67,0.19;...
            0,1,0;...
            0.87,0.49,0;...
            0.75,0,0.75;...
            1,0,0;...
            ];
        
        o = [1,1];
        
        for n1 = [1,2]
            subplot(1,2,n1); hold on;
            for n2 = [2:8]
                x = dat{n1}{n2}(startind:end-endind,1);
                y = dat{n1}{n2}(startind:end-endind,2);
                
                ind = find(x == min(x));
                freq0{n1}(n2) = y(ind);
                
                plot( x, (y-freq0{n1}(n2))/(1e3*amp/w0(n2)) + (n2-2)*o(n1) ,'o','color',colz(n2-1,:),...
                    'displayname', [num2str(T_RTM(n2)) ' K Data'],...
                    'markerfacecolor',[1,1,1],'linewidth',1.5,...
                    'markersize',8);
            end
        end
        
        H = .25:.25:60; t = 4;
        
        HCEF = HCEFf(paramsB);
        JJ = JJf(paramsJJ);
        
        n2 = 1;
        for t = [4,6,12,20, 30,50,70]
            subplot(1,2,2); hold on;
            Y = [0 k_vH(HCEF,JJ,t,H,0   ,.001)] + (n2-1)*o(2) ;
            plot([0 H],Y,'k:','displayname' ,[num2str(t) ' K model'],...
                'linewidth',2.5)
            %display(sum(isnan(Y)));
            subplot(1,2,1); hold on;
            Y = [0 k_vH(HCEF,JJ,t,H,pi/2,.001)] + (n2-1)*o(1);
            plot([0 H],Y,'k:','displayname' ,[num2str(t) ' K  model'],...
                'linewidth',2.5)
            %display(sum(isnan(Y)));
            n2 = n2+1;
        end
        
        subplot(1,2,2); set(gca,'fontsize',24); title('H || c');
        %        legend('location','northwest'); grid on;
        box on;
        xlabel('\mu_0H [T]','fontsize',26.4)
        ylabel('k [meV per Yb^{3+}]','fontsize',26.4); ylim([-9,7.5])
        
        subplot(1,2,1); set(gca,'fontsize',24); title('H || ab','fontsize',26.4)
        box on;
        xlabel('\mu_0H [T]')
        ylabel('k [meV per Yb^{3+}]','fontsize',26.4); ylim([-2.5,6.5])
        
    case 'H-Dept. k alt phi'
        
        
        
        dat{1}{1} = load('p048_031020_TDO004.dat');
        dat{1}{2} = load('p043_031020_TDO003.dat');
        dat{1}{3} = load('p050_031020_TDO002.dat');
        dat{1}{4} = load('p052_031020_TDO001.dat');
        dat{1}{5} = load('p058_031020_TDO001.dat');
        dat{1}{6} = load('p067_031020_TDO001.dat');
        dat{1}{7} = load('p079_031020_TDO001.dat');
        dat{1}{8} = load('p092_031020_TDO001.dat');
        
        % 0 deg
        %cd /Users/iCAPOS/Google Drive/CsYbSe2_RTM
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
        
        T_RTM = [1.5,4,6,12,20, 30,50,70];
        
        for n1 = [1,2]
            subplot(1,2,n1); hold on;
            for n2 = [2:8]
                x = dat{n1}{n2}(startind:end-endind,1);
                y = dat{n1}{n2}(startind:end-endind,2);
                
                ind = find(x == min(x));
                freq0{n1}(n2) = y(ind);
                
                plot( x, y-freq0{n1}(n2) ,'sq','color',[.5,.5,.5], 'displayname', [num2str(T_RTM(n2)) ' K Data']);
            end
        end
        
        H = 4:8:60; t = 4;
        
        HCEF = HCEFf(paramsB);
        JJ = JJf(paramsJJ);
        
        n2 = 1; ph = pi/2;
        for t = [1.5,4,6,12,20, 30,50,70]
            subplot(1,2,2); hold on;
            Y = [0 1e3*amp*k_vH_aph(HCEF,JJ,t,H,0   , ph, .001)]/w0(n2);
            plot([0 H],Y,'displayname' ,[num2str(t) ' K model'])
            %display(sum(isnan(Y)));
            subplot(1,2,1); hold on;
            Y = [0 1e3*amp*k_vH_aph(HCEF,JJ,t,H,pi/2, ph, .001)]/w0(n2);
            plot([0 H],Y,'displayname' ,[num2str(t) ' K  model'])
            %display(sum(isnan(Y)));
            n2 = n2+1;
        end
        
        subplot(1,2,2); title('\theta = 0, H || c')
        legend('location','northwest'); grid on; box on;
        xlabel('\mu_0H [T]')
        ylabel('\Deltaf [Hz]')
        
        subplot(1,2,1); title('\theta = \pi/2, H || ab')
        grid on; box on;
        xlabel('\mu_0H [T]')
        ylabel('\Deltaf [Hz]')
        
    case 'H-Dept. k SHF'
        
        
        dat{1}{1} = load('p048_031020_TDO004.dat');
        dat{1}{2} = load('p043_031020_TDO003.dat');
        dat{1}{3} = load('p050_031020_TDO002.dat');
        dat{1}{4} = load('p052_031020_TDO001.dat');
        dat{1}{5} = load('p058_031020_TDO001.dat');
        dat{1}{6} = load('p067_031020_TDO001.dat');
        dat{1}{7} = load('p079_031020_TDO001.dat');
        dat{1}{8} = load('p092_031020_TDO001.dat');
        
        % 0 deg
        %cd /Users/iCAPOS/Google Drive/CsYbSe2_RTM
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
        
        T_RTM = [1.5,4,6,12,20, 30,50,70];
        
        for n1 = [1,2]
            subplot(1,2,n1); hold on;
            for n2 = [2:8]
                x = dat{n1}{n2}(startind:end-endind,1);
                y = dat{n1}{n2}(startind:end-endind,2);
                
                ind = find(x == min(x));
                freq0{n1}(n2) = y(ind);
                
                plot( x, y-freq0{n1}(n2) ,'sq','color',[.5,.5,.5], 'displayname', [num2str(T_RTM(n2)) ' K Data']);
            end
        end
        
        H = [10:10:50 100:50:200] ; t = 4;
        
        HCEF = HCEFf(paramsB);
        JJ = JJf(paramsJJ);
        
        n2 = 1;
        for t = [1.5,4,6,12,20, 30,50,70]
            subplot(1,2,2); hold on;
            Y = [0 1e3*amp*k_vH(HCEF,JJ,t,H,0   ,.001)]/w0(n2);
            plot([0 H],Y,'displayname' ,[num2str(t) ' K model'])
            %display(sum(isnan(Y)));
            subplot(1,2,1); hold on;
            Y = [0 1e3*amp*k_vH(HCEF,JJ,t,H,pi/2,.001)]/w0(n2);
            plot([0 H],Y,'displayname' ,[num2str(t) ' K  model'])
            %display(sum(isnan(Y)));
            n2 = n2+1;
        end
        
        subplot(1,2,2); title('\theta = 0, H || c')
        legend('location','northwest'); grid on; box on;
        xlabel('\mu_0H [T]')
        ylabel('\Deltaf [Hz]')
        
        subplot(1,2,1); title('\theta = \pi/2, H || ab')
        grid on; box on;
        xlabel('\mu_0H [T]')
        ylabel('\Deltaf [Hz]')
        
    case 'H-Dept. k alt'
        H = 4:4:60; t = 4;
        
        HCEF = HCEFf(paramsB);
        JJ = JJf(paramsJJ);
        
        figure; n2 = 1;
        for t = [4, 70]
            subplot(1,2,1); hold on;
            plot([0 H],[0 1e3*amp*k_vH(HCEF,JJ,t,H,0   ,.001)]/w0(n2),'displayname' ,[num2str(t) ' K'])
            subplot(1,2,2); hold on;
            plot([0 H],[0 1e3*amp*k_vH(HCEF,JJ,t,H,pi/2,.001)]/w0(n2),'displayname' ,[num2str(t) ' K'])
            n2 = n2+1;
        end
        legend show;
        
    case 'Angle Dept. k'
        
        HCEF = HCEFf(paramsB); JJ = JJf(paramsJJ);
        
        T1 = 4;
        N_T = 1;
        
        N_Th = 50; dth = pi/N_Th;
        Th = linspace(0 - 2*dth, pi + 2*dth, N_Th+4);
        
        H = [4:4:28];
        N_H = length(H);
        
        n1 = 1; Z = ones(N_H,N_Th);
        for t = T1
            n2 = 1;
            b = 1./(kB*t);
            for h = H
                k(n2,:) = k_vTh_equispc(HCEF,JJ,t,h,Th,dth);
                n2=n2+1;
            end
            n1=n1+1;
        end
        
        % Plotting
        colz = varycolor(N_H);
        figure; hold on;
        
        th_p  = Th(3:end-2);
        
        for n2 = 1:N_H
            plot(th_p, amp*(1e3)*k(n2,:)/w0(2) , 'color', colz(n2,:), 'displayname', num2str(H(n2)) );
        end
        
        title({['B_n^m = [' num2str(paramsB) '],']
            ['J_{xx} = ' num2str(paramsJJ(1)) ',  J_{zz} = ' num2str(paramsJJ(2)) ]});
        
    case 'Angle Dept. k fringing field'
        
        HCEF = HCEFf(paramsB); JJ = JJf(paramsJJ);
        
        T1 = 4;
        N_T = 1;
        
        N_Th = 30; dth = pi/N_Th;
        Th = linspace(0 - 2*dth, pi + 2*dth, N_Th+4);
        
        H = 28;
        N_H = length(H);
        
        offset_case = 3;
        
        switch offset_case
            
            case 1
                
                A0 =  -0.0004243; % (-0.0004244, -0.0004242)
                A2 =    0.001273; % (0.001273, 0.001273)
                B1 =  -0.0009163; %  (-0.0009165, -0.0009161)
                B3 =    0.000531; % (0.0005308, 0.0005312)
                
                A1 =  -0.0007039; %  (-0.0007041, -0.0007037)
                A3 =   0.0005309; %  (0.0005307, 0.0005311)
                B2 =   -0.001273; %  (-0.001273, -0.001272)
                
                
                fa = @(x) H*(A0 + A2*cos(2*x) + (B1+1)*sin(x) + B3*sin(3*x));
                fc = @(x) H*(B2*sin(2*x) + (A1+1)*cos(x) + A3*cos(3*x));
                
            case 2
                
                A0 =   -0.002153; %  (-0.002153, -0.002152)
                A2 =    0.006458; %  (0.006457, 0.006458)
                B1 =    -0.02144; %  (-0.02144, -0.02144)
                B3 =   0.0005531; %  (0.0005524, 0.0005538)
                
                A1 =    -0.02122;%  (-0.02122, -0.02122)
                A3 =   0.0005531;%  (0.0005524, 0.0005537)
                B2 =   -0.006457;%  (-0.006458, -0.006456)
                
                fa = @(x) H*(A0 + A2*cos(2*x) + (B1+1)*sin(x) + B3*sin(3*x));
                fc = @(x) H*(B2*sin(2*x) + (A1+1)*cos(x) + A3*cos(3*x));
                
            case 3
                
                B1 =     0.01046; %  (0.01046, 0.01046)
                B2 =    0.004183; %  (0.004182, 0.004184)
                B3 =   0.0005108; %   (0.00051, 0.0005117)
                
                A1 =     0.01067; %  (0.01067, 0.01067)
                A2 =    0.004182; %  (0.004181, 0.004183)
                A3 =   0.0005107; %  (0.0005098, 0.0005116)
                
                fa = @(x) H*((B1+1)*sin(x) + B2*sin(2*x) + B3*sin(3*x));
                fc = @(x) H*((A1+1)*cos(x) + A2*cos(2*x) + A3*cos(3*x));
                
        end
        
        n1 = 1; Z = ones(N_H,N_Th);
        for t = T1
            n2 = 1;
            b = 1./(kB*t);
            for th = Th
                
                hx0 = H*sin(th);
                hz0 = H*cos(th);
                
                hx = fa(th);
                hz = fc(th);
                
                hv0 = [hx0,0,hz0];
                FF0(n2) = F_1pt_hv(HCEF,JJ,t,hv0);
                
                hv1 = [hx,0,hz];
                FF1(n2) = F_1pt_hv(HCEF,JJ,t,hv1);
                
                
                n2=n2+1;
            end
            
            k0 = (C2(1)*FF0(1:end-4) + C2(2)*FF0(2:end-3) + C2(3)*FF0(3:end-2)...
                + C2(4)*FF0(4:end-1) + C2(5)*FF0(5:end))/dth/dth;
            
            k1 = (C2(1)*FF1(1:end-4) + C2(2)*FF1(2:end-3) + C2(3)*FF1(3:end-2)...
                + C2(4)*FF1(4:end-1) + C2(5)*FF1(5:end))/dth/dth;
            
            n1=n1+1;
        end
        
        % Plotting
        colz = varycolor(N_H);
        figure; hold on;
        
        th_p  = Th(3:end-2);
        
        plot(th_p, k0, 'k--');
        plot(th_p, k1, 'b');
        
        figure; hold on; plot(th_p, k1-k0)
        
        %         for n2 = 1:N_H
        %             plot(th_p, amp*(1e3)*k(n2,:)/w0(2) , 'color', colz(n2,:), 'displayname', num2str(H(n2)) );
        %         end
        
        title({['B_n^m = [' num2str(paramsB) '],']
            ['J_{xx} = ' num2str(paramsJJ(1)) ',  J_{zz} = ' num2str(paramsJJ(2)) ]});
        
    case 'Th-Dept. k corrected 4K'
        
        Hset = 4:4:28;
        
        colz = varycolor(7);
        
        clear Th1; clear Fv1;
        
        % 4 K Data
        figure; set(gca,'fontsize',16);
        hold on; title('4 K'); xlabel('\theta'); ylabel('\Delta f [Hz]');
        grid on; box on;
        Dat = load('C:\Users\Chris\Documents\Lab Work\CsYbSe2\CsYbSe2_FvTh.mat');
        for n = 1:7
            Th1{n} = Dat.CsYbSe2_ThDeptCell.Th{n};
            Fv1{n}  = Dat.CsYbSe2_ThDeptCell.F{n};
            
            plot(Th1{n},Fv1{n},'sq','color',colz(n,:))
            
        end
        
        HCEF = HCEFf(paramsB); JJ = JJf(paramsJJ);
        
        T1 = 4;
        N_T = 1;
        
        th_i = -pi/4;
        th_f = pi+pi/4;
        
        N_Th = 50; dth = (th_f - th_i)/(N_Th);
        Th = linspace( th_i - 2*dth  , th_f + 2*dth, N_Th+4+1);
        
        
        H = Hset;
        N_H = length(H);
        
        n1 = 1; Z = ones(N_H,N_Th);
        for t = T1
            n2 = 1;
            b = 1./(kB*t);
            for h = H
                F(n2,:) = F_vTh_equispc(HCEF,JJ,t,h,Th);
                n2=n2+1;
            end
            n1=n1+1;
        end
        
        % Plotting
        colz = varycolor(N_H);
        
        th_p  = Th(3:end-2);
        K0 = 123.2481/(2*pi);
        
        
        DatZFw0 = [-49.0214	49526.47947
            -21.0736	49525.782526
            -12.2291	49524.032328
            -9.3972	49524.28404
            0.1	49521.00881
            10.5889	49517.251865
            20.6854	49524.882365
            30.6758	49527.947547
            40.8927	49523.060674
            50.8253	49525.337303
            60.7131	49524.740536
            70.6209	49526.060354
            80.4684	49526.518149
            90.0276	49527.012619
            99.4588	49527.092806
            108.9467	49526.214482
            118.5334	49525.950081
            127.8561	49525.038617
            137.0115	49524.842469
            145.2634	49524.219632
            153.0882	49524.296639
            159.2969	49523.795006];
        
        x = DatZFw0(:,1)*pi/180;
        y = DatZFw0(:,2);
        
        x1 = [(x-pi) x (x+pi)];
        y1 = [y,y,y];
        
        [x1, ind] = sort(x1);
        y1 = y1(ind);
        
        w = Interp1NonUnique(x1,y1, th_p);
        
        for n2 = 1:N_H
            
            tau   = FiniteDiff5Point(C1,F(n2,:))/dth;
            kv    = FiniteDiff5Point(C2,F(n2,:))/dth/dth;
            dkdth = FiniteDiff5Point(C3,F(n2,:))/dth/dth/dth;
            
            %dkdth = (C1(1)*kv(1:end-4) + C1(2)*kv(2:end-3) + C1(3)*kv(3:end-2) + C1(4)*kv(4:end-1) + C1(5)*kv(5:end))/dth ;
            
            
            
            %Df0 = (3/K0)*(kv/2)*1e3*w0(2);
            %Df1 = (3/K0)*( -(7/6)*3*tau.*dkdth./(K0 + 3.*kv) )*1e3*w0(2);
            
            %Df0 = (1./K0)*(.5*kv)*1e3*w0(2);
            %Df1 = (1./K0)*( -3*tau.*dkdth./(K0 + 3.*kv) )*1e3*w0(2);
            
            Df0 = (1./K0)*(kv/2)*1e3.*w0(2);
            Df1 = (1/(2*K0)).*(kv - dkdth.*tau./(K0 + kv)).*1e3*w0(2);
            
            plot(th_p, Df0/(2*pi) , '--', 'color', colz(n2,:), 'displayname', num2str(H(n2)) );
            plot(th_p, Df1/(2*pi) , '-', 'color', colz(n2,:), 'displayname', num2str(H(n2)) );
            
            X1 = linspace(0,pi,1000);
            Y1 = Interp1NonUnique(th_p,Df1,X1);
            
            weight(n2) = sum((Y1(1:end-1) +  Y1(2:end)).*diff(X1)/2);
            
            %plot(th_p, Df1 , '-', 'color', colz(n2,:), 'displayname', num2str(H(n2)) );
        end
        
        xlim([-pi/4,pi+pi/4])
        xticks([-pi/2,0,pi/2,pi,3*pi/2])
        xticklabels({'-\pi/2','0','\pi/2','\pi','3\pi/2'})
        
        
        figure; hold on;
        plot(Hset,weight)
        
    case 'Calculate Weight'
        T = 4;
        N_T = 1;
        
        N_Th = 50; dth = pi/N_Th;
        Th = linspace(0 - 2*dth, pi + 2*dth, N_Th+4);
        
        HCEF = HCEFf(paramsB);
        JJ = JJf(paramsJJ);
        
        H = 8:4:32;
        N_H = length(H);
        
        n1 = 1; Z = ones(N_H,N_Th);
        for t = T
            
            n2 = 1;
            b = 1./(kB*t);
            for h = H
                
                n3 = 1;
                for th = Th
                    
                    hx = h*sin(th);
                    hz = h*cos(th);
                    hV = [hx,0,hz];
                    
                    mV = M(HCEF, JJ, hV, t);
                    HMF0 = HMF(HCEF, JJ, hV, mV);
                    
                    Z(n2,n3) = trace(expm(-b*HMF0));
                    
                    n3=n3+1;
                end
                
                n2=n2+1;
            end
            n1=n1+1;
        end
        
        F = -kB*T*log(Z);
        
        % Plotting
        colz = varycolor(N_H);
        figure; hold on;
        
        % 2nd Derivative -> k
        C2 = [-1/12 	4/3 	-5/2 	4/3 	-1/12];
        
        for n2 = 1:N_H
            
            k = (C2(1)*F(n2,1:end-4) + C2(2)*F(n2,2:end-3) + C2(3)*F(n2,3:end-2)...
                + C2(4)*F(n2,4:end-1) + C2(5)*F(n2,5:end))/dth/dth;
            
            %plot(Th(3:end-2), k , 'color', colz(n2,:), 'displayname', num2str(H(n2)) );
            
            w(n2) = sum(k)*dth;
        end
        
        %         title({['B_n^m = [' num2str(paramsB) '],']
        %                ['J_{xx} = ' num2str(JJxx) ',  J_{zz} = ' num2str(JJzz) ]});
        
        %figure; hold on;
        plot(H,w);
        
    case 'replot k(th)'
        
        Hset = 4:4:28;
        
        colz = varycolor(7);
        
        clear Th1; clear Fv1;
        
        % 4 K Data
        figure; set(gca,'fontsize',16);
        hold on; title('4 K'); xlabel('\theta'); ylabel('\Delta f [Hz]');
        grid on; box on;
        Dat = load('C:\Users\Chris\Documents\Lab Work\CsYbSe2\CsYbSe2_FvTh.mat');
        for n = 1:7
            Th1{n} = Dat.CsYbSe2_ThDeptCell.Th{n};
            Fv1{n}  = Dat.CsYbSe2_ThDeptCell.F{n};
            %             F_int = Interp1NonUnique(Th1{n},Fv1{n},Th_int1);
            %
            %             W(n) = sum((F_int(1:end-1) + F_int(2:end)).*diff(Th_int1))/2;
            %            Fv1{n} = Fv1{n}-W(n)/pi;
            X1 = [(Th1{n}(2:end)-pi) Th1{n}(2:end) (Th1{n}(2:end)+pi)];
            Y1 = [Fv1{n}(2:end) Fv1{n}(2:end) Fv1{n}(2:end)]; [X1 ind] = sort(X1);
            plot(X1,Y1,'sq-','color',colz(n,:))
            
            Th_int1 = linspace(-pi/2,pi/2,22);
            F_int = Interp1NonUnique(X1, Y1, Th_int1);
            
            W(n) = sum((F_int(1:end-1) + F_int(2:end)).*diff(Th_int1))/2;
        end
        
        xlim([-pi/4,pi+pi/4])
        xticks([0,pi/2,pi])
        xticklabels({'-\pi/2','0','\pi/2','\pi','3\pi/2'})
        
        
        if 1
            figure; hold on;
            
            plot(Hset, W ,'sq-')
        end
        
    case 'Fmin Th-dept Data 4K'
        
        %% Data Sets
        Hset = 4:4:28;
        
        colz = varycolor(7);
        
        clear Th1; clear Fv1;
        
        % 4 K Data
        figure;
        hold on; title('4 K'); xlabel('\theta'); ylabel('\Delta f [Hz]');
        grid on; box on;
        Dat = load('/Users/iCAPOS/Desktop/ML/CsYbSe2_FvTh.mat');
        for n = 1:7
            Th1{n} = Dat.CsYbSe2_ThDeptCell.Th{n};
            Fv1{n}  = Dat.CsYbSe2_ThDeptCell.F{n};
            %             F_int = Interp1NonUnique(Th1{n},Fv1{n},Th_int1);
            %
            %             W(n) = sum((F_int(1:end-1) + F_int(2:end)).*diff(Th_int1))/2;
            %            Fv1{n} = Fv1{n}-W(n)/pi;
            plot(Th1{n},Fv1{n},'sq','color',colz(n,:))
        end
        Thset = Th1{1};
        
        x0 = [amp,paramsB(1),paramsB(2),paramsB(3),paramsB(4),paramsB(5),paramsB(6),paramsJJ(1),paramsJJ(2)];
        f = @(x) CummErr_ThDept0(x(1),HCEFf([x(2),x(3),x(4),x(5),x(6),x(7)]),JJf([x(8),x(9)]));
        opts = optimset('MaxIter',1);
        %         %p = fminsearchbnd(f,x0,[0,0,0,-Inf,-Inf],[Inf,2,2,0,0],opts);
        p = fminsearchbnd(f,x0,[250, -1, -1, -1, -1, -1, -1, 0, 0],[250, 1, 1, 1, 1, 1, 1, Inf, Inf],opts);
        display(p)
        
        HCEF = HCEFf(p(2:7));
        JJ = JJf(p(8:9));
        
        Th_start = -pi/4;
        Th_end   =  pi;
        N_Th = 70; dth = (Th_end - Th_start)/N_Th;
        Th = linspace(Th_start - 2*dth, Th_end + 2*dth, N_Th+4);
        
        t = 4;
        th_int = linspace(-pi,pi,200);%)(pi/180)*linspace(-51,180,1000);
        for n = 1:7
            y = p(1)*k_vTh(HCEF,JJ,t,Hset(n),Th,dth);
            plot(Th(3:end-2),y,'-','color',colz(n,:));
            
        end
        grid on; box on; xlabel('\Theta'); ylabel('\Delta f [Hz]')
        
    case 'Nonlinlsq Th-dept Data 4K'
        
        Th_int1 = (pi/180)*[-45:10:135];
        
        %% Data Sets
        Hset = 4:4:28;
        
        colz = varycolor(7);
        
        clear Th1; clear Fv1;
        
        % 4 K Data
        figure;
        hold on; title('4 K'); xlabel('\theta'); ylabel('\Delta f [Hz]');
        grid on; box on;
        Dat = load('/Users/iCAPOS/Desktop/ML/CsYbSe2_FvTh.mat');
        for n = 1:7
            Th1{n} = Dat.CsYbSe2_ThDeptCell.Th{n};
            Fv1{n}  = Dat.CsYbSe2_ThDeptCell.F{n};
            F_int = Interp1NonUnique(Th1{n},Fv1{n},Th_int1);
            %
            W(n) = sum((F_int(1:end-1) + F_int(2:end)).*diff(Th_int1))/2;
            Fv1{n} = Fv1{n}-W(n)/pi;
            plot(Th1{n},Fv1{n},'sq','color',colz(n,:))
        end
        
        Xdata = [Th1{1} Th1{2} Th1{3} Th1{4} Th1{5} Th1{6} Th1{7}];
        Ydata = [Fv1{1} Fv1{2} Fv1{3} Fv1{4} Fv1{5} Fv1{6} Fv1{7}];
        
        Thset = Th1{1};
        
        
        x0 = [amp,paramsB(1),paramsB(2),paramsB(3),paramsB(4),paramsB(5),paramsB(6),paramsJJ(1),paramsJJ(2)];
        
        f = @(x,X) Func1(x,X);
        p = lsqcurvefit(f,x0,Xdata,Ydata,[14.5587, -1, -1, -1, -1, -1, -1, 0, 0],[14.5587, 1, 1, 1, 1, 1, 1, Inf, Inf]);
        
        display(p)
        
        output = p;
        
        HCEF = HCEFf(p(2:7));
        JJ = JJf(p(8:9));
        
        Th_start = -pi/4;
        Th_end   =  pi;
        N_Th = 70; dth = (Th_end - Th_start)/N_Th;
        Th = linspace(Th_start - 2*dth, Th_end + 2*dth, N_Th+4);
        
        t = 4;
        th_int = linspace(-pi,pi,200);%)(pi/180)*linspace(-51,180,1000);
        for n = 1:7
            y = p(1)*k_vTh_equispc(HCEF,JJ,t,Hset(n),Th,dth);
            plot(Th(3:end-2),y,'-','color',colz(n,:));
            
        end
        grid on; box on; xlabel('\Theta'); ylabel('\Delta f [Hz]')
        
    case 'Replot (H,T)-data'
        
        Hset1 = [2:3:29];
        Hset2 = [4:3:58];
        
        Hd_int{1} = {Hset2, Hset2, Hset1, Hset1, Hset1, Hset2, Hset1, Hset2};
        Hd_int{2} = {Hset2, Hset2, Hset1, Hset2, Hset1, Hset2, Hset1, Hset2};
        
        T= [1.6, 4, 6, 12, 20, 30, 50, 70];
        
        % 90 deg
        dat{1}{1} = load('p048_031020_TDO004.dat');
        dat{1}{2} = load('p043_031020_TDO003.dat');
        dat{1}{3} = load('p050_031020_TDO002.dat');
        dat{1}{4} = load('p052_031020_TDO001.dat');
        dat{1}{5} = load('p058_031020_TDO001.dat');
        dat{1}{6} = load('p067_031020_TDO001.dat');
        dat{1}{7} = load('p079_031020_TDO001.dat');
        dat{1}{8} = load('p092_031020_TDO001.dat');
        
        % 0 deg
        %cd /Users/iCAPOS/Google Drive/CsYbSe2_RTM
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
        
        for n1 = [1:2]
            subplot(1,2,n1); hold on;
            for n2 = [1:8]
                x = dat{n1}{n2}(startind:end-endind,1);
                y = dat{n1}{n2}(startind:end-endind,2);
                
                Hd{n1}{n2} = x;
                ind = find(x == min(x));
                f0{n1}(n2) = y(ind);
                kd{n1}{n2} = y-f0{n1}(n2);
                plot( Hd{n1}{n2}, kd{n1}{n2})
                
                kd_int{n1}{n2} = Interp1NonUnique(x,kd{n1}{n2}, Hd_int{n1}{n2} );
                plot(Hd_int{n1}{n2},kd_int{n1}{n2} ,'ksq')
                display(sum(isnan(kd_int{n1}{n2})));
            end
        end
        
        subplot(1,2,1); title('H||ab');
        box on; xlabel('\mu_0H [T]'); ylabel('\Delta f [Hz]')
        
        subplot(1,2,2); title('H||c');
        box on; xlabel('\mu_0H [T]'); ylabel('\Delta f [Hz]')
        
    case 'Nonlinlsq (H,T) Data'
        
        Hset1 = [2:3:29];
        Hset2 = [4:3:58];
        
        Hd_int{1} = {Hset2, Hset2, Hset1, Hset1, Hset1, Hset2, Hset1, Hset2};
        Hd_int{2} = {Hset2, Hset2, Hset1, Hset2, Hset1, Hset2, Hset1, Hset2};
        
        T= [1.6, 4, 6, 12, 20, 30, 50, 70];
        
        % 90 deg
        dat{1}{1} = load('p048_031020_TDO004.dat');
        dat{1}{2} = load('p043_031020_TDO003.dat');
        dat{1}{3} = load('p050_031020_TDO002.dat');
        dat{1}{4} = load('p052_031020_TDO001.dat');
        dat{1}{5} = load('p058_031020_TDO001.dat');
        dat{1}{6} = load('p067_031020_TDO001.dat');
        dat{1}{7} = load('p079_031020_TDO001.dat');
        dat{1}{8} = load('p092_031020_TDO001.dat');
        
        % 0 deg
        %cd /Users/iCAPOS/Google Drive/CsYbSe2_RTM
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
        
        for n1 = [1:2]
            subplot(1,2,n1); hold on;
            for n2 = [1:8]
                x = dat{n1}{n2}(startind:end-endind,1);
                y = dat{n1}{n2}(startind:end-endind,2);
                
                Hd{n1}{n2} = x;
                ind = find(x == min(x));
                f0{n1}(n2) = y(ind);
                kd{n1}{n2} = y-f0{n1}(n2);
                plot( Hd{n1}{n2}, kd{n1}{n2} ,'sq')
                
                kd_int{n1}{n2} = Interp1NonUnique(x,kd{n1}{n2}, Hd_int{n1}{n2} );
                %plot(Hd_int{n1}{n2},kd_int{n1}{n2} ,'ksq')
                %display(sum(isnan(kd_int{n1}{n2})));
            end
        end
        
        if 1
            Xdata = [Hd_int{1}{2} Hd_int{1}{3} Hd_int{1}{4} Hd_int{1}{5} Hd_int{1}{6} Hd_int{1}{7} Hd_int{1}{8}...
                Hd_int{2}{2} Hd_int{2}{3} Hd_int{2}{4} Hd_int{2}{5} Hd_int{2}{6} Hd_int{2}{7} Hd_int{2}{8}];
            %];
            Ydata = [kd_int{1}{2} kd_int{1}{3} kd_int{1}{4} kd_int{1}{5} kd_int{1}{6} kd_int{1}{7} kd_int{1}{8}...
                kd_int{2}{2} kd_int{2}{3} kd_int{2}{4} kd_int{2}{5} kd_int{2}{6} kd_int{2}{7} kd_int{2}{8}];
            % ];
            f = @(x,X) Func2(x,X);
            
        elseif 1
            
            Xdata = [Hd_int{1}{2} Hd_int{1}{3} Hd_int{1}{4} Hd_int{1}{5} Hd_int{1}{6} Hd_int{1}{7} Hd_int{1}{8}];%;...
            %Hd_int{2}{2} Hd_int{2}{3} Hd_int{2}{4} Hd_int{2}{5} Hd_int{2}{6} Hd_int{2}{7} Hd_int{2}{8}];
            %];
            Ydata = [kd_int{1}{2} kd_int{1}{3} kd_int{1}{4} kd_int{1}{5} kd_int{1}{6} kd_int{1}{7} kd_int{1}{8}];%...
            %;kd_int{2}{2} kd_int{2}{3} kd_int{2}{4} kd_int{2}{5} kd_int{2}{6} kd_int{2}{7} kd_int{2}{8}];
            % ];
            f = @(x,X) Func3(x,X);
            
        else
            
            Xdata = [...%Hd_int{1}{2} Hd_int{1}{3} Hd_int{1}{4} Hd_int{1}{5} Hd_int{1}{6} Hd_int{1}{7} Hd_int{1}{8}];%;...
                Hd_int{2}{2} Hd_int{2}{3} Hd_int{2}{4} Hd_int{2}{5} Hd_int{2}{6} Hd_int{2}{7} Hd_int{2}{8}];
            %];
            Ydata = [...%kd_int{1}{2} kd_int{1}{3} kd_int{1}{4} kd_int{1}{5} kd_int{1}{6} kd_int{1}{7} kd_int{1}{8}];%...
                kd_int{2}{2} kd_int{2}{3} kd_int{2}{4} kd_int{2}{5} kd_int{2}{6} kd_int{2}{7} kd_int{2}{8}];
            % ];
            f = @(x,X) Func4(x,X);
            
            
        end
        
        
        subplot(1,2,1); title('H||ab');
        box on; xlabel('\mu_0H [T]'); ylabel('\Delta f [Hz]')
        
        subplot(1,2,2); title('H||c');
        box on; xlabel('\mu_0H [T]'); ylabel('\Delta f [Hz]')
        
        
        
        x0 = [amp,paramsB(1),paramsB(2),paramsB(3),paramsB(4),paramsB(5),paramsB(6),paramsJJ(1),paramsJJ(2)];
        
        
        %plot(Xdata,f(x0,Xdata),'k-')
        
        p = lsqcurvefit(f,x0,Xdata,Ydata,[9.95, -1, -1, -1, -1, -1, -1, 0, 0],[148.5, 1, 1, 1, 1, 1, 1, Inf, Inf]);
        %
        display(p)
        %
        output = p;
        %
        HCEF = HCEFf(p(2:7));
        JJ = JJf(p(8:9));
        
        Hset3 = 1:60;
        Thset1 = [pi/2,0];
        
        for n1 = [1,2]
            subplot(1,2,n1)
            for n2 = [2:8]
                plot(Hset3, 1e3*p(1)*k_vH(HCEF,JJ,T(n2),Hset3,Thset1(n1),.001),'k--')
            end
            
        end
        
    case 'Nonlinlsq (H,T) Data vary single param'
        
        Hset1 = [2:3:29];
        Hset2 = [4:3:58];
        
        Hd_int{1} = {Hset2, Hset2, Hset1, Hset1, Hset1, Hset2, Hset1, Hset2};
        Hd_int{2} = {Hset2, Hset2, Hset1, Hset2, Hset1, Hset2, Hset1, Hset2};
        
        T= [1.6, 4, 6, 12, 20, 30, 50, 70];
        
        % 90 deg
        dat{1}{1} = load('p048_031020_TDO004.dat');
        dat{1}{2} = load('p043_031020_TDO003.dat');
        dat{1}{3} = load('p050_031020_TDO002.dat');
        dat{1}{4} = load('p052_031020_TDO001.dat');
        dat{1}{5} = load('p058_031020_TDO001.dat');
        dat{1}{6} = load('p067_031020_TDO001.dat');
        dat{1}{7} = load('p079_031020_TDO001.dat');
        dat{1}{8} = load('p092_031020_TDO001.dat');
        
        % 0 deg
        %cd /Users/iCAPOS/Google Drive/CsYbSe2_RTM
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
        
        %figure; hold on;
        
        for n1 = [1:2]
            subplot(1,2,n1); hold on;
            for n2 = [1:8]
                x = dat{n1}{n2}(startind:end-endind,1);
                y = dat{n1}{n2}(startind:end-endind,2);
                
                Hd{n1}{n2} = x;
                ind = find(x == min(x));
                f0{n1}(n2) = y(ind);
                kd{n1}{n2} = y-f0{n1}(n2);
                %plot( Hd{n1}{n2}, kd{n1}{n2} ,'sq')
                
                kd_int{n1}{n2} = Interp1NonUnique(x,kd{n1}{n2}, Hd_int{n1}{n2} );
                %plot(Hd_int{n1}{n2},kd_int{n1}{n2} ,'ksq')
                %display(sum(isnan(kd_int{n1}{n2})));
            end
        end
        
        if 1
            Xdata = [Hd_int{1}{2} Hd_int{1}{3} Hd_int{1}{4} Hd_int{1}{5} Hd_int{1}{6} Hd_int{1}{7} Hd_int{1}{8}...
                Hd_int{2}{2} Hd_int{2}{3} Hd_int{2}{4} Hd_int{2}{5} Hd_int{2}{6} Hd_int{2}{7} Hd_int{2}{8}];
            %];
            Ydata = [kd_int{1}{2} kd_int{1}{3} kd_int{1}{4} kd_int{1}{5} kd_int{1}{6} kd_int{1}{7} kd_int{1}{8}...
                kd_int{2}{2} kd_int{2}{3} kd_int{2}{4} kd_int{2}{5} kd_int{2}{6} kd_int{2}{7} kd_int{2}{8}];
            % ];
            %f = @(x,X) Func2(x,X);
            
        else
            Xdata = [Hd_int{1}{4} Hd_int{1}{5} Hd_int{1}{6} Hd_int{1}{7} Hd_int{1}{8}...
                Hd_int{2}{4} Hd_int{2}{5} Hd_int{2}{6} Hd_int{2}{7} Hd_int{2}{8}];
            %];
            Ydata = [kd_int{1}{4} kd_int{1}{5} kd_int{1}{6} kd_int{1}{7} kd_int{1}{8}...
                kd_int{2}{4} kd_int{2}{5} kd_int{2}{6} kd_int{2}{7} kd_int{2}{8}];
            % ];
            %f = @(x,X) Func3(x,X);
            
            
        end
        
        %subplot(1,2,1); title('H||ab');
        %box on; xlabel('\mu_0H [T]'); ylabel('\Delta f [Hz]')
        
        %subplot(1,2,2); title('H||c');
        %box on; xlabel('\mu_0H [T]'); ylabel('\Delta f [Hz]')
        
        %plot(Xdata,f(x0,Xdata),'k-')
        
        %         xLB = x0;
        %         xUB = x0;
        %
        %         p_var_ind = opt2;
        %
        %         for n =  p_var_ind
        %             if n <= 7
        %                xLB(n) = -10;
        %                xUB(n) = 10;
        %             else
        %                 xLB(n) = 0;
        %                 xUB(n) = 100;
        %             end
        %         end
        %
        %         if p_var_ind  <= 6
        %             xLB(p_var_ind) = -10;
        %             xUB(p_var_ind) = 10;
        %         else
        %
        %         end
        %
        %          p = lsqcurvefit(f,x0,Xdata,Ydata,xLB,xUB);
        
        x0 = [amp,paramsB(1),paramsB(2),paramsB(3),paramsB(4),paramsB(5),paramsB(6),paramsJJ(1),paramsJJ(2)];
        p_var_ind = opt2;
        
        x0(p_var_ind) = x0(p_var_ind).*(1 + .05*(rand(1,length(p_var_ind)) - .5));
        
        pstr = '['; n = 1;
        for i = 1:9
            if ismember(i,p_var_ind)
                if i == 1
                    pstr = [pstr 'x(' num2str(n) ')' ];
                else
                    pstr = [pstr ', x(' num2str(n) ')' ];
                end
                n = n+1;
            else
                if i == 1
                    pstr = [pstr num2str(x0(i))];
                else
                    pstr = [pstr ', ' num2str(x0(i))];
                end
            end
        end
        pstr = [pstr ']'];
        
        n = 1;
        for i =  p_var_ind
            if  i <= 7
                xLB(n) =  min([.8 1.2]*x0(i));
                xUB(n) =  max([.8 1.2]*x0(i));
            else
                xLB(n) = 0;
                xUB(n) = 100;
            end
            n = n+1;
        end
        
        eval(['f = @(x,X) Func2(' pstr ',X);']);
        %x = fminsearchbnd(f,x0(p_var_ind),xLB,xUB);
        x = lsqcurvefit(f,x0(p_var_ind),Xdata,Ydata,xLB,xUB);
        eval(['p = ' pstr]);
        
        display(p)
        %
        
        output.params = p;
        output.res = sum(power(f(x,Xdata) - Ydata,2));
        %
        HCEF = HCEFf(p(2:7));
        JJ = JJf(p(8:9));
        
        Hset3 = 1:60;
        Thset1 = [pi/2,0];
        
        clf('reset')
        for n1 = [1,2]
            subplot(1,2,n1); hold on;
            for n2 = [2:8]
                plot([0, Hset3], [0, 1e3*amp*k_vH(HCEF,JJ,T(n2),Hset3,Thset1(n1),.001)/w0(n2)],'k--')
            end
            drawnow
        end
        
    case 'Nonlinlsq (H,T,Th) Data'
        
        Th_int1 = (pi/180)*[-45:10:135];
        
        %% Data Sets
        Hset = 4:4:28;
        
        colz = varycolor(7);
        
        clear Th1; clear Fv1;
        
        
        Hset1 = [4:2:28];
        Hset2 = [4:2:58];
        
        Hd_int{1} = {Hset2, Hset2, Hset1, Hset1, Hset1, Hset2, Hset1, Hset2};
        Hd_int{2} = {Hset2, Hset2, Hset1, Hset2, Hset1, Hset2, Hset1, Hset2};
        
        T= [1.6, 4, 6, 12, 20, 30, 50, 70];
        
        % 90 deg
        dat{1}{1} = load('p048_031020_TDO004.dat');
        dat{1}{2} = load('p043_031020_TDO003.dat');
        dat{1}{3} = load('p050_031020_TDO002.dat');
        dat{1}{4} = load('p052_031020_TDO001.dat');
        dat{1}{5} = load('p058_031020_TDO001.dat');
        dat{1}{6} = load('p067_031020_TDO001.dat');
        dat{1}{7} = load('p079_031020_TDO001.dat');
        dat{1}{8} = load('p092_031020_TDO001.dat');
        
        % 0 deg
        %cd /Users/iCAPOS/Google Drive/CsYbSe2_RTM
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
        
        for n1 = [1:2]
            subplot(2,2,n1); hold on;
            for n2 = [1:8]
                x = dat{n1}{n2}(startind:end-endind,1);
                y = dat{n1}{n2}(startind:end-endind,2);
                
                Hd{n1}{n2} = x;
                ind = find(x == min(x));
                f0{n1}(n2) = y(ind);
                kd{n1}{n2} = y-f0{n1}(n2);
                plot( Hd{n1}{n2}, kd{n1}{n2} ,'sq')
                
                kd_int{n1}{n2} = Interp1NonUnique(x,kd{n1}{n2}, Hd_int{n1}{n2} );
                %plot(Hd_int{n1}{n2},kd_int{n1}{n2} ,'ksq')
                display(sum(isnan(kd_int{n1}{n2})));
            end
        end
        
        subplot(2,2,3); hold on;
        title('4 K'); xlabel('\theta'); ylabel('\Delta f [Hz]');
        grid on; box on;
        Dat = load('C:\Users\Chris\Documents\Lab Work\CsYbSe2_FvTh.mat');
        for n = 1:7
            Th1{n} = Dat.CsYbSe2_ThDeptCell.Th{n};
            Fv1{n}  = Dat.CsYbSe2_ThDeptCell.F{n};
            F_int = Interp1NonUnique(Th1{n},Fv1{n},Th_int1);
            %
            W(n) = sum((F_int(1:end-1) + F_int(2:end)).*diff(Th_int1))/2;
            Fv1{n} = Fv1{n};%-W(n)/pi;
            plot(Th1{n},Fv1{n},'sq','color',colz(n,:))
        end
        
        %Xdata = [Th1{1} Th1{2} Th1{3} Th1{4} Th1{5} Th1{6} Th1{7}];
        %Ydata = [Fv1{1} Fv1{2} Fv1{3} Fv1{4} Fv1{5} Fv1{6} Fv1{7}];
        
        Thset = Th1{1};
        
        
        if 0
            Xdata = [Hd_int{1}{2} Hd_int{1}{3} Hd_int{1}{4} Hd_int{1}{5} Hd_int{1}{6} Hd_int{1}{7} Hd_int{1}{8}...
                Hd_int{2}{2} Hd_int{2}{3} Hd_int{2}{4} Hd_int{2}{5} Hd_int{2}{6} Hd_int{2}{7} Hd_int{2}{8}...
                Th1{1} Th1{2} Th1{3} Th1{4} Th1{5} Th1{6} Th1{7}...
                ];
            %];
            Ydata = [kd_int{1}{2} kd_int{1}{3} kd_int{1}{4} kd_int{1}{5} kd_int{1}{6} kd_int{1}{7} kd_int{1}{8}...
                kd_int{2}{2} kd_int{2}{3} kd_int{2}{4} kd_int{2}{5} kd_int{2}{6} kd_int{2}{7} kd_int{2}{8}...
                Fv1{1} Fv1{2} Fv1{3} Fv1{4} Fv1{5} Fv1{6} Fv1{7}
                ];
            % ];
            f = @(x,X) Func5(x,X);
            
        elseif 1
            Xdata = [Hd_int{1}{4} Hd_int{1}{5} Hd_int{1}{6} Hd_int{1}{7} Hd_int{1}{8}];
            %];
            Ydata = [kd_int{1}{4} kd_int{1}{5} kd_int{1}{6} kd_int{1}{7} kd_int{1}{8}];
            % ];
            f = @(x,X) Func3([9.95,x(1),x(2),x(3),x(4),x(5),x(6),params(8),params(9)],X);
            
            
        elseif 0
            Xdata = [Hd_int{1}{5} Hd_int{1}{6} Hd_int{1}{7} Hd_int{1}{8}...
                Hd_int{2}{5} Hd_int{2}{6} Hd_int{2}{7} Hd_int{2}{8}];
            %];
            Ydata = [kd_int{1}{5} kd_int{1}{6} kd_int{1}{7} kd_int{1}{8}...
                kd_int{2}{5} kd_int{2}{6} kd_int{2}{7} kd_int{2}{8}];
            % ];
            f = @(x,X) Func10(x,X);
            
        elseif 0
            Xdata = [Hd_int{1}{2} Hd_int{1}{8} Hd_int{2}{2} Hd_int{2}{8} Th1{7}];
            %];
            Ydata = [kd_int{1}{2} kd_int{1}{8} kd_int{2}{2} kd_int{2}{8} Fv1{7}];
            % ];
            f = @(x,X) Func8(x,X);
            
        elseif 0
            Xdata = [Hd_int{1}{8} Hd_int{2}{8}];
            %];
            Ydata = [kd_int{1}{8} kd_int{2}{8}];
            % ];
            f = @(x,X) Func9(x,X);
            
            
        elseif 0
            Xdata = [Hd_int{1}{2} Hd_int{1}{3} Hd_int{1}{4} Hd_int{1}{5} Hd_int{1}{6} Hd_int{1}{7} Hd_int{1}{8}...
                ...%Hd_int{2}{2} Hd_int{2}{3} Hd_int{2}{4} Hd_int{2}{5} Hd_int{2}{6} Hd_int{2}{7} Hd_int{2}{8}...
                Th1{1} Th1{2} Th1{3} Th1{4} Th1{5} Th1{6} Th1{7}...
                ];
            %];
            Ydata = [kd_int{1}{2} kd_int{1}{3} kd_int{1}{4} kd_int{1}{5} kd_int{1}{6} kd_int{1}{7} kd_int{1}{8}...
                ...%kd_int{2}{2} kd_int{2}{3} kd_int{2}{4} kd_int{2}{5} kd_int{2}{6} kd_int{2}{7} kd_int{2}{8}...
                Fv1{1} Fv1{2} Fv1{3} Fv1{4} Fv1{5} Fv1{6} Fv1{7}
                ];
            % ];
            f = @(x,X) Func6(x,X);
            
        else
            Xdata = [...%Hd_int{1}{2} Hd_int{1}{3} Hd_int{1}{4} Hd_int{1}{5} Hd_int{1}{6} Hd_int{1}{7} Hd_int{1}{8}...
                ...%Hd_int{2}{2} Hd_int{2}{3} Hd_int{2}{4} Hd_int{2}{5} Hd_int{2}{6} Hd_int{2}{7} Hd_int{2}{8}...
                Th1{1} Th1{2} Th1{3} Th1{4} Th1{5} Th1{6} Th1{7}...
                ];
            %];
            Ydata = [...%kd_int{1}{2} kd_int{1}{3} kd_int{1}{4} kd_int{1}{5} kd_int{1}{6} kd_int{1}{7} kd_int{1}{8}...
                ...%kd_int{2}{2} kd_int{2}{3} kd_int{2}{4} kd_int{2}{5} kd_int{2}{6} kd_int{2}{7} kd_int{2}{8}...
                Fv1{1} Fv1{2} Fv1{3} Fv1{4} Fv1{5} Fv1{6} Fv1{7}
                ];
            % ];
            f = @(x,X) Func7(x,X);
            
        end
        
        
        subplot(2,2,1); title('H||ab');
        box on; xlabel('\mu_0H [T]'); ylabel('\Delta f [Hz]')
        
        subplot(2,2,2); title('H||c');
        box on; xlabel('\mu_0H [T]'); ylabel('\Delta f [Hz]')
        
        
        
        %x0 = [amp,paramsB(1),paramsB(2),paramsB(3),paramsB(4),paramsB(5),paramsB(6),paramsJJ(1),paramsJJ(2)];
        x0 = params(2:7);
        
        %plot(Xdata,f(x0,Xdata),'k-')
        
        %p = lsqcurvefit(f,x0,Xdata,Ydata,[amp,0,0,0,0,1000,1e-10,1e-10,1e-10,1e-10,1e-10,1e-10,1e-10,1e-10],[amp,Inf,Inf,Inf,Inf,1000,7/2,7/2,7/2,7/2,7/2,7/2,7/2,7/2]);
        %         p = lsqcurvefit(f,x0,Xdata,Ydata,...
        %         [amp,  1,  2,  0,  0,  0,1e-10,1e-10,1e-10,1e-10,1e-10,1e-10,1e-10,1e-10],... %lb
        %         [amp, 1, 2,Inf,Inf,Inf,  7/2,  7/2,  7/2,  7/2,  7/2,  7/2,  7/2,  7/2]); %ub
        %
        % Highly constrained
        p = lsqcurvefit(f,x0,Xdata,Ydata,...
            [-10, -10, -10, -10, -10, -10],...
            [10, 10, 10, 10, 10, 10]);
        %
        %p = x0;
        display(p)
        %
        output = p;
        %
        HCEF = HCEFf(p(1:6));
        JJ = JJf(params(8:9));
        
        Hset3 = 1:60;
        Thset1 = [pi/2,0];
        
        for n1 = [1,2]
            subplot(2,2,n1)
            for n2 = [2:8]
                plot(Hset3, 1e3*p(1)*k_vH(HCEF,JJ,T(n2),Hset3,Thset1(n1),.001)/w0(n2),'k--')
            end
            
        end
        
        thstart = -.86;
        thend = 2.78;
        N_Th = 100; dth = (thend-thstart)/N_Th;
        Th = linspace(thstart - 2*dth, thend + 2*dth, N_Th+4);
        
        H = [4:4:28];
        N_H = length(H);
        
        colz = varycolor(N_H);
        
        n1 = 1; Z = ones(N_H,N_Th);
        
        subplot(2,2,3)
        for t = 4
            
            n2 = 1;
            for h = H
                plot(Th(3:end-2), 1e3*p(1)*k_vTh_equispc(HCEF,JJ,t,h,Th,dth)/w0(2) , 'color', colz(n2,:), 'displayname', num2str(H(n2)) );
                n2=n2+1;
            end
            n1=n1+1;
        end
        
    case 'replot fit'
        
        Th_int1 = (pi/180)*[-45:10:135];
        
        %% Data Sets
        Hset = 4:4:28;
        
        colz = varycolor(7);
        
        clear Th1; clear Fv1;
        
        
        Hset1 = [4:2:28];
        Hset2 = [4:2:58];
        
        Hd_int{1} = {Hset2, Hset2, Hset1, Hset1, Hset1, Hset2, Hset1, Hset2};
        Hd_int{2} = {Hset2, Hset2, Hset1, Hset2, Hset1, Hset2, Hset1, Hset2};
        
        T= [1.6, 4, 6, 12, 20, 30, 50, 70];
        
        % 90 deg
        dat{1}{1} = load('p048_031020_TDO004.dat');
        dat{1}{2} = load('p043_031020_TDO003.dat');
        dat{1}{3} = load('p050_031020_TDO002.dat');
        dat{1}{4} = load('p052_031020_TDO001.dat');
        dat{1}{5} = load('p058_031020_TDO001.dat');
        dat{1}{6} = load('p067_031020_TDO001.dat');
        dat{1}{7} = load('p079_031020_TDO001.dat');
        dat{1}{8} = load('p092_031020_TDO001.dat');
        
        % 0 deg
        %cd /Users/iCAPOS/Google Drive/CsYbSe2_RTM
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
        
        for n1 = [1:2]
            subplot(2,2,n1); hold on;
            for n2 = [1:8]
                x = dat{n1}{n2}(startind:end-endind,1);
                y = dat{n1}{n2}(startind:end-endind,2);
                
                Hd{n1}{n2} = x;
                ind = find(x == min(x));
                f0{n1}(n2) = y(ind);
                kd{n1}{n2} = y-f0{n1}(n2);
                plot( Hd{n1}{n2}, kd{n1}{n2} ,'sq')
                
                kd_int{n1}{n2} = Interp1NonUnique(x,kd{n1}{n2}, Hd_int{n1}{n2} );
                %plot(Hd_int{n1}{n2},kd_int{n1}{n2} ,'ksq')
                %display(sum(isnan(kd_int{n1}{n2})));
            end
        end
        
        subplot(2,2,3); hold on;
        title('4 K'); xlabel('\theta'); ylabel('\Delta f [Hz]');
        grid on; box on;
        Dat = load('C:\Users\Chris\Documents\Lab Work\CsYbSe2_FvTh.mat');
        for n = 1:7
            Th1{n} = Dat.CsYbSe2_ThDeptCell.Th{n};
            Fv1{n}  = Dat.CsYbSe2_ThDeptCell.F{n};
            F_int = Interp1NonUnique(Th1{n},Fv1{n},Th_int1);
            %
            W(n) = sum((F_int(1:end-1) + F_int(2:end)).*diff(Th_int1))/2;
            Fv1{n} = Fv1{n};%-W(n)/pi;
            plot(Th1{n},Fv1{n},'sq','color',colz(n,:))
        end
        
        subplot(2,2,1); title('H||ab');
        box on; xlabel('\mu_0H [T]'); ylabel('\Delta f [Hz]')
        
        subplot(2,2,2); title('H||c');
        box on; xlabel('\mu_0H [T]'); ylabel('\Delta f [Hz]')
        
        
        
        %x0 = [amp,paramsB(1),paramsB(2),paramsB(3),paramsB(4),paramsB(5),paramsB(6),paramsJJ(1),paramsJJ(2)];
        p = params;
        
        HCEF = HCEFf(p(2:7));
        JJ = JJf(p(8:9));
        
        Hset3 = 1:60;
        Thset1 = [pi/2,0];
        
        for n1 = [1,2]
            subplot(2,2,n1)
            for n2 = [2:8]
                plot(Hset3, 1e3*p(1)*k_vH(HCEF,JJ,T(n2),Hset3,Thset1(n1),.001)/w0(n2),'k--')
            end
            
        end
        
        thstart = -.86;
        thend = 2.78;
        N_Th = 100; dth = (thend-thstart)/N_Th;
        Th = linspace(thstart - 2*dth, thend + 2*dth, N_Th+4);
        
        H = [4:4:28];
        N_H = length(H);
        
        colz = varycolor(N_H);
        
        n1 = 1; Z = ones(N_H,N_Th);
        
        subplot(2,2,3)
        for t = 4
            
            n2 = 1;
            for h = H
                plot(Th(3:end-2), 1e3*p(1)*k_vTh_equispc(HCEF,JJ,t,h,Th,dth)/w0(2) , 'color', colz(n2,:), 'displayname', num2str(H(n2)) );
                n2=n2+1;
            end
            n1=n1+1;
        end
        
    case 'High Field 70K model test'
        
        a1 = [-335.2,   -629,  82.92,  237.1];
        a2 = [  -209, -993.5, -12.01, -128.8];
        
        X = linspace(0,pi,30);
        
        figure; hold on;
        
        for n = 1:4
            Y{n} = a1(n)*cos(2*X)+a2(n)*cos(4*X);
            plot(X,Y{n})
        end
        
        Xdata = [X X X X];
        Ydata = [Y{1} Y{2} Y{3} Y{4}];
        
        
        
        f = @(x,X) Func_test2(x,X);
        %p = params;
        p = lsqcurvefit(f,params,Xdata,Ydata,...
            [amp  , -1, -1, -1, -1, -1, -1,   0,   0],...
            [amp,  1,  1,  1,  1,  1,  1, Inf, Inf]);
        
        HCEF = HCEFf(p(2:7));
        JJ = JJf(p(8:9));
        display(p);
        plot(X,1e3*p(1)*k_vTh(HCEF,JJ,4,28,X,.001)/w0(2),'k--')
        plot(X,1e3*p(1)*k_vTh(HCEF,JJ,4,56,X,.001)/w0(2),'k--')
        plot(X,1e3*p(1)*k_vTh(HCEF,JJ,70,28,X,.001)/w0(8),'k--')
        plot(X,1e3*p(1)*k_vTh(HCEF,JJ,70,56,X,.001)/w0(8),'k--')
        
    case 'leading order 2theta 4theta'
        Th_int1 = (pi/180)*[-45:10:135];
        
        %% Data Sets
        Hset = 4:4:28;
        
        colz = varycolor(7);
        
        clear Th1; clear Fv1;
        
        
        Hset1 = [4:4:28];
        Hset2 = [4:4:58];
        
        Hd_int{1} = {Hset2, Hset2, Hset1, Hset1, Hset1, Hset2, Hset1, Hset2};
        Hd_int{2} = {Hset2, Hset2, Hset1, Hset1, Hset1, Hset2, Hset1, Hset2};
        
        T= [1.6, 4, 6, 12, 20, 30, 50, 70];
        
        % 90 deg
        dat{1}{1} = load('p048_031020_TDO004.dat');
        dat{1}{2} = load('p043_031020_TDO003.dat');
        dat{1}{3} = load('p050_031020_TDO002.dat');
        dat{1}{4} = load('p052_031020_TDO001.dat');
        dat{1}{5} = load('p058_031020_TDO001.dat');
        dat{1}{6} = load('p067_031020_TDO001.dat');
        dat{1}{7} = load('p079_031020_TDO001.dat');
        dat{1}{8} = load('p092_031020_TDO001.dat');
        
        % 0 deg
        %cd /Users/iCAPOS/Google Drive/CsYbSe2_RTM
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
        
        %figure; hold on;
        
        for n1 = [1:2]
            
            for n2 = [1:8]
                x = dat{n1}{n2}(startind:end-endind,1);
                y = dat{n1}{n2}(startind:end-endind,2);
                
                Hd{n1}{n2} = x;
                ind = find(x == min(x));
                f0{n1}(n2) = y(ind);
                kd{n1}{n2} = y-f0{n1}(n2);
                %plot( Hd{n1}{n2}, kd{n1}{n2} ,'sq')
                
                kd_int{n1}{n2} = Interp1NonUnique(x,kd{n1}{n2}, Hset2);
                %plot(Hd_int{n1}{n2},kd_int{n1}{n2} ,'ksq')
                display(sum(isnan(kd_int{n1}{n2})));
            end
        end
        
        figure; hold on;
        
        %ft = fittype( 'a*power(x,2) + b*power(x,4)', 'independent', 'x', 'dependent', 'y' );
        %opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
        %opts.Display = 'Off';
        
        
        for n2 = [1:8]
            a1(n2,:) = (kd_int{2}{n2} - kd_int{1}{n2})./2;
            a2(n2,:) = (kd_int{2}{n2} + kd_int{1}{n2})./2;
            
            subplot(1,2,1);  hold on;
            plot(Hset2, a1(n2,:), 'displayname', [num2str(T(n2)) ' K'])
            
            %[fitresult, gof] = fit( Hd_int{1}{n2}', a1', ft, opts );
            %plot(fitresult,'k--')
            
            subplot(1,2,2);  hold on;
            plot(Hset2, a2(n2,:), 'displayname', [num2str(T(n2)) ' K'])
            
            %[fitresult, gof] = fit( Hd_int{1}{n2}', a2', ft, opts );
            %plot(fitresult,'k--')
            
        end
        
        figure; hold on;
        colz = varycolor(14);
        X = linspace(0,pi,100);
        for n1 = 2:14
            subplot(1,2,1);  hold on;
            plot(X, a1(2,n1)*cos(2*X) + a2(2,n1)*cos(4*X) , 'color', colz(n1,:))
            
            subplot(1,2,2);  hold on;
            plot(X, a1(8,n1)*cos(2*X) + a2(8,n1)*cos(4*X) , 'color', colz(n1,:))
        end
        
        %         subplot(1,2,1)
        %         legend show;
        %
        %         subplot(1,2,2)
        %         legend show;
        
    case 'HC 0,9 T'
        
        HCEF = HCEFf(paramsB);
        JJ = JJf(paramsJJ);
        
        h = 9;
        
        N = 1e3;
        B2 = linspace(1./(kB*2), 1./(kB*330) , N);
        T2 = 1./(kB*B2);
        
        HZeeman_x = -muB*h*Jx;
        HZeeman_z = -muB*h*Jz;
        
        % Calculate Partion Functions : expm (not exp) is crucial for proper
        % exponentiation of matrices.
        n = 1;
        for t = T2
            b = 1./(kB*t);
            
            mv2 = M(HCEF, JJ, [0,0,0], t);
            HMF2 = HMF(HCEF, JJ, [0,0,0], mv2);
            [P2,D2] = eig(HMF2);
            Ev2 = sort(real([D2(1,1) D2(2,2) D2(3,3) D2(4,4) D2(5,5) D2(6,6) D2(7,7) D2(8,8)]));
            E02(n) = Ev2(1);
            
            Z0(n) = trace(expm(-b*(HMF2 - E02(n)*ID)));
            
            mv2 = M(HCEF, JJ, [h,0,0], t);
            HMF2 = HMF(HCEF, JJ, [h,0,0], mv2);
            
            Z1(n) = trace(expm(-b*(HMF2 - E02(n)*ID)));
            
            mv2 = M(HCEF, JJ, [0,0,h], t);
            HMF2 = HMF(HCEF, JJ, [0,0,h], mv2);
            
            Z2(n) = trace(expm(-b*(HMF2 - E02(n)*ID)));
            n = n + 1;
        end
        
        dB = (max(B2)-min(B2))/(N-1);
        
        
        FF0 = log(Z0);
        FF1 = log(Z1);
        FF2 = log(Z2);
        
        
        %         C1 = [-1/2, 0, 1/2];
        %
        %         dZ0 =  (C1(1)*Z0(1:end-2) + C1(2)*Z0(2:end-1) + C1(3)*Z0(3:end))/dB;
        %         dZ1 =  (C1(1)*Z1(1:end-2) + C1(2)*Z1(2:end-1) + C1(3)*Z1(3:end))/dB;
        %         dZ2 =  (C1(1)*Z2(1:end-2) + C1(2)*Z2(2:end-1) + C1(3)*Z2(3:end))/dB;
        %
        %         A0 = dZ0./Z0(2:end-1) ;
        %         A1 = dZ1./Z1(2:end-1) ;
        %         A2 = dZ2./Z2(2:end-1) ;
        %
        %         dZ0_2 = (C1(1)*A0(1:end-2) + C1(2)*A0(2:end-1) + C1(3)*A0(3:end))/dB;
        %         dZ1_2 = (C1(1)*A1(1:end-2) + C1(2)*A1(2:end-1) + C1(3)*A1(3:end))/dB;
        %         dZ2_2 = (C1(1)*A2(1:end-2) + C1(2)*A2(2:end-1) + C1(3)*A2(3:end))/dB;
        
        dZ0_2 = (C2(1)*FF0(1:end-4) + C2(2)*FF0(2:end-3) + C2(3)*FF0(3:end-2) + C2(4)*FF0(4:end-1) + C2(5)*FF0(5:end))/dB/dB;
        dZ1_2 = (C2(1)*FF1(1:end-4) + C2(2)*FF1(2:end-3) + C2(3)*FF1(3:end-2) + C2(4)*FF1(4:end-1) + C2(5)*FF1(5:end))/dB/dB;
        dZ2_2 = (C2(1)*FF2(1:end-4) + C2(2)*FF2(2:end-3) + C2(3)*FF2(3:end-2) + C2(4)*FF2(4:end-1) + C2(5)*FF2(5:end))/dB/dB;
        
        figure; hold on;
        
        R = 8.3145; % [J/mol K]
        Y0 = power(B2(3:end-2),2).*dZ0_2*R;
        Y1 = power(B2(3:end-2),2).*dZ1_2*R;
        Y2 = power(B2(3:end-2),2).*dZ2_2*R;
        
        plot(T2(3:end-2), Y0, 'k', 'displayname', '0 T')
        plot(T2(3:end-2), Y1, 'r', 'displayname', '9 T || ab')
        plot(T2(3:end-2), Y2, 'b', 'displayname', '9 T || c')
        
        grid on; box on; set(gca, 'fontsize', 20);
        xlabel('T [K]'); ylabel('C [J/K mol]')
        
    case 'LT HC'
        
        HCEF = HCEFf(paramsB);
        JJ = JJf(paramsJJ);
        
        h = 9;
        
        N = 1e2;
        B2 = linspace(1./(kB*.4), 1./(kB*10) , N); dB = (max(B2)-min(B2))/(N-1);
        T2 = 1./(kB*B2);
        
        HZeeman_x = -muB*h*Jx;
        HZeeman_z = -muB*h*Jz;
        
        
        fields = [0:2:10];
        % Calculate Partion Functions : expm (not exp) is crucial for proper
        % exponentiation of matrices.
        
        figure; hold on; R = 8.3145; % [J/mol K]
        
        q = 1;
        for h = fields
            
            n = 1;
            for t = T2
                b = 1./(kB*t);
                
                mv2 = M(HCEF, JJ, [h,0,0], t);
                HMF2 = HMF(HCEF, JJ, [h,0,0], mv2);
                [P2,D2] = eig(HMF2);
                Ev2 = sort(real([D2(1,1) D2(2,2) D2(3,3) D2(4,4) D2(5,5) D2(6,6) D2(7,7) D2(8,8)]));
                E02(n) = Ev2(1);
                
                Za{q}(n) = trace(expm(-b*(HMF2 - E02(n)*ID)));
                
                
                mv2 = M(HCEF, JJ, [0,0,h], t);
                HMF2 = HMF(HCEF, JJ, [0,0,h], mv2);
                [P2,D2] = eig(HMF2);
                Ev2 = sort(real([D2(1,1) D2(2,2) D2(3,3) D2(4,4) D2(5,5) D2(6,6) D2(7,7) D2(8,8)]));
                E02(n) = Ev2(1);
                
                Zc{q}(n) = trace(expm(-b*(HMF2 - E02(n)*ID)));
                
                n = n + 1;
                
            end
            
            FFa{q} = log(Za{q});
            FFc{q} = log(Zc{q});
            
            
            dZa{q} = (C2(1)*FFa{q}(1:end-4) + C2(2)*FFa{q}(2:end-3) + C2(3)*FFa{q}(3:end-2) + C2(4)*FFa{q}(4:end-1) + C2(5)*FFa{q}(5:end))/dB/dB;
            dZc{q} = (C2(1)*FFc{q}(1:end-4) + C2(2)*FFc{q}(2:end-3) + C2(3)*FFc{q}(3:end-2) + C2(4)*FFc{q}(4:end-1) + C2(5)*FFc{q}(5:end))/dB/dB;
            
            Ya = power(B2(3:end-2),2).*dZa{q}*R;
            Yc = power(B2(3:end-2),2).*dZc{q}*R;
            
            subplot(1,2,1); hold on
            plot(T2(3:end-2), Ya, 'sq-', 'displayname', [num2str(h) ' T'])
            
            subplot(1,2,2); hold on
            plot(T2(3:end-2), Yc, 'sq-', 'displayname', [num2str(h) ' T'])
            
            q = q + 1;
            
        end
        
        subplot(1,2,1);
        grid on; box on; set(gca, 'fontsize', 20);
        xlabel('T [K]'); ylabel('C [J/K mol]')
        
        subplot(1,2,2);
        grid on; box on; set(gca, 'fontsize', 20);
        xlabel('T [K]'); ylabel('C [J/K mol]')
        
        
        
        
end





end