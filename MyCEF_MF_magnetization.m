%% Constants
% % fundamentals
muB = 5.78828e-2; % [meV/T];
kB  = 8.617e-2  ; % [meV/K];
% gj = 1+(J*(J+1)+S*(S+1)-L*(L+1))/(2*J*(J+1))
gJ = 6/5; %  Er3+ : L=6, S=3/2, J=15/2 g-lande factor; 
z= 6; % nearest octahedron environment
A1 = 6*kB/(muB*gJ);

% % 2nd Deriv coefficients
C1 = [1/12, 	-2/3, 	0,      2/3, 	-1/12];
C2 = [-1/12 	4/3 	-5/2 	4/3 	-1/12];

% Allen's PRB parameter
 B20 = -3.559e-2 %[meV]; err1 = 0.64e-2
 B40 = -3.849e-4 %[meV]; err2 = 0.11e-4
 B43 = -1.393e-2 %[meV]; err3 = 0.03e-2
 B60 = 3.154e-6 %[meV]; err4 = 0.04e-6
 B63 = -4.695e-6 %[meV];  err5 = 3.56e-6
 B66 = 3.381e-5 %[meV]; err6 = 0.37e-5
% 
% % % Hope 's parameter
%  B20 = -4.910422e-2 %[meV]; err1 = 0.64e-2
%  B40 = -3.64830e-4 %[meV]; err2 = 0.11e-4
%  B43 = -1.474077e-2 %[meV]; err3 = 0.03e-2
%  B60 = 3.1547e-6 %[meV]; err4 = 0.04e-6
%  B63 = 3.2378e-6 %[meV];  err5 = 3.56e-6
%  B66 = 4.2797e-5 %[meV]; err6 = 0.37e-5

 %% Spin 15/2 matrices, {|15/2,15/2>.,...,|15/2,-15/2>| basis

ID = diag(ones(1,16));

Jz = diag([15/2:-1:-15/2]);

% Jp =  [[zeros(15,1) , diag(sqrt([15, 28, 39, 48, 55, 60, 63, 64, 63, 60, 55, 48, 39, 28, 15]))] ; zeros(1,16) ];
% Jm =  [zeros(1,16) ; [diag(sqrt([15, 28, 39, 48, 55, 60, 63, 64, 63, 60, 55, 48, 39, 28, 15])) , zeros(15,1) ]];

Jp=spinOp(15/2,'p');
Jm=spinOp(15/2,'m');

Jx = (Jp+Jm)/2;

Jy = (Jp-Jm)/2i;

%% Steven's Operators
J=15/2; X = J*(J+1); A = Jp*Jp*Jp + Jm*Jm*Jm;

O20 = 3*Jz*Jz - X*ID;
O40 = 35*power(Jz,4) - (30*X - 25)*Jz*Jz + (3*X*X - 6*X)*ID;
O60 = 231*power(Jz,6) - (315*X-735)*power(Jz,4) + (105*X*X - 525*X +294)*power(Jz,2) - (5*X*X*X + 40*X*X -60*X)*ID;
O43 = (1/4)*( (A)*Jz + Jz*(A) );
O63 = (1/4)*( A*(11*power(Jz,3) - (3*X + 59)*Jz ) + (11*power(Jz,3) -(3*X + 59)*Jz)*A );
O66 = (1/2)*(Jp*Jp*Jp*Jp*Jp*Jp + Jm*Jm*Jm*Jm*Jm*Jm);



Bnm = [B20 B40 B43 B60 B63 B66];
HCEFf =  @(Bnm) Bnm(1)*O20 + Bnm(2)*O40 + Bnm(3)*O43 + Bnm(4)*O60 + Bnm(5)*O63 + Bnm(6)*O66;


%% Exchange interactions in order to avoid confusion with spin operator Jx, Jy, Jz
%  exchange interaction is JJab, JJc in the unit of   temperature in meV 
%  exchang interaction matrix, JJ = [JJab,0,0; 0,JJab,0; 0,0,JJc]; 
JJf = @(JJ) [JJ(1),0,0;0,JJ(1),0;0,0,JJ(2)]*kB;  

%%  Define mean field MF Hamiltonian
% given a HCEF matrix (HCEFm) with Bnm
% given a 3x3 J-J coupling constant matrix JJ =[JJ(1),0,0;0,JJ(1),0;0,0,JJ(2)]
% at finite field hv = (hx,hy,hz)
% in terms of MF magnetizations mv = (mx,my,mz)
HMF = @(HCEFm, JJm, hv, mv) HCEFm - (z/2)*(mv*JJm*mv')*ID...
    - (muB*gJ*hv(1) - z*JJm(1,:)*mv')*Jx...
    - (muB*gJ*hv(2) - z*JJm(2,:)*mv')*Jy...
    - (muB*gJ*hv(3) - z*JJm(3,:)*mv')*Jz;


%% Magnetization calc
%function out = M(HCEFm, JJab, hv, t)
    t = 2; HCEF = HCEFf(Bnm);    
    JJ =JJf([-0.0002 0.0024]); %JJx = 5K, JJz = 1K; 
    B =  1./(kB*t);
    field = 0:0.1:8;
    dh = 0.1;
    
   
    HCEF = HCEFf(Bnm);
        
        [P,D] = eig(HCEF + Jz*1e-10);        
    
        En= real(diag(D));
        Ev = sort(En);
        
        E0 = Ev(1);
        
        ind1 = find(En==Ev(1));
        ind2 = find(En==Ev(2));
        
        ev1 = P(:,ind1);
        ev2 = P(:,ind2);
        
        
        
        gzz0 = 2*gJ*abs(ev1'*Jz*ev1);
        gxx0 = 2*gJ*abs(ev1'*Jx*ev2);
        gyy0 = 2*gJ*abs(ev1'*Jy*ev2);


    for i=1:length(field)
        hv1 = [field(i) 0 0];
        hv3= [0,0,field(i)];
        f1 = @(mx,my,mz) trace( Jx*expm(-B*(HMF(HCEF, JJ, hv1, [mx,my,mz]) - E0*ID)) )/trace( expm(-B*(HMF(HCEF, JJ, hv1, [mx,my,mz]) - E0*ID)) ) - mx;
        %f2 = @(mx,my,mz) trace( Jy*expm(-B*(HMF(HCEF, JJab, hv, [mx,my,mz]) - E0*ID)) )/trace( expm(-B*(HMF(HCEFm, JJab, hv, [mx,my,mz]) - E0*ID)) ) - my;
        f3 = @(mx,my,mz) trace( Jz*expm(-B*(HMF(HCEF, JJ, hv3, [mx,my,mz]) - E0*ID)) )/trace( expm(-B*(HMF(HCEF, JJ, hv3, [mx,my,mz]) - E0*ID)) ) - mz;
        
        mx_g = 0; my_g = 0; mz_g = 0;
        
        % iteration to ontain converging  value 
        itrnum = 10; 
        for m = 1:itrnum
          
            f4 = @(x) real(f1(x,my_g,mz_g)); mx_g = fzero(f4,mx_g);
            %f5 = @(x) real(f2(mx_g,  x ,mz_g)); my_g = fzero(f5,my_g);
            f6 = @(x) real(f3(mx_g,my_g,  x )); mz_g = fzero(f6,mz_g); 
            
        end
        mxarr(i)= mx_g; 
%        myarr(n)=my_g;
        mzarr(i) =mz_g;
    end
        %out = [mx_g, my_g, mz_g];
    figure(2); clf;
    plot( field, mxarr,'-o','DisplayName','Mx'); hold on; 
     plot( field, mzarr,'-o','DisplayName','Mz');hold off; legend;
 %end
%% Field dependence  spectrum calculation only with HCEG + Zeeman
        HCEF = HCEFf(Bnm);
        
        [P,D] = eig(HCEF + Jz*1e-10);        
    
        En= real(diag(D));
        Ev = sort(En);
        
        E0 = Ev(1);
        
        ind1 = find(En==Ev(1));
        ind2 = find(En==Ev(2));
        
        ev1 = P(:,ind1);
        ev2 = P(:,ind2);
        
        
        
        gzz0 = 2*gJ*abs(ev1'*Jz*ev1);
        gxx0 = 2*gJ*abs(ev1'*Jx*ev2);
        gyy0 = 2*gJ*abs(ev1'*Jy*ev2);
        
        H = 0:1:10;
        
        % [P,D] = eig(HCEF);
        % E0 = min(min(D));
        
        n = 1;

        for h = H
            % H || x
            HZeeman = -gJ*muB*h*Jx;
            [P,D] = eig(HCEF+HZeeman);
            
            for m = 1:16
                Eigenval(m) = real(D(m,m));
            end
            Eigenval = sort(Eigenval);
            
            for m = 1:16
                E1{m}(n) = Eigenval(m);
            end
            
            % H || y
            HZeeman = -gJ*muB*h*Jy;
            [P,D] = eig(HCEF+HZeeman);
            
            for m = 1:16
                Eigenval(m) = real(D(m,m));
            end
            Eigenval = sort(Eigenval);
            for m = 1:16
                E2{m}(n) = Eigenval(m);
            end
            
            % H || z
            HZeeman = -gJ*muB*h*Jz;
            [P,D] = eig(HCEF+HZeeman);
            
            for m = 1:16
                Eigenval(m) = real(D(m,m));
            end
            Eigenval = sort(Eigenval);
            for m = 1:16
                E3{m}(n) = Eigenval(m);
            end
            
            n = n + 1;
        end
        
        figure(10); clf ; hold on
        colz = {'k','k','r','r','b','b','c','c','g','g','m','m', [0.87 0.49 0],[0.87 0.49 0],[0.5 0.5 0.5],[0.5 0.5 0.5]};
        for m = 1:16
            subplot(1,3,1); hold on;
            plot(H,E1{m}-E0,'color',colz{m},'Linewidth',2)
            
            subplot(1,3,2); hold on;
            plot(H,E2{m}-E0,'color',colz{m},'Linewidth',2)
            
            subplot(1,3,3); hold on;
            plot(H,E3{m}-E0,'color',colz{m},'Linewidth',2)
        end
        
        subplot(1,3,1);
        xlabel('\mu_0H [T]'); ylabel('\omega [meV]'); grid on; box on;
        %title(['H||a,   g_{\perp} = ' num2str(gxx0)])
        title(['H||a'])%
        
        subplot(1,3,2);
        xlabel('\mu_0H [T]'); ylabel('\omega [meV]'); grid on; box on;
        title(['H||b'])%   g_{yy} = ' num2str(gyy0)])
        
        subplot(1,3,3);
        xlabel('\mu_0H [T]'); ylabel('\omega [meV]'); grid on; box on;
        %title(['H||c,   g_{||} = ' num2str(gzz0)])
        title(['H||c'])%
%% Full spectra with meanfield 
        temp = 2; 
        JJ =JJf ([-0.0002 0.0024]); % exchange energy JJab= 0.0002 meV
        
        H = 0:10;
        n= 1 
        for h=H
            hv1 = [h, 0, 0];
            mv0 = M(HCEF, JJ, hv1, temp);
            HMFa = HMF(HCEF, JJ, hv1, mv0);

            [P1, D1] = eig(HMFa)
            for q = 1:16
                Eigenval(q) = real(D1(q,q));
            end
            Eigenval = sort(Eigenval);
            for q = 1:16
                E1{q}(n) = Eigenval(q);
            end

            hv1 = [0,0,h];
            mv0 = M(HCEF, JJ, hv1, t);
            HMFc = HMF(HCEF, JJ, hv1, mv0);
            
            [P1,D1] = eig(HMFc);
            
            for q = 1:16
                Eigenval(q) = real(D1(q,q));
            end
            Eigenval = sort(Eigenval);
            for q = 1:16
                E2{q}(n) = Eigenval(q);
            end
            
            n = n + 1;
        end

        figure(9); hold on
        
        for q = 1:16
            subplot(1,3,1); hold on;
            plot(H,E1{q}-E_0,'--','color',colz{q},'Linewidth', 2)
            
            subplot(1,3,3); hold on;
            plot(H,E2{q}-E_0,'--','color',colz{q},'Linewidth',2)
        end
        
        % subplot(1,2,1);
        % xlabel('\mu_0H [T]'); ylabel('\omega [meV]'); grid on; box on;
        % title(['H||ab,   g_{\perp} = ' num2str(gxx0)])
        % 
        % subplot(1,2,2);
        % xlabel('\mu_0H [T]'); ylabel('\omega [meV]'); grid on; box on;
        % title(['H||c,   g_{||} = ' num2str(gzz0)])
 %% test 

  mv2 = M(HCEF, JJ, [0,0,1], 2);
        %% Self consistent determination of magnetization 
function out = M(HCEFm, JJab, hv, t)
        
        B =  1./(kB*t);
        
        [P,D] = eig(HCEFm);
        Ev = sort(real(diag(D)));
        E0 = Ev(1);
        
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
