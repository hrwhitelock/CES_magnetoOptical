function output = CEF_Spectrum_D3d_S15o2(option,params)

%% 2nd Deriv coefficients
C1 = [1/12, 	-2/3, 	0,      2/3, 	-1/12];
C2 = [-1/12 	4/3 	-5/2 	4/3 	-1/12];

%% Spin 15/2 matrices, {|15/2,15/2>.,...,|15/2,-15/2>| basis


ID = diag(ones(1,16));

Jz = diag([15/2:-1:-15/2]);

Jp =  [[zeros(15,1) , diag(sqrt([15, 28, 39, 48, 55, 60, 63, 64, 63, 60, 55, 48, 39, 28, 15]))] ; zeros(1,16) ];

Jm =  [zeros(1,16) ; [diag(sqrt([15, 28, 39, 48, 55, 60, 63, 64, 63, 60, 55, 48, 39, 28, 15])) , zeros(15,1) ]];

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

B20 = params(1);
B40 = params(2);
B43 = params(3);
B60 = params(4);
B63 = params(5);
B66 = params(6);

HCEF = B20*O20 + B40*O40 + B43*O43 + B60*O60 + B63*O63 + B66*O66;

muB = 5.78828e-2; % [meV/T];
kB  = 8.617e-2  ; % [meV/K];
gJ = 6/5; % L=3, S=1/2, J=7/2 g-lande factor;

%% Options
Bmn_CsErSe2 = [-3.559e-2, -3.849e-4,  -1.393e-2,  3.154e-6,  -4.695e-6,  3.381e-5];

switch option
    
    
    
    case 'ZF wavefunctions'
        
        [P,D] = eig(HCEF + Jz*1e-6);
        display(round(P,4));
        
    case 'eig'
        
        [P,D] = eig(HCEF);
        
        En = real([D(1,1) D(2,2) D(3,3) D(4,4) D(5,5) D(6,6) D(7,7) D(8,8) ]);
        Ev = sort(En);
        
        E0 = Ev(1);
        E1 = Ev(3);
        E2 = Ev(5);
        E3 = Ev(7);
        
        %         E2 = 1000;
        %         E3 = 1000;
        
        display([E0,E1,E2,E3])
        
    case 'Schottky'
        
        [P,D] = eig(HCEF);
        
        En = real([D(1,1) D(2,2) D(3,3) D(4,4) D(5,5) D(6,6) D(7,7) D(8,8)]);
        Ev = sort(En);
        
        E0 = Ev(1);
        E1 = Ev(3) - E0;
        E2 = Ev(5) - E0;
        E3 = Ev(7) - E0;
        
        %         E2 = 1000;
        %         E3 = 1000;
        
        display([E1,E2,E3])
        
        %E1 = 30;
        %E2 = 1000;
        %E3 = 1000;
        
        N  = 1e4;
        B  = linspace(1/(kB*301),.9/kB,N); T = 1./(kB*B);
        dB = (max(B)-min(B))/(N-1);
        B_ext = [(B(1)-2*dB) (B(1)-dB) B (B(end)+dB) (B(end)+2*dB)];
        
        logZ = log(2*(1 + exp(-B_ext*E1) + exp(-B_ext*E2) + exp(-B_ext*E3)));
        
        Csch = power(B,2).*(C2(1)*logZ(1:end-4) + C2(2)*logZ(2:end-3)...
            + C2(3)*logZ(3:end-2) + C2(4)*logZ(4:end-1)...
            + C2(5)*logZ(5:end))/dB/dB;
        
        figure; subplot(3,1,1); hold on;
        R = 8.3145; % [J/mol K]
        plot(T, R*Csch); xlim([0,300]); ylim([0,6])
        xlabel('T [K]'); ylabel('C [J/mol K]')
        
        %         C2lvl1 = power(B*E1,2).*exp(B*E1)./power(1 + exp(B*E1),2);
        %         C2lvl2 = power(B*E2,2).*exp(B*E2)./power(1 + exp(B*E2),2);
        %         plot(T, C2lvl1 ,'r--')
        %         plot(T, C2lvl2 ,'g--')
        %         plot(T, C2lvl1 + C2lvl2,'c--')
        
        plot([E1 E1]/kB,[0,6],'k--')
        plot([E2 E2]/kB,[0,6],'k--')
        plot([E3 E3]/kB,[0,6],'k--')
        
        
        dx = .001;
        X = 1-2*dx:dx:300+2*dx;
        Y = Interp1NonUnique(T, R*Csch, X);
        
        dYdX = (C1(1)*Y(1:end-4) + C1(2)*Y(2:end-3) + C1(3)*Y(3:end-2) + C1(4)*Y(4:end-1) + C1(5)*Y(5:end))/dx;
        
        subplot(3,1,2); hold on;
        plot(X(3:end-2), dYdX);
        xlim([0,300]);
        xlabel('T [K]'); ylabel('dC/dT [J/mol K^2]')
        
        
        plot([E1 E1]/kB,[-.05,.2],'k--')
        plot([E2 E2]/kB,[-.05,.2],'k--')
        plot([E3 E3]/kB,[-.05,.2],'k--')
        
        dYdX2 = (C2(1)*Y(1:end-4) + C2(2)*Y(2:end-3) + C2(3)*Y(3:end-2) + C2(4)*Y(4:end-1) + C2(5)*Y(5:end))/dx/dx;
        
        subplot(3,1,3); hold on;
        plot(X(3:end-2), dYdX2);
        xlim([0,300]);
        
    case 'all eigenvectors'
        
        [P,D] = eig(HCEF + Jz*1e-10);
        
        En = real([D(1,1) D(2,2) D(3,3) D(4,4) D(5,5) D(6,6) D(7,7) D(8,8)]);
        Ev = sort(En);
        
        for n = 1:4
            %display([2*(n-1) + 1, 2*n]);
            ind1 = find(En==Ev(2*(n-1) + 1));
            ind2 = find(En==Ev(2*n));
            
            ev1 = P(:,ind1);
            ev2 = P(:,ind2);
            
            display(ev1')
            display(ev2')
        end
          
    case 'all g tensors'
        
        colz = {'k','r','b','c'};
        
        [P,D] = eig(HCEF + Jz*1e-10);
        
        En = real([D(1,1) D(2,2) D(3,3) D(4,4) D(5,5) D(6,6) D(7,7) D(8,8)]);
        Ev = sort(En);
        E0 = Ev(1);
        
        
        H = 0:60;
        
        for n = 1:4
            %display([2*(n-1) + 1, 2*n]);
            ind1 = find(En==Ev(2*(n-1) + 1));
            ind2 = find(En==Ev(2*n));
            
            ev1 = P(:,ind1);
            ev2 = P(:,ind2);
            
            %gzz = 2*gJ*abs(ev1'*Jz*ev1);
            %gxx = 2*gJ*abs(ev1'*Jx*ev2);
            
            gzz = 2*gJ*(ev1'*Jz*ev1);
            gxx = 2*gJ*(ev1'*Jx*ev2);
            
            display(['gxx' num2str(n-1) ' = ' num2str(gxx)]);
            display(['gzz' num2str(n-1) ' = ' num2str(gzz)]);
            
            subplot(1,2,1); hold on;
            plot(H, Ev(2*n) - E0 + muB*H*gxx/2, '--' ,'color',colz{n})
            plot(H, Ev(2*n) - E0 - muB*H*gxx/2, '--' ,'color',colz{n})
            
            subplot(1,2,2); hold on;
            plot(H, Ev(2*n) - E0 + muB*H*gzz/2, '--' ,'color',colz{n})
            plot(H, Ev(2*n) - E0 - muB*H*gzz/2, '--' ,'color',colz{n})
            
            
            
        end
        
    case 'Secondary Params'
        
        %colz = {'k','r','b','c'};
        
        [P,D] = eig(HCEF + Jz*1e-10);
        
        En = real([D(1,1) D(2,2) D(3,3) D(4,4) D(5,5) D(6,6) D(7,7) D(8,8)]);
        Ev = sort(En);
        E0 = Ev(1);
        
        
        H = 0:60;
        
        for n = 1:4
            %display([2*(n-1) + 1, 2*n]);
            ind1 = find(En==Ev(2*(n-1) + 1));
            ind2 = find(En==Ev(2*n));
            
            ev1{n} = P(:,ind1(1));
            ev2{n} = P(:,ind2(end));
            
            gzz{n} = abs(2*gJ*(ev1{n}'*Jz*ev1{n}));
            gxx{n} = abs(2*gJ*(ev1{n}'*Jx*ev2{n}));
            
            %display(['gxx' num2str(n-1) ' = ' num2str(gxx)]);
            %display(['gzz' num2str(n-1) ' = ' num2str(gzz)]);
            
        end
        
        n = 1;
        az0{n} = -power(muB*gJ,2).*((power(ev1{n}'*Jz*ev1{2},2) + power(ev1{n}'*Jz*ev2{2},2) )./(Ev(3) - E0)...
            +(power(ev1{n}'*Jz*ev1{3},2) + power(ev1{n}'*Jz*ev2{3},2) )./(Ev(5) - E0)...
            +(power(ev1{n}'*Jz*ev1{4},2) + power(ev1{n}'*Jz*ev2{4},2) )./(Ev(7) - E0));
        %
        ax0{n} = -power(muB*gJ,2).*((power(ev1{n}'*Jx*ev1{2},2) + power(ev1{n}'*Jx*ev2{2},2) )./(Ev(3) - E0)...
            +(power(ev1{n}'*Jx*ev1{3},2) + power(ev1{n}'*Jx*ev2{3},2) )./(Ev(5) - E0)...
            +(power(ev1{n}'*Jx*ev1{4},2) + power(ev1{n}'*Jx*ev2{4},2) )./(Ev(7) - E0));
        
        n = 2;
        az0{n} = -power(muB*gJ,2).*((power(ev1{n}'*Jz*ev1{1},2) + power(ev1{n}'*Jz*ev2{1},2) )./(Ev(1) - Ev(3))...
            +(power(ev1{n}'*Jz*ev1{3},2) + power(ev1{n}'*Jz*ev2{3},2) )./(Ev(5) - Ev(3))...
            +(power(ev1{n}'*Jz*ev1{4},2) + power(ev1{n}'*Jz*ev2{4},2) )./(Ev(7) - Ev(3)));
        %
        ax0{n} = -power(muB*gJ,2).*((power(ev1{n}'*Jx*ev1{1},1) + power(ev1{n}'*Jx*ev2{1},2) )./(Ev(1) - Ev(3))...
            +(power(ev1{n}'*Jx*ev1{3},2) + power(ev1{n}'*Jx*ev2{3},2) )./(Ev(5) - Ev(3))...
            +(power(ev1{n}'*Jx*ev1{4},2) + power(ev1{n}'*Jx*ev2{4},2) )./(Ev(7) - Ev(3)));
        
        n = 3;
        az0{n} = -power(muB*gJ,2).*((power(ev1{n}'*Jz*ev1{1},2) + power(ev1{n}'*Jz*ev2{1},2) )./(Ev(1) - Ev(5))...
            +(power(ev1{n}'*Jz*ev1{2},2) + power(ev1{n}'*Jz*ev2{2},2) )./(Ev(3) - Ev(5))...
            +(power(ev1{n}'*Jz*ev1{4},2) + power(ev1{n}'*Jz*ev2{4},2) )./(Ev(7) - Ev(5)));
        %
        ax0{n} = -power(muB*gJ,2).*((power(ev1{n}'*Jx*ev1{1},1) + power(ev1{n}'*Jx*ev2{1},2) )./(Ev(1) - Ev(5))...
            +(power(ev1{n}'*Jx*ev1{2},2) + power(ev1{n}'*Jx*ev2{2},2) )./(Ev(3) - Ev(5))...
            +(power(ev1{n}'*Jx*ev1{4},2) + power(ev1{n}'*Jx*ev2{4},2) )./(Ev(7) - Ev(5)));
        
        n = 4;
        az0{n} = -power(muB*gJ,2).*((power(ev1{n}'*Jz*ev1{1},2) + power(ev1{n}'*Jz*ev2{1},2) )./(Ev(1) - Ev(7))...
            +(power(ev1{n}'*Jz*ev1{2},2) + power(ev1{n}'*Jz*ev2{2},2) )./(Ev(3) - Ev(7))...
            +(power(ev1{n}'*Jz*ev1{3},2) + power(ev1{n}'*Jz*ev2{3},2) )./(Ev(5) - Ev(7)));
        %
        ax0{n} = -power(muB*gJ,2).*((power(ev1{n}'*Jx*ev1{1},1) + power(ev1{n}'*Jx*ev2{1},2) )./(Ev(1) - Ev(7))...
            +(power(ev1{n}'*Jx*ev1{2},2) + power(ev1{n}'*Jx*ev2{2},2) )./(Ev(3) - Ev(7))...
            +(power(ev1{n}'*Jx*ev1{3},2) + power(ev1{n}'*Jx*ev2{3},2) )./(Ev(5) - Ev(7)));
        output.ax = ax0;
        output.az = az0;
        output.gxx = gxx;
        output.gzz = gzz;
        output.D10 = Ev(3)-E0;
        output.D20 = Ev(5)-E0;
        output.D30 = Ev(7)-E0;
        
    case 'Perturbative expression'
        
        colz = {'k','r','b','c'};
        
        [P,D] = eig(HCEF + Jz*1e-10);
        
        En = real([D(1,1) D(2,2) D(3,3) D(4,4) D(5,5) D(6,6) D(7,7) D(8,8)]);
        Ev = sort(En);
        E0 = Ev(1);
        
        
        H = 0:60;
        phi = pi/7; display(phi); Jp = Jx*cos(phi)+Jy*sin(phi);
        for n = 1:4
            %display([2*(n-1) + 1, 2*n]);
            ind1 = find(En==Ev(2*(n-1) + 1));
            ind2 = find(En==Ev(2*n));
            
            ev1{n} = P(:,ind1(1));
            ev2{n} = P(:,ind2(end));
            
            
            
            
            gzz{n} = abs(2*gJ*(ev1{n}'*Jz*ev1{n}));
            gxx{n} = abs(2*gJ*(ev1{n}'*Jp*ev2{n}));
            
            %display(['gxx' num2str(n-1) ' = ' num2str(gxx)]);
            %display(['gzz' num2str(n-1) ' = ' num2str(gzz)]);
            
        end
        
        n = 1;
        az0{n} = -power(muB*gJ,2).*((power(ev1{n}'*Jz*ev1{2},2) + power(ev1{n}'*Jz*ev2{2},2) )./(Ev(3) - E0)...
            +(power(ev1{n}'*Jz*ev1{3},2) + power(ev1{n}'*Jz*ev2{3},2) )./(Ev(5) - E0)...
            +(power(ev1{n}'*Jz*ev1{4},2) + power(ev1{n}'*Jz*ev2{4},2) )./(Ev(7) - E0));
        %
        ax0{n} = -power(muB*gJ,2).*((power(ev1{n}'*Jp*ev1{2},2) + power(ev1{n}'*Jp*ev2{2},2) )./(Ev(3) - E0)...
            +(power(ev1{n}'*Jp*ev1{3},2) + power(ev1{n}'*Jp*ev2{3},2) )./(Ev(5) - E0)...
            +(power(ev1{n}'*Jp*ev1{4},2) + power(ev1{n}'*Jp*ev2{4},2) )./(Ev(7) - E0));
        
        n = 2;
        az0{n} = -power(muB*gJ,2).*((power(ev1{n}'*Jz*ev1{1},2) + power(ev1{n}'*Jz*ev2{1},2) )./(Ev(1) - Ev(3))...
            +(power(ev1{n}'*Jz*ev1{3},2) + power(ev1{n}'*Jz*ev2{3},2) )./(Ev(5) - Ev(3))...
            +(power(ev1{n}'*Jz*ev1{4},2) + power(ev1{n}'*Jz*ev2{4},2) )./(Ev(7) - Ev(3)));
        %
        ax0{n} = -power(muB*gJ,2).*((power(ev1{n}'*Jp*ev1{1},2) + power(ev1{n}'*Jp*ev2{1},2) )./(Ev(1) - Ev(3))...
            +(power(ev1{n}'*Jp*ev1{3},2) + power(ev1{n}'*Jp*ev2{3},2) )./(Ev(5) - Ev(3))...
            +(power(ev1{n}'*Jp*ev1{4},2) + power(ev1{n}'*Jp*ev2{4},2) )./(Ev(7) - Ev(3)));
        
        n = 3;
        az0{n} = -power(muB*gJ,2).*((power(ev1{n}'*Jz*ev1{1},2) + power(ev1{n}'*Jz*ev2{1},2) )./(Ev(1) - Ev(5))...
            +(power(ev1{n}'*Jz*ev1{2},2) + power(ev1{n}'*Jz*ev2{2},2) )./(Ev(3) - Ev(5))...
            +(power(ev1{n}'*Jz*ev1{4},2) + power(ev1{n}'*Jz*ev2{4},2) )./(Ev(7) - Ev(5)));
        %
        ax0{n} = -power(muB*gJ,2).*((power(ev1{n}'*Jp*ev1{1},2) + power(ev1{n}'*Jp*ev2{1},2) )./(Ev(1) - Ev(5))...
            +(power(ev1{n}'*Jp*ev1{2},2) + power(ev1{n}'*Jp*ev2{2},2) )./(Ev(3) - Ev(5))...
            +(power(ev1{n}'*Jp*ev1{4},2) + power(ev1{n}'*Jp*ev2{4},2) )./(Ev(7) - Ev(5)));
        
        n = 4;
        az0{n} = -power(muB*gJ,2).*((power(ev1{n}'*Jz*ev1{1},2) + power(ev1{n}'*Jz*ev2{1},2) )./(Ev(1) - Ev(7))...
            +(power(ev1{n}'*Jz*ev1{2},2) + power(ev1{n}'*Jz*ev2{2},2) )./(Ev(3) - Ev(7))...
            +(power(ev1{n}'*Jz*ev1{3},2) + power(ev1{n}'*Jz*ev2{3},2) )./(Ev(5) - Ev(7)));
        %
        ax0{n} = -power(muB*gJ,2).*((power(ev1{n}'*Jp*ev1{1},2) + power(ev1{n}'*Jp*ev2{1},2) )./(Ev(1) - Ev(7))...
            +(power(ev1{n}'*Jp*ev1{2},2) + power(ev1{n}'*Jp*ev2{2},2) )./(Ev(3) - Ev(7))...
            +(power(ev1{n}'*Jp*ev1{3},2) + power(ev1{n}'*Jp*ev2{3},2) )./(Ev(5) - Ev(7)));
        
        
        for n = 1:4
            
            subplot(1,2,1); hold on;
            plot(H, Ev(2*n) - E0 + muB*H*gxx{n}/2 + ax0{n}*power(H,2), '--' ,'color',colz{n})
            plot(H, Ev(2*n) - E0 - muB*H*gxx{n}/2 + ax0{n}*power(H,2), '--' ,'color',colz{n})
            
            subplot(1,2,2); hold on;
            plot(H, Ev(2*n) - E0 + muB*H*gzz{n}/2 + az0{n}*power(H,2), '--' ,'color',colz{n})
            plot(H, Ev(2*n) - E0 - muB*H*gzz{n}/2 + az0{n}*power(H,2), '--' ,'color',colz{n})
            
        end
        
        output.D = {Ev(3) - E0, Ev(5) - E0, Ev(7) - E0};
        output.ax0 = ax0;
        output.az0 = az0;
        output.gxx = gxx;
        output.gzz = gzz;
        
        
        %display([num2str(gxx{1}) ', ' num2str(gzz{1}) ', ' num2str(ax0{1}) ', ' num2str(az0{1})])
        
        %display(Ev(3) - E0);% gxx{2}, gzz{2}, ax0{2}, az0{2}])
        
    case 'get quartic term'
        
        colz = {'k','r','b','c'};
        
        [P,D] = eig(HCEF + Jz*1e-10);
        
        En = real([D(1,1) D(2,2) D(3,3) D(4,4) D(5,5) D(6,6) D(7,7) D(8,8)]);
        Ev = sort(En);
        E0 = Ev(1);
        
        
        H = 0:60;
        
        for n = 1:4
            %display([2*(n-1) + 1, 2*n]);
            ind1 = find(En==Ev(2*(n-1) + 1));
            ind2 = find(En==Ev(2*n));
            
            ev1{n} = P(:,ind1(1));
            ev2{n} = P(:,ind2(end));
            
            gzz{n} = abs(2*gJ*(ev1{n}'*Jz*ev1{n}));
            gxx{n} = abs(2*gJ*(ev1{n}'*Jx*ev2{n}));
            
            %display(['gxx' num2str(n-1) ' = ' num2str(gxx)]);
            %display(['gzz' num2str(n-1) ' = ' num2str(gzz)]);
            
        end
        
        n = 1;
        az0{n} = -power(muB*gJ,2).*((power(ev1{n}'*Jz*ev1{2},2) + power(ev1{n}'*Jz*ev2{2},2) )./(Ev(3) - E0)...
            +(power(ev1{n}'*Jz*ev1{3},2) + power(ev1{n}'*Jz*ev2{3},2) )./(Ev(5) - E0)...
            +(power(ev1{n}'*Jz*ev1{4},2) + power(ev1{n}'*Jz*ev2{4},2) )./(Ev(7) - E0));
        %
        ax0{n} = -power(muB*gJ,2).*((power(ev1{n}'*Jx*ev1{2},2) + power(ev1{n}'*Jx*ev2{2},2) )./(Ev(3) - E0)...
            +(power(ev1{n}'*Jx*ev1{3},2) + power(ev1{n}'*Jx*ev2{3},2) )./(Ev(5) - E0)...
            +(power(ev1{n}'*Jx*ev1{4},2) + power(ev1{n}'*Jx*ev2{4},2) )./(Ev(7) - E0));
        
        n = 2;
        az0{n} = -power(muB*gJ,2).*((power(ev1{n}'*Jz*ev1{1},2) + power(ev1{n}'*Jz*ev2{1},2) )./(Ev(1) - Ev(3))...
            +(power(ev1{n}'*Jz*ev1{3},2) + power(ev1{n}'*Jz*ev2{3},2) )./(Ev(5) - Ev(3))...
            +(power(ev1{n}'*Jz*ev1{4},2) + power(ev1{n}'*Jz*ev2{4},2) )./(Ev(7) - Ev(3)));
        %
        ax0{n} = -power(muB*gJ,2).*((power(ev1{n}'*Jx*ev1{1},1) + power(ev1{n}'*Jx*ev2{1},2) )./(Ev(1) - Ev(3))...
            +(power(ev1{n}'*Jx*ev1{3},2) + power(ev1{n}'*Jx*ev2{3},2) )./(Ev(5) - Ev(3))...
            +(power(ev1{n}'*Jx*ev1{4},2) + power(ev1{n}'*Jx*ev2{4},2) )./(Ev(7) - Ev(3)));
        
        n = 3;
        az0{n} = -power(muB*gJ,2).*((power(ev1{n}'*Jz*ev1{1},2) + power(ev1{n}'*Jz*ev2{1},2) )./(Ev(1) - Ev(5))...
            +(power(ev1{n}'*Jz*ev1{2},2) + power(ev1{n}'*Jz*ev2{2},2) )./(Ev(3) - Ev(5))...
            +(power(ev1{n}'*Jz*ev1{4},2) + power(ev1{n}'*Jz*ev2{4},2) )./(Ev(7) - Ev(5)));
        %
        ax0{n} = -power(muB*gJ,2).*((power(ev1{n}'*Jx*ev1{1},1) + power(ev1{n}'*Jx*ev2{1},2) )./(Ev(1) - Ev(5))...
            +(power(ev1{n}'*Jx*ev1{2},2) + power(ev1{n}'*Jx*ev2{2},2) )./(Ev(3) - Ev(5))...
            +(power(ev1{n}'*Jx*ev1{4},2) + power(ev1{n}'*Jx*ev2{4},2) )./(Ev(7) - Ev(5)));
        
        n = 4;
        az0{n} = -power(muB*gJ,2).*((power(ev1{n}'*Jz*ev1{1},2) + power(ev1{n}'*Jz*ev2{1},2) )./(Ev(1) - Ev(7))...
            +(power(ev1{n}'*Jz*ev1{2},2) + power(ev1{n}'*Jz*ev2{2},2) )./(Ev(3) - Ev(7))...
            +(power(ev1{n}'*Jz*ev1{3},2) + power(ev1{n}'*Jz*ev2{3},2) )./(Ev(5) - Ev(7)));
        %
        ax0{n} = -power(muB*gJ,2).*((power(ev1{n}'*Jx*ev1{1},1) + power(ev1{n}'*Jx*ev2{1},2) )./(Ev(1) - Ev(7))...
            +(power(ev1{n}'*Jx*ev1{2},2) + power(ev1{n}'*Jx*ev2{2},2) )./(Ev(3) - Ev(7))...
            +(power(ev1{n}'*Jx*ev1{3},2) + power(ev1{n}'*Jx*ev2{3},2) )./(Ev(5) - Ev(7)));
        
        h1 = figure; subplot(1,2,1); hold on;
        xlabel('\mu_0H [T]'); ylabel('\omega [meV]'); grid on; box on;
        title(['H||ab,   g_{\perp} = ' num2str(gxx{1})])
        
        subplot(1,2,2); hold on;
        xlabel('\mu_0H [T]'); ylabel('\omega [meV]'); grid on; box on;
        title(['H||c,   g_{||} = ' num2str(gzz{1})])
        
        n = 1;
        for h = H
            HZeeman = -gJ*muB*h*Jx;
            [P,D] = eig(HCEF+HZeeman);
            
            for m = 1:8
                Eigenval(m) = real(D(m,m));
            end
            Eigenval = sort(Eigenval);
            for m = 1:8
                E1{m}(n) = Eigenval(m);
            end
            
            HZeeman = -gJ*muB*h*Jz;
            [P,D] = eig(HCEF+HZeeman);
            
            for m = 1:8
                Eigenval(m) = real(D(m,m));
            end
            Eigenval = sort(Eigenval);
            for m = 1:8
                E2{m}(n) = Eigenval(m);
            end
            
            n = n + 1;
        end
        
        colz = {'k','k','r','r','b','b','c','c'};
        for m = 1:8
            subplot(1,2,1);
            plot(H,E1{m}-E0,'color',colz{m})
            
            subplot(1,2,2);
            plot(H,E2{m}-E0,'color',colz{m})
        end
        
        colz = {'k','r','b','c'};
        
        
        for n = 1:4
            subplot(1,2,1);
            plot(H, Ev(2*n) - E0 + muB*H*gxx{n}/2 + ax0{n}*power(H,2), '--' ,'color',colz{n})
            plot(H, Ev(2*n) - E0 - muB*H*gxx{n}/2 + ax0{n}*power(H,2), '--' ,'color',colz{n})
            
            subplot(1,2,2);
            plot(H, Ev(2*n) - E0 + muB*H*gzz{n}/2 + az0{n}*power(H,2), '--' ,'color',colz{n})
            plot(H, Ev(2*n) - E0 - muB*H*gzz{n}/2 + az0{n}*power(H,2), '--' ,'color',colz{n})
            
        end
        
        plot2 = 0;
        if plot2
            h2 = figure; hold on;
        end
        for n = 1:4
            
            Y1p{n} = Ev(2*n) + muB*H*gxx{n}/2 + ax0{n}*power(H,2) - E1{2*n};
            Y1m{n} = Ev(2*n) - muB*H*gxx{n}/2 + ax0{n}*power(H,2) - E1{2*n-1};
            %plot(H, abs(Y1p), '-' ,'color',colz{n})
            %plot(H, abs(Y1m), '-' ,'color',colz{n})
            if plot2
                subplot(1,2,1); hold on;
                plot(H, abs(Y1p{n} - Y1m{n})./2 , '-' ,'color',colz{n})
                plot(H, abs(Y1p{n} + Y1m{n})./2, '--' ,'color',colz{n})
            end
            
            Y2p{n} = Ev(2*n) + muB*H*gzz{n}/2 + az0{n}*power(H,2) - E2{2*n};
            Y2m{n} = Ev(2*n) - muB*H*gzz{n}/2 + az0{n}*power(H,2) - E2{2*n-1};
            %plot(H, abs(Y2p), '-' ,'color',colz{n})
            %plot(H, abs(Y2m), '-' ,'color',colz{n})
            if plot2
                subplot(1,2,2); hold on;
                plot(H, abs(Y2p{n} - Y2m{n})./2, '-' ,'color',colz{n})
                plot(H, abs(Y2p{n} + Y2m{n})./2, '--' ,'color',colz{n})
            end
        end
        
        if plot2
            subplot(1,2,1); set(gca,'xscale','log','yscale','log');grid on; box on;
            subplot(1,2,2); set(gca,'xscale','log','yscale','log');grid on; box on;
        end
        
        ft1 = fittype('b*power(x,3)', 'independent', 'x', 'dependent', 'y');
        ft2 = fittype('c*power(x,4)', 'independent', 'x', 'dependent', 'y');
        
        opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
        opts.Display = 'Off';
        opts.StartPoint = [0]; ind = 1:60;
        [fitresult{1}, gof(1)] = fit( H(ind)', (Y1p{1}(ind) - Y1m{1}(ind))'/2, ft1, opts );
        [fitresult{2}, gof(1)] = fit( H(ind)', (Y1p{1}(ind) + Y1m{1}(ind))'/2, ft2, opts );
        
        [fitresult{3}, gof(1)] = fit( H(ind)', (Y2p{1}(ind) - Y2m{1}(ind))'/2, ft1, opts );
        [fitresult{4}, gof(1)] = fit( H(ind)', (Y2p{1}(ind) + Y2m{1}(ind))'/2, ft2, opts );
        
        
        bx = -fitresult{1}.b;
        cx = -fitresult{2}.c;
        
        bz = -fitresult{3}.b;
        cz = -fitresult{4}.c;
        
        if plot2
            subplot(1,2,1)
            plot(H,abs(bx).*power(H,3),'y')
            plot(H,abs(cx).*power(H,4),'y')
            
            subplot(1,2,2)
            plot(H,abs(bz).*power(H,3),'y')
            plot(H,abs(cz).*power(H,4),'y')
        end
        
        figure(h1);
        
        subplot(1,2,1); n = 1;
        plot(H, Ev(2*n) - E0 + muB*H*gxx{n}/2 + ax0{n}*power(H,2) + bx*power(H,3) + cx*power(H,4) , '--' ,'color','y')
        plot(H, Ev(2*n) - E0 - muB*H*gxx{n}/2 + ax0{n}*power(H,2) - bx*power(H,3) + cx*power(H,4), '--' ,'color','y')
        
        subplot(1,2,2);
        plot(H, Ev(2*n) - E0 + muB*H*gzz{n}/2 + az0{n}*power(H,2) + bz*power(H,3) + cz*power(H,4), '--' ,'color','y')
        plot(H, Ev(2*n) - E0 - muB*H*gzz{n}/2 + az0{n}*power(H,2) - bz*power(H,3) + cz*power(H,4), '--' ,'color','y')
        
        display(['bx0, bz0, cx0, cz0 = ' num2str(bx) ', ' num2str(bz) ', ' num2str(cx) ', ' num2str(cz)])
              
    case 'arb angle'
        th = pi/4;
        
        [P,D] = eig(HCEF + Jz*1e-10);
        
        En = real([D(1,1) D(2,2) D(3,3) D(4,4) D(5,5) D(6,6) D(7,7) D(8,8)]);
        Ev = sort(En);
        E0 = Ev(1);
        H = 0:60;
        
        for n = 1:4
            ind1 = find(En==Ev(2*(n-1) + 1));
            ind2 = find(En==Ev(2*n));
            
            ev1{n} = P(:,ind1(1));
            ev2{n} = P(:,ind2(end));
            
            gzz{n} = abs(2*gJ*(ev1{n}'*Jz*ev1{n}));
            gxx{n} = abs(2*gJ*(ev1{n}'*Jx*ev2{n}));   
        end
       
        n = 1;
        a{n} =  @(th) -power(muB*gJ,2).*(...
             (power(ev1{n}'*(Jx*sin(th) + Jz*cos(th))*ev1{2},2) + power(ev1{n}'*(Jx*sin(th) + Jz*cos(th))*ev2{2},2) )./(Ev(3) - E0)...
            +(power(ev1{n}'*(Jx*sin(th) + Jz*cos(th))*ev1{3},2) + power(ev1{n}'*(Jx*sin(th) + Jz*cos(th))*ev2{3},2) )./(Ev(5) - E0)...
            +(power(ev1{n}'*(Jx*sin(th) + Jz*cos(th))*ev1{4},2) + power(ev1{n}'*(Jx*sin(th) + Jz*cos(th))*ev2{4},2) )./(Ev(7) - E0));
        
        n = 2;
        a{n} =  @(th) -power(muB*gJ,2).*(...
             (power(ev1{n}'*(Jx*sin(th) + Jz*cos(th))*ev1{1},2) + power(ev1{n}'*(Jx*sin(th) + Jz*cos(th))*ev2{1},2) )./(Ev(1) - Ev(3))...
            +(power(ev1{n}'*(Jx*sin(th) + Jz*cos(th))*ev1{3},2) + power(ev1{n}'*(Jx*sin(th) + Jz*cos(th))*ev2{3},2) )./(Ev(5) - Ev(3))...
            +(power(ev1{n}'*(Jx*sin(th) + Jz*cos(th))*ev1{4},2) + power(ev1{n}'*(Jx*sin(th) + Jz*cos(th))*ev2{4},2) )./(Ev(7) - Ev(3)));
        
        n = 3;
        a{n} =  @(th) -power(muB*gJ,2).*(...
             (power(ev1{n}'*(Jx*sin(th) + Jz*cos(th))*ev1{1},2) + power(ev1{n}'*(Jx*sin(th) + Jz*cos(th))*ev2{1},2) )./(Ev(1) - Ev(5))...
            +(power(ev1{n}'*(Jx*sin(th) + Jz*cos(th))*ev1{2},2) + power(ev1{n}'*(Jx*sin(th) + Jz*cos(th))*ev2{2},2) )./(Ev(3) - Ev(5))...
            +(power(ev1{n}'*(Jx*sin(th) + Jz*cos(th))*ev1{4},2) + power(ev1{n}'*(Jx*sin(th) + Jz*cos(th))*ev2{4},2) )./(Ev(7) - Ev(5)));
        
        n = 4;
        a{n} =  @(th) -power(muB*gJ,2).*(...
             (power(ev1{n}'*(Jx*sin(th) + Jz*cos(th))*ev1{1},2) + power(ev1{n}'*(Jx*sin(th) + Jz*cos(th))*ev2{1},2) )./(Ev(1) - Ev(7))...
            +(power(ev1{n}'*(Jx*sin(th) + Jz*cos(th))*ev1{3},2) + power(ev1{n}'*(Jx*sin(th) + Jz*cos(th))*ev2{3},2) )./(Ev(5) - Ev(7))...
            +(power(ev1{n}'*(Jx*sin(th) + Jz*cos(th))*ev1{2},2) + power(ev1{n}'*(Jx*sin(th) + Jz*cos(th))*ev2{2},2) )./(Ev(3) - Ev(7)));
        
%         
%         for n = 1:4
%            ev{2*n-1} = ev1{n};
%            ev{2*n}   = ev2{n};
%         end
%         
%         function out = b_cubic(n,th) 
%                 nset = setdiff(1:8,[2*n-1,2*n]);
%                 
%                 out = 0;
%                 for m1 = nset 
%                     for m2 = nset
%                        out = out + (ev{n}'*(Jx*sin(th) + Jz*cos(th))*ev{m1} ).*(...
%                                     ev{m1}'*(Jx*sin(th) + Jz*cos(th))*ev{m2}).*(...
%                                     ev{m2}'*(Jx*sin(th) + Jz*cos(th))*ev{n} )./((Ev(2*n)-Ev(m1))*(Ev(2*n)-Ev(m2)));
%                     end
%                     
%                     
%                 end
%              
%         end
        
        
        n = 1;
        for h = H
            HZeeman = -gJ*muB*h*Jx*sin(th) -gJ*muB*h*Jz*cos(th) ;
            [P,D] = eig(HCEF+HZeeman);
            
            for m = 1:8
                Eigenval(m) = real(D(m,m));
            end
            Eigenval = sort(Eigenval);
            for m = 1:8
                E3{m}(n) = Eigenval(m);
            end
            
            n = n + 1;
        end
        
        figure; hold on;
        colz = {'k','k','r','r','b','b','c','c'};
        for m = 1:8
            plot(H,E3{m}-E0,'color',colz{m})
        end
        
        
        H1 = H*sin(th); H3 = H*cos(th);
        for n = 1:4
            plot(H, Ev(2*n) + (muB/2)*sqrt(power(gxx{n}*H1,2) + power(gzz{n}*H3,2)) + a{n}(th)*power(H,2)...
            - E0,'g--');
            plot(H, Ev(2*n) - (muB/2)*sqrt(power(gxx{n}*H1,2) + power(gzz{n}*H3,2)) + a{n}(th)*power(H,2)...
            - E0,'g--');
        end
        
        Th = linspace(-pi/2,pi/2,1000);
        figure; hold on;
        for n = 1:4
            m = 1;
            for th = Th
               aa(m) = a{n}(th); m = m+1; 
            end
            plot(Th, aa ,'color',colz{2*n})
        end
             
    case 'arb angle perturbative modelling'
        
        [P,D] = eig(HCEF + Jz*1e-10);
        
        En = real([D(1,1) D(2,2) D(3,3) D(4,4) D(5,5) D(6,6) D(7,7) D(8,8)]);
        Ev = sort(En);
        E0 = Ev(1);
        
        H = 0:60;
        
        for n = 1:4
            ind1 = find(En==Ev(2*(n-1) + 1));
            ind2 = find(En==Ev(2*n));
            
            ev1{n} = P(:,ind1(1));
            ev2{n} = P(:,ind2(end));
            
            gzz{n} = abs(2*gJ*(ev1{n}'*Jz*ev1{n}));
            gxx{n} = abs(2*gJ*(ev1{n}'*Jx*ev2{n}));
            
            
        end
        
        n = 1;
        az0{n} = -power(muB*gJ,2).*((power(ev1{n}'*Jz*ev1{2},2) + power(ev1{n}'*Jz*ev2{2},2) )./(Ev(3) - E0)...
            +(power(ev1{n}'*Jz*ev1{3},2) + power(ev1{n}'*Jz*ev2{3},2) )./(Ev(5) - E0)...
            +(power(ev1{n}'*Jz*ev1{4},2) + power(ev1{n}'*Jz*ev2{4},2) )./(Ev(7) - E0));
        %
        ax0{n} = -power(muB*gJ,2).*((power(ev1{n}'*Jx*ev1{2},2) + power(ev1{n}'*Jx*ev2{2},2) )./(Ev(3) - E0)...
            +(power(ev1{n}'*Jx*ev1{3},2) + power(ev1{n}'*Jx*ev2{3},2) )./(Ev(5) - E0)...
            +(power(ev1{n}'*Jx*ev1{4},2) + power(ev1{n}'*Jx*ev2{4},2) )./(Ev(7) - E0));
        
        n = 2;
        az0{n} = -power(muB*gJ,2).*((power(ev1{n}'*Jz*ev1{1},2) + power(ev1{n}'*Jz*ev2{1},2) )./(Ev(1) - Ev(3))...
            +(power(ev1{n}'*Jz*ev1{3},2) + power(ev1{n}'*Jz*ev2{3},2) )./(Ev(5) - Ev(3))...
            +(power(ev1{n}'*Jz*ev1{4},2) + power(ev1{n}'*Jz*ev2{4},2) )./(Ev(7) - Ev(3)));
        %
        ax0{n} = -power(muB*gJ,2).*((power(ev1{n}'*Jx*ev1{1},1) + power(ev1{n}'*Jx*ev2{1},2) )./(Ev(1) - Ev(3))...
            +(power(ev1{n}'*Jx*ev1{3},2) + power(ev1{n}'*Jx*ev2{3},2) )./(Ev(5) - Ev(3))...
            +(power(ev1{n}'*Jx*ev1{4},2) + power(ev1{n}'*Jx*ev2{4},2) )./(Ev(7) - Ev(3)));
        
        n = 3;
        az0{n} = -power(muB*gJ,2).*((power(ev1{n}'*Jz*ev1{1},2) + power(ev1{n}'*Jz*ev2{1},2) )./(Ev(1) - Ev(5))...
            +(power(ev1{n}'*Jz*ev1{2},2) + power(ev1{n}'*Jz*ev2{2},2) )./(Ev(3) - Ev(5))...
            +(power(ev1{n}'*Jz*ev1{4},2) + power(ev1{n}'*Jz*ev2{4},2) )./(Ev(7) - Ev(5)));
        %
        ax0{n} = -power(muB*gJ,2).*((power(ev1{n}'*Jx*ev1{1},1) + power(ev1{n}'*Jx*ev2{1},2) )./(Ev(1) - Ev(5))...
            +(power(ev1{n}'*Jx*ev1{2},2) + power(ev1{n}'*Jx*ev2{2},2) )./(Ev(3) - Ev(5))...
            +(power(ev1{n}'*Jx*ev1{4},2) + power(ev1{n}'*Jx*ev2{4},2) )./(Ev(7) - Ev(5)));
        
        n = 4;
        az0{n} = -power(muB*gJ,2).*((power(ev1{n}'*Jz*ev1{1},2) + power(ev1{n}'*Jz*ev2{1},2) )./(Ev(1) - Ev(7))...
            +(power(ev1{n}'*Jz*ev1{2},2) + power(ev1{n}'*Jz*ev2{2},2) )./(Ev(3) - Ev(7))...
            +(power(ev1{n}'*Jz*ev1{3},2) + power(ev1{n}'*Jz*ev2{3},2) )./(Ev(5) - Ev(7)));
        %
        ax0{n} = -power(muB*gJ,2).*((power(ev1{n}'*Jx*ev1{1},1) + power(ev1{n}'*Jx*ev2{1},2) )./(Ev(1) - Ev(7))...
            +(power(ev1{n}'*Jx*ev1{2},2) + power(ev1{n}'*Jx*ev2{2},2) )./(Ev(3) - Ev(7))...
            +(power(ev1{n}'*Jx*ev1{3},2) + power(ev1{n}'*Jx*ev2{3},2) )./(Ev(5) - Ev(7)));
        
        n = 1;
        for h = H
            HZeeman = -gJ*muB*h*Jx;
            [P,D] = eig(HCEF+HZeeman);
            
            for m = 1:8
                Eigenval(m) = real(D(m,m));
            end
            Eigenval = sort(Eigenval);
            for m = 1:8
                E1{m}(n) = Eigenval(m);
            end
            
            HZeeman = -gJ*muB*h*Jz;
            [P,D] = eig(HCEF+HZeeman);
            
            for m = 1:8
                Eigenval(m) = real(D(m,m));
            end
            Eigenval = sort(Eigenval);
            for m = 1:8
                E2{m}(n) = Eigenval(m);
            end
            
            n = n + 1;
        end
        
        
        for n = 1:4
            Y1p{n} = Ev(2*n) + muB*H*gxx{n}/2 + ax0{n}*power(H,2) - E1{2*n};
            Y1m{n} = Ev(2*n) - muB*H*gxx{n}/2 + ax0{n}*power(H,2) - E1{2*n-1};
            
            Y2p{n} = Ev(2*n) + muB*H*gzz{n}/2 + az0{n}*power(H,2) - E2{2*n};
            Y2m{n} = Ev(2*n) - muB*H*gzz{n}/2 + az0{n}*power(H,2) - E2{2*n-1};
            
        end
        
        
        
        ft1 = fittype('b*power(x,3)', 'independent', 'x', 'dependent', 'y');
        ft2 = fittype('c*power(x,4)', 'independent', 'x', 'dependent', 'y');
        
        opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
        opts.Display = 'Off';
        opts.StartPoint = [0]; ind = 1:60;
        
        for n = 1:4
            [fitresult1, gof] = fit( H(ind)', (Y1p{n}(ind) - Y1m{n}(ind))'/2, ft1, opts );
            [fitresult2, gof] = fit( H(ind)', (Y1p{n}(ind) + Y1m{n}(ind))'/2, ft2, opts );
            
            bx{n} =  fitresult1.b;
            cx{n} = -fitresult2.c;
        
            [fitresult1, gof] = fit( H(ind)', (Y2p{n}(ind) - Y2m{n}(ind))'/2, ft1, opts );
            [fitresult2, gof] = fit( H(ind)', (Y2p{n}(ind) + Y2m{n}(ind))'/2, ft2, opts );
            
            bz{n} =  fitresult1.b;
            cz{n} = -fitresult2.c;
        
        end

        th = 0;
        
        n = 1;
        for h = H
            HZeeman = -gJ*muB*h*Jx*sin(th) -gJ*muB*h*Jz*cos(th) ;
            [P,D] = eig(HCEF+HZeeman);
            
            for m = 1:8
                Eigenval(m) = real(D(m,m));
            end
            Eigenval = sort(Eigenval);
            for m = 1:8
                E3{m}(n) = Eigenval(m);
            end
            
            n = n + 1;
        end
        
        figure; hold on; grid on; box on;
        colz = {'k','k','r','r','b','b','c','c'};
        for m = 1:8
            plot(H,E3{m}-E0,'color',colz{m})
        end
        
        H1 = H*sin(th); H3 = H*cos(th);
        for n = 1:4
            
            plot(H, Ev(2*n) - E0 + (muB/2)*sqrt(power(gxx{n}*H1,2) + power(gzz{n}*H3,2))...
                + ax0{n}*power(H1,2) + az0{n}*power(H3,2)...
                , '--' ,'color','m')
            
            plot(H, Ev(2*n) - E0 -(muB/2)*sqrt(power(gxx{n}*H1,2) + power(gzz{n}*H3,2))...
                + ax0{n}*power(H1,2) + az0{n}*power(H3,2)...
                , '--' ,'color','m')
            
            plot(H, Ev(2*n) - E0 + (muB/2)*sqrt(power(gxx{n}*H1,2) + power(gzz{n}*H3,2))...
                + ax0{n}*power(H1,2) + az0{n}*power(H3,2)...
                + cx{n}*power(H1,4) +  cz{n}*power(H3,4)...
                ...- power(power(bx{n},2/3)*power(H1,2) + power(bz{n},2/3)*power(H3,2),3/2)...
                - bx{n}*power(H1,3) - bz{n}*power(H3,3) ...
                , '--' ,'color','g')
            
            plot(H, Ev(2*n) - E0 -(muB/2)*sqrt(power(gxx{n}*H1,2) + power(gzz{n}*H3,2))...
                + ax0{n}*power(H1,2) + az0{n}*power(H3,2)...
                + cx{n}*power(H1,4) +  cz{n}*power(H3,4)...
                ... + power(power(bx{n},2/3)*power(H1,2) + power(bz{n},2/3)*power(H3,2),3/2)...
                - bx{n}*power(H1,3) - bz{n}*power(H3,3)...
                , '--' ,'color','g')
            
            plot(H, Ev(2*n) - E0 + (muB/2)*sqrt(power(gxx{n}*H1,2) + power(gzz{n}*H3,2))...
                + ax0{n}*power(H1,2) + az0{n}*power(H3,2)...
                + cx{n}*power(H1,4) +  cz{n}*power(H3,4)...
                - power(power(bx{n},2/3)*power(H1,2) + power(bz{n},2/3)*power(H3,2),3/2)...
                , '--' ,'color','y')
            
            plot(H, Ev(2*n) - E0 -(muB/2)*sqrt(power(gxx{n}*H1,2) + power(gzz{n}*H3,2))...
                + ax0{n}*power(H1,2) + az0{n}*power(H3,2)...
                + cx{n}*power(H1,4) +  cz{n}*power(H3,4)...
                + power(power(bx{n},2/3)*power(H1,2) + power(bz{n},2/3)*power(H3,2),3/2)...
                , '--' ,'color','y')
            
        end
        
        
        %         plot(H, Ev(2*n) - E0 + muB*H*gzz{n}/2 + az0{n}*power(H,2) + bz*power(H,3) + cz*power(H,4), '--' ,'color','y')
        %         plot(H, Ev(2*n) - E0 - muB*H*gzz{n}/2 + az0{n}*power(H,2) - bz*power(H,3) + cz*power(H,4), '--' ,'color','y')
          
    case 'Transition MEz'
        
        [P,D] = eig(HCEF + Jz*1e-10);
        
        En = real([D(1,1) D(2,2) D(3,3) D(4,4) D(5,5) D(6,6) D(7,7) D(8,8)]);
        Ev = sort(En);
        
        for n = 1:4
            ind1 = find(En==Ev(2*(n-1) + 1));
            ind2 = find(En==Ev(2*n));
            
            ev1{n} = P(:,ind1);
            ev2{n} = P(:,ind2);
        end
        
        for n = 2:4
            MEz1 = abs(ev1{1}'*Jz*ev1{n});
            MEz2 = abs(ev1{1}'*Jz*ev2{n});
            MEz3 = abs(ev2{1}'*Jz*ev1{n});
            MEz4 = abs(ev2{1}'*Jz*ev2{n});
            
            display(MEz1);
            display(MEz2);
            display(MEz3);
            display(MEz4);
        end
        
    case 'Transition MEz2'
        
        [P,D] = eig(HCEF + Jz*1e-10);
        
        En = real([D(1,1) D(2,2) D(3,3) D(4,4) D(5,5) D(6,6) D(7,7) D(8,8)]);
        Ev = sort(En);
        
        for n = 1:4
            ind1 = find(En==Ev(2*(n-1) + 1));
            ind2 = find(En==Ev(2*n));
            
            ev1{n} = P(:,ind1);
            ev2{n} = P(:,ind2);
        end
        
        for n = 1:4
            MEz1 = abs(ev1{1}'*Jz*Jz*ev1{n});
            MEz2 = abs(ev1{1}'*Jz*Jz*ev2{n});
            MEz3 = abs(ev2{1}'*Jz*Jz*ev1{n});
            MEz4 = abs(ev2{1}'*Jz*Jz*ev2{n});
            
            display(MEz1);
            display(MEz2);
            display(MEz3);
            display(MEz4);
        end
        
    case 'Transition MEx'
        
        [P,D] = eig(HCEF + Jz*1e-10);
        
        En = real([D(1,1) D(2,2) D(3,3) D(4,4) D(5,5) D(6,6) D(7,7) D(8,8)]);
        Ev = sort(En);
        
        for n = 1:4
            ind1 = find(En==Ev(2*(n-1) + 1));
            ind2 = find(En==Ev(2*n));
            
            ev1{n} = P(:,ind1);
            ev2{n} = P(:,ind2);
        end
        
        for n = 2:4
            MEx1 = abs(ev1{1}'*Jx*ev1{n});
            MEx2 = abs(ev1{1}'*Jx*ev2{n});
            MEx3 = abs(ev2{1}'*Jx*ev1{n});
            MEx4 = abs(ev2{1}'*Jx*ev2{n});
            
            display(MEx1);
            display(MEx2);
            display(MEx3);
            display(MEx4);
        end
        
    case 'Transition MEx2'
        
        [P,D] = eig(HCEF + Jz*1e-10);
        
        En = real([D(1,1) D(2,2) D(3,3) D(4,4) D(5,5) D(6,6) D(7,7) D(8,8)]);
        Ev = sort(En);
        
        for n = 1:4
            ind1 = find(En==Ev(2*(n-1) + 1));
            ind2 = find(En==Ev(2*n));
            
            ev1{n} = P(:,ind1);
            ev2{n} = P(:,ind2);
        end
        
        for n = 2:4
            MEz1 = abs(ev1{1}'*Jx*Jx*ev1{n});
            MEz2 = abs(ev1{1}'*Jx*Jx*ev2{n});
            MEz3 = abs(ev2{1}'*Jx*Jx*ev1{n});
            MEz4 = abs(ev2{1}'*Jx*Jx*ev2{n});
            
            display(MEz1);
            display(MEz2);
            display(MEz3);
            display(MEz4);
        end
    
    case 'full spectrum ff ellipse ab plane'
        
        
        
        [P,D] = eig(HCEF + Jz*1e-10);
        
        En = real([D(1,1) D(2,2) D(3,3) D(4,4) D(5,5) D(6,6) D(7,7) D(8,8)...
            D(9,9) D(10,10) D(11,11) D(12,12) D(13,13) D(14,14) D(15,15) D(16,16)]);
        Ev = sort(En);
        
        E0 = Ev(1);
        
        ind1 = find(En==Ev(1));
        ind2 = find(En==Ev(2));
        
        ev1 = P(:,ind1);
        ev2 = P(:,ind2);
        
        
        
        
        
        gzz0 = 2*gJ*abs(ev1'*Jz*ev1);
        gxx0 = 2*gJ*abs(ev1'*Jx*ev2);
        
        %H = 0:1:20;
        
        [P,D] = eig(HCEF);
        E0 = min(min(D));
        
        n = 1; Phi = linspace(0,2*pi,1e3); h = 18;
        
        for ph = Phi
            HZeeman = -gJ*muB*h*(Jx*cos(ph) + Jy*sin(ph) ) ;
            [P,D] = eig(HCEF+HZeeman);
            
            for m = 1:16
                Eigenval(m) = real(D(m,m));
            end
            Eigenval = sort(Eigenval);
            for m = 1:16
                E1{m}(n) = Eigenval(m);
            end
            
            
            n = n + 1;
        end
        
        E0 = min( E1{1});
        
        figure; hold on
        colz = {'k','k','r','r','b','b','c','c','g','g','m','m','y','y','k','k'};
        for m = 1:16
            
            R = E1{m} - E0;
            X = R.*cos(Phi);
            Y = R.*sin(Phi);
            
            plot(X,Y,'color',colz{m})

        end
        
         R0 = 45;
         r = -R0:R0; for th = 0:pi/6:5*pi/6; plot(r*cos(th), r*sin(th),'--', 'color', [.7,.7,.7], 'displayname', ' '); end
         r = R0; th = linspace(0,2*pi,1e3); plot(r*cos(th), r*sin(th),'-', 'color', [.7,.7,.7], 'displayname', ' ');
         xlabel('X = E_k cos(\phi) [meV]')
         ylabel('Y = E_k sin(\phi) [meV]')
         box on; title(['H = ' num2str(h) ' T'])
    
    case 'full spectrum ff ellipse ac plane'
        
        
        
        [P,D] = eig(HCEF + Jz*1e-10);
        
        En = real([D(1,1) D(2,2) D(3,3) D(4,4) D(5,5) D(6,6) D(7,7) D(8,8)...
            D(9,9) D(10,10) D(11,11) D(12,12) D(13,13) D(14,14) D(15,15) D(16,16)]);
        Ev = sort(En);
        
        E0 = Ev(1);
        
        ind1 = find(En==Ev(1));
        ind2 = find(En==Ev(2));
        
        ev1 = P(:,ind1);
        ev2 = P(:,ind2);
        
        
        
        
        
        gzz0 = 2*gJ*abs(ev1'*Jz*ev1);
        gxx0 = 2*gJ*abs(ev1'*Jx*ev2);
        
        H = 0:1:20;
        
        [P,D] = eig(HCEF);
        E0 = min(min(D));
        
        n = 1; Theta = linspace(0,2*pi,1e3); 
        
        h = 18; ph = pi/6;
        
        for th = Theta
            HZeeman = -gJ*muB*h*(Jz*cos(th) + Jx*sin(th)*cos(ph) + Jy*sin(th)*sin(ph) ) ;
            [P,D] = eig(HCEF+HZeeman);
            
            for m = 1:16
                Eigenval(m) = real(D(m,m));
            end
            Eigenval = sort(Eigenval);
            for m = 1:16
                E1{m}(n) = Eigenval(m);
            end
            
            
            n = n + 1;
        end
        
        E0 = min( E1{1});
        
        figure; hold on
        colz = {'k','k','r','r','b','b','c','c','g','g','m','m','y','y','k','k'};
        for m = 1:16
            
            R = E1{m} - E0;
            X = R.*cos(Theta);
            Y = R.*sin(Theta);
            
            plot(X,Y,'color',colz{m})

        end
        
        R0 = 45;
         r = -R0:R0; for th = 0:pi/6:5*pi/6; plot(r*cos(th), r*sin(th),'--', 'color', [.7,.7,.7], 'displayname', ' '); end
         r = R0; th = linspace(0,2*pi,1e3); plot(r*cos(th), r*sin(th),'-', 'color', [.7,.7,.7], 'displayname', ' ');
         xlabel('X = E_k cos(\theta) [meV]')
         ylabel('Y = E_k sin(\theta) [meV]')
         box on; title(['H = ' num2str(h) ' T,  \phi = ' num2str(ph)])
           
    case 'full spectrum'
        
        
        
        [P,D] = eig(HCEF + Jz*1e-10);
        
        En = real([D(1,1) D(2,2) D(3,3) D(4,4) D(5,5) D(6,6) D(7,7) D(8,8)...
            D(9,9) D(10,10) D(11,11) D(12,12) D(13,13) D(14,14) D(15,15) D(16,16)]);
        Ev = sort(En);
        
        E0 = Ev(1);
        
        ind1 = find(En==Ev(1));
        ind2 = find(En==Ev(2));
        
        ev1 = P(:,ind1);
        ev2 = P(:,ind2);
        
        
        
        
        
        gzz0 = 2*gJ*abs(ev1'*Jz*ev1);
        gxx0 = 2*gJ*abs(ev1'*Jx*ev2);
        gyy0 = 2*gJ*abs(ev1'*Jy*ev2);
        
        H = 0:1:20;
        
        [P,D] = eig(HCEF);
        E0 = min(min(D));
        
        n = 1;
        for h = H
            HZeeman = -gJ*muB*h*Jx;
            [P,D] = eig(HCEF+HZeeman);
            
            for m = 1:16
                Eigenval(m) = real(D(m,m));
            end
            Eigenval = sort(Eigenval);
            for m = 1:16
                E1{m}(n) = Eigenval(m);
            end
            
            HZeeman = -gJ*muB*h*Jy;
            [P,D] = eig(HCEF+HZeeman);
            
            for m = 1:16
                Eigenval(m) = real(D(m,m));
            end
            Eigenval = sort(Eigenval);
            for m = 1:16
                E2{m}(n) = Eigenval(m);
            end
            
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
        
        figure; hold on
        colz = {'k','k','r','r','b','b','c','c','g','g','m','m','y','y','k','k'};
        for m = 1:16
            subplot(1,3,1); hold on;
            plot(H,E1{m}-E0,'color',colz{m})
            
            subplot(1,3,2); hold on;
            plot(H,E2{m}-E0,'color',colz{m})
            
            subplot(1,3,3); hold on;
            plot(H,E3{m}-E0,'color',colz{m})
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
        %         figure; hold on;
        %         plot(H,E1{2}-E1{1},'b');
        %         plot(H,E2{2}-E2{1},'r');
        %         legend({'H||ab', 'H||c'})
        %         xlim([0,20]);
        %         xlabel('\mu_0H [T]'); ylabel('\Delta_{0,+-}');
        
    case 'gap H-dept'
        
        
        
        [P,D] = eig(HCEF + Jz*1e-10);
        
        En = real([D(1,1) D(2,2) D(3,3) D(4,4) D(5,5) D(6,6) D(7,7) D(8,8)]);
        Ev = sort(En);
        
        E0 = Ev(1);
        
        ind1 = find(En==Ev(1));
        ind2 = find(En==Ev(2));
        
        ev1 = P(:,ind1);
        ev2 = P(:,ind2);
        
        
        
        
        
        gzz0 = 2*gJ*abs(ev1'*Jz*ev1);
        gxx0 = 2*gJ*abs(ev1'*Jx*ev2);
        
        H = 0:1:60;
        
        [P,D] = eig(HCEF);
        E0 = min(min(D));
        
        n = 1;
        for h = H
            HZeeman = -gJ*muB*h*Jx;
            [P,D] = eig(HCEF+HZeeman);
            
            for m = 1:8
                Eigenval(m) = real(D(m,m));
            end
            Eigenval = sort(Eigenval);
            for m = 1:8
                E1{m}(n) = Eigenval(m);
            end
            
            HZeeman = -gJ*muB*h*Jz;
            [P,D] = eig(HCEF+HZeeman);
            
            for m = 1:8
                Eigenval(m) = real(D(m,m));
            end
            Eigenval = sort(Eigenval);
            for m = 1:8
                E2{m}(n) = Eigenval(m);
            end
            
            n = n + 1;
        end
        
        
        colz = {'k','k','r','r','b','b','c','c'};
        
        for m = 1:4
            Em1{m} = (E1{2*(m-1)+1} + E1{2*m})./2;
            Em2{m} = (E2{2*(m-1)+1} + E2{2*m})./2;
            
        end
        
        
        for m = 1:4
            %subplot(1,2,1);
            Y1 = Em1{m};
            %plot(H,Y1,'color',colz{2*m})
            p = polyfit(H(1:3),Y1(1:3),2);
            %plot(H,polyval(p,H),'g--')
            display(p(1))
            
            %subplot(1,2,2);
            Y2 = Em2{m};
            %plot(H,Y2,'color',colz{2*m})
            p = polyfit(H(1:3),Y2(1:3),2);
            %plot(H,polyval(p,H),'g--')
            display(p(1))
        end
        
    case 'display gxx gzz'
        [P,D] = eig(HCEF);
        E0 = min(min(D));
        
        ind = find(min(D)==E0);
        
        ev1 = P(:,ind(1));
        ev2 = P(:,ind(2));
        
        gzz0 = 2*gJ*abs(ev1'*Jz*ev1);
        gxx0 = gJ*abs(ev1'*Jp*ev2);
        
        display(['g_{\perp} = ' num2str(gxx0) ',    g_{||} = ' num2str(gzz0) ])
        
    case 'ZF eigenvalues'
        [P,D] = eig(HCEF);
        display(D);
        
    case 'Field Dept H||c'
        
        H = 0:1:60;
        
        [~,D] = eig(HCEF);
        E0 = min(min(D));
        
        n = 1;
        for h = H
            HZeeman = -muB*h*Jz;
            [P,D] = eig(HCEF+HZeeman);
            
            for m = 1:8
                Eigenval(m) = D(m,m);
            end
            Eigenval = sort(Eigenval);
            for m = 1:8
                E{m}(n) = Eigenval(m);
            end
            
            n = n + 1;
        end
        
        figure; hold on
        colz = {'k','k','r','r','b','b','c','c'};
        for m = 1:8
            plot(H,E{m}-E0,'color',colz{m})
        end
        
        xlabel('\mu_0H [T]'); ylabel('\omega [meV]'); grid on; box on;
        title('H||c')
        
    case 'Field Dept H||ab'
        
        H = 0:1:60;
        
        [P,D] = eig(HCEF);
        E0 = min(min(D));
        
        n = 1;
        for h = H
            HZeeman = -muB*h*Jx;
            [P,D] = eig(HCEF+HZeeman);
            
            for m = 1:8
                Eigenval(m) = real(D(m,m));
            end
            Eigenval = sort(Eigenval);
            for m = 1:8
                E{m}(n) = Eigenval(m);
            end
            
            n = n + 1;
        end
        
        figure; hold on
        colz = {'k','k','r','r','b','b','c','c'};
        for m = 1:8
            plot(H,E{m}-E0,'color',colz{m})
        end
        
        xlabel('\mu_0H [T]'); ylabel('\omega [meV]'); grid on; box on;
        title('H||ab')
        
    case 'Chi v T'
        
        h = .1;
        T = 2:300;
        
        phi = pi/6;
        
        HZeeman_x = -muB*h*Jx;
        HZeeman_y = -muB*h*(cos(phi)*Jx + sin(phi)*Jy);
        HZeeman_z = -muB*h*Jz;
        
        % Calculate Partion Functions : expm (not exp) is crucial for proper
        % exponentiation of matrices.
        n = 1;
        for t = T
            b = 1./(kB*t);
            Z0(n) = trace(expm(-b*(HCEF)));
            Z1(n) = trace(expm(-b*(HCEF+HZeeman_x)));
            Z2(n) = trace(expm(-b*(HCEF+HZeeman_y)));
            Z3(n) = trace(expm(-b*(HCEF+HZeeman_z)));
            n = n + 1;
        end
        
        
        F0 = (kB*T).*log(Z0);
        F1 = (kB*T).*log(Z1);
        F2 = (kB*T).*log(Z2);
        F3 = (kB*T).*log(Z3);
        
        Xa = (F1-F0)/h/h;
        Xb = (F2-F0)/h/h;
        Xc = (F3-F0)/h/h;
        
        C = 0.0422*108.3/88.14;
        
        figure; hold on;
        plot(T,C./Xa,'r')
        plot(T,abs(C./Xb),'g')
        plot(T,C./Xc,'b')
        
        xlabel('T [K]'); ylabel('\chi^{-1} [mol Oe/emu]')
        grid on; box on;
        
    case 'Angle Dept. k'
        
        T = 4;
        N_T = 1;
        
        N_Th = 100; dth = pi/N_Th;
        Th = linspace(0 - 2*dth, pi + 2*dth, N_Th+4);
        
        
        H = 4:4:32;
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
                    HZeeman = -muB*hx*Jx-muB*hz*Jz;
                    
                    Z(n2,n3) = trace(expm(-b*(HCEF+HZeeman)));
                    
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
            
            plot(Th(3:end-2), k , 'color', colz(n2,:), 'displayname', num2str(H(n2)) );
        end
        
end



end