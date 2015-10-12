% Calculates the displacements for layered anisotropic materials using
% the Legendre and Laguerre polynomial approach
% Originally supplied by O. Bou-Matar, Lille
clear all

%% LAYUP, PLY PROPERTIES AND ORIENTATION
NomMat = {'Plexiglass'};
Phi = [0];                  % First Euler angle
Theta = [0];                % Second Euler angle
Psi = [0];                  % Third Euler angle
h = 8e-3;                   % thickness
psip = 0;                   % propagation direction

% NomMat = {'Alamsa' 'Alamsa' 'Alamsa' 'Alamsa' 'Alamsa' 'Alamsa' 'Alamsa' 'Alamsa' ...
%     'Alamsa' 'Alamsa' 'Alamsa' 'Alamsa'};            
% Phi = [0 -pi/4 0 -pi/4 0 -pi/4 -pi/4 0 -pi/4 0 -pi/4 0]; % First Euler angle
% Theta = [0 0 0 0 0 0 0 0 0 0 0 0];                       % second euler angle
% Psi = [0 0 0 0 0 0 0 0 0 0 0 0];                         % third euler angle
% h=0.23e-3;                              % Thickness of the ply in [m]
% psip = 0;                               % propagation direction with respect to the main in-plane coordinate axis 
               
%% COMPUTATIONAL PARAMETERS
freq = 50e3;                            % Frequency to be investigated
mMode = 3;                              % Index of the mode to be inspected
legDeg=10;                              % Legendre polynomial degree      
nPlies = size(NomMat,2);                % number of plies
Nodes = ones(nPlies,1)*legDeg;          % vector of degrees of polynomial expansions in different plies
H = ones(size(NomMat,2),1)*h;           % vector of ply thicknesses
w = 2*pi*freq;                          % Angular frequency to be inspected
nTot = 3*sum(Nodes);
nMax = 3*sum(Nodes);

%% NORMALIZATION PARAMETERS
Ca = 1e11;          % Pa = N/m^2
rhoa = 1e3;         % kg/m^3         
ka = w*sqrt(rhoa/Ca);
dr = 0.01;         % spacing in the normalized ply z coordinates
r = -1:dr:1;        % normalized z coordinates for legendre polynomials
PLeg=nan(length(Nodes(1)),length(r));
for ii=0:(Nodes(1)-1)
    Pn = legendre(ii,r);        % legendre associated functions evaluated at r
    PLeg(ii+1,:) = Pn(1,:);     % take only mth-order
end

%% PREALLOCATION OF THE VARIABLES FOR CALCULATION
Un = zeros(nPlies,length(r),3,nMax);
rho0=nan(nPlies,1);     % Density vector
F11=nan(3,3);
F12=nan(size(F11));
F21=nan(size(F11));
F22=nan(size(F11));
F33=nan(3,3,nPlies);
F31=nan(size(F11));
F13=nan(size(F11));
F32=nan(size(F11));
F23=nan(size(F11));
A1=nan(3,3,nPlies);
BB=nan(3,3,nPlies);
CC=nan(3,3,nPlies);
ABC=nan(3,3,nPlies);
A2=nan(3,3,nPlies);
F1 = zeros(nTot,nTot);
G1 = zeros(nTot,nTot);
H1 = zeros(nTot,nTot);

%% LOAD THE MATERIAL PROPERTIES
for ply = 1:nPlies
    Matp = LoadElasticConstants(fullfile('./Materials',strcat(NomMat{ply},'.dat')));
    C = RotateElasticConstants(Matp.C,Phi(ply),Theta(ply),Psi(ply)); % stiffness tensor rotated to principal axis (of anisotropy)
    C = C./Ca;                  % normalization of stiffness tensor 
    rho0(ply) = Matp.rho/rhoa;  % convert kg/m^3 to g/cm^3
    F11(1,1) = C(1,1);
    F11(1,2) = C(1,6);
    F11(1,3) = C(1,5);
    F11(2,1) = C(1,6);
    F11(2,2) = C(6,6);
    F11(2,3) = C(5,6);
    F11(3,1) = C(1,5);
    F11(3,2) = C(5,6);
    F11(3,3) = C(5,5);

    F12(1,1) = C(1,6);
    F12(1,2) = C(1,2);
    F12(1,3) = C(1,4);
    F12(2,1) = C(6,6);
    F12(2,2) = C(2,6);
    F12(2,3) = C(4,6);
    F12(3,1) = C(5,6);
    F12(3,2) = C(2,5);
    F12(3,3) = C(4,5);

    F22(1,1) = C(6,6);
    F22(1,2) = C(2,6);
    F22(1,3) = C(4,6);
    F22(2,1) = C(2,6);
    F22(2,2) = C(2,2);
    F22(2,3) = C(2,4);
    F22(3,1) = C(4,6);
    F22(3,2) = C(2,4);
    F22(3,3) = C(4,4);

    F33(1,1,ply) = C(5,5);
    F33(1,2,ply) = C(4,5);
    F33(1,3,ply) = C(3,5);
    F33(2,1,ply) = C(4,5);
    F33(2,2,ply) = C(4,4);
    F33(2,3,ply) = C(3,4);
    F33(3,1,ply) = C(3,5);
    F33(3,2,ply) = C(3,4);
    F33(3,3,ply) = C(3,3);

    F31(1,1) = C(1,5);
    F31(1,2) = C(5,6);
    F31(1,3) = C(5,5);
    F31(2,1) = C(1,4);
    F31(2,2) = C(4,6);
    F31(2,3) = C(4,5);
    F31(3,1) = C(1,3);
    F31(3,2) = C(3,6);
    F31(3,3) = C(3,5);

    F32(1,1) = C(5,6);
    F32(1,2) = C(2,5);
    F32(1,3) = C(4,5);
    F32(2,1) = C(4,6);
    F32(2,2) = C(2,4);
    F32(2,3) = C(4,4);
    F32(3,1) = C(3,6);
    F32(3,2) = C(2,3);
    F32(3,3) = C(3,4);

    F21 = F12';
    F13 = F31';
    F23 = F32';

    A1(:,:,ply) = cos(psip)^2*F11+cos(psip)*sin(psip)*(F12+F21)+sin(psip)^2*F22;
    BB(:,:,ply) = cos(psip)*(F13+F31)+sin(psip)*(F23+F32);
    CC(:,:,ply) = -F33(:,:,ply);
    ABC(:,:,ply) = cos(psip)*F31+sin(psip)*F32;
    A2(:,:,ply) = -rho0(ply)*eye(3);
end

%% PREPARATION OF THE LEGENDRE POLYNOMIAL (degree N(ply)) IN ALL PLIES
Ntemp = 1;
for ply = 1:nPlies
    As = zeros(3*(Nodes(ply)-2),3*Nodes(ply));
    Bs = zeros(3*(Nodes(ply)-2),3*Nodes(ply));
    Cs = zeros(3*(Nodes(ply)-2),3*Nodes(ply));
    for m=0:Nodes(ply)-3
        for n=0:Nodes(ply)-1
            Bs(3*m+1:3*(m+1),3*n+1:3*(n+1)) = 2/H(ply)*BB(:,:,ply)/ka*PmdPn(m,n); 
            Cs(3*m+1:3*(m+1),3*n+1:3*(n+1)) = 4/(H(ply)^2)*CC(:,:,ply)/ka^2*Pmd2Pn(m,n);
            if (m == n)
                As(3*m+1:3*(m+1),3*n+1:3*(n+1)) = 2*A1(:,:,ply)/(2*n+1);
                Cs(3*m+1:3*(m+1),3*n+1:3*(n+1)) = Cs(3*m+1:3*(m+1),3*n+1:3*(n+1))+2*A2(:,:,ply)/(2*n+1);
            end
        end
    end          
    H1(Ntemp:Ntemp+3*(Nodes(ply)-2)-1,3*sum(Nodes(1:ply-1))+1:3*sum(Nodes(1:ply))) = As;
    F1(Ntemp:Ntemp+3*(Nodes(ply)-2)-1,3*sum(Nodes(1:ply-1))+1:3*sum(Nodes(1:ply))) = Bs;
    G1(Ntemp:Ntemp+3*(Nodes(ply)-2)-1,3*sum(Nodes(1:ply-1))+1:3*sum(Nodes(1:ply))) = Cs;
    Ntemp = Ntemp+3*(Nodes(ply)-2);
end

%% BOUNDARY CONDITIONS - CONTINUITY ON THE INTERFACE ply and ply+1
for ply = 1:nPlies-1           
    % Continuity of displacemt between ply and ply+1
    Ds = zeros(3,3*Nodes(ply+1));
    Es = zeros(3,3*Nodes(ply));
    for n=0:Nodes(ply)-1
        Es(:,3*n+1:3*(n+1)) = eye(3);
    end
    for n=0:Nodes(ply+1)-1
        Ds(:,3*n+1:3*(n+1)) = -(-1)^n*eye(3);
    end       
    G1(Ntemp:Ntemp+2,3*sum(Nodes(1:ply-1))+1:3*sum(Nodes(1:ply))) = Es;
    G1(Ntemp:Ntemp+2,3*sum(Nodes(1:ply))+1:3*sum(Nodes(1:ply+1))) = Ds;
    Ntemp = Ntemp+3;     
    % Continuity of normal stress between ply and ply+1
    Ds = zeros(3,3*Nodes(ply));
    Es = zeros(3,3*Nodes(ply));
    Dp = zeros(3,3*Nodes(ply+1));
    Ep = zeros(3,3*Nodes(ply+1));
    for n=0:Nodes(ply)-1
        Ds(:,3*n+1:3*(n+1)) = -ABC(:,:,ply);
        Es(:,3*n+1:3*(n+1)) = 2/H(ply)*F33(:,:,ply)/ka*n*(n+1)/2;
    end
    for n=0:Nodes(ply+1)-1
        Dp(:,3*n+1:3*(n+1)) = ABC(:,:,ply+1)*(-1)^n;
        Ep(:,3*n+1:3*(n+1)) = -2/H(ply+1)*F33(:,:,ply+1)/ka*(-1)^(n+1)*n*(n+1)/2;
    end
    F1(Ntemp:Ntemp+2,3*sum(Nodes(1:ply-1))+1:3*sum(Nodes(1:ply))) = Ds;
    F1(Ntemp:Ntemp+2,3*sum(Nodes(1:ply))+1:3*sum(Nodes(1:ply+1))) = Dp;
    G1(Ntemp:Ntemp+2,3*sum(Nodes(1:ply-1))+1:3*sum(Nodes(1:ply))) = Es;
    G1(Ntemp:Ntemp+2,3*sum(Nodes(1:ply))+1:3*sum(Nodes(1:ply+1))) = Ep;
   
    Ntemp = Ntemp+3;       
end

%% BOUNDARY CONDITIONS AT THE LOWER INTERFACE
Ds = zeros(3,3*Nodes(1));
Es = zeros(3,3*Nodes(1));
for n=0:Nodes(1)-1
   Ds(:,3*n+1:3*(n+1)) = ABC(:,:,1)*(-1)^n;
   Es(:,3*n+1:3*(n+1)) = -2/H(1)*F33(:,:,1)/ka*(-1)^(n+1)*n*(n+1)/2;
end
F1(Ntemp:Ntemp+2,1:3*Nodes(1)) = Ds;
G1(Ntemp:Ntemp+2,1:3*Nodes(1)) = Es;
Ntemp = Ntemp+3;
      
%% BOUNDARY CONDITIONS AT THE UPPER INTERFACE
Dp = zeros(3,3*Nodes(nPlies));
Ep = zeros(3,3*Nodes(nPlies));
for n=0:Nodes(nPlies)-1
    Dp(:,3*n+1:3*(n+1)) = -ABC(:,:,nPlies);
    Ep(:,3*n+1:3*(n+1)) = 2/H(nPlies)*F33(:,:,nPlies)/ka*n*(n+1)/2;
end
F1(Ntemp:Ntemp+2,3*sum(Nodes(1:nPlies-1))+1:3*sum(Nodes(1:nPlies))) = Dp;
G1(Ntemp:Ntemp+2,3*sum(Nodes(1:nPlies-1))+1:3*sum(Nodes(1:nPlies))) = Ep;
    
%% CALCULATION OF THE SYSTEM EIGENVALUES
M1 = [F1 -eye(nTot);-H1 zeros(nTot)];
M2 = [G1 zeros(nTot);zeros(nTot) eye(nTot)];
[Z1,K] = eig(M1,M2);

kp = zeros(2*nTot,1);
for ii=1:2*nTot
    if (K(ii,ii) == 0)
        kp(ii) = NaN;
    else
        kp(ii) = 1i/K(ii,ii)*ka;
    end
end

for ii=1:2*nTot
    if (real(kp(ii)) == 0)
        kp(ii) = NaN;
    else
        if (abs(imag(kp(ii)))/abs(real(kp(ii))) > 1e-8)
            kp(ii) = NaN;        
        end
    end
end

%% POST-PROCESSING OF THE RESULTS
% interm = zeros(1,2*Ntot);
[interm, Ind] = sort(kp);
k = interm(1:nMax);                 % complex wavenumber 
ZUn = Z1(1:nTot,Ind(1:nMax));       % coefficients of the Legednre ploynamials 

%% RECONSTRUNCTION OF THE U VECTOR USING LEGENDRE POLYNOMIALS
for kk=1:nMax
    for ll=1:3
        for ply=1:nPlies
            Un(ply,:,ll,kk) = zeros(1,length(r));
            for jj=1:Nodes(ply)
                Un(ply,:,ll,kk) = Un(ply,:,ll,kk) + ZUn(ll+(jj-1)*3 + sum(Nodes(1:ply-1))*3,kk)*PLeg(jj,:);
            end
        end
    end
end

%% VELOCITY CALCULATIONS
wavenumber = abs(real(k));          % wavenumber
vit = 2*pi*freq./wavenumber;%*sqrt(max(rho0)*rhoa/(max(max(max(C)))*Ca))
l = 1;
Velocity(l) = vit(l);
for i = 2:length(vit)
   if (vit(i) < 0.99999*vit(i-1) && isfinite(vit(i)))
    l = l+1;
    Velocity(l) = vit(i);
    Displacement(:,:,:,l) = Un(:,:,:,i);
   end
end
fprintf('\nPlotting mode with velocity\t v = %.2f m/s\n',Velocity(mMode));

%% VISUALIZATION
figure
ht = 0;
hTot = sum(H(1:nPlies));
for j=1:nPlies
    ht = H(j)+ht;
    zm = ht - H(j)/2;
    z = H(j)*r/2+zm;
    plot(abs(Displacement(j,:,1,mMode)),(hTot-z)*1e3,'b','LineWidth',2)
    hold on
    plot(abs(Displacement(j,:,2,mMode)),(hTot-z)*1e3,'r--','LineWidth',2)
    plot(abs(Displacement(j,:,3,mMode)),(hTot-z)*1e3,'g-.','LineWidth',2)
end
xlabel('Displacement [-]');
ylabel('Z-position [mm]');
legend('u_x','u_y','u_z')