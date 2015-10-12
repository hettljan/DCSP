% Calculates the displacements for layered anisotropic materials using
% the Legendre and Laguerre polynomial approach
% Originally supplied by O. Bou-Matar, Lille
clear all

%% LAYUP, PLY PROPERTIES AND ORIENTATION
NomMat = {'Alamsa' 'Alamsa' 'Alamsa' 'Alamsa' 'Alamsa' 'Alamsa' 'Alamsa' 'Alamsa' ...
    'Alamsa' 'Alamsa' 'Alamsa' 'Alamsa'};            
Phi = [0 -pi/4 0 -pi/4 0 -pi/4 -pi/4 0 -pi/4 0 -pi/4 0]; % First Euler angle
Theta = [0 0 0 0 0 0 0 0 0 0 0 0];                       % second euler angle
Psi = [0 0 0 0 0 0 0 0 0 0 0 0];                         % third euler angle
h=0.23e-3;                              % Thickness of the ply in [m]
psip = 0;                               % propagation direction with respect to the main in-plane coordinate axis 
freq = 300e3;                           % Frequency to be investigated
mMode = 2;                              % ?
Npts = 100;                             % ?
legDeg=10;                              % Legendre polynomial degree                      

%% ADDITIONAL VARIABLES
Nodes = ones(nPlies,1)*legDeg;          % vector of degrees of polynomial expansions in different plies
nPlies = size(NomMat,2);                % number of plies
H = ones(size(NomMat,2),1)*h;           % vector of ply thicknesses
w = 2*pi*freq;                          % Angular frequency to be inspected
nTot = 3*sum(Nodes);
nMax = 3*sum(Nodes);

%% NORMALIZATION PARAMETERS
Ca = 1e11;          % Pa = N/m^2
rhoa = 1e3;         % kg/m^3
% htot = sum(H);                
ka = w*sqrt(rhoa/Ca);
dr = 0.001;         % spacing in the normalized ply z coordinates
r = -1:dr:1;        % normalized z coordinates for legendre polynomials
PLeg=nan(Nodes(1));
for ii=0:(Nodes(1)-1)
    Pn = legendre(ii,r);        % legendre associated functions evaluated at r
    PLeg(ii+1,:) = Pn(1,:);     % take only mth-order
end

%% LOAD THE MATERIAL PROPERTIES
Un = zeros(nPlies,length(r),3,nMax);
for ply = 1:nPlies
 Matp = LoadElasticConstants(fullfile('./Materials',strcat(NomMat{ply},'.dat')));
    C = RotateElasticConstants(Matp.C,Phi(ply),Theta(ply),Psi(ply)); % stiffness tensor rotated to principal axis (of anisotropy)
    C = C./Ca;                  % normalization of stiffness tensor 
    rho0(ply) = Matp.rho/rhoa;  % convert kg/m^3 to g/cm^3
    C(1,1,ii) = Cm(1,1)/Ca;
    C(1,2,ii) = Cm(1,2)/Ca;
    C(1,3,ii) = Cm(1,3)/Ca;
    C(1,4,ii) = Cm(1,4)/Ca;
    C(1,5,ii) = Cm(1,5)/Ca;
    C(1,6,ii) = Cm(1,6)/Ca;
    C(2,2,ii) = Cm(2,2)/Ca;
    C(2,3,ii) = Cm(2,3)/Ca;
    C(2,4,ii) = Cm(2,4)/Ca;
    C(2,5,ii) = Cm(2,5)/Ca;
    C(2,6,ii) = Cm(2,6)/Ca;
    C(3,3,ii) = Cm(3,3)/Ca;
    C(3,4,ii) = Cm(3,4)/Ca;
    C(3,5,ii) = Cm(3,5)/Ca;
    C(3,6,ii) = Cm(3,6)/Ca;
    C(4,4,ii) = Cm(4,4)/Ca;
    C(4,5,ii) = Cm(4,5)/Ca;
    C(4,6,ii) = Cm(4,6)/Ca;
    C(5,5,ii) = Cm(5,5)/Ca;
    C(5,6,ii) = Cm(5,6)/Ca;
    C(6,6,ii) = Cm(6,6)/Ca;
    rho0(ii) = rho0m/rhoa;
end


for ii = 1:nPlies
    A11(1,1,ii) = C(1,1,ii);
    A11(1,2,ii) = C(1,6,ii);
    A11(1,3,ii) = C(1,5,ii);
    A11(2,1,ii) = C(1,6,ii);
    A11(2,2,ii) = C(6,6,ii);
    A11(2,3,ii) = C(5,6,ii);
    A11(3,1,ii) = C(1,5,ii);
    A11(3,2,ii) = C(5,6,ii);
    A11(3,3,ii) = C(5,5,ii);

    A12(1,1,ii) = C(1,6,ii);
    A12(1,2,ii) = C(1,2,ii);
    A12(1,3,ii) = C(1,4,ii);
    A12(2,1,ii) = C(6,6,ii);
    A12(2,2,ii) = C(2,6,ii);
    A12(2,3,ii) = C(4,6,ii);
    A12(3,1,ii) = C(5,6,ii);
    A12(3,2,ii) = C(2,5,ii);
    A12(3,3,ii) = C(4,5,ii);

    A22(1,1,ii) = C(6,6,ii);
    A22(1,2,ii) = C(2,6,ii);
    A22(1,3,ii) = C(4,6,ii);
    A22(2,1,ii) = C(2,6,ii);
    A22(2,2,ii) = C(2,2,ii);
    A22(2,3,ii) = C(2,4,ii);
    A22(3,1,ii) = C(4,6,ii);
    A22(3,2,ii) = C(2,4,ii);
    A22(3,3,ii) = C(4,4,ii);

    A33(1,1,ii) = C(5,5,ii);
    A33(1,2,ii) = C(4,5,ii);
    A33(1,3,ii) = C(3,5,ii);
    A33(2,1,ii) = C(4,5,ii);
    A33(2,2,ii) = C(4,4,ii);
    A33(2,3,ii) = C(3,4,ii);
    A33(3,1,ii) = C(3,5,ii);
    A33(3,2,ii) = C(3,4,ii);
    A33(3,3,ii) = C(3,3,ii);

    A31(1,1,ii) = C(1,5,ii);
    A31(1,2,ii) = C(5,6,ii);
    A31(1,3,ii) = C(5,5,ii);
    A31(2,1,ii) = C(1,4,ii);
    A31(2,2,ii) = C(4,6,ii);
    A31(2,3,ii) = C(4,5,ii);
    A31(3,1,ii) = C(1,3,ii);
    A31(3,2,ii) = C(3,6,ii);
    A31(3,3,ii) = C(3,5,ii);

    A32(1,1,ii) = C(5,6,ii);
    A32(1,2,ii) = C(2,5,ii);
    A32(1,3,ii) = C(4,5,ii);
    A32(2,1,ii) = C(4,6,ii);
    A32(2,2,ii) = C(2,4,ii);
    A32(2,3,ii) = C(4,4,ii);
    A32(3,1,ii) = C(3,6,ii);
    A32(3,2,ii) = C(2,3,ii);
    A32(3,3,ii) = C(3,4,ii);

    A12p = A12(:,:,ii);
    A21(:,:,ii) = A12p';
    A31p = A31(:,:,ii);
    A13(:,:,ii) = A31p';
    A32p = A32(:,:,ii);
    A23(:,:,ii) = A32p';

    A1(:,:,ii) = cos(psip)^2*A11(:,:,ii)+cos(psip)*sin(psip)*(A12(:,:,ii)+A21(:,:,ii))+sin(psip)^2*A22(:,:,ii);

    BB(:,:,ii) = cos(psip)*(A13(:,:,ii)+A31(:,:,ii))+sin(psip)*(A23(:,:,ii)+A32(:,:,ii));

    CC(:,:,ii) = -A33(:,:,ii);

    ABC(:,:,ii) = cos(psip)*A31(:,:,ii)+sin(psip)*A32(:,:,ii);
end

% Debut du calcul
ZUn = zeros(nTot,nMax);

F1 = zeros(nTot,nTot);
G1 = zeros(nTot,nTot);
H1 = zeros(nTot,nTot);

Ntemp = 1;

%% Développement sur une base de polynômes de Legendre d'ordre N(ii) pour la
% couche ii = 1 à NPlies
for ii = 1:nPlies
    As = zeros(3*(Nodes(ii)-2),3*Nodes(ii));
    Bs = zeros(3*(Nodes(ii)-2),3*Nodes(ii));
    Cs = zeros(3*(Nodes(ii)-2),3*Nodes(ii));
    
    A2(:,:,ii) = -rho0(ii)*eye(3);
    
    for mm=0:Nodes(ii)-3
        for nn=0:Nodes(ii)-1
            Bs(3*mm+1:3*(mm+1),3*nn+1:3*(nn+1)) = 2/H(ii)*BB(:,:,ii)/ka*PmdPn(mm,nn); 
            Cs(3*mm+1:3*(mm+1),3*nn+1:3*(nn+1)) = 4/(H(ii)^2)*CC(:,:,ii)/ka^2*Pmd2Pn(mm,nn);
            if (mm == nn)
                As(3*mm+1:3*(mm+1),3*nn+1:3*(nn+1)) = 2*A1(:,:,ii)/(2*nn+1);
                Cs(3*mm+1:3*(mm+1),3*nn+1:3*(nn+1)) = Cs(3*mm+1:3*(mm+1),3*nn+1:3*(nn+1))+2*A2(:,:,ii)/(2*nn+1);
            end
        end
    end          
    H1(Ntemp:Ntemp+3*(Nodes(ii)-2)-1,3*sum(Nodes(1:ii-1))+1:3*sum(Nodes(1:ii))) = As;
    F1(Ntemp:Ntemp+3*(Nodes(ii)-2)-1,3*sum(Nodes(1:ii-1))+1:3*sum(Nodes(1:ii))) = Bs;
    G1(Ntemp:Ntemp+3*(Nodes(ii)-2)-1,3*sum(Nodes(1:ii-1))+1:3*sum(Nodes(1:ii))) = Cs;
    
    Ntemp = Ntemp+3*(Nodes(ii)-2);
end

%% Conditions de continuité entre les couches ii et ii+1
for ii = 1:nPlies-1
            
    % Continuité des déplacements entre les couches ii et ii+1
    Ds = zeros(3,3*Nodes(ii+1));
    Es = zeros(3,3*Nodes(ii));
    for nn=0:Nodes(ii)-1
        Es(:,3*nn+1:3*(nn+1)) = eye(3);
    end
    for nn=0:Nodes(ii+1)-1
        Ds(:,3*nn+1:3*(nn+1)) = -(-1)^nn*eye(3);
    end
        
    G1(Ntemp:Ntemp+2,3*sum(Nodes(1:ii-1))+1:3*sum(Nodes(1:ii))) = Es;
    G1(Ntemp:Ntemp+2,3*sum(Nodes(1:ii))+1:3*sum(Nodes(1:ii+1))) = Ds;
    
    Ntemp = Ntemp+3; 
            
    % Continuité des contraintes normales entre les couches ii et
    % ii+1
    Ds = zeros(3,3*Nodes(ii));
    Es = zeros(3,3*Nodes(ii));
    Dp = zeros(3,3*Nodes(ii+1));
    Ep = zeros(3,3*Nodes(ii+1));
    for nn=0:Nodes(ii)-1
        Ds(:,3*nn+1:3*(nn+1)) = -ABC(:,:,ii);
        Es(:,3*nn+1:3*(nn+1)) = 2/H(ii)*A33(:,:,ii)/ka*nn*(nn+1)/2;
    end
    for nn=0:Nodes(ii+1)-1
        Dp(:,3*nn+1:3*(nn+1)) = ABC(:,:,ii+1)*(-1)^nn;
        Ep(:,3*nn+1:3*(nn+1)) = -2/H(ii+1)*A33(:,:,ii+1)/ka*(-1)^(nn+1)*nn*(nn+1)/2;
    end
    F1(Ntemp:Ntemp+2,3*sum(Nodes(1:ii-1))+1:3*sum(Nodes(1:ii))) = Ds;
    F1(Ntemp:Ntemp+2,3*sum(Nodes(1:ii))+1:3*sum(Nodes(1:ii+1))) = Dp;
    G1(Ntemp:Ntemp+2,3*sum(Nodes(1:ii-1))+1:3*sum(Nodes(1:ii))) = Es;
    G1(Ntemp:Ntemp+2,3*sum(Nodes(1:ii))+1:3*sum(Nodes(1:ii+1))) = Ep;
   
    Ntemp = Ntemp+3;
           
end

%% Conditions aux limites sur la surface basse
Ds = zeros(3,3*Nodes(1));
Es = zeros(3,3*Nodes(1));
for nn=0:Nodes(1)-1
   Ds(:,3*nn+1:3*(nn+1)) = ABC(:,:,1)*(-1)^nn;
   Es(:,3*nn+1:3*(nn+1)) = -2/H(1)*A33(:,:,1)/ka*(-1)^(nn+1)*nn*(nn+1)/2;
end
F1(Ntemp:Ntemp+2,1:3*Nodes(1)) = Ds;
G1(Ntemp:Ntemp+2,1:3*Nodes(1)) = Es;
    
Ntemp = Ntemp+3;
      
%% Conditions aux limites sur la surface haute
Dp = zeros(3,3*Nodes(nPlies));
Ep = zeros(3,3*Nodes(nPlies));
for nn=0:Nodes(nPlies)-1
    Dp(:,3*nn+1:3*(nn+1)) = -ABC(:,:,nPlies);
    Ep(:,3*nn+1:3*(nn+1)) = 2/H(nPlies)*A33(:,:,nPlies)/ka*nn*(nn+1)/2;
end
F1(Ntemp:Ntemp+2,3*sum(Nodes(1:nPlies-1))+1:3*sum(Nodes(1:nPlies))) = Dp;
G1(Ntemp:Ntemp+2,3*sum(Nodes(1:nPlies-1))+1:3*sum(Nodes(1:nPlies))) = Ep;
    
%% Calcul des valeurs propres (modes) de la structure
M1 = [F1 -eye(nTot);-H1 zeros(nTot)];
M2 = [G1 zeros(nTot);zeros(nTot) eye(nTot)];

[Z1,K] = eig(M1,M2);

kp = zeros(2*nTot,1);
for ii=1:2*nTot
    if (K(ii,ii) == 0)
        kp(ii) = NaN;
    else
        kp(ii) = i/K(ii,ii)*ka;
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

% interm = zeros(1,2*Ntot);
[interm, Ind] = sort(kp);
k = interm(1:nMax);
ZUn = Z1(1:nTot,Ind(1:nMax));


% Reconstruction du vecteur U
for kk=1:nMax
    for ll=1:3
        for ii=1:nPlies
            Un(ii,:,ll,kk) = zeros(1,length(r));
            for jj=1:Nodes(ii)
                Un(ii,:,ll,kk) = Un(ii,:,ll,kk) + ZUn(ll+(jj-1)*3+sum(Nodes(1:ii-1))*3,kk)*PLeg(jj,:);
            end
        end
    end
end

wavenumber = abs(real(k));
vit = 2*pi*freq./wavenumber;%*sqrt(max(rho0)*rhoa/(max(max(max(C)))*Ca))

vitesse(1) = vit(1);
ll = 1;
for ii = 2:length(vit)
   if (vit(ii) < 0.99999*vit(ii-1) && isfinite(vit(ii)))
    ll = ll+1;
    vitesse(ll) = vit(ii);
    Dep(:,:,:,ll) = Un(:,:,:,ii);
   end
end

vitesse
vitesse(mMode)

%% VISUALIZATION
figure
ht = 0;
hdec = sum(H(1:nPlies));
for ii=1:nPlies
    ht = H(ii)+ht;
    zm = ht - H(ii)/2;
    z = H(ii)*r/2+zm;
    plot((hdec-z)*1e3,abs(Dep(ii,:,1,mMode)))
    hold on
    plot((hdec-z)*1e3,abs(Dep(ii,:,2,mMode)),'k')
    plot((hdec-z)*1e3,abs(Dep(ii,:,3,mMode)),'r')
end