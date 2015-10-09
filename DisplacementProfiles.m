% LambAnisotropic3DLegendre_U
% Calcul analytique des courbes de dispersion
% dans des matériaux multicouches avec prise en compte
% des effets piézoélectriques et piézomagnétiques
%
% Distribution spatiale pour un mode
%
% Méthode de Legendre pour les NPlies couches
% Conditions aux limites "approchées" sur les deux surfaces
%
clear all

%% LAYUP, PLY PROPERTIES AND ORIENTATION
NomMat = {'Alamsa' 'Alamsa' 'Alamsa' 'Alamsa' 'Alamsa' 'Alamsa' 'Alamsa' 'Alamsa' ...
    'Alamsa' 'Alamsa' 'Alamsa' 'Alamsa'};
h = ones(size(NomMat,2),1)*0.23e-3;        % Thickness of the ply in [m]
NPlies = length(h);                     % Number of plies
freq = 300e3;                           % Frequency to be investigated
w = 2*pi*freq;                          % Angular frequency
NMode = 2;
N = [10 10 10 10 10 10 10 10 10 10 10 10];
phi = [0 -pi/4 0 -pi/4 0 -pi/4 -pi/4 0 -pi/4 0 -pi/4 0]; % First Euler angle
theta = [0 0 0 0 0 0 0 0 0 0 0 0];
psi = [0 0 0 0 0 0 0 0 0 0 0 0];
Ntot = 3*sum(N);
Nmax = 3*sum(N);
psip = 0;
Npts = 100;

%% Parametres de normalisation
Ca = 1e11; % Pa = N/m^2
rhoa = 1e3; % kg/m^3
htot = sum(h);
ka = w*sqrt(rhoa/Ca);
dr = 0.001;
r = -1:dr:1;
for ii=0:(N-1)
    Pn = legendre(ii,r);
    PLeg(ii+1,:) = Pn(1,:);
end

%%
Un = zeros(NPlies,length(r),3,Nmax);

for ii = 1:NPlies
    matp = load_material_elastic(strcat('./',NomMat{ii},'.dat'));
    mat  = RotationConstantesMateriaux_elastic(matp,phi(ii),theta(ii),psi(ii));
    Cm = mat.C;
    rho0m = mat.rho;
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


for ii = 1:NPlies
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
ZUn = zeros(Ntot,Nmax);

F1 = zeros(Ntot,Ntot);
G1 = zeros(Ntot,Ntot);
H1 = zeros(Ntot,Ntot);

Ntemp = 1;

%% Développement sur une base de polynômes de Legendre d'ordre N(ii) pour la
% couche ii = 1 à NPlies
for ii = 1:NPlies
    As = zeros(3*(N(ii)-2),3*N(ii));
    Bs = zeros(3*(N(ii)-2),3*N(ii));
    Cs = zeros(3*(N(ii)-2),3*N(ii));
    
    A2(:,:,ii) = -rho0(ii)*eye(3);
    
    for mm=0:N(ii)-3
        for nn=0:N(ii)-1
            Bs(3*mm+1:3*(mm+1),3*nn+1:3*(nn+1)) = 2/h(ii)*BB(:,:,ii)/ka*PmdPn(mm,nn); 
            Cs(3*mm+1:3*(mm+1),3*nn+1:3*(nn+1)) = 4/(h(ii)^2)*CC(:,:,ii)/ka^2*Pmd2Pn(mm,nn);
            if (mm == nn)
                As(3*mm+1:3*(mm+1),3*nn+1:3*(nn+1)) = 2*A1(:,:,ii)/(2*nn+1);
                Cs(3*mm+1:3*(mm+1),3*nn+1:3*(nn+1)) = Cs(3*mm+1:3*(mm+1),3*nn+1:3*(nn+1))+2*A2(:,:,ii)/(2*nn+1);
            end
        end
    end          
    H1(Ntemp:Ntemp+3*(N(ii)-2)-1,3*sum(N(1:ii-1))+1:3*sum(N(1:ii))) = As;
    F1(Ntemp:Ntemp+3*(N(ii)-2)-1,3*sum(N(1:ii-1))+1:3*sum(N(1:ii))) = Bs;
    G1(Ntemp:Ntemp+3*(N(ii)-2)-1,3*sum(N(1:ii-1))+1:3*sum(N(1:ii))) = Cs;
    
    Ntemp = Ntemp+3*(N(ii)-2);
end

%% Conditions de continuité entre les couches ii et ii+1
for ii = 1:NPlies-1
            
    % Continuité des déplacements entre les couches ii et ii+1
    Ds = zeros(3,3*N(ii+1));
    Es = zeros(3,3*N(ii));
    for nn=0:N(ii)-1
        Es(:,3*nn+1:3*(nn+1)) = eye(3);
    end
    for nn=0:N(ii+1)-1
        Ds(:,3*nn+1:3*(nn+1)) = -(-1)^nn*eye(3);
    end
        
    G1(Ntemp:Ntemp+2,3*sum(N(1:ii-1))+1:3*sum(N(1:ii))) = Es;
    G1(Ntemp:Ntemp+2,3*sum(N(1:ii))+1:3*sum(N(1:ii+1))) = Ds;
    
    Ntemp = Ntemp+3; 
            
    % Continuité des contraintes normales entre les couches ii et
    % ii+1
    Ds = zeros(3,3*N(ii));
    Es = zeros(3,3*N(ii));
    Dp = zeros(3,3*N(ii+1));
    Ep = zeros(3,3*N(ii+1));
    for nn=0:N(ii)-1
        Ds(:,3*nn+1:3*(nn+1)) = -ABC(:,:,ii);
        Es(:,3*nn+1:3*(nn+1)) = 2/h(ii)*A33(:,:,ii)/ka*nn*(nn+1)/2;
    end
    for nn=0:N(ii+1)-1
        Dp(:,3*nn+1:3*(nn+1)) = ABC(:,:,ii+1)*(-1)^nn;
        Ep(:,3*nn+1:3*(nn+1)) = -2/h(ii+1)*A33(:,:,ii+1)/ka*(-1)^(nn+1)*nn*(nn+1)/2;
    end
    F1(Ntemp:Ntemp+2,3*sum(N(1:ii-1))+1:3*sum(N(1:ii))) = Ds;
    F1(Ntemp:Ntemp+2,3*sum(N(1:ii))+1:3*sum(N(1:ii+1))) = Dp;
    G1(Ntemp:Ntemp+2,3*sum(N(1:ii-1))+1:3*sum(N(1:ii))) = Es;
    G1(Ntemp:Ntemp+2,3*sum(N(1:ii))+1:3*sum(N(1:ii+1))) = Ep;
   
    Ntemp = Ntemp+3;
           
end

%% Conditions aux limites sur la surface basse
Ds = zeros(3,3*N(1));
Es = zeros(3,3*N(1));
for nn=0:N(1)-1
   Ds(:,3*nn+1:3*(nn+1)) = ABC(:,:,1)*(-1)^nn;
   Es(:,3*nn+1:3*(nn+1)) = -2/h(1)*A33(:,:,1)/ka*(-1)^(nn+1)*nn*(nn+1)/2;
end
F1(Ntemp:Ntemp+2,1:3*N(1)) = Ds;
G1(Ntemp:Ntemp+2,1:3*N(1)) = Es;
    
Ntemp = Ntemp+3;
      
%% Conditions aux limites sur la surface haute
Dp = zeros(3,3*N(NPlies));
Ep = zeros(3,3*N(NPlies));
for nn=0:N(NPlies)-1
    Dp(:,3*nn+1:3*(nn+1)) = -ABC(:,:,NPlies);
    Ep(:,3*nn+1:3*(nn+1)) = 2/h(NPlies)*A33(:,:,NPlies)/ka*nn*(nn+1)/2;
end
F1(Ntemp:Ntemp+2,3*sum(N(1:NPlies-1))+1:3*sum(N(1:NPlies))) = Dp;
G1(Ntemp:Ntemp+2,3*sum(N(1:NPlies-1))+1:3*sum(N(1:NPlies))) = Ep;
    
%% Calcul des valeurs propres (modes) de la structure
M1 = [F1 -eye(Ntot);-H1 zeros(Ntot)];
M2 = [G1 zeros(Ntot);zeros(Ntot) eye(Ntot)];

[Z1,K] = eig(M1,M2);

kp = zeros(2*Ntot,1);
for ii=1:2*Ntot
    if (K(ii,ii) == 0)
        kp(ii) = NaN;
    else
        kp(ii) = i/K(ii,ii)*ka;
    end
end

for ii=1:2*Ntot
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
k = interm(1:Nmax);
ZUn = Z1(1:Ntot,Ind(1:Nmax));


% Reconstruction du vecteur U
for kk=1:Nmax
    for ll=1:3
        for ii=1:NPlies
            Un(ii,:,ll,kk) = zeros(1,length(r));
            for jj=1:N(ii)
                Un(ii,:,ll,kk) = Un(ii,:,ll,kk) + ZUn(ll+(jj-1)*3+sum(N(1:ii-1))*3,kk)*PLeg(jj,:);
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
vitesse(NMode)

%% VISUALIZATION
figure
ht = 0;
hdec = sum(h(1:NPlies));
for ii=1:NPlies
    ht = h(ii)+ht;
    zm = ht - h(ii)/2;
    z = h(ii)*r/2+zm;
    plot((hdec-z)*1e3,abs(Dep(ii,:,1,NMode)))
    hold on
    plot((hdec-z)*1e3,abs(Dep(ii,:,2,NMode)),'k')
    plot((hdec-z)*1e3,abs(Dep(ii,:,3,NMode)),'r')
end