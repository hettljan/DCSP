% Calculates the dispersion curves for layered anisotropic materials using
% the Legendre and Laguerre polynomial approach
% Originally supplied by O. Bou-Matar, Lille
clear all

%%
% NomMat = {'Plexiglass'};
% h = 8e-3*ones(size(NomMat,2),1);     % Thickness of the ply in [m]
% NPlies = length(h);                     % Number of plies
% phi = [0]; % First Euler angle
% theta = [0];                       % Second Euler angle
% psi = [0];                 % Third Euler angle
% N = [10];
% Ntot = 3*sum(N);           % Total number of modes
% Nmax = 3*sum(N);           % Max number of modes????        
% psip = 0;
% Npts = 655;                % Number of frequency steps
% df=1.2207e3;                   %  frequency step
% nModes=15;

%% LAYUP, PLY PROPERTIES AND ORIENTATION
NomMat = {'Alamsa' 'Alamsa' 'Alamsa' 'Alamsa' 'Alamsa' 'Alamsa'...
    'Alamsa' 'Alamsa' 'Alamsa' 'Alamsa' 'Alamsa' 'Alamsa'};
h=0.23e-3;                                             % Thickness of the ply in [m]
Phi=[0 -pi/4 0 -pi/4 0 -pi/4 -pi/4 0 -pi/4 0 -pi/4 0]; % First Euler angle
Theta=[0 0 0 0 0 0 0 0 0 0 0 0];                       % Second Euler angle
Psi = [0 0 0 0 0 0 0 0 0 0 0 0];                       % Third Euler angle                
psip = 0;
Npts = 135;                     % number of frequency steps
df=6.1035e3;                    % frequency step size [Hz]
nModes=12;

%%
nPlies = size(NomMat,2);        % Number of plies
N=ones(nPlies,1)*10;            % vector with number of nodes per layer 
Ntot=3*sum(N);                  % total number of modes
Nmax=3*sum(N);                  % max number of modes????   
H=ones(nPlies,1)*h;             % vector of ply thicknesses

%% NORMALIZATION PARAMETERS
Ca = 1e11;              % Pa = N/m^2
rhoa = 1e3;             % kg/m^3
dw = 2*pi*df;           % Angular frequency step
htot=sum(H);            % Total thickness of the plate
k = zeros(Nmax,Npts);
Un = zeros(Ntot,Nmax,Npts);

%% PREPARE PARTIAL MATRICES
freq=nan(Npts,1);       % Initialize the frequency vector
rho0=nan(nPlies,1);     % Density vector
A11=nan(3,3);
A12=nan(size(A11));
A21=nan(size(A11));
A22=nan(size(A11));
A33=nan(3,3,nPlies);
A31=nan(size(A11));
A13=nan(size(A11));
A32=nan(size(A11));
A23=nan(size(A11));
A1=nan(3,3,nPlies);
BB=nan(3,3,nPlies);
CC=nan(3,3,nPlies);
ABC=nan(3,3,nPlies);
A2=nan(3,3,nPlies);

%% LOAD AND STACK THE MATERIAL PROPERTIES FOR PLIES
for ply = 1:nPlies
    matp = LoadElasticConstants(fullfile('./Materials',strcat(NomMat{ply},'.dat')));
    mat  = RotateElasticConstants(matp,Phi(ply),Theta(ply),Psi(ply));
    C = mat.C./Ca;
    rho0(ply) = mat.rho/rhoa;
    A11(1,1) = C(1,1);
    A11(1,2) = C(1,6);
    A11(1,3) = C(1,5);
    A11(2,1) = C(1,6);
    A11(2,2) = C(6,6);
    A11(2,3) = C(5,6);
    A11(3,1) = C(1,5);
    A11(3,2) = C(5,6);
    A11(3,3) = C(5,5);

    A12(1,1) = C(1,6);
    A12(1,2) = C(1,2);
    A12(1,3) = C(1,4);
    A12(2,1) = C(6,6);
    A12(2,2) = C(2,6);
    A12(2,3) = C(4,6);
    A12(3,1) = C(5,6);
    A12(3,2) = C(2,5);
    A12(3,3) = C(4,5);

    A22(1,1) = C(6,6);
    A22(1,2) = C(2,6);
    A22(1,3) = C(4,6);
    A22(2,1) = C(2,6);
    A22(2,2) = C(2,2);
    A22(2,3) = C(2,4);
    A22(3,1) = C(4,6);
    A22(3,2) = C(2,4);
    A22(3,3) = C(4,4);

    A33(1,1,ply) = C(5,5);
    A33(1,2,ply) = C(4,5);
    A33(1,3,ply) = C(3,5);
    A33(2,1,ply) = C(4,5);
    A33(2,2,ply) = C(4,4);
    A33(2,3,ply) = C(3,4);
    A33(3,1,ply) = C(3,5);
    A33(3,2,ply) = C(3,4);
    A33(3,3,ply) = C(3,3);

    A31(1,1) = C(1,5);
    A31(1,2) = C(5,6);
    A31(1,3) = C(5,5);
    A31(2,1) = C(1,4);
    A31(2,2) = C(4,6);
    A31(2,3) = C(4,5);
    A31(3,1) = C(1,3);
    A31(3,2) = C(3,6);
    A31(3,3) = C(3,5);

    A32(1,1) = C(5,6);
    A32(1,2) = C(2,5);
    A32(1,3) = C(4,5);
    A32(2,1) = C(4,6);
    A32(2,2) = C(2,4);
    A32(2,3) = C(4,4);
    A32(3,1) = C(3,6);
    A32(3,2) = C(2,3);
    A32(3,3) = C(3,4);

    A21 = A12';
    A13 = A31';
    A23 = A32';

    A1(:,:,ply) = cos(psip)^2*A11+cos(psip)*sin(psip)*(A12+A21)+sin(psip)^2*A22;
    BB(:,:,ply) = cos(psip)*(A13+A31)+sin(psip)*(A23+A32);
    CC(:,:,ply) = -A33(:,:,ply);
    ABC(:,:,ply) = cos(psip)*A31+sin(psip)*A32;
    A2(:,:,ply) = -rho0(ply)*eye(3);
end

%% CALCULATION LOOP
tic;
parfor_progress(Npts); % Initialize 
parfor kk = 0:Npts-1
%     fprintf('Step %d\n',kk+1);
    w = dw+dw*kk;
    freq(kk+1) = w/(2*pi);
    ka = w*sqrt(rhoa/Ca);
    F1 = zeros(Ntot,Ntot);
    G1 = zeros(Ntot,Ntot);
    H1 = zeros(Ntot,Ntot);
    Ntemp = 1;
    % Développement sur une base de polynômes de Legendre d'ordre N(ii) pour la
    % couche ii = 1 à NCouche
    for ply = 1:nPlies
        As = zeros(3*(N(ply)-2),3*N(ply));
        Bs = zeros(3*(N(ply)-2),3*N(ply));
        Cs = zeros(3*(N(ply)-2),3*N(ply));
        for mm=0:N(ply)-3
            for nn=0:N(ply)-1
                Bs(3*mm+1:3*(mm+1),3*nn+1:3*(nn+1)) = 2/H(ply)*BB(:,:,ply)/ka*PmdPn(mm,nn); 
                Cs(3*mm+1:3*(mm+1),3*nn+1:3*(nn+1)) = 4/(H(ply)^2)*CC(:,:,ply)/ka^2*Pmd2Pn(mm,nn);
                if (mm == nn)
                    As(3*mm+1:3*(mm+1),3*nn+1:3*(nn+1)) = 2*A1(:,:,ply)/(2*nn+1);
                    Cs(3*mm+1:3*(mm+1),3*nn+1:3*(nn+1)) = Cs(3*mm+1:3*(mm+1),3*nn+1:3*(nn+1))+2*A2(:,:,ply)/(2*nn+1);
                end
            end
        end          
        H1(Ntemp:Ntemp+3*(N(ply)-2)-1,3*sum(N(1:ply-1))+1:3*sum(N(1:ply))) = As;
        F1(Ntemp:Ntemp+3*(N(ply)-2)-1,3*sum(N(1:ply-1))+1:3*sum(N(1:ply))) = Bs;
        G1(Ntemp:Ntemp+3*(N(ply)-2)-1,3*sum(N(1:ply-1))+1:3*sum(N(1:ply))) = Cs;
        Ntemp = Ntemp+3*(N(ply)-2);
    end

% Conditions de continuité entre les couches ii et ii+1
    for ply = 1:nPlies-1
        % Continuité des déplacements entre les couches ii et ii+1
        Ds = zeros(3,3*N(ply+1));
        Es = zeros(3,3*N(ply));
        for nn=0:N(ply)-1
            Es(:,3*nn+1:3*(nn+1)) = eye(3);
        end
        for nn=0:N(ply+1)-1
            Ds(:,3*nn+1:3*(nn+1)) = -(-1)^nn*eye(3);
        end

        G1(Ntemp:Ntemp+2,3*sum(N(1:ply-1))+1:3*sum(N(1:ply))) = Es;
        G1(Ntemp:Ntemp+2,3*sum(N(1:ply))+1:3*sum(N(1:ply+1))) = Ds;

        Ntemp = Ntemp+3; 

        % Continuité des contraintes normales entre les couches ii et
        % ii+1
        Ds = zeros(3,3*N(ply));
        Es = zeros(3,3*N(ply));
        Dp = zeros(3,3*N(ply+1));
        Ep = zeros(3,3*N(ply+1));
        for nn=0:N(ply)-1
            Ds(:,3*nn+1:3*(nn+1)) = -ABC(:,:,ply);
            Es(:,3*nn+1:3*(nn+1)) = 2/H(ply)*A33(:,:,ply)/ka*nn*(nn+1)/2;
        end
        for nn=0:N(ply+1)-1
            Dp(:,3*nn+1:3*(nn+1)) = ABC(:,:,ply+1)*(-1)^nn;
            Ep(:,3*nn+1:3*(nn+1)) = -2/H(ply+1)*A33(:,:,ply+1)/ka*(-1)^(nn+1)*nn*(nn+1)/2;
        end
        F1(Ntemp:Ntemp+2,3*sum(N(1:ply-1))+1:3*sum(N(1:ply))) = Ds;
        F1(Ntemp:Ntemp+2,3*sum(N(1:ply))+1:3*sum(N(1:ply+1))) = Dp;
        G1(Ntemp:Ntemp+2,3*sum(N(1:ply-1))+1:3*sum(N(1:ply))) = Es;
        G1(Ntemp:Ntemp+2,3*sum(N(1:ply))+1:3*sum(N(1:ply+1))) = Ep;

        Ntemp = Ntemp+3;            
    end

    % Conditions aux limites sur la surface basse
    Ds = zeros(3,3*N(1));
    Es = zeros(3,3*N(1));
    for nn=0:N(1)-1
        Ds(:,3*nn+1:3*(nn+1)) = ABC(:,:,1)*(-1)^nn;
        Es(:,3*nn+1:3*(nn+1)) = -2/h(1)*A33(:,:,1)/ka*(-1)^(nn+1)*nn*(nn+1)/2;
    end
    F1(Ntemp:Ntemp+2,1:3*N(1)) = Ds;
    G1(Ntemp:Ntemp+2,1:3*N(1)) = Es;        

    Ntemp = Ntemp+3;

    % Conditions aux limites sur la surface haute
    Dp = zeros(3,3*N(nPlies));
    Ep = zeros(3,3*N(nPlies));
    for nn=0:N(nPlies)-1
        Dp(:,3*nn+1:3*(nn+1)) = -ABC(:,:,nPlies);
        Ep(:,3*nn+1:3*(nn+1)) = 2/H(nPlies)*A33(:,:,nPlies)/ka*nn*(nn+1)/2;
    end
    F1(Ntemp:Ntemp+2,3*sum(N(1:nPlies-1))+1:3*sum(N(1:nPlies))) = Dp;
    G1(Ntemp:Ntemp+2,3*sum(N(1:nPlies-1))+1:3*sum(N(1:nPlies))) = Ep;

    % Calcul des valeurs propres (modes) de la structure
    M1 = [F1 -eye(Ntot);-H1 zeros(Ntot)];
    M2 = [G1 zeros(Ntot);zeros(Ntot) eye(Ntot)];

    [Z1,K] = eig(M1,M2);

    kp = zeros(2*Ntot,1);
    for ply=1:2*Ntot
        if (K(ply,ply) == 0)
            kp(ply) = NaN;
        else
            kp(ply) = i/K(ply,ply)*ka;
        end
    end

    for ply=1:2*Ntot
        if (real(kp(ply)) == 0)
            kp(ply) = NaN;
        else
            if (abs(imag(kp(ply)))/abs(real(kp(ply))) > 1e-8)
                kp(ply) = NaN;        
            end
        end
    end

    [interm, Ind] = sort(kp);
    k(:,kk+1) = interm(1:Nmax);
    Un(:,:,kk+1) = Z1(1:Ntot,Ind(1:Nmax));
    parfor_progress;
end
parfor_progress(0);
toc;

%% ADJUST THE UNITS
Wavenumber = abs(real(k))./(2*pi);
Wavenumber(Wavenumber>2000)=nan;    % Delete the insanely high wavenumbers
Wavenumber=Wavenumber(1:2:end,:);
Wavenumber=DispersionCurveSorting(freq,Wavenumber,nModes);

%% CALCULATE THE VELOCITY
Velocity=nan(size(Wavenumber));
for i=1:size(Wavenumber,1)
    Velocity(i,:)=freq'./squeeze(Wavenumber(i,:));
end 

%% VISUALIZATION
figure
subplot(1,2,1)
plot(freq*1e-3,Wavenumber,'*')
xlim([freq(1),freq(end)]*1e-3)
ylim([0,700]);
xlabel(strcat('Frequency [kHz]'),'FontSize',14)
ylabel(strcat('Wavenumber [m^{-1}]'),'FontSize',14)

subplot(1,2,2)
hold on
plot(freq*1e-3,Velocity,'*')
xlim([freq(1) freq(end)]*1e-3)
ylim([0 1e4]);
xlabel(strcat('Frequency [kHz]'),'FontSize',14)
ylabel(strcat('Phase velocity [ms^{-1}]'),'FontSize',14)
