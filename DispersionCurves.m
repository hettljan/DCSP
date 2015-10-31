function [Freq,Wavenumber,Velocity]=DispersionCurves(NomMat,EulerAngles,H,...
    psip,df,nFreqs,legDeg,nModes2Track,varargin)
% Calculates the dispersion curves for layered anisotropic materials using
% the Legendre and Laguerre polynomial approach
% Originally supplied by O. Bou-Matar, Lille
%
% INPUT: 
% NomMat        - Cell array with names of materials (link to text files .dat)
% EulerAngles   - 1st, 2nd, 3rd Euler angle in [3xnPlies] matrix, phi is in first
%                 row, theta in 2nd row, psi in 3rd row, all in [rad]
% H             - Vector with ply/layer thicknesses [m]
% psip          - Angle of wave propagation with respect to the main
%                 in-plane coordinate axis [rad]
% df            - Frequency step in [Hz]
% nFreqs        - Number of frequency steps
% legDeg        - Degree of Legendre polynomial expansion - determines the maximum number of modes 3/2*legDeg
% nModes2Track  - Number of modes to be tracked
% OPTIONAL:
%   saveOn      - Enable saving of the dispersion data
%   figOn       - Enable plotting

%% INPUT PARSING
numvarargs = length(varargin);
if numvarargs > 2
    error('myfuns:somefun2Alt:TooManyInputs', ...
        'requires at most 2 optional inputs');
end
optargs = {0,0};
optargs(1:numvarargs) = varargin;
[saveOn,figOn] = optargs{:};

%% PREPARE THE OTHER PARAMETERS
folder='\\plato.kulak.be\groupdata\acoustics\Jan_Hettler\MATLAB\Simulation\DCSP\Materials';
Phi=EulerAngles(1,:);
Theta=EulerAngles(2,:);
Psi=EulerAngles(3,:);
nPlies = size(NomMat,2);        % Number of plies
Nodes=ones(nPlies,1)*legDeg;    % vector with the number of nodes/order of polynomial expansion  per layer 
nNodes=3*sum(Nodes);            % total number of nodes in laminate x 3 components of displacement
nMax=3*sum(Nodes);              % max number of modes????   

%% NORMALIZATION PARAMETERS
Ca = 1e11;                      % Pa = N/m^2 normalization parameter
rhoa = 1e3;                     % kg/m^3 second normalization parameter 
dw = 2*pi*df;                   % Angular frequency step
k = zeros(nMax,nFreqs);
Un = zeros(nNodes,nMax,nFreqs);

%% PREPARE PARTIAL MATRICES
Freq=nan(nFreqs,1);     % Initialize the frequency vector
Rho0=nan(nPlies,1);     % Density vector
F11=nan(3,3);
F12=nan(size(F11));
F22=nan(size(F11));
F33=nan(3,3,nPlies);
F31=nan(size(F11));
F32=nan(size(F11));
A1=nan(3,3,nPlies);
BB=nan(3,3,nPlies);
CC=nan(3,3,nPlies);
ABC=nan(3,3,nPlies);
A2=nan(3,3,nPlies);

%% LOAD AND STACK THE MATERIAL PROPERTIES FOR PLIES AND STROH MATRIX
for ply = 1:nPlies
    Matp = LoadElasticConstants(fullfile(folder,strcat(NomMat{ply},'.dat')));
    C = RotateElasticConstants(Matp.C,Phi(ply),Theta(ply),Psi(ply)); % stiffness tensor rotated to principal axis (of anisotropy)
    C = C./Ca;                  % normalization of stiffness tensor 
    Rho0(ply) = Matp.rho/rhoa;  % convert kg/m^3 to g/cm^3
    F11(1,1) = C(1,1);          % F11 matrix component for building Stroh matrix
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

    A1(:,:,ply) = F11*cos(psip)^2+(F12+F21)*cos(psip)*sin(psip)+F22*sin(psip)^2;
    BB(:,:,ply) = (F13+F31)*cos(psip)+(F23+F32)*sin(psip);
    CC(:,:,ply) = -F33(:,:,ply);
    ABC(:,:,ply) = F31*cos(psip)+F32*sin(psip);
    A2(:,:,ply) = -Rho0(ply)*eye(3);
end

%% CALCULATION LOOP
parfor_progress(nFreqs);
for kk=0:nFreqs-1
    w = dw+dw*kk;
    Freq(kk+1) = w/(2*pi);
    ka = w*sqrt(rhoa/Ca);
    F1 = zeros(nNodes,nNodes);
    G1 = zeros(nNodes,nNodes);
    H1 = zeros(nNodes,nNodes);
    Ntemp = 1;
    for ply = 1:nPlies      % Preparation of the computational values from Legendre polynomials
        As = zeros(3*(Nodes(ply)-2),3*Nodes(ply));
        Bs = zeros(3*(Nodes(ply)-2),3*Nodes(ply));
        Cs = zeros(3*(Nodes(ply)-2),3*Nodes(ply));
        for mm=0:Nodes(ply)-3
            for nn=0:Nodes(ply)-1
                Bs(3*mm+1:3*(mm+1),3*nn+1:3*(nn+1)) = 2/H(ply)*BB(:,:,ply)/ka*PmdPn(mm,nn); 
                Cs(3*mm+1:3*(mm+1),3*nn+1:3*(nn+1)) = 4/(H(ply)^2)*CC(:,:,ply)/ka^2*Pmd2Pn(mm,nn);
                if (mm == nn)
                    As(3*mm+1:3*(mm+1),3*nn+1:3*(nn+1)) = 2*A1(:,:,ply)/(2*nn+1);
                    Cs(3*mm+1:3*(mm+1),3*nn+1:3*(nn+1)) = Cs(3*mm+1:3*(mm+1),3*nn+1:3*(nn+1))+2*A2(:,:,ply)/(2*nn+1);
                end
            end
        end          
        H1(Ntemp:Ntemp+3*(Nodes(ply)-2)-1,3*sum(Nodes(1:ply-1))+1:3*sum(Nodes(1:ply))) = As;
        F1(Ntemp:Ntemp+3*(Nodes(ply)-2)-1,3*sum(Nodes(1:ply-1))+1:3*sum(Nodes(1:ply))) = Bs;
        G1(Ntemp:Ntemp+3*(Nodes(ply)-2)-1,3*sum(Nodes(1:ply-1))+1:3*sum(Nodes(1:ply))) = Cs;
        Ntemp = Ntemp+3*(Nodes(ply)-2);
    end

    % Conditions of continuity between plies - ply and ply+1
    for ply = 1:nPlies-1
        % Continuity of displacement on the interface - ply and ply+1
        Ds = zeros(3,3*Nodes(ply+1));
        Es = zeros(3,3*Nodes(ply));
        for nn=0:Nodes(ply)-1
            Es(:,3*nn+1:3*(nn+1)) = eye(3);
        end
        for nn=0:Nodes(ply+1)-1
            Ds(:,3*nn+1:3*(nn+1)) = -(-1)^nn*eye(3);
        end

        G1(Ntemp:Ntemp+2,3*sum(Nodes(1:ply-1))+1:3*sum(Nodes(1:ply))) = Es;
        G1(Ntemp:Ntemp+2,3*sum(Nodes(1:ply))+1:3*sum(Nodes(1:ply+1))) = Ds;

        Ntemp = Ntemp+3; 

        % Continuity of normal stress on the interface - ply and ply+1
        % ii+1
        Ds = zeros(3,3*Nodes(ply));
        Es = zeros(3,3*Nodes(ply));
        Dp = zeros(3,3*Nodes(ply+1));
        Ep = zeros(3,3*Nodes(ply+1));
        for nn=0:Nodes(ply)-1
            Ds(:,3*nn+1:3*(nn+1)) = -ABC(:,:,ply);
            Es(:,3*nn+1:3*(nn+1)) = 2/H(ply)*F33(:,:,ply)/ka*nn*(nn+1)/2;
        end
        for nn=0:Nodes(ply+1)-1
            Dp(:,3*nn+1:3*(nn+1)) = ABC(:,:,ply+1)*(-1)^nn;
            Ep(:,3*nn+1:3*(nn+1)) = -2/H(ply+1)*F33(:,:,ply+1)/ka*(-1)^(nn+1)*nn*(nn+1)/2;
        end
        F1(Ntemp:Ntemp+2,3*sum(Nodes(1:ply-1))+1:3*sum(Nodes(1:ply))) = Ds;
        F1(Ntemp:Ntemp+2,3*sum(Nodes(1:ply))+1:3*sum(Nodes(1:ply+1))) = Dp;
        G1(Ntemp:Ntemp+2,3*sum(Nodes(1:ply-1))+1:3*sum(Nodes(1:ply))) = Es;
        G1(Ntemp:Ntemp+2,3*sum(Nodes(1:ply))+1:3*sum(Nodes(1:ply+1))) = Ep;

        Ntemp = Ntemp+3;            
    end

    % Boundary conditions at the lower interface
    Ds = zeros(3,3*Nodes(1));
    Es = zeros(3,3*Nodes(1));
    for nn=0:Nodes(1)-1
        Ds(:,3*nn+1:3*(nn+1)) = ABC(:,:,1)*(-1)^nn;
        Es(:,3*nn+1:3*(nn+1)) = -2/H(1)*F33(:,:,1)/ka*(-1)^(nn+1)*nn*(nn+1)/2;
    end
    F1(Ntemp:Ntemp+2,1:3*Nodes(1)) = Ds;
    G1(Ntemp:Ntemp+2,1:3*Nodes(1)) = Es;        

    Ntemp = Ntemp+3;

    % Boundary conditions at the upper interface
    Dp = zeros(3,3*Nodes(nPlies));
    Ep = zeros(3,3*Nodes(nPlies));
    for nn=0:Nodes(nPlies)-1
        Dp(:,3*nn+1:3*(nn+1)) = -ABC(:,:,nPlies);
        Ep(:,3*nn+1:3*(nn+1)) = 2/H(nPlies)*F33(:,:,nPlies)/ka*nn*(nn+1)/2;
    end
    F1(Ntemp:Ntemp+2,3*sum(Nodes(1:nPlies-1))+1:3*sum(Nodes(1:nPlies))) = Dp;
    G1(Ntemp:Ntemp+2,3*sum(Nodes(1:nPlies-1))+1:3*sum(Nodes(1:nPlies))) = Ep;

    % Calculation of the eigenvalue problem
    M1 = [F1 -eye(nNodes);-H1 zeros(nNodes)];
    M2 = [G1 zeros(nNodes);zeros(nNodes) eye(nNodes)];
    [Z1,K] = eig(M1,M2);        % bottleneck, calculates the mode shapes and "wavenumbers"
    kp = zeros(2*nNodes,1);
    
    % clean-up the unreasonable wavenumbers
    for ply=1:2*nNodes
        if (K(ply,ply) == 0)
            kp(ply) = NaN;
        else
            kp(ply) = 1i/K(ply,ply)*ka;
        end
    end
    for ply=1:2*nNodes
        if (real(kp(ply)) == 0)
            kp(ply) = NaN;
        else
            if (abs(imag(kp(ply)))/abs(real(kp(ply))) > 1e-8)
                kp(ply) = NaN;        
            end
        end
    end

    [interm, Ind] = sort(kp);
    k(:,kk+1) = interm(1:nMax);
    Un(:,:,kk+1) = Z1(1:nNodes,Ind(1:nMax));
    parfor_progress;
end
parfor_progress(0);

%% ADJUST THE UNITS
Wavenumber = abs(real(k))./(2*pi);
Wavenumber(Wavenumber>2000)=nan;      % Delete the insanely high wavenumbers
Wavenumber=Wavenumber(1:2:end,:);     % two subsequent modes are duplicates, so take just one of them
Wavenumber=DispersionCurveSorting(Freq,Wavenumber,nModes2Track);

%% CALCULATE THE VELOCITY
Velocity=nan(size(Wavenumber));
for mode=1:size(Wavenumber,1)
    Velocity(mode,:)=Freq'./squeeze(Wavenumber(mode,:));
end 

%% VISUALIZATION
if figOn == 1
    fUnits='f';
    switch fUnits
        case 'f'
            Freq=Freq*1e-3; % frequency in kHz
            xLab='Frequency [kHz]';
        case 'fd'
            Freq=Freq*sum(H)*1e-3;   % frequency-thickness product in MHz*mm
            xLab='Frequency-thickness [MHz \cdot mm]';
    end

    figure
    subplot(2,1,1)
    plot(Freq,Wavenumber,'*')
    xlim([Freq(1),Freq(end)])
    ylim([0,700]);
    xlabel(xLab,'FontSize',14)
    ylabel(strcat('Wavenumber [m^{-1}]'),'FontSize',14)

    subplot(2,1,2)
    hold on
    plot(Freq,Velocity,'*')
    xlim([Freq(1) Freq(end)])
    ylim([0 1e4]);
    xlabel(xLab,'FontSize',14)
    ylabel(strcat('Phase velocity [ms^{-1}]'),'FontSize',14)
end

%% SAVING
if saveOn == 1
    save 'DispData' 'Freq' 'Wavenumber' 'Velocity' 
end