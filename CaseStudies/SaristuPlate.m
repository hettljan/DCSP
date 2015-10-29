% This script calculates the dispersion curves for a thin CFRP plate for
% SARISTU project. There are 11 layers with different elastic properties
clear all
psip = 0;           % angle of wave propagation with respect to the main in-plane coordinate axis 
nFreqs = 100;       % number of frequency steps
df=5e3;             % frequency step size [Hz]
legDeg=10;          % degree of Legendre polynomial expansion - determines the maximum number of modes 3/2*legDeg
nModes2Track=10;    % number of modes to be tracked
saveOn=1;
figOn=1;
NomMat = {'5Harness' 'NCFBiaxial' 'Uniweave' 'NCFBiaxial' 'NCFBiaxial' 'NCFBiaxial'...
    'NCFBiaxial' 'Uniweave' 'NCFBiaxial' '5Harness' 'Glass'};
Phi=[0 pi/4 0 pi/4 0 pi/2 -pi/4 0 -pi/4 pi/2 0]; % First Euler angle
Theta=[0 0 0 0 0 0 0 0 0 0 0];                       % Second Euler angle
Psi = [0 0 0 0 0 0 0 0 0 0 0];                       % Third Euler angle                
h=0.4801e-3;                      	% Thickness of the ply in [m]
nPlies = size(NomMat,2);         	% Number of plies
H=[ones(nPlies-1,1)*h; 0.127e-3] ;  % vector of ply thicknesses

%% CALCULATE THE DISPERSION CURVES
[~,~,~]=DispersionCurves(NomMat,[Phi;Theta;Psi],H,psip,df,nFreqs,legDeg,nModes2Track,saveOn,figOn);

%% CALCULATE THE DISPLACEMENT PROFILE FOR GIVEN MODE AND FREQUENCY
freq=50e3;
nMode=1;
DisplacementProfiles(NomMat,[Phi;Theta;Psi],H,psip,freq,nMode,legDeg)

