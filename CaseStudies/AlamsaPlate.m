% This script calculates the dispersion curves for a thin CFRP plate for
% ALAMSA project. There are 12 layers of the same material with
% quasiisotropic properties.

psip = 0;           % angle of wave propagation with respect to the main in-plane coordinate axis 
nFreqs = 400;       % number of frequency steps
df=2.5e3;             % frequency step size [Hz]
legDeg=10;          % degree of Legendre polynomial expansion - determines the maximum number of modes 3/2*legDeg
nModes2Track=5;     % number of modes to be tracked
NomMat = {'Alamsa' 'Alamsa' 'Alamsa' 'Alamsa' 'Alamsa' 'Alamsa'...
    'Alamsa' 'Alamsa' 'Alamsa' 'Alamsa' 'Alamsa' 'Alamsa'};
Phi=[0 -pi/4 0 -pi/4 0 -pi/4 -pi/4 0 -pi/4 0 -pi/4 0]; % First Euler angle
Theta=[0 0 0 0 0 0 0 0 0 0 0 0];                       % Second Euler angle
Psi = [0 0 0 0 0 0 0 0 0 0 0 0];                       % Third Euler angle                
h=0.23e-3;                      % Thickness of the ply in [m]
nPlies = size(NomMat,2);        % Number of plies
H=ones(nPlies,1)*h;             % vector of ply thicknesses

DispersionCurves(NomMat,[Phi;Theta;Psi],H,psip,df,nFreqs,legDeg,nModes2Track)