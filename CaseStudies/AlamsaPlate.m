% This script calculates the dispersion curves for a thin CFRP plate for
% ALAMSA project. There are 12 layers of the same material with
% quasiisotropic properties.
clearvars

%% DEFINE LAYERING [Material,phi [rad],theta [rad],psi [rad],thickness [m]]
Layers{1}=Layer('Alamsa',0,0,0,0.23e-3);
Layers{2}=Layer('Alamsa',-pi/4,0,0,0.23e-3);
Layers{3}=Layer('Alamsa',0,0,0,0.23e-3);
Layers{4}=Layer('Alamsa',-pi/4,0,0,0.23e-3);
Layers{5}=Layer('Alamsa',0,0,0,0.23e-3);
Layers{6}=Layer('Alamsa',-pi/4,0,0,0.23e-3);
Layers{7}=Layer('Alamsa',-pi/4,0,0,0.23e-3);
Layers{8}=Layer('Alamsa',0,0,0,0.23e-3);
Layers{9}=Layer('Alamsa',-pi/4,0,0,0.23e-3);
Layers{10}=Layer('Alamsa',0,0,0,0.23e-3);
Layers{11}=Layer('Alamsa',-pi/4,0,0,0.23e-3);
Layers{12}=Layer('Alamsa',0,0,0,0.23e-3);
AlamsaPlt=Sample(Layers);  % create the object representing plate

%% DEFINE OTHER PARAMETERS
psip = 0;          % angle of wave propagation with respect to the main in-plane coordinate axis 
nFreqs = 120;       % number of frequency steps
df=5e3;           % frequency step size [Hz]
legDeg=10;         % degree of Legendre polynomial expansion - determines the maximum number of modes 3/2*legDeg
nModes2Track=10;   % number of modes to be tracked
saveOn=1;
figOn=1;         

%% CALCULATE THE DISPERSION CURVES
[~,~,~]=DispersionCurves(AlamsaPlt,psip,df,nFreqs,legDeg,nModes2Track,saveOn,figOn);

%% CALCULATE THE DISPLACEMENT PROFILE FOR GIVEN MODE AND FREQUENCY
% freq=150e3;
% nMode=3;
% DisplacementProfiles(NomMat,[Phi;Theta;Psi],H,psip,freq,nMode,legDeg)