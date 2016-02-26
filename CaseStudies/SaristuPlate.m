% This script calculates the dispersion curves for a thin CFRP plate for
% SARISTU project. There are 11 layers with different elastic properties
clearvars

%% DEFINE LAYERING [Material,phi [rad],theta [rad],psi [rad],thickness [m]]
Layers{1}=Layer('5Harness',0,0,0,0.4801e-3);
Layers{2}=Layer('NCFBiaxial',pi/4,0,0,0.4801e-3);
Layers{3}=Layer('Uniweave',0,0,0,0.4801e-3);
Layers{4}=Layer('NCFBiaxial',pi/4,0,0,0.4801e-3);
Layers{5}=Layer('NCFBiaxial',0,0,0,0.4801e-3);
Layers{6}=Layer('NCFBiaxial',pi/2,0,0,0.4801e-3);
Layers{7}=Layer('NCFBiaxial',-pi/4,0,0,0.4801e-3);
Layers{8}=Layer('Uniweave',0,0,0,0.4801e-3);
Layers{9}=Layer('NCFBiaxial',-pi/4,0,0,0.4801e-3);
Layers{10}=Layer('5Harness',pi/2,0,0,0.4801e-3);
Layers{11}=Layer('Glass',0,0,0,0.127e-3);
SaristuPlt=Sample(Layers);  % create the object representing plate

%% DEFINE OTHER PARAMETERS
psip = 0;          % angle of wave propagation with respect to the main in-plane coordinate axis 
nFreqs = 300;       % number of frequency steps
df=2e3;           % frequency step size [Hz]
legDeg=10;         % degree of Legendre polynomial expansion - determines the maximum number of modes 3/2*legDeg
nModes2Track=15;   % number of modes to be tracked
saveOn=0;
figOn=1;         

%% CALCULATE THE DISPERSION CURVES
[~,~,~]=DispersionCurves(SaristuPlt,psip,df,nFreqs,legDeg,nModes2Track,...
    saveOn,figOn);

%% CALCULATE THE DISPLACEMENT PROFILE FOR GIVEN MODE AND FREQUENCY
% freq=50e3;
% nMode=3;
% DisplacementProfiles(SaristuPlt,psip,freq,nMode,legDeg)

