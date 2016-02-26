% This script calculates the dispersion curves for a thin CFRP plate for
% SARISTU project. There are 11 layers with different elastic properties
clearvars

%% DEFINE LAYERING [Material,phi [rad],theta [rad],psi [rad],thickness [m]]
Layers{1}=Layer('Materials/Uniweave.dat',0,0,0,2e-3);
SaristuPlt=Sample(Layers);  % create the object representing plate

%% DEFINE OTHER PARAMETERS
psip = 0;         % angle of wave propagation with respect to the main in-plane coordinate axis 
nFreqs = 500;     % number of frequency steps
df=1e3;           % frequency step size [Hz]
legDeg=10;        % degree of Legendre polynomial expansion - determines the maximum number of modes 3/2*legDeg
nModes2Track=15;  % number of modes to be tracked
saveOn=0;
figOn=1; 
parOn=0;

%% CALCULATE THE DISPERSION CURVES
[~,~,~]=DispersionCurves(SaristuPlt,psip,df,nFreqs,legDeg,nModes2Track,...
    saveOn,figOn,parOn);

%% GEOMETRICAL DISPERSION
freq=20e3;
AngleRange=linspace(-pi,pi,360)'; % propagation angle from -pi to pi (360 degrees);
nModes2Track=3;     % Number of modes to be tracked
legDeg=10;          % Degree of Legendre polynomial expansion - 
saveOn=0;           % Enable saving of the dispersion data
figOn=1; 
[PhaseVelocity]=GeomDispersion(SaristuPlt,AngleRange,freq,nModes2Track,...
    legDeg,saveOn,figOn);

