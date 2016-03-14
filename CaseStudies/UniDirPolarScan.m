% This script calculates the dispersion curves for a thin CFRP plate for
% SARISTU project. There are 11 layers with different elastic properties
clearvars

%% DEFINE LAYERING [Material,phi [rad],theta [rad],psi [rad],thickness [m]]
folder='C:\Users\u0088749\Box Sync\MATLAB\Simulation\DCSP\Materials';
Layers{1}=Layer(fullfile(folder,'UniDirPolarScan.dat'),0,0,0,1e-3);
SaristuPlt=Sample(Layers);  % create the object representing plate

%% PARAMETERS
freq=5e6;
AngleRange=linspace(-pi,pi,360)'; % propagation angle from -pi to pi (360 degrees);
nModes2Track=8;     % Number of modes to be tracked
legDeg=15;          % Degree of Legendre polynomial expansion - 
saveOn=0;           % Enable saving of the dispersion data
figOn=1; 

%% GEOMETRICAL DISPERSION
[PhaseVelocity]=GeomDispersion(SaristuPlt,AngleRange,freq,nModes2Track,...
    legDeg,saveOn,figOn);
