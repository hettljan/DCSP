% This script calculates the dispersion curves for a thin CFRP plate for
% SARISTU project. There are 11 layers with different elastic properties
clearvars

%% DEFINE LAYERING [Material,phi [rad],theta [rad],psi [rad],thickness [m]]
folder='C:\Users\u0088749\Box Sync\MATLAB\Simulation\DCSP\Materials';
Layers{1}=Layer(fullfile(folder,'Aluminum.dat'),0,0,0,3e-3);
UniweavePlate=Sample(Layers);  % create the object representing plate

%% CALCULATE THE DISPERSION CURVES FOR GIVEN ANGLE
propAngle = 0;         % angle of wave propagation with respect to the main in-plane coordinate axis 
nFreqs = 500;     % number of frequency steps
df=1e3;           % frequency step size [Hz]
legDeg=10;        % degree of Legendre polynomial expansion - determines the maximum number of modes 3/2*legDeg
nModes2Track=15;  % number of modes to be tracked
saveOn=0;
figOn=1; 
parOn=0;

% calculate the dispersion curves for the given orientation
[~,~,~]=DispersionCurves(UniweavePlate,propAngle,df,nFreqs,legDeg,...
    nModes2Track,saveOn,figOn,parOn);

%% GEOMETRICAL DISPERSION
freq=220e3;
AngleRange=linspace(-pi,pi,360)'; % propagation angle from -pi to pi (360 degrees);
nModes2Track=5;     % Number of modes to be tracked
figOn=1; 

% calculate the geometrical dispersion with given parameters
[PhaseVelocity]=GeomDispersion(UniweavePlate,AngleRange,freq,nModes2Track,...
    legDeg,saveOn,figOn);

%% DISPLACEMENT PROFILES
nMode=3;
DisplacementProfiles(UniweavePlate,propAngle,freq,nMode,legDeg)