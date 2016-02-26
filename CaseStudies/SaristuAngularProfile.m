% This script calculates the angular profile of the dispersion curves for a 
% CFRP plate for SARISTU project. There are 11 layers with different elastic properties
clearvars

%% DEFINE LAYERING [Material,phi [rad],theta [rad],psi [rad],thickness [m]]
folder='C:\Users\u0088749\Box Sync\MATLAB\Simulation\DCSP\Materials';
Layers{1}=Layer(fullfile(folder,'5Harness.dat'),0,0,0,0.4801e-3);
Layers{2}=Layer(fullfile(folder,'NCFBiaxial.dat'),pi/4,0,0,0.4801e-3);
Layers{3}=Layer(fullfile(folder,'Uniweave.dat'),0,0,0,0.4801e-3);
Layers{4}=Layer(fullfile(folder,'NCFBiaxial.dat'),pi/4,0,0,0.4801e-3);
Layers{5}=Layer(fullfile(folder,'NCFBiaxial.dat'),0,0,0,0.4801e-3);
Layers{6}=Layer(fullfile(folder,'NCFBiaxial.dat'),pi/2,0,0,0.4801e-3);
Layers{7}=Layer(fullfile(folder,'NCFBiaxial.dat'),-pi/4,0,0,0.4801e-3);
Layers{8}=Layer(fullfile(folder,'Uniweave.dat'),0,0,0,0.4801e-3);
Layers{9}=Layer(fullfile(folder,'NCFBiaxial.dat'),-pi/4,0,0,0.4801e-3);
Layers{10}=Layer(fullfile(folder,'5Harness.dat'),pi/2,0,0,0.4801e-3);
Layers{11}=Layer(fullfile(folder,'Glass.dat'),0,0,0,0.127e-3);
SaristuPlt=Sample(Layers);  % create the object representing plate

%% DEFINE THE OTHER PARAMETERS
AngleRange=linspace(-pi,pi,360)'; % propagation angle from -pi to pi (360 degrees)
nFreqs = 1;                 % number of frequency steps
df=350e3;                   % frequency step size [Hz]
legDeg=10;                  % degree of Legendre polynomial expansion - determines the maximum number of modes 3/2*legDeg
nModes2Track=8;             % number of modes to be tracked
saveOn=0;
figOn=1;

%%
tic
PhaseVelocity=GeomDispersion(SaristuPlt,AngleRange,nFreqs*df,nModes2Track,...
    legDeg,saveOn,figOn);
toc
