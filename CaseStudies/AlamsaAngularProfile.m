% This script calculates the angular profile of the dispersion curves for a 
% CFRP plate for SARISTU project. There are 11 layers with different elastic properties
clearvars

%% DEFINE LAYERING [Material,phi [rad],theta [rad],psi [rad],thickness [m]]
folder='C:\Users\u0088749\Box Sync\MATLAB\Simulation\DCSP\Materials';
Layers{1}=Layer(fullfile(folder,'Alamsa.dat'),0,0,0,0.23e-3);
Layers{2}=Layer(fullfile(folder,'Alamsa.dat'),-pi/4,0,0,0.23e-3);
Layers{3}=Layer(fullfile(folder,'Alamsa.dat'),0,0,0,0.23e-3);
Layers{4}=Layer(fullfile(folder,'Alamsa.dat'),-pi/4,0,0,0.23e-3);
Layers{5}=Layer(fullfile(folder,'Alamsa.dat'),0,0,0,0.23e-3);
Layers{6}=Layer(fullfile(folder,'Alamsa.dat'),-pi/4,0,0,0.23e-3);
Layers{7}=Layer(fullfile(folder,'Alamsa.dat'),-pi/4,0,0,0.23e-3);
Layers{8}=Layer(fullfile(folder,'Alamsa.dat'),0,0,0,0.23e-3);
Layers{9}=Layer(fullfile(folder,'Alamsa.dat'),-pi/4,0,0,0.23e-3);
Layers{10}=Layer(fullfile(folder,'Alamsa.dat'),0,0,0,0.23e-3);
Layers{11}=Layer(fullfile(folder,'Alamsa.dat'),-pi/4,0,0,0.23e-3);
Layers{12}=Layer(fullfile(folder,'Alamsa.dat'),0,0,0,0.23e-3);
AlamsaPlt=Sample(Layers);  % create the object representing plate

%% DEFINE THE OTHER PARAMETERS
AngleRange=linspace(-pi,pi,360)'; % propagation angle from -pi to pi (360 degrees)
freq=283e3;                 % frequency step size [Hz]
legDeg=5;                  % degree of Legendre polynomial expansion - determines the maximum number of modes 3/2*legDeg
nModes2Track=6;             % number of modes to be tracked
saveOn=0;
figOn=1;

%% CALCULATION LOOP
[PhaseVelocity]=GeomDispersion(AlamsaPlt,AngleRange,freq,nModes2Track,...
    legDeg,saveOn,figOn);

%% SORT THE ANGULAR PROFILES
PhaseVelocity=DispersionCurveSorting(AngleRange,PhaseVelocity',nModes2Track);
PhaseVelocity=PhaseVelocity';

%% PLOTTING
if figOn == 1 
    figure
    h=polar(AngleRange,PhaseVelocity(:,2));
    set(h,'linewidth',2)
    hold on
    h=polar(AngleRange,PhaseVelocity(:,3));
    set(h,'linewidth',2)
    h=polar(AngleRange,PhaseVelocity(:,5));
    set(h,'linewidth',2)
    legend('S0','SH0','A0')
end