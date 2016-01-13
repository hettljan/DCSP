% This script calculates the angular profile of the dispersion curves for a 
% CFRP plate for SARISTU project. There are 11 layers with different elastic properties
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

%% DEFINE THE OTHER PARAMETERS
Psip=linspace(-pi,pi,180)'; % propagation angle from -pi to pi (360 degrees)
nFreqs = 1;                 % number of frequency steps
df=283e3;                   % frequency step size [Hz]
legDeg=10;                  % degree of Legendre polynomial expansion - determines the maximum number of modes 3/2*legDeg
nModes2Track=5;             % number of modes to be tracked
saveOn=1;
figOn=1;
parOn=0;

%% CALCULATION LOOP
Velocity=nan(length(Psip),nModes2Track);
parfor_progress(length(Psip));
parfor i=1:length(Psip)
    [Freq,Wvn,Vel]=DispersionCurves(AlamsaPlt,Psip(i),df,nFreqs,legDeg,nModes2Track,0,0,parOn);
    Velocity(i,:)=Vel(1:nModes2Track,end);
    parfor_progress;
end
parfor_progress(0);

%% SORT THE ANGULAR PROFILES
Velocity=DispersionCurveSorting(Psip,Velocity',nModes2Track);
Velocity=Velocity';

%% PLOTTING
if figOn == 1 
    figure
    h=polar(Psip,Velocity(:,2));
    set(h,'linewidth',2)
    hold on
    h=polar(Psip,Velocity(:,3));
    set(h,'linewidth',2)
    h=polar(Psip,Velocity(:,5));
    set(h,'linewidth',2)
    legend('S0','SH0','A0')
%     figure
%     h=polar(Psip,Velocity(:,3));
%     set(h,'linewidth',2)
end

%% SAVING
if saveOn ==1
    save Data Psip Velocity
end