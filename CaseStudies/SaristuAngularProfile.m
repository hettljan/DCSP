% This script calculates the angular profile of the dispersion curves for a 
% CFRP plate for SARISTU project. There are 11 layers with different elastic properties
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

%% DEFINE THE OTHER PARAMETERS
Psip=linspace(-pi,pi,360)'; % propagation angle from -pi to pi (360 degrees)
nFreqs = 1;                 % number of frequency steps
df=350e3;                   % frequency step size [Hz]
legDeg=10;                  % degree of Legendre polynomial expansion - determines the maximum number of modes 3/2*legDeg
nModes2Track=8;             % number of modes to be tracked
saveOn=1;
figOn=1;

%% CALCULATION LOOP
Velocity=nan(length(Psip),nModes2Track);
parfor_progress(length(Psip));
parfor i=1:length(Psip)
    [Freq,Wvn,Vel]=DispersionCurves(SaristuPlt,Psip(i),df,nFreqs,legDeg,nModes2Track,0,0);
    Velocity(i,:)=Vel(1:nModes2Track,end);
    parfor_progress;
end
parfor_progress(0);

%% PLOTTING
if figOn == 1 
    figure
    h=polar(Psip,Velocity(:,1));
    set(h,'linewidth',2)
    hold on
    h=polar(Psip,Velocity(:,2));
    set(h,'linewidth',2)
    h=polar(Psip,Velocity(:,3));
    set(h,'linewidth',2)
    legend('S0','SH0','A0')
    figure
    h=polar(Psip,Velocity(:,3));
    set(h,'linewidth',2)
end

%% SAVING
if saveOn ==1
    save Data Psip Velocity
end