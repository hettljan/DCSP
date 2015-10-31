% This script calculates the angular profile of the dispersion curves for a 
% CFRP plate for SARISTU project. There are 11 layers with different elastic properties
clearvars
Psip=linspace(-pi,pi,360)';
nFreqs = 1;       % number of frequency steps
df=350e3;          % frequency step size [Hz]
legDeg=10;        % degree of Legendre polynomial expansion - determines the maximum number of modes 3/2*legDeg
nModes2Track=5;   % number of modes to be tracked
saveOn=0;
figOn=0;
NomMat = {'5Harness' 'NCFBiaxial' 'Uniweave' 'NCFBiaxial' 'NCFBiaxial' 'NCFBiaxial'...
    'NCFBiaxial' 'Uniweave' 'NCFBiaxial' '5Harness' 'Glass'};
Phi=[0 pi/4 0 pi/4 0 pi/2 -pi/4 0 -pi/4 pi/2 0]; % First Euler angle
Theta=[0 0 0 0 0 0 0 0 0 0 0];                       % Second Euler angle
Psi = [0 0 0 0 0 0 0 0 0 0 0];                       % Third Euler angle                
h=0.4801e-3;                      	% Thickness of the ply in [m]
nPlies = size(NomMat,2);         	% Number of plies
H=[ones(nPlies-1,1)*h; 0.127e-3] ;  % vector of ply thicknesses

%% CALCULATION LOOP
Velocity=nan(length(Psip),3);
parfor_progress(length(Psip));
parfor i=1:length(Psip)
    [Freq,Wvn,Vel]=DispersionCurves(NomMat,[Phi;Theta;Psi],H,Psip(i),df,nFreqs,legDeg,nModes2Track,saveOn,figOn);
    Velocity(i,:)=Vel(1:3,end);
    parfor_progress;
end
parfor_progress(0);

%%
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
h=polar(Psip,Velocity(:,1));
set(h,'linewidth',2)

%% 
save Data Psip Velocity