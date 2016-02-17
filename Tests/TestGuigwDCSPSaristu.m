% This script tests whether the GUIGW and DCSP give the same results when using
% the input from GUGIGW to generate the input for DCSP

guigwFile='C:\Users\u0088749\Desktop\GUIGUW\application\CaseStudies\Saristu19Layers.mat';
Layers=Guigw2DCSP(guigwFile);
SartistuPlt=Sample(Layers);  % create the object representing plate

%% DEFINE OTHER PARAMETERS
psip = 0;           % angle [rad] of wave propagation with respect to the main in-plane coordinate axis 
nFreqs = 100;        % number of frequency steps
df=5e3;            % frequency step size [Hz]
legDeg=10;          % degree of Legendre polynomial expansion - determines the maximum number of modes 3/2*legDeg
nModes2Track=8;     % number of modes to be tracked
saveOn=0;
figOn=0;         

%% CALCULATE THE DISPERSION CURVES
[Freq,Wavenumber,PhaseVelocity]=DispersionCurves(SartistuPlt,psip,df,nFreqs,legDeg,nModes2Track,saveOn,figOn,1);

%%
folder='C:\Users\u0088749\KULeuven\PROJECTS\SARISTU\Results\DispersionCurves';
GUIGWFileName='GuiGugw0-500kHz19Layers.mat';
GUIGUWData=load(fullfile(folder,GUIGWFileName));

%% PLOT GUIGW DATA AS BACKGROUND
figure(1)
plot(GUIGUWData.Frequency_Hz(:,1:end-1),GUIGUWData.Phase_Velocity_m_s(:,1:end-1),'k*')
hold on
plot(Freq,PhaseVelocity,'-','Linewidth',2);
ylim([1 8000])
xlabel('Frequency [kHz]');
ylabel('Phase Velocity [m/s]');