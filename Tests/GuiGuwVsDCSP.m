% Comparison of the results of GUIGUW (SAFE) code and the code based on 
% Legendre poylnomial expansion. The dispersion is measured in the 0 deg
% propagation direction.

clearvars

%% LOAD DATA
folder='C:\Users\u0088749\Box Sync\Projects\ALAMSA\Results\DispersionCurves\IAI-REF-A';
GUIGWFileName='Guigw0-500kHz.mat';
DCSPFileName='DCSP0-500kHz.mat';
GUIGUWData=load(fullfile(folder,GUIGWFileName));
DCSPData=load(fullfile(folder,DCSPFileName));

%% PLOT GUIGW DATA AS BACKGROUND
figure(1)
plot(GUIGUWData.Frequency_Hz(:,1:end-1),GUIGUWData.Phase_Velocity_m_s(:,1:end-1),'k*')
hold on
plot(DCSPData.Freq*1e3,DCSPData.Velocity,'-','Linewidth',2);
ylim([1 8000])
xlabel('Frequency [kHz]');
ylabel('Phase Velocity [m/s]');
