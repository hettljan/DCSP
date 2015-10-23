% Sample is defined a array of structs that contain the following fields:
% material  - Name of the material
% phi       - First Euler angle [rad]
% theta     - Second Euler angle [rad]
% pis       - Third Euler angle [rad]
% h         - Thickness of the layer [m]

%% BACKUP DATA
% NomMat = {'Alamsa' 'Alamsa' 'Alamsa' 'Alamsa' 'Alamsa' 'Alamsa'...
%     'Alamsa' 'Alamsa' 'Alamsa' 'Alamsa' 'Alamsa' 'Alamsa'};
% Phi=[0 -pi/4 0 -pi/4 0 -pi/4 -pi/4 0 -pi/4 0 -pi/4 0]; % First Euler angle
% Theta=[0 0 0 0 0 0 0 0 0 0 0 0];                       % Second Euler angle
% Psi = [0 0 0 0 0 0 0 0 0 0 0 0];                       % Third Euler angle                
% h=0.23e-3;                          % Thickness of the ply in [m]

clear all

%% DEFINE LYAERING
Layers{1}=struct('material','Alamsa','phi',0,'theta',0,'psi',0,'h',0.23e-3);
Layers{2}=struct('material','Alamsa','phi',-pi/4,'theta',0,'psi',0,'h',0.23e-3);
Layers{3}=struct('material','Alamsa','phi',0,'theta',0,'psi',0,'h',0.23e-3);
Layers{4}=struct('material','Alamsa','phi',-pi/4,'theta',0,'psi',0,'h',0.23e-3);
Layers{5}=struct('material','Alamsa','phi',0,'theta',0,'psi',0,'h',0.23e-3);
Layers{6}=struct('material','Alamsa','phi',-pi/4,'theta',0,'psi',0,'h',0.23e-3);
Layers{7}=struct('material','Alamsa','phi',-pi/4,'theta',0,'psi',0,'h',0.23e-3);
Layers{8}=struct('material','Alamsa','phi',0,'theta',0,'psi',0,'h',0.23e-3);
Layers{9}=struct('material','Alamsa','phi',-pi/4,'theta',0,'psi',0,'h',0.23e-3);
Layers{10}=struct('material','Alamsa','phi',0,'theta',0,'psi',0,'h',0.23e-3);
Layers{11}=struct('material','Alamsa','phi',-pi/4,'theta',0,'psi',0,'h',0.23e-3);
Layers{12}=struct('material','Alamsa','phi',0,'theta',0,'psi',0,'h',0.23e-3);

%%
psip = pi/6;            % angle of wave propagation with respect to the main in-plane coordinate axis 
nFreqs = 250;           % number of frequency steps
df=2e3;                 % frequency step size [Hz]
legDeg=10;              % degree of Legendre polynomial expansion - determines the maximum number of modes 3/2*legDeg
nModes2Track=15;        % number of modes to be tracked

%% 
DispersionCurves(Layers,psip,df,nFreqs,legDeg,nModes2Track) 
