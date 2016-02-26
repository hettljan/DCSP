function PhaseVelocity=GeomDispersion(TestSample,AngleRange,freq,varargin)
% Calculates the geometrical dispersion for a layered anisotropic materials
% using the Legendre and Laguerre polynomial approach
% Originally supplied by O. Bou-Matar, Lille
%
% INPUT: 
%   Sample        - Object with description of the layering and layer
%                   properties, see Sample.mat and Layer.mat for details 
%   AngleRange    - Vector with propagation ngles of wave propagation with
%                   respect to the main in-plane coordinate axis [rad]
%   freq          - Frequency to analyze
% OPTIONAL:
%   nModes2Track    - Number of modes to be tracked
%   legDeg          - Degree of Legendre polynomial expansion - 
%                     determines the maximum number of modes 3/2*legDeg
%   saveOn          - Enable saving of the dispersion data
%   figOn           - Enable plotting
% OUTPUT:
%   PhaseVelocity   - Matrix with phase velocities for identified modes 
%                     in [m/s]

%% INPUTS PARSING
numvarargs = length(varargin);
if numvarargs > 4
    error('myfuns:somefun2Alt:TooManyInputs', ...
        'requires at most 4 optional inputs');
end
optargs = {3,10,0,0};
optargs(1:numvarargs) = varargin;
[nModes2Track,legDeg,saveOn,figOn] = optargs{:};

%% CALCULATION
PhaseVelocity=nan(length(AngleRange),nModes2Track);
parfor_progress(size(AngleRange,1));
parfor i=1:length(AngleRange)
    [~,~,Vel]=DispersionCurves(TestSample,AngleRange(i),freq,1,legDeg,...
        nModes2Track,saveOn,0,0);
    PhaseVelocity(i,:)=Vel(1:nModes2Track,end);
    parfor_progress;
end
parfor_progress(0);

%% PLOTTING
if figOn == 1 
    figure
    for i=1:size(PhaseVelocity,2)
        polar(AngleRange,PhaseVelocity(:,i));
        hold on
    end
end

%% SAVING
if saveOn ==1
    save Data AngleRange PhaseVelocity
end
