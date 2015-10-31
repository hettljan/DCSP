function SortedWvns=DispersionCurveSorting(Freq,Wavenumbers,nModes,varargin)
% This function sorts the dispersion curves that are generated using the
% LambAnisotropic3DLegendre.m. The output is the sorted matrix in which
% each row correspond to one mode
% INPUT:
%   Freq        -   Frequency vector in [Hz]
%   Wavenumbers -   Matrix with dispersion curves ordered row wise
%   nModes      -   Number of modes to be tracked
% OTPIONAL:
%   figOn       -   Enable plotting
% OUTPUT:
%   SortedWvns  -   Matrix of sorted modes ordered row-wise

% fprintf('\n### STARTING DISPERSION CURVES SORTING PROGRAM ###\n\n');

%% INPUT PARSING
numvarargs = length(varargin);
if numvarargs > 1
    error('myfuns:somefun2Alt:TooManyInputs', ...
        'requires at most 1 optional inputs');
end
optargs = {0};
optargs(1:numvarargs) = varargin;
[figOn] = optargs{:};

%% CALCULATION
for idx=1:nModes;
    for i=2:size(Wavenumbers,2)-1
        if isnan(Wavenumbers(idx,i))~=1
            y0=Wavenumbers(idx,i-1);            % previous wavenumber    
            y1=Wavenumbers(idx,i);              % current wavenumber
            p=polyfit(Freq(i-1:i)',[y0 y1],1);  % linear fit to get extrapolation
            y2t=polyval(p,Freq(i+1));           % extrapolated value of y2
            y2=Wavenumbers(idx:end,i+1);        % real available wavenumbers
            y2=y2(isnan(y2)==0);                % real available non-nan wavenumbers
            Diff=abs(y2-y2t);                   % calculate the absolute difference
            [~,minIdx]=min(Diff);               % select the minimum difference between extrpolated and real
            minIdx=minIdx+idx-1;
            if minIdx ~= idx
                Wavenumbers([idx minIdx],i+1)=Wavenumbers([minIdx idx],i+1);    % flip the values
            end
        end
    end
end
SortedWvns=Wavenumbers(1:nModes,:);     % rearrange to the output array

%% PLOTTING
if figOn == 1 
    figure
    plot(SortedWvns','*');
end