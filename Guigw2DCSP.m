function Layers=Guigw2DCSP(guigwFile)
% This function converts the plate/layer generated and stored by GUIGUW
% program to the input for the DCSP program. 
% INPUT:
%   guigwFile   -   String containing absolute path the GUIGUW file
% OPTIONAL:
% OUTPUT:
%   Layers      -   Cell strcuture representing the layup of the sample,                                
%                   each element of the cell represents one layer in the 
%                   real structure. This layer contains the following fields:
%                   Material - string with material name
%                   phi,psi,theta - Euler angles
%                   h   - layer thickness in [m]
%                   rho - layer density in [kg/m^3]
%                   C - 6x6 elasticity C matrix in [Pa]

%% LOAD THE GUIGUW FILE
GuigwData=load(guigwFile);

%% CREATE THE OUTPUT STRUCTURE
nLayers=size(GuigwData.layers,1);
Layers=cell(1,nLayers);

%% PARSING OF THE ORIGINAL DATA
folder='\\plato.kulak.be\groupdata\acoustics\Jan_Hettler\MATLAB\Simulation\DCSP\Materials';
for i=1:nLayers
  C=[GuigwData.C11(i),GuigwData.C12(i),GuigwData.C13(i),0,0,0;...
    GuigwData.C12(i),GuigwData.C22(i),GuigwData.C23(i),0,0,0;...
    GuigwData.C13(i),GuigwData.C23(i),GuigwData.C33(i),0,0,0;
    0,0,0,GuigwData.C44(i),0,0;
    0,0,0,0,GuigwData.C55(i),0;
    0,0,0,0,0,GuigwData.C66(i)];
  
    fileName=fullfile(folder,strcat(GuigwData.layers{i},'.dat'));
    file=fopen(fileName,'w');
    fprintf(file,'$Density\n%d\n$EndDensity\n\n$LinearElasticConstants\n',GuigwData.rho(i));
    fclose(file);
    dlmwrite(fileName,C,'delimiter','\t','-append','precision','%.4e')
    file=fopen(fileName,'a');
    fprintf(file,'$EndLinearElasticConstants');
    fclose(file);
    
    currLayer=Layer(GuigwData.layers{i},deg2rad(GuigwData.teta(i)),0,0,GuigwData.thick(i));        
    Layers{i}=currLayer;
end
% delete(fileName)
                    