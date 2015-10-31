classdef Sample
    properties
        Layers      % cell structure containing Layer objects, each representing one factual layer
        Phi         % vector with 1st Euler angle [rad]
        Theta       % vector with second Euler angle [rad]
        Psi         % vector with 3rd Euler angle [rad]
        C           % 3d array with unrotated C matrix in Voigt notation
        Rho         % vector with material densities in [kg/m^3]
        H           % vector will layer thicknesses [m]
        hTot        % total thickness of the sample (sum of layer thicknesses)
        nLayers     % the total number of layers in the lamina
        LayerNames  % cell with the layer names
        
    end
    methods
        function smp = Sample(layers)
            smp.Layers=layers;
            smp.nLayers=length(layers);
            smp.H=nan(length(layers),1);
            smp.Phi=nan(length(layers),1);
            smp.Theta=nan(length(layers),1);
            smp.Psi=nan(length(layers),1);
            smp.Rho=nan(length(layers),1);
            smp.C=nan(6,6,length(layers));
            smp.LayerNames=cell(1,length(layers));
            for i=1:length(layers)
                smp.H(i)=layers{i}.h;
                smp.Phi(i)=layers{i}.phi;
                smp.Theta(i)=layers{i}.theta;
                smp.Psi(i)=layers{i}.psi;
                smp.Rho(i)=layers{i}.rho;
                smp.C(:,:,i)=layers{i}.C;
                smp.LayerNames{i}=layers{i}.Material;
            end
            smp.hTot=sum(smp.H); 
        end
    end
end