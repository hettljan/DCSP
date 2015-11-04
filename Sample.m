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
        function smp = Sample(layers)   % class constructor, the input is a cell containing layers
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
        function plotLayers(smp)
            figure('units','normalized','outerposition',[0 0 1/3 1]) 
            h1=subplot(2,1,1);
            hold on
            ylim([0 smp.hTot]*1e3);
            ylabel('z [mm]')
            h2=subplot(2,1,2);
            hold on
            LayerCols = containers.Map(smp.LayerNames{1},'b'); % constructs an empty Map container for storing colors
            Colors=['r';'g';'m';'c';'y'];     % vector with cute colors
            z=0;
            colInd=1;                           % index of the Colors vector 
            for i=1:smp.nLayers
                if LayerCols.isKey(smp.LayerNames{i})~=1
                    LayerCols(smp.LayerNames{i})=Colors(colInd);
                    colInd=colInd+1;
                end
                axes(h1)
                rectangle('Position',[0,z,1,smp.H(i)*1e3],'FaceColor',LayerCols(smp.LayerNames{i}));
                text(0.425,z+smp.H(i)/2*1e3,smp.LayerNames{i})
                quiver3(h2,0,0,z,cos(smp.Phi(i)),sin(smp.Phi(i)),0,LayerCols(smp.LayerNames{i}),'LineWidth',2)
                z=z+smp.H(i)*1e3;
            end
            axes(h2)
            legend(smp.LayerNames)
            view([50 50])
        end          
    end
    
end