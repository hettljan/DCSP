classdef Layer
    properties
        Material  % string with the name of the material
        phi       % 1st Euler angle [rad]
        theta     % second Euler angle [rad]
        psi       % 3rd Euler angle [rad]
        h         % layer thickness [m]
        C         % C matrix in Voigt notation
        rho       % material density in [kg/m^3]
    end
    methods
        function lyr = Layer(material,phi,theta,psi,h)
            lyr.Material=material;
            lyr.phi=phi;
            lyr.theta=theta;
            lyr.psi=psi;
            lyr.h=h;
            % READ THE MASS DENSITY
            folder='\\plato.kulak.be\groupdata\acoustics\Jan_Hettler\MATLAB\Simulation\DCSP\Materials';
            fid = fopen(fullfile(folder,strcat(material,'.dat')), 'rt');
            myline = fgetl(fid);
            while ( ~strcmp(myline, '$Density') && ischar(myline) )
                myline = fgetl(fid);
            end
            if ischar(myline)
                lyr.rho = fscanf(fid, '%g');
                if ( ~strcmp(fgetl(fid), '$EndDensity') )
                    error('Material file has a corrupt $Density section.')
                end
            end
            % READ THE (SYMMETRIC) STIFFNESS TENSOR
            frewind(fid);
            myline = fgetl(fid);
            while (~strcmp(myline, '$LinearElasticConstants') && ischar(myline))
                myline = fgetl(fid);
            end
            if ischar(myline)
                Cpom = zeros(6,6);
                for row = 1:6
                    fscanf(fid, '%s', row-1);                                  % Read the stars.
                    Cpom(row,row:end) = fscanf(fid, '%g', 6-row+1);      % Read the non-stars.
                end
                lyr.C=Cpom+triu(Cpom,1)';             % populate the stiffness tensor with symmetric values
                fgetl(fid);                 % dummy to get newline
                myline = fgetl(fid);
                if strcmp(myline, '$EndLinearElasticConstants')~=1
                    error('Material file has a corrupt $LinearElasticConstants section.')
                end
            end
        fclose(fid);
      end
   end
end