function [material] = LoadElasticConstants(filename)
% Read the material properties from the specified file.  Currently, we fill
% in the material structure as follows:
%   m.rho               -> scalar
%   m.C                 -> 6x6 matrix
% INPUT:
%   filename -  Path to the file with elastic properties
% OUTPUT:
%   material -  Structure that contains the material density and elastic
%               matrix

%% OPEN THE FILE
fid = fopen(filename, 'rt');

%% READ THE MASS DENSITY
myline = fgetl(fid);
while ( ~strcmp(myline, '$Density') && ischar(myline) )
  myline = fgetl(fid);
end
if ischar(myline)
  material.rho = fscanf(fid, '%g');
  if ( ~strcmp(fgetl(fid), '$EndDensity') )
    error('Material file has a corrupt $Density section.')
  end
end

%% READ THE (SYMMETRIC) STIFFNESS TENSOR
frewind(fid);
myline = fgetl(fid);
while (~strcmp(myline, '$LinearElasticConstants') && ischar(myline))
  myline = fgetl(fid);
end
if ischar(myline)
  material.C = zeros(6,6);
  for row = 1:6
    fscanf(fid, '%s', row-1);                           % Read the stars.
    material.C(row,row:end) = fscanf(fid, '%g', 6-row+1);      % Read the non-stars.
  end
  material.C = material.C + triu(material.C,1)';     % populate the stiffness tensor with symmetric values
  fgetl(fid);                   % dummy to get newline
  myline = fgetl(fid);
  if strcmp(myline, '$EndLinearElasticConstants')~=1
    error('Material file has a corrupt $LinearElasticConstants section.')
  end
end
fclose(fid);
