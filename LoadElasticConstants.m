function [m] = load_material_elastic(filename)
% Read the material properties from the specified file.  Currently, we fill
% in the material structure as follows:
%
%   m.rho               -> scalar
%   m.C                 -> 6x6 matrix

fid = fopen(filename, 'rt');


% Read the mass density.
myline = fgetl(fid);
while ( ~strcmp(myline, '$Density') && ischar(myline) )
  myline = fgetl(fid);
end
if ischar(myline)
  m.rho = fscanf(fid, '%g');
  if ( ~strcmp(fgetl(fid), '$EndDensity') )
    error('Material file has a corrupt $Density section.')
  end
end

% Read the (symmetric) stiffness tensor.
frewind(fid);
myline = fgetl(fid);
while ( ~strcmp(myline, '$LinearElasticConstants') && ischar(myline) )
  myline = fgetl(fid);
end
if ischar(myline)
  m.C = zeros(6,6);
  for row = 1:6
  
    % Read the stars.
    fscanf(fid, '%s', row-1);
  
    % Read the non-stars.
    m.C(row,row:end) = fscanf(fid, '%g', 6-row+1);
  
  end

  % Create symmetric stiffness tensor.
  m.C = m.C + triu(m.C,1)';

  fgetl(fid); % dummy to get newline
  myline = fgetl(fid);
  if ~strcmp(myline, '$EndLinearElasticConstants')
    error('Material file has a corrupt $LinearElasticConstants section.')
  end
end

fclose(fid);
