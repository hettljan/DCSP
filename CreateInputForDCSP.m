E=70.75e9;
nu=0.34;
rho=2700;
v1=6320;
v2=3130;
fileName='Aluminum.dat';

% convert E and nu to Cmatrix
C=Eng2CmatrixIsotropic(E,nu);

% save the Cmatrix to file
file=fopen(fileName,'w');
fprintf(file,'$Density\n%.2f\n$EndDensity\n\n$LinearElasticConstants\n',rho);
fclose(file);
dlmwrite(fileName,C,'delimiter','\t','-append','precision','%.4e')
file=fopen(fileName,'a');
fprintf(file,'$EndLinearElasticConstants');
fclose(file);

