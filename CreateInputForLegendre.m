E=5e9;      % Young's modulus in GPa
mu=0.35;    % Poisson ratio    

C=Eng2CmatrixIsotropic(E,mu);
fileName='Plexiglass.dat';
file=fopen(fileName,'w');
fprintf(file,'$Density\n1210\n$EndDensity\n\n$LinearElasticConstants\n');
fclose(file);
dlmwrite(fileName,C,'delimiter','\t','-append','precision','%.4e')
file=fopen(fileName,'a');
fprintf(file,'$EndLinearElasticConstants');
fclose(file);

