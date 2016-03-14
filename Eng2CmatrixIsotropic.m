function C=Eng2CmatrixIsotropic(E,mu)
% converts engineering elastic properties Young modulus and Poissons ratio
%to the full elastic matrix of the isotropic material
% INPUT:
%   E   -   Young's modulus [Pa]
%   mu  -   Poisson's ratio [-]
% OUTPUT:
%   C   -   Stiffness matrix

%% CALCULATE REMAINING LIN. DEP. ELASTIC  CONSTANTS
a=E/((1+mu)*(1-2*mu));

%%
C11=1-mu;
C22=1-mu;
C33=1-mu;
C44=(1-2*mu)/2;
C55=(1-2*mu)/2;
C66=(1-2*mu)/2;
C12=mu;
C13=mu;
C21=mu;
C23=mu;
C31=mu;
C32=mu;

%% CREATE THE MATRIX
C=a*[C11 C12 C13 0 0 0;...
    C21 C22 C23 0 0 0;...
    C31 C32 C33 0 0 0;...
    0 0 0 C44 0 0;...
    0 0 0 0 C55 0;...
    0 0 0 0 0 C66]