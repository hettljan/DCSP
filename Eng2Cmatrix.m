function CMatrix=Eng2Cmatrix(E1,E2,E3,G12,G23,G31,v21,v31,v32)
% converts engineering elastic properties to the full 
% elastic matrix of the orthotropic material (plate)
% expects the 9 independent elastic constants in the following order
% (E1,E2,E3,G12,G23,G13,v12,v23,v13) 
% [Eij]=[Gij] = GPa
% returns the C matrix

%% INPUT CHECKING
if nargin < 9
    disp('WRONG NUMBER OF ELASTIC CONSTANTS')
end

%% CALCULATE REMAINING LIN. DEP. ELASTIC  CONSTANTS
v12=v21*E1/E2;
v13=v31*E1/E3;
v23=v32*E2/E3;

%% CONVERSION OF THE VALUES
Delta=(1-v12*v21-v23*v32-v31*v13-2*v12*v23*v31)/(E1*E2*E3);
C11=(1-v23*v32)/(E2*E3*Delta);
C22=(1-v31*v13)/(E1*E3*Delta);
C33=(1-v12*v21)/(E1*E2*Delta);
C44=G23;
C55=G31;
C66=G12;
C12=(v21+v31*v23)/(E2*E3*Delta);
C13=(v31+v21*v32)/(E2*E3*Delta);
C21=(v12+v13*v32)/(E1*E3*Delta);
C23=(v32+v31*v12)/(E1*E3*Delta);
C31=(v13+v12*v23)/(E1*E2*Delta);
C32=(v23+v13*v21)/(E1*E2*Delta);
CMatrix=[C11 C12 C13 0 0 0;...
    C21 C22 C23 0 0 0;...
    C31 C32 C33 0 0 0;...
    0 0 0 C44 0 0;...
    0 0 0 0 C55 0;...
    0 0 0 0 0 C66];