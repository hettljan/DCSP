function [RotMat] = RotateElasticConstants(C,phi,theta,psi)
% This function rotates the symmetric stiffness matrix to the principal
% axis directions using the 3 given Euler's angles phi, theta, psi
% References:
% Auld, Bertram Alexander. Acoustic fields and waves in solids. ????? ???????, 1973.
% Sun, Miao. 2002. Optimal Recovery of Elastic Properties for Anisotropic Materials through Ultrasonic Measurements.
%
% INPUT:
%   C -  Symmetric 6x6 stiffness matrix [GPa]
%   phi  -  First Euler angle [rad]
%   theta-  Second Euler angle [rad]
%   psi  -  Third Euler angle [rad]
% OUTPUT:
%   RotMat - Rotated stiffness matrix [GPa]

%% TRANSFORMATION MATRIX FOR 3D ROTATION OF COORDINATES
R(1,1) = cos(psi)*cos(phi)-cos(theta)*sin(phi)*sin(psi);
R(1,2) = cos(psi)*sin(phi)+cos(theta)*cos(phi)*sin(psi);
R(1,3) = sin(psi)*sin(theta);
R(2,1) = -sin(psi)*cos(phi)-cos(theta)*sin(phi)*cos(psi);
R(2,2) = -sin(psi)*sin(phi)+cos(theta)*cos(phi)*cos(psi);
R(2,3) = cos(psi)*sin(theta);
R(3,1) = sin(theta)*sin(phi);
R(3,2) = -sin(theta)*cos(phi);
R(3,3) = cos(theta);

%% CREATE THE BOND TRANSFORMATION MATRIX M
M(1,1:3) = R(1,1:3).^2;
M(1,2) = R(1,2)^2;
M(1,3) = R(1,3)^2;
M(1,4) = 2*R(1,2)*R(1,3);
M(1,5) = 2*R(1,3)*R(1,1);
M(1,6) = 2*R(1,1)*R(1,2);
M(2,1:3) = R(2,1:3).^2;
M(2,4) = 2*R(2,2)*R(2,3);
M(2,5) = 2*R(2,3)*R(2,1);
M(2,6) = 2*R(2,1)*R(2,2);
M(3,1:3) = R(3,1:3).^2;
M(3,4) = 2*R(3,2)*R(3,3);
M(3,5) = 2*R(3,3)*R(3,1);
M(3,6) = 2*R(3,1)*R(3,2);
M(4,1) = R(2,1)*R(3,1);
M(4,2) = R(2,2)*R(3,2);
M(4,3) = R(2,3)*R(3,3);
M(4,4) = R(2,2)*R(3,3)+R(2,3)*R(3,2);
M(4,5) = R(2,1)*R(3,3)+R(2,3)*R(3,1);
M(4,6) = R(2,2)*R(3,1)+R(2,1)*R(3,2);
M(5,1) = R(3,1)*R(1,1);
M(5,2) = R(3,2)*R(1,2);
M(5,3) = R(3,3)*R(1,3);
M(5,4) = R(1,2)*R(3,3)+R(1,3)*R(3,2);
M(5,5) = R(1,3)*R(3,1)+R(1,1)*R(3,3);
M(5,6) = R(1,1)*R(3,2)+R(1,2)*R(3,1);
M(6,1) = R(1,1)*R(2,1);
M(6,2) = R(1,2)*R(2,2);
M(6,3) = R(1,3)*R(2,3);
M(6,4) = R(1,2)*R(2,3)+R(1,3)*R(2,2);
M(6,5) = R(1,3)*R(2,1)+R(1,1)*R(2,3);
M(6,6) = R(1,1)*R(2,2)+R(1,2)*R(2,1);

%% MULTIPLY THE ORIGINAL STIFFNESS MATRIX BY TRANSFORMATION MATRIX
RotMat = M*C*M.';
% RotMat = Cp.*(abs(Cp)>1e5);
RotMat(abs(RotMat)<1e4)=0;      % replace the very small elements with 0
