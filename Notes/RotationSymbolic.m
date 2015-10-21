syms a b c
D=[cos(a) sin(a) 0;... 
    -sin(a) cos(a) 0;...
    0 0 1];
C=[1 0 0;...
    0 cos(b) sin(b);...
    0 -sin(b) cos(b)];
B=[cos(c) sin(c) 0;...
    -sin(c) cos(c) 0;...
    0 0 1];

%%
Rsym=B*C*D

%%
R(1,1) = cos(c)*cos(a)-cos(b)*sin(a)*sin(c);
R(1,2) = cos(c)*sin(a)+cos(b)*cos(a)*sin(c);
R(1,3) = sin(c)*sin(b);
R(2,1) = -sin(c)*cos(a)-cos(b)*sin(a)*cos(c);
R(2,2) = -sin(c)*sin(a)+cos(b)*cos(a)*cos(c);
R(2,3) = cos(c)*sin(b);
R(3,1) = sin(b)*sin(a);
R(3,2) = -sin(b)*cos(a);
R(3,3) = cos(b);
R
%%
isequal(Rsym,R)