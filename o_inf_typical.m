m=10; %mass
A=[0 1; 0 0];
B=[0;1/m];
C=[1 0 ;0 0];
D=[0;1];
dt=0.03;
mass_force= ss(A,B,C,D); %continuous 
Dis=c2d(mass_force,dt);
Gu=0;
Gx=[1;0];
%K=[-1 -5];
r=10;
Acl=A-B*K;
Bcl=B*(Gu+K*Gx);
Ccl=(C-D*K);
Dcl=D*(Gu+K*Gx);
xmin=0;
xmax=12;
umax=10;
umin=-1;
H=[1 0;-1 0;0 1;0 -1];
h=[xmax;-xmin;umax;-umin];
e=0.05;
%steady state costraints 
s_s=Polyhedron([zeros(4,1) H*(C*Gx+D*Gu)],(1-e)*h);
s_s.plot
gamma0=H*Ccl;
delta0=H*Dcl;
 