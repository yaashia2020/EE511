%%initialisation file 
global A B C D K P yub ylb x0 r uub ulb
A=[1 .1; 0 1];
B=[0;.1];
C=[1 0 ;
   0 1];
D=[0;0];
yub = [5 1]';
ylb = [-5 -1]';
uub = 2;
ulb = -2;
Q = diag([1 1]);
R = 1;
[K,P,~] = dlqr(A,B,Q,R);
x0=[0;0];
r=10;
