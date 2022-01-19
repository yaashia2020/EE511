global final_poly
A=[1 .1; 0 1];
B=[0;.1];
C=[1 0 ;
   0 1;
   0 0];
D=[0;
   0; 
   1];

Gu=0;
Gx=[1;0];

Q = diag([1 1]);
R = 1;
[K,P,~] = dlqr(A,B,Q,R);

Acl=A-B*K;
Bcl=B*(Gu+K*Gx);
Ccl=(C-D*K);
Dcl=D*(Gu+K*Gx);

yub = [5 1 2]';
ylb = [-5 -1 -2]';
H=[eye(length(D));-eye(length(D))];
h=[yub;-ylb];

gamma0=H*Ccl;
delta0=H*Dcl;

e=0.05;
tolerance=10^-3;

%steady state costraints + constraints on x0
s_s=Polyhedron([H*[zeros(size(C)) (C*Gx+D*Gu)]; gamma0 delta0],[(1-e)*h; h]);
s_s.plot;
axis equal
pause(1)
[new]= next_constraint(gamma0, delta0,Acl,Bcl);
[final_poly,A_new,b_new] = elim_linprog(s_s,new,h,tolerance);
final_poly.plot;
axis equal
i=0;
while 1
    old=new;
    initial_poly=final_poly;
    [new]=next_constraint(old(:,1:2),old(:,3),Acl,Bcl);
    [final_poly,A_new,b_new] = elim_linprog(initial_poly,new,h,tolerance);
    final_poly.plot;
    axis equal
    pause(1)
    i=i+1
    if isempty(A_new)==1
        break
    end 
    
end
 %% function to compute the next set of constraints 
function [new]= next_constraint(gamma0, delta0,Acl,Bcl)
[new]=[gamma0 delta0] *[Acl Bcl;zeros(1,2) 1];
end
