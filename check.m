%% check 
%fix v
global final_poly
v=0;
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

x2_sweep_in=[];x2_sweep_out=[];

eq_pt=[v;0;0];

y=[];x=[];
%% set x2 and sweep
for x2=-1.5:0.1:1.5
   
    [inside,outside]= check_poly_x2sweep(final_poly,x2,v);
    x2_sweep_in=[x2_sweep_in [inside]];
    x2_sweep_out=[x2_sweep_out [outside]];
end

%% checking with the system dynamics
i=1;
x(:,1)=x2_sweep_in(1:2,1);
while 1
    x(:,i+1)=Acl*x(:,i)+Bcl*v; 
    out=Ccl*(x(:,i))+Dcl*v;
    y=[y [out]];
        
    if norm(y(:,i)-eq_pt)<10^-1
            break 
    end
        
    i=i+1;
    if H*y(i)>h
        disp('error')
        break
    end
end

%% 

function [inside,outside]= check_poly_x2sweep(final_poly,x2,v)
inside=[];outside=[];
    for x1=-6:0.01:6
         if final_poly.contains([x1;x2;v])==1
             inside=[x1;x2;v];
             outside=[x1-0.01;x2;v];
             break
         end    
    end 
    if x1<6
        for x1_new=x1:0.01:6 
            if final_poly.contains([x1_new;x2;v])==0
                inside=[inside [x1_new-0.01;x2;v]];
                outside=[outside [x1_new;x2;v]];
                break
            end
        end
    end 
% % i=final_poly.contains([0;0;0]);
end