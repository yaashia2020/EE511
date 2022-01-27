%% function to solve Oinf 
%region 1 
function [final_poly,Gamma, Delta] = set_computation_1(H,h,h_e,Acl,Bcl,C,D,Gu,Gx,uub,ulb,gamma0,K,delta0,poly_1)
tolerance=10^-3;


% the initial polyhedron - steady state and region constraints
% initialising the process
final_poly=poly_1;
A_new=final_poly.A;b_new=final_poly.b;
[new]= [gamma0 delta0];


i=0; % no of steps

Gamma = [];
Delta = [];
b_nr=h;
while 1
    % initialising the loop
    old=new;
    initial_poly=final_poly;
    A_old=A_new;b_old=b_new;
    
    
    
    %getting the final poly and the new constraints
    [new]=next_constraint(old(:,1:2),old(:,3),Acl,Bcl);
    
%   [~,A_nr,b_nr] = elim_linprog(Polyhedron(new,h),old,h,tolerance);
    [~,A_nr,b_nr] = elim_linprog(Polyhedron(new,b_nr),old,b_nr,tolerance);
        if length(h)-length(b_nr)>0
            [new]=next_constraint(A_nr(:,1:2),A_nr(:,3),Acl,Bcl);
        end
    [final_poly,A_new,b_new] = elim_linprog(initial_poly,new,h,tolerance);
    
    
    final_poly.plot;
    axis equal
    pause(1)
    i=i+1
    
    if isempty(A_new)==1
        break
    end
    Gamma = [Gamma; A_new];
    Delta = [Delta; b_new];
    
    A_new=[A_new;A_old];b_new=[b_new;b_old];
end



end
 %% function to compute the next set of constraints 
 
function [new]= next_constraint(gamma0, delta0,Acl,Bcl)
[new]=[gamma0 delta0] *[Acl Bcl;zeros(1,2) 1];
end
