%% function to solve Oinf 
%region 1 
function final_poly = set_computation_1(H,h,h_e,Acl,Bcl,C,D,Gu,Gx,uub,ulb,gamma0,K,delta0)
tolerance=10^-3;

s_s=Polyhedron([-K (K*Gx+Gu); K -(K*Gx+Gu);H*[zeros(size(C)) (C*Gx+D*Gu)]; gamma0 delta0],[uub;-ulb;h_e; h]);
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



end
 %% function to compute the next set of constraints 
 
function [new]= next_constraint(gamma0, delta0,Acl,Bcl)
[new]=[gamma0 delta0] *[Acl Bcl;zeros(1,2) 1];
end