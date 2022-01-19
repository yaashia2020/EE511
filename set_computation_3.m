%% function to solve Oinf 
%region 3
function final_poly = set_computation_3(H,h,h_e,A,B,C,D,Gu,Gx,uub,gamma0,K,delta0,oinf_set_1,oinf_set_2)
tolerance=10^-3;
s_s=Polyhedron([ K -(K*Gx+Gu);H*[zeros(size(C)) (C*Gx+D*Gu)]],[-uub;h_e]);
% s_s.plot;
% axis equal
% pause(1)
add_A=[oinf_set_1.A(3:end,:);oinf_set_2.A(2:end,:)];
add_b=[oinf_set_1.b(3:end,:);oinf_set_2.b(2:end,:)];
% [new]= next_constraint(gamma0,delta0,h,A,B,uub);
[final_poly,A_new,b_new] = elim_linprog(s_s,add_A,add_b,tolerance);
final_poly.plot;
axis equal
new=[A_new b_new];
i=0;
while 1
    old=new;
    initial_poly=final_poly;
    [new]=next_constraint(old(:,1:2),old(:,3),old(:,4),A,B,uub);
    [final_poly,A_new,b_new] = elim_linprog(initial_poly,new(:,1:3),new(:,4),tolerance);
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
 %region 2 
function [new]= next_constraint(gamma0, delta0, h, A,B,uub)
[new]=[gamma0 delta0 h] *[A zeros(2,1) -B*uub;zeros(1,2) 1 0;zeros(1,3) 1];
end
