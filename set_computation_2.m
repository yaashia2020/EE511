%% function to solve Oinf 
%region 2
function final_poly = set_computation_2(H,h,h_e,A,B,C,D,Gu,Gx,ulb,gamma0,K,delta0,oinf_set_1)
tolerance=10^-3;
s_s=Polyhedron([-K +(K*Gx+Gu);H*[zeros(size(C)) (C*Gx+D*Gu)]],[ulb;h_e]);
% s_s.plot;
% axis equal
% pause(1)

add_A=oinf_set_1.A(3:end,:);
add_b=oinf_set_1.b(3:end,:);


[final_poly,A_new,b_new] = elim_linprog(s_s,add_A,add_b,tolerance);
final_poly.plot;
axis equal
new=[A_new b_new];
i=0;
while 1
    old=new;
    initial_poly=final_poly;
    [new]=next_constraint(old(:,1:2),old(:,3),old(:,4),A,B,ulb);
    [final_poly,A_new,b_new] = elim_linprog(initial_poly,new(:,1:3),new(:,4),tolerance);
    final_poly.plot;
    axis equal
    pause(1)
    i=i+1
    
    if i==1
        [~,~,err] = elim_linprog(Polyhedron(A_new,b_new),initial_poly.A,initial_poly.b,tolerance);
        if isempty(err)
        else
            break;
        end
    
    end
    
    if isempty(A_new)==1
        break
    end 
    
end

end
 %% function to compute the next set of constraints 
 %region 2 
function [new]= next_constraint(gamma0, delta0, h, A,B,ulb)
[new]=[gamma0 delta0 h] *[A zeros(2,1) -B*ulb;zeros(1,2) 1 0;zeros(1,3) 1];
end
