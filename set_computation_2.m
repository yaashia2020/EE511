%% function to solve Oinf 
%region 2
function [final_poly,Gamma,Delta] = set_computation_2(H,h,h_e,A,B,C,D,Gu,Gx,ulb,gamma0,K,delta0,poly_1)
tolerance=10^-3;

final_poly=poly_1;
A_new=final_poly.A;b_new=final_poly.b;
[new]= [gamma0 delta0 h];
Gamma = [];
Delta = [];
i=0;
while 1
    old=new;
    initial_poly=final_poly;
    A_old=A_new;b_old=b_new;
    
    [new]=next_constraint(old(:,1:2),old(:,3),old(:,4),A,B,ulb);
    [~,A_nr,b_nr] = elim_linprog(Polyhedron(new(:,1:3),new(:,4)),old(:,1:3),old(:,4),tolerance);
        if length(b_old)-length(b_nr)>0
            [new]=next_constraint(A_nr(:,1:2),A_nr(:,3),b_nr,A,B,ulb);
        end
    [final_poly,A_new,b_new] = elim_linprog(initial_poly,new(:,1:3),new(:,4),tolerance);
    final_poly.plot;
    axis equal
    pause(1)
    i=i+1
    
%     if i==2
%         [~,~,err] = elim_linprog(Polyhedron(A_new,b_new),initial_poly.A,initial_poly.b,tolerance);
%         if length(initial_poly.b)-length(err)>0
%    A_new=gamma0;b_new=delta0;final_poly=initial_poly;            break;
%         end
%     
%     end
    
    if isempty(A_new)==1
        break
    end 
    
    Gamma = [Gamma; A_new];
    Delta = [Delta; b_new];
    
    A_new=[A_new;A_old];b_new=[b_new;b_old];
end

end
 %% function to compute the next set of constraints 
 %region 2 
function [new]= next_constraint(gamma0, delta0, h, A,B,ulb)
[new]=[gamma0 delta0 h] *[A zeros(2,1) -B*ulb;zeros(1,2) 1 0;zeros(1,3) 1];
end
