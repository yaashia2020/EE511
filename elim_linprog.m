
%% 
function [final_poly,A_new,b_new] = elim_linprog(poly_initial,A,b,tolerance)
A_new=[]; b_new=[];

for i = 1:size(A,1)
    f=-A(i,:);
    C=[poly_initial.A 
        A((1:size(A,1))~=i,:)];
    d=[poly_initial.b 
        b((1:size(A,1))~=i,:)];
    [~,fval,eflag] = linprog(f,C,d);         
        if eflag == -2
            break
        end
        
        if eflag == -3 || -fval > b(i)+tolerance
            A_new = [A_new;A(i,:)];
            b_new = [b_new;b(i)];
        end
end 

final_poly= Polyhedron([poly_initial.A; A_new],[poly_initial.b; b_new]);

end


