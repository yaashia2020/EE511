%% command governor
function v = reference_gov(final_poly,Xk,r)
constraints=final_poly.A;

Gamma=constraints(:,1:length(Xk));
Delta=constraints(:,length(Xk)+1:end);
H=eye(size(r,1));
f=-2*r;
b=h-Gamma*Xk;
v=quadprog(H,f,Delta,b);
h=final_poly.b;
end