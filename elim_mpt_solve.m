A=[0 1; -9 0; 0 -4; 1 0];
b=[6;5;1;0];
initial=Polyhedron(A,b);

initial.plot
C=[0 1];d=5;
lp.f=-C;lp.A=A;lp.b=b;
result=mpt_solve(lp);
tolerance=10^-6;
if -lp.f*result.xopt>d + tolerance
    A=[A;C];b=[b;d];
end 
final=Polyhedron(A,b)
final.plot
