global MPTOPTIONS

if isempty(MPTOPTIONS)
	MPTOPTIONS = mptopt('lpsolver', 'lcp', 'abs_tol', 1e-8);
end

%infeasible set
% A=[0 1; -9 0; 0 -4; 1 0];
% b=[6;5;1;0];
% 
% C=[0 1;1 0];d=[5;4];
A=[0 1; 0 -4; -1 0];
b=[6;1;0];

C=[0 1;1 0];d=[5;1];
tolerance=10^-6;
for i=1: size(C,1)
    lp.f=-C(i,:);
    lp.A=A;
    lp.b=b;
    result=mpt_solve(lp);
    [X,~,flag] = linprog(lp.f,lp.A,lp.b);
    if result.exitflag==MPTOPTIONS.OK 
    
        if -lp.f*result.xopt>d(i) + tolerance
            A=[A;C(i,:)];b=[b;d(i)];
        end 
        elseif result.exitflag==MPTOPTIONS.UNBOUNDED
            A=[A;C(i,:)];b=[b;d(i)];
     end 
end


final=Polyhedron(A,b);
final.plot
A=[A;C];b=[b;d];
poly=Polyhedron(A,b)
E=poly.minHRep()
E.plot