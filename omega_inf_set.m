%% the max admissible set
close all
global A B C D K x0 r N Gu Gx H h_1 h_e Acl Bcl uub ulb gamma0_1 delta0_1 h_2 gamma0_2 delta0_2 h_3
% oinf_set_1= set_computation_1(H,h_1,h_e,Acl,Bcl,C,D,Gu,Gx,uub,ulb,gamma0_1,K,delta0_1);

oinf_set_1= Polyhedron([-K (K*Gx+Gu); K -(K*Gx+Gu);H*[zeros(size(C)) (C*Gx+D*Gu)]; gamma0_1 delta0_1],[uub;-ulb;h_e; h_1]);
figure
oinf_set_2 = set_computation_2(H,h_2,h_e,A,B,C,D,Gu,Gx,ulb,gamma0_2,K,delta0_2,oinf_set_1);
figure 
oinf_set_3 = set_computation_3(H,h_3,h_e,A,B,C,D,Gu,Gx,uub,gamma0_2,K,delta0_2,oinf_set_1,oinf_set_2);
figure
plot(oinf_set_1, 'color', 'r', oinf_set_2, 'color', 'b',oinf_set_3, 'color', 'b')
axis equal
