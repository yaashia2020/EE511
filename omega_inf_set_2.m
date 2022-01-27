close all
% run(initialise)
% run(system_parameters)
global A B C D K x0 r N Gu Gx H h_1 h_e Acl Bcl 
global uub ulb gamma0_1 delta0_1 h_2 gamma0_2 delta0_2 h_3

i=0;
while 1
    if i==0
        s_s1=Polyhedron([-K (K*Gx+Gu); K -(K*Gx+Gu);H*[zeros(size(C)) (C*Gx+D*Gu)]],[uub;-ulb;h_e]); 
        s_s1=s_s1.minHRep();
        [oinf_set_1,Gamma1_1,Delta1_1]= set_computation_1(H,h_1,h_e,Acl,Bcl,C,D,Gu,Gx,uub,ulb,gamma0_1,K,delta0_1,s_s1);
        figure
        s_s2=Polyhedron([-K +(K*Gx+Gu);H*[zeros(size(C)) (C*Gx+D*Gu)]],[ulb;h_e]);
        s_s2=s_s2.minHRep();
        [oinf_set_2,Gamma2_1,Delta2_1] = set_computation_2(H,h_2,h_e,A,B,C,D,Gu,Gx,ulb,gamma0_2,K,delta0_2,s_s2);
        figure
        s_s3=Polyhedron([ K -(K*Gx+Gu);H*[zeros(size(C)) (C*Gx+D*Gu)]],[-uub;h_e]);
        s_s3=s_s3.minHRep();
        [oinf_set_3,Gamma3_1, Delta3_1] = set_computation_3(H,h_3,h_e,A,B,C,D,Gu,Gx,uub,gamma0_2,K,delta0_2,s_s3);
        figure
        plot(oinf_set_1, 'color', 'r', oinf_set_2, 'color', 'b',oinf_set_3, 'color', 'b')
        axis equal
        figure 
    end
    [oinf_set1_2,Gamma1_2,Delta1_2]= set_computation_1(H,[Delta2_1;Delta3_1],h_e,Acl,Bcl,C,D,Gu,Gx,uub,ulb,[Gamma2_1(:,1:2);Gamma3_1(:,1:2)],K,[Gamma2_1(:,3);Gamma3_1(:,3)],oinf_set_1);
    axis equal
    figure
    [oinf_set2_2,Gamma2_2,Delta2_2] = set_computation_2(H,[Delta1_1;Delta3_1],h_e,A,B,C,D,Gu,Gx,ulb,[Gamma1_1(:,1:2);Gamma3_1(:,1:2)],K,[Gamma1_1(:,3);Gamma3_1(:,3)],oinf_set_2);
    axis equal
    figure
    [oinf_set3_2,Gamma3_2, Delta3_2] = set_computation_3(H,[Delta1_1;Delta2_1],h_e,A,B,C,D,Gu,Gx,uub,[Gamma1_1(:,1:2);Gamma2_1(:,1:2)],K,[Gamma1_1(:,3);Gamma2_1(:,3)],oinf_set_3);
    axis equal
    figure
     plot(oinf_set1_2, 'color', 'r', oinf_set2_2, 'color', 'b',oinf_set3_2, 'color', 'b')
     axis equal
     figure
     [oinf_set1_3,Gamma1_2,Delta1_2]= set_computation_1(H,[Delta2_2;Delta3_2],h_e,Acl,Bcl,C,D,Gu,Gx,uub,ulb,[Gamma2_2(:,1:2);Gamma3_2(:,1:2)],K,[Gamma2_2(:,3);Gamma3_2(:,3)],oinf_set1_2);
    axis equal
    figure
    [oinf_set2_3,Gamma2_3,Delta2_3] = set_computation_2(H,[Delta1_2;Delta3_2],h_e,A,B,C,D,Gu,Gx,ulb,[Gamma1_2(:,1:2);Gamma3_2(:,1:2)],K,[Gamma1_2(:,3);Gamma3_2(:,3)],oinf_set2_2);
    axis equal
    figure
    [oinf_set3_3,Gamma3_3, Delta3_3] = set_computation_3(H,[Delta1_2;Delta2_2],h_e,A,B,C,D,Gu,Gx,uub,[Gamma1_2(:,1:2);Gamma2_2(:,1:2)],K,[Gamma1_2(:,3);Gamma2_2(:,3)],oinf_set3_2);
    axis equal
    figure
     plot(oinf_set1_3, 'color', 'r', oinf_set2_3, 'color', 'b',oinf_set3_3, 'color', 'b')
     axis equal
     
    i=i+1;
        
end
