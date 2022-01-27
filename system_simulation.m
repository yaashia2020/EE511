%%main loop for the system
global A B C D K x0 r N Gu Gx H h_1 h_e Acl Bcl uub ulb gamma0_1 delta0_1
xk=[];
yk=[];
uk=[];
oinf_set_1= set_computation_1(H,h_1,h_e,Acl,Bcl,C,D,Gu,Gx,uub,ulb,gamma0_1,K,delta0_1)
for k=1:N
    if k==1
        xk(:,k)=x0;
    end
    vk(:,k) = reference_gov(oinf_set_1,xk(:,k),r);
    uk(:,k) = Gu*vk(:,k)-K*(xk(:,k)-Gx*vk(:,k));
    xk(:,k+1)=A*xk(:,k)+B*uk(:,k);
    yk(:,k)=C*xk(:,k)+D*uk(:,k);
end
figure 
plot(xk(1,:));
hold on
plot(Gx(1,:)*vk(1,:));
hold all
plot(xk(2,:));
hold all
plot(Gx(2,:)*vk(1,:));
grid on
figure 
plot(uk(1,:));
figure 
plot(H(1,:)*yk-h_1(1))
hold all
plot(H(2,:)*yk-h_1(2))
hold all
plot(H(3,:)*yk-h_1(3))
hold all
plot(H(4,:)*yk-h_1(4))
hold all
plot(H(5,:)*yk-h_1(5))
hold all
plot(H(6,:)*yk-h_1(6))
hold all
