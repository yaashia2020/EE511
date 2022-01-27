%%closed loop and system parameters
global A B C D K P yub ylb Acl Bcl Ccl Dcl H h_1 h_e N Gu Gx gamma0_1 delta0_1 gamma0_2 delta0_2
global ulb uub h_2 h_3
%%generalise Gu and Gx
Gu=0;
Gx=[1;0];
e=0.05;

Acl=A-B*K;
Bcl=B*(Gu+K*Gx);
Ccl=(C-D*K);
Dcl=D*(Gu+K*Gx);

H=[eye(length(D));-eye(length(D))];
h_1=[yub;-ylb];
h_2=h_1-ulb*H*D;
h_3=h_1-uub*H*D;
h_e=[yub-e*(yub-ylb);(-ylb -e*(yub-ylb))];
gamma0_1=H*Ccl;
delta0_1=H*Dcl;
gamma0_2=H*C;
delta0_2=zeros(4,1);


N=100;