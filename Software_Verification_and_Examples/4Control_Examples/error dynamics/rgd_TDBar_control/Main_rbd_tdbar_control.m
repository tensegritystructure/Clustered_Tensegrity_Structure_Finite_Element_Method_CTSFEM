%% Main_rbd_tdbar_control
% THis defines all points for morphing.
clear all; clc ; warning off
N=[0 0 0;1 1 0;2 0 0;1 -1 0;0 2 0;0 -2 0]';    
% Manually specify connectivity indices.
C_s_in = [1 3;2 4;2 5;4 6]; C_b_in = [1 2;2 3;3 4;4 1;1 5;1 6];
% Convert the above matrices into full connectivity matrices.
CBT = tenseg_ind2C(C_b_in,N); CST = tenseg_ind2C(C_s_in,N);
% tenseg_plot(N,C_b,C_s); 
%%
% figure(1)
dt = 0.01; n_tar = 2.6; N_cell = cell(1,1/dt);
for  i = 1:1/dt
    x_tem = 0;  y_tem = 0;  N_tem = [];
    x_tem = 1+(n_tar/2-1)*i*dt; y_tem = sqrt(2 - x_tem^2);
    N_tem=[0 0 0;x_tem y_tem 0;2*x_tem 0 0;x_tem -y_tem 0;0 2 0;0 -2 0]';    
    N_cell{:,i} = N_tem;
% tenseg_plot(N_tem,CBT,CST,1); hold off
end
% tenseg_plot(N,CBT,CST);
%%
S_len = []; S_0 = N*CST'; S_len0 = [];
for ss = 1:length(CST(:,1))
    S_len0 = [S_len0;norm(S_0(:,ss),2)];
end
B_len = []; B_0 = N*CBT'; B_len0 = [];
for bb = 1:length(CBT(:,1))
    B_len0 = [B_len0;norm(B_0(:,bb),2)];
end
%%
N_simple = N; pinned_nodes = [1 5 6]; C_b=CBT; C_s=CST;
[N_new,Cb,Cs,P,~] = tenseg_class_k_convert(N_simple,C_b,C_s,pinned_nodes);
D = zeros(3,size(P,2));
D(:,end-length(pinned_nodes)+1:end) = N_simple(:,pinned_nodes);
control_point=1:length(N_simple(1,:)); control_point(pinned_nodes) = [];
L = [1 0 0;0 1 0;0 0 0]; R = zeros(length(N_new(1,:)),length(control_point));
for i =1:length(control_point)
    R(control_point(i),i) = 1;
end
% R = zeros(size(N_new,2),1);
length(N_cell);
tenseg_plot_target_3D(N_cell{1,end},C_b,C_s,N_cell{1,end}')
%%
for i =1:length(N_cell)
    N_cell = cell_part_value(N_cell,[1,i],pinned_nodes,2);
end
N = N_new; B = N*Cb';S=N*Cs';b0=sqrt(diag(diag(B'*B)));
s0=sqrt(diag(diag(S'*S))); W = zeros(3,size(N,2));

[U,Ez,V] = svd(P);
U1 = U(:,1:size(P,2));U2 = U(:,size(P,2)+1:end); E0 = Ez(1:size(P,2),:);
Nd = 0*N; x0 = [N*U2 , Nd*U2]; dt = 0.001; tf =10;
t = 0:dt:tf;
a=2.5; b=4;
Inp.Yt = N_cell; Inp.a=a;Inp.b=b;Inp.Cs = Cs; Inp.Cb=Cb; Inp.W=W; Inp.s0=s0;
Inp.b0=b0;Inp.tf=tf;Inp.P=P;Inp.D=D;Inp.B=B;Inp.S=S;Inp.L=L;Inp.R=R; Inp.tf = tf;
%%
y = ode41(@classkcontrol_ljc_stable,t,x0,Inp);
%%
X1 = D*V*E0^-1;
figure(2)
E=[];ble=zeros(size(Cb,1),size(y,1));
for i = 1:size(y,1)
    X2 = (y(i,1:size(y,2)/2));
    X2 = reshape(X2,3,size(U2,2));
    N = [X1 X2]*U^-1;
    Ntrace(:,:,i) = N;
    Btrace(:,:,i) = N*Cb';
    Strace(:,:,i) = N*Cs';
    RBtrace(:,:,i) = 0.5*N*abs(Cb');
    RStrace(:,:,i) = 0.5*N*abs(Cs');
    ble(:,i)=diag(B'*B)-diag(Btrace(:,:,i)'*Btrace(:,:,i));
    E = [E;(L*N*R-N_cell{1,end})];
end
span = [1:2:(length(t)-1)*2]'+1;
% for i =1:4
% %     figure(i+1)
%     holld = [1 3 4 6];
% %     err(:,i) = E(span,holld(i));
% %     err(:,i) = E(span,i);
%
% %     plot(t(1:length(span))/10,err(:,i)),axis([0 tf/10 -0.05 0.05]),hold on
% end
% str=['a=' mat2str(a) '_b=' mat2str(b) '_TIME=' mat2str(STIME)  '_STEP=' mat2str(SSTEP) '_AMPT='  mat2str(AMPT) ];
% save ([str '.mat']);
% xlabel('Time (sec)')
% ylabel('Error (m)')
% tenseg_plot(N,Cb,Cs)
%%
% cartoon(Ntrace,Cb,Cs,['exw' '.avi'],N_cell{1,end}');    %cartoon of control

%%
% tenseg_cartoon_gif(Ntrace,N,Cb,Cs,['exw' '.avi'],N_cell{1,end}');
tenseg_cartoon_gif(Ntrace,N,Cb,Cs,['exw' '.avi'],N_cell{1,end}');

%%
dt = 0.01;  S_error = []; B_error = [];
for i = 1:size(y,1)
    N_new2 = Ntrace(:,:,i);
    S_len = []; B_len = []; S = N_new2*Cs';
    for ss = 1: length(S(1,:))
        S_len = [S_len;norm(S(:,ss),2)];
    end
    S_len_error = abs(S_len - S_len0);
    %     S_len_error = S_len;
    S_error = [S_error S_len_error];
    B = N_new2*Cb';
    for bb = 1: length(B(1,:))
        B_len = [B_len;norm(B(:,bb),2)];
    end
    B_len_error = (B_len - B_len0); B_error = [B_error B_len_error];
end
%%
figure()
for  i = 1:length(B(1,:))
    plot([1:size(y,1)]/100,B_error(i,:),'-','linewidth',2);
    set(gca,'fontsize', 14,'linewidth',1.2);
    legend_str{i} = ['Bar ' num2str(i)]; hold on;
    %     title('Bar Length Error','Interpreter','latex');
    xlabel('Time t(s)','Interpreter','latex')
    ylabel('Bar Length Error (m)','Interpreter','latex')
        set(gca,'ticklength',1*get(gca,'ticklength'));

end
legend(legend_str);
%%
figure()
for  i = 1:length(S(1,:))
    plot([1:size(y,1)]/100,S_error(i,:),'-','linewidth',2);
    set(gca,'fontsize',14,'linewidth',1.2);
    legend_str{i} = ['String ' num2str(i)]; hold on;
    %     title('String Length Change','Interpreter','latex');
    xlabel('Time t(s)','Interpreter','latex')
    ylabel('String Length Change (m)','Interpreter','latex')
    set(gca,'ticklength',1*get(gca,'ticklength'));

    %     set(gca,'ticklength',1.2*get(gca,'ticklength'));
end
legend(legend_str);
%%
% save Ntrace.mat